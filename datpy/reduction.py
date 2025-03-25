import numpy as np
from copy import deepcopy
from .helpers import fort_write
from .constants import MAXF, ULI
from .propagation import propagate_prior_to_dataset


def reduce_dataset(
    dataset, E11, E22, EQ, Q, mxm1, file_IO2
):
    dataset = deepcopy(dataset)
    new_dataset = deepcopy(dataset)
    #  Interpolation type (this is very system specific -
    # does not apply for other simultaneous evaluations
    interp_type = 'lin-lin'
    if (dataset.quantity_type in (1, 2, 5, 8) and
            dataset.reaction_ids[0] not in (2, 5, 10) and
            dataset.reaction_ids[0] <= 10):
        interp_type = 'log-log'

    # GET GRID VALUES  - try at all apriori energies to find data
    new_dataset.num_values = 0
    new_dataset.energies = np.zeros(mxm1-1, dtype=float)
    new_dataset.measured_values = np.zeros(mxm1-1, dtype=float)

    for L in range(1, mxm1):  # 40
        E1 = (EQ[L-1] + EQ[L]) / 2.
        E2 = (EQ[L] + EQ[L+1]) / 2.
        E11 = 0.6 * EQ[L]
        E22 = 1.4 * EQ[L]
        if E1 < E11:
            E1 = E11
        if E2 > E22:
            E2 = E22

        AV = 0.
        WTS = 0.
        NKOT = 0
        for N in range(12):  # 133
            dataset.uncertainties[N, MAXF-1] = 0.

        if E1 > .03:
            interp_type = 'lin-lin'

        # INTERPOLATION CONST.
        if interp_type == 'lin-lin':
            # LIN LIN
            AL = (Q[L-1]-Q[L])/(EQ[L-1]-EQ[L])
            BL = Q[L]-AL*EQ[L]
            AR = (Q[L]-Q[L+1])/(EQ[L]-EQ[L+1])
            BR = Q[L]-AR*EQ[L]
        if interp_type == 'log-log':
            # LOG LOG
            QBL = (np.log(Q[L-1])-np.log(Q[L]))/(np.log(EQ[L])-np.log(EQ[L-1]))
            QAL = Q[L]*(EQ[L]**QBL)
            QBR = (np.log(Q[L])-np.log(Q[L+1]))/(np.log(EQ[L+1])-np.log(EQ[L]))
            QAR = Q[L]*(EQ[L]**QBR)

        # GRID VALUES
        num_values = len(dataset.measured_values)
        for K in range(num_values):  # 35
            if dataset.energies[K] < E1*(1.-1e-5) or dataset.energies[K] >= E2*(1.+1e-5):
                continue

            WT = 1./dataset.uncertainties[11, K]
            WT = WT*WT
            if dataset.energies[K] > EQ[L]:
                # right of energy grid point
                if interp_type == 'lin-lin':
                    ADD = AR*dataset.energies[K] + BR
                    AD = dataset.measured_values[K] + Q[L] - ADD
                elif interp_type == 'log-log':
                    ADD = QAR / (dataset.energies[K]**QBR)
                    AD = dataset.measured_values[K] * Q[L] / ADD

            elif dataset.energies[K] < EQ[L]:
                # left o energy grid point
                if interp_type == 'lin-lin':
                    ADD = AL * dataset.energies[K] + BL
                    AD = dataset.measured_values[K] + Q[L] - ADD
                elif interp_type == 'log-log':
                    ADD = QAL / (dataset.energies[K]**QBL)
                    AD = dataset.measured_values[K] * Q[L] / ADD
            else:
                # same energy as grid point
                AD = dataset.measured_values[K]

            # check if difference is within requested limit of ULI*sigma
            if ULI != 0:
                T1X = 100.*(AD-Q[L])/Q[L]
                T2X = T1X*T1X
                TEST = np.sqrt(WT*T2X)
                if TEST >= ULI:
                    F33 = dataset.uncertainties[2, K] * dataset.uncertainties[2, K]
                    F44 = 1./WT - F33
                    FNEW = T2X / (ULI*ULI)
                    F33N = FNEW - F44
                    dataset.uncertainties[11, K] = np.sqrt(FNEW)
                    dataset.uncertainties[2, K] = np.sqrt(F33N)
                    WT = 1./FNEW

                    format511 = "(20X,' VALUE OUTSIDE ',F5.2,' SIGMA BY ',F10.2)"
                    fort_write(file_IO2, format511, [ULI, TEST])

            AV = AV + AD*WT
            WTS = WTS + WT

            # statistical uncertainty reduces if more than one value contributes,
            # average for all other uncertainties
            for M in range(11):  # 38
                if dataset.NETG[M] != 9:
                    dataset.uncertainties[M, MAXF-1] = dataset.uncertainties[M, MAXF-1] + dataset.uncertainties[M, K]
                elif dataset.uncertainties[M, K] != 0.0:
                    dataset.uncertainties[M, MAXF-1] = dataset.uncertainties[M, MAXF-1] + (1./dataset.uncertainties[M, K])**2

            NKOT = NKOT + 1

        AKOT = NKOT
        if AV != 0.0:
            # GRID VALUE AND OUT
            EEE = EQ[L]
            QQQ = AV / WTS
            DIF = QQQ / Q[L]
            for N in range(11):  # 39
                if dataset.NETG[N] != 9:
                    dataset.uncertainties[N, MAXF-1] = dataset.uncertainties[N, MAXF-1] / AKOT
                elif dataset.uncertainties[N, MAXF-1] > 0.0:
                    dataset.uncertainties[N, MAXF-1] = 1. / np.sqrt(dataset.uncertainties[N, MAXF-1])
                else:
                    dataset.uncertainties[N, MAXF-1] = 0.

            # OUTPUT
            new_dataset.energies[new_dataset.num_values] = EEE
            new_dataset.measured_values[new_dataset.num_values] = QQQ
            new_dataset.uncertainties[0:12, new_dataset.num_values] = dataset.uncertainties[0:12, MAXF-1]
            new_dataset.num_values += 1

            if not hasattr(new_dataset, 'DIF'):
                new_dataset.DIF = []
            new_dataset.DIF.append(DIF)

    new_dataset.uncertainties[0:12, MAXF-1] = dataset.uncertainties[0:12, MAXF-1]
    return new_dataset


def reduce_datablocks(datablocks, reaction_prior, file_IO2):

    reduced_datablocks = []
    for datablock in datablocks:
        reduced_datasets = []
        for dataset in datablock:

            num_values = len(dataset.measured_values)
            E11, E22, EQ, Q, mxm1 = propagate_prior_to_dataset(dataset, reaction_prior)

            if dataset.quantity_type != 6:
                # skip datasets whose energies are byeond limits
                if dataset.energies[0] > E22 or dataset.energies[num_values-1] < E11:
                    continue

            if dataset.quantity_type != 6:
                reduced_dataset = reduce_dataset(
                    dataset, E11, E22, EQ, Q, mxm1, file_IO2
                )
            else:
                # no reduction necessary for fission spectrum average dataset
                reduced_dataset = deepcopy(dataset)

            reduced_datasets.append(reduced_dataset)

        if len(reduced_datasets) > 0:
            reduced_datablocks.append(reduced_datasets)

    return reduced_datablocks


# for tracking function calls
import atexit
from debug import must_be_called

# other python packages
import os
import numpy as np

# input/output
from input_output import (
    copy_gma_controls,
    read_apriori,
    transfer_apriori_to_output_file,
    DATRCL,
)

# helper functions
from helpers import (
    fort_write,
    Bunch,
    find_indices_with_tol,
)
from constants import (
    NQND,
    NQST,
    NOM,
    MAXF,
    ULI,
    SHOULD_TEST_OUTPUT,
)


@must_be_called
def deal_with_CS_and_CS_SHAPE(
    xp, prior_number_points, prior_energy_mesh, prior_cross_section
):
    EQ = np.empty((200,), dtype=float)
    Q = np.zeros((NOM,), dtype=float)
    # CS + CS SHAPE
    M1 = xp.NID[0] - 1
    NON = prior_number_points[M1]
    NON1 = NON-1
    E11 = (prior_energy_mesh[M1, 0] + prior_energy_mesh[M1, 1]) / 2.
    E22 = (prior_energy_mesh[M1, NON-1] + prior_energy_mesh[M1, NON1-1]) / 2.
    for K in range(NON):
        EQ[K] = prior_energy_mesh[M1, K]
        Q[K] = prior_cross_section[M1, K]
    format4611 = "(i6,'  cr. sec. apriori for interp. ')"
    fort_write(None, format4611, [NON-1])
    return E11, E22, EQ, Q, NON-1


@must_be_called
def deal_with_RATIO_and_RATIO_SHAPE(
    xp, prior_number_points, prior_energy_mesh, prior_cross_section
):
    EQ = np.empty((200,), dtype=float)
    Q = np.zeros((NOM,), dtype=float)
    # RATIO + RATIO SHAPE
    M1 = xp.NID[0] - 1
    M2 = xp.NID[1] - 1
    NO1 = prior_number_points[M1]
    NON = prior_number_points[M2]
    mxm = 0

    for K in range(NO1):
        # find matching energies
        L = find_indices_with_tol(prior_energy_mesh[M2,:NON], prior_energy_mesh[M1, K:K+1], atol=1e-8, rtol=1e-4).item()
        if L == -1:  # not found
            continue
        EQ[mxm] = prior_energy_mesh[M1, K]
        if prior_cross_section[M2, L] == 0.:
            format4618 = ":(' apriori constr. ',4i5,f10.7,'*****************')"
            fort_write(None, format4618, [M1-1, K+1, M2-1, L+1, EQ[mxm]])
            exit()

        Q[mxm] = prior_cross_section[M1, K] / prior_cross_section[M2, L]
        mxm += 1

    E11 = (EQ[0] + EQ[1]) / 2.
    E22 = (EQ[mxm-1] + EQ[mxm-2]) / 2.

    format4612 = "(i6,'  ratio apriori for interp. ')"
    fort_write(None, format4612, [mxm])
    return E11, E22, EQ, Q, mxm-1


@must_be_called
def deal_with_SUM_and_SHAPE_OF_SUM(
    xp, prior_number_points, prior_energy_mesh, prior_cross_section
):
    EQ = np.empty((200,), dtype=float)
    Q = np.zeros((NOM,), dtype=float)
    # SUM AND SHAPE OF SUM
    M1 = xp.NID[0] - 1
    M2 = xp.NID[1] - 1
    M3 = xp.NID[2] - 1
    NO1 = prior_number_points[M1]
    NON = prior_number_points[M2]
    if M3 != -1:
        NO3 = prior_number_points[M3]

    mxm = 0
    for K in range(NO1):  # 23

        L = find_indices_with_tol(
            prior_energy_mesh[M2,:NON], prior_energy_mesh[M1, K:K+1], atol=1e-8, rtol=1e-4
        ).item()
        if L == -1:  # not found
            continue

        if M3 != -1:
            J = find_indices_with_tol(
                prior_energy_mesh[M3,:NO3], prior_energy_mesh[M1, K:K+1], atol=1e-8, rtol=1e-4
            ).item()
            if J == -1:  # not found
                continue

        EQ[mxm] = prior_energy_mesh[M1, K]
        Q[mxm] = prior_cross_section[M1, K] + prior_cross_section[M2, L]
        if M3 != -1:
            Q[mxm] = Q[mxm] + prior_cross_section[M3, J]

        mxm += 1 

    E11 = (EQ[0] + EQ[1]) / 2.
    E22 = (EQ[mxm-1] + EQ[mxm-2]) / 2.

    format4613 = "(i6,'  sum apriori for interp. ')"
    fort_write(None, format4613, [mxm])
    return E11, E22, EQ, Q, mxm-1


@must_be_called
def deal_with_CS_VS_SUM_PLUS_SHAPE(
    xp, prior_number_points, prior_energy_mesh, prior_cross_section
):
    EQ = np.empty((200,), dtype=float)
    Q = np.zeros((NOM,), dtype=float)
    #  RATIO OF CS VS. SUM + SHAPE
    M1 = xp.NID[0] - 1
    M2 = xp.NID[1] - 1
    M3 = xp.NID[2] - 1
    NO1 = prior_number_points[M1]
    NON = prior_number_points[M2]
    NO3 = prior_number_points[M3]
    mxm = 0
    pymxm = mxm - 1

    for K in range(NO1):

        L = find_indices_with_tol(
            prior_energy_mesh[M2,:NON], prior_energy_mesh[M1, K:K+1], atol=1e-8, rtol=1e-4
        ).item()
        if L == -1:  # not found
            continue

        if M3 != -1:
            J = find_indices_with_tol(
                prior_energy_mesh[M3,:NO3], prior_energy_mesh[M1, K:K+1], atol=1e-8, rtol=1e-4
            ).item()
            if J == -1:  # not found
                continue

        EQ[mxm] = prior_energy_mesh[M1, K]
        Q[mxm] = (
            prior_cross_section[M1, K] /
            (prior_cross_section[M2, L] + prior_cross_section[M3, J])
        )
        mxm += 1

    mxm1 = mxm - 1
    E11 = (EQ[0] + EQ[1]) / 2.
    E22 = (EQ[pymxm] + EQ[mxm1-1]) / 2.

    format4614 = "(i6,'  cr. sec. vs sum apriori for interp. ')"
    fort_write(None, format4614, [mxm])
    return E11, E22, EQ, Q, mxm1


@must_be_called
def reduce_data():

    NQQA = NQST

    basedir = '.'
    # OPEN(14,FILE='DAT.INP')
    prior_file_handle = open(os.path.join(basedir, 'DAT.INP'), 'r')
    # OPEN(15,FILE='DAT.LST')
    file_IO2 = open(os.path.join(basedir, 'DAT.LST'), 'w')
    # OPEN(12,FILE='GMDATA.CRD')
    expdata_file_handle = open(os.path.join(basedir, 'GMDATA.CRD'), 'r')
    # OPEN(13,FILE='DAT.RES')
    gma_file_handle = open(os.path.join(basedir, 'DAT.RES'), 'w')

    format250 = '(4HEDBL,1X,2I5)'
    # NOTE: format200 and format290 will be used
    #       only much later
    if not SHOULD_TEST_OUTPUT:
        format200 = '(2E10.4,12F5.1)'
        format290 = '(2E10.4,12F5.1,F7.3)'
    else:
        format200 = '(2E14.6,12F12.7)'
        format290 = '(2E14.6,12F12.7,F9.5)'

    copy_gma_controls(prior_file_handle, file_IO2, gma_file_handle)
    prior_number_points, LAB, prior_energy_mesh, prior_cross_section = read_apriori(prior_file_handle)
    transfer_apriori_to_output_file(
        file_IO2, gma_file_handle, prior_number_points, LAB, prior_energy_mesh, prior_cross_section
    )

    # START OF REDUCTION AND TRANSFER
    while True:
        # Bunch allows to access the dictionary elements
        # returned by DATRCL using the syntax expdata.varname
        xp = Bunch(DATRCL(expdata_file_handle, 1, 1))
        if xp.NR == 9999:
            break
        format3733 = "(' read data set  ',i7)"
        fort_write(None, format3733, [xp.NR])

        #  Interpolation type (this is very system specific -
        # does not apply for other simultaneous evaluations
        INT = 1
        if (xp.NT in (1, 2, 5, 8) and
                xp.NID[0] not in (2, 5, 10) and
                xp.NID[0] <= 10):
            INT = 2

        # CONSTRUCT APRIORI

        # NOTE: computed goto of fortran replaced
        #       by if-else statements
        if xp.NT in (1, 2):
            E11, E22, EQ, Q, mxm1 = deal_with_CS_and_CS_SHAPE(
                xp, prior_number_points, prior_energy_mesh, prior_cross_section
            )
        elif xp.NT in (3, 4):
            E11, E22, EQ, Q, mxm1 = deal_with_RATIO_and_RATIO_SHAPE(
                xp, prior_number_points, prior_energy_mesh, prior_cross_section
            )
        elif xp.NT in (5, 8):
            E11, E22, EQ, Q, mxm1 = deal_with_SUM_and_SHAPE_OF_SUM(
                xp, prior_number_points, prior_energy_mesh, prior_cross_section
            )
        elif xp.NT in (7, 9):
            E11, E22, EQ, Q, mxm1 = deal_with_CS_VS_SUM_PLUS_SHAPE(
                xp, prior_number_points, prior_energy_mesh, prior_cross_section
            )
        assert xp.NT >= 1 and xp.NT <= 9

        # REDUCTION

        if xp.NT != 6:
            # FIND USEFUL DATA RANGE
            if xp.E[0] > E22 or xp.E[xp.NO-1] < E11:
                # out of range
                if NQQA != NQND and xp.NQQ == NQND:
                    NQQA = xp.NQQ
                continue

        NQMM = 0
        if NQQA == NQND:
            fort_write(gma_file_handle, format250, [NQMM, NQMM])
            fort_write(file_IO2, format250, [NQMM, NQMM])

        if NQQA == NQND or NQQA == NQST:
            format251 = '(4HBLCK,1X,2I5)'
            fort_write(gma_file_handle, format251, [NQMM, NQMM])
            fort_write(file_IO2, format251, [NQMM, NQMM])

        NQQA = xp.NQQ

        # find nr of CS involved
        for MN in range(1, 5):  # 66
            pyMN = MN - 1
            if xp.NID[pyMN] == 0:
                break

        NNN = MN - 1
        xp.NID[3] = 0

        # special marker for thermal constants
        if xp.NR >= 910 and xp.NR <= 934:
            xp.NID[3] = 1

        format253 = '(5HDATA ,9I5)'
        fort_write(gma_file_handle, format253, [xp.NR, xp.NT, xp.NCO, NNN, xp.NID[0:4]])
        fort_write(file_IO2, format253, [xp.NR, xp.NT, xp.NCO, NNN, xp.NID[0:4]])

        format254 = '(3I5,A28,8X,A20)'
        fort_write(gma_file_handle, format254, [xp.NY, xp.NQ, xp.NCS, xp.NAU, xp.NREF])
        fort_write(file_IO2, format254, [xp.NY, xp.NQ, xp.NCS, xp.NAU, xp.NREF])

        if xp.NT not in (2, 4, 8, 9):
            # normalization uncertainties
            format261 = '(10F5.1,10I3)'
            fort_write(gma_file_handle, format261, [xp.ENF[0:10], xp.NENF[0:10]])
            fort_write(file_IO2, format261, [xp.ENF[0:10], xp.NENF[0:10]])

        # energy dep. unc. parameters
        format262 = '(3F5.2,I3)'
        # NOTE: this loop is implicit in the write statement
        #       of the fortran code
        for K in range(11):
            fort_write(file_IO2, format262, [xp.EPA[0:3, K], xp.NETG[K]])
            fort_write(gma_file_handle, format262, [xp.EPA[0:3, K], xp.NETG[K]])
        if xp.NCS != 0:
            # cross correlations
            format263 = '(I5,20I3)'
            format293 = '(10F5.1)'
            for K in range(xp.NCS):  # 83
                # NOTE: during flattening in fort_write first index should
                #       change fastest
                fort_write(gma_file_handle, format263, [xp.NCST[K], xp.NEC[:, :, K]])
                fort_write(file_IO2, format263, [xp.NCST[K], xp.NEC[:, :, K]])
                fort_write(gma_file_handle, format293, [xp.FCFC[0:10, K]])
                fort_write(file_IO2, format293, [xp.FCFC[0:10, K]])

        if xp.NT == 6:
            # fission spectrum average data set
            fort_write(file_IO2, format200, [xp.E[0], xp.S[0], xp.F[0:12, 0]])
            fort_write(gma_file_handle, format200, [xp.E[0], xp.S[0], xp.F[0:12, 0]])
            fort_write(file_IO2, format200, [0, 0, xp.F[0:12, MAXF-1]])
            fort_write(gma_file_handle, format200, [0, 0, xp.F[0:12, MAXF-1]])
            continue

        # GET GRID VALUES  - try at all apriori energies to find data

        format5173 = "(/' ENERGY/MEV  VALUE       UNCERTAINTIES                     RATIO TO APRIORI'/)"
        fort_write(file_IO2, format5173, [None])

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
                xp.F[N, MAXF-1] = 0.

            if E1 > .03:
                INT = 1

            # INTERPOLATION CONST.
            if INT == 1:
                # LIN LIN
                AL = (Q[L-1]-Q[L])/(EQ[L-1]-EQ[L])
                BL = Q[L]-AL*EQ[L]
                AR = (Q[L]-Q[L+1])/(EQ[L]-EQ[L+1])
                BR = Q[L]-AR*EQ[L]
            if INT == 2:
                # LOG LOG
                QBL = (np.log(Q[L-1])-np.log(Q[L]))/(np.log(EQ[L])-np.log(EQ[L-1]))
                QAL = Q[L]*(EQ[L]**QBL)
                QBR = (np.log(Q[L])-np.log(Q[L+1]))/(np.log(EQ[L+1])-np.log(EQ[L]))
                QAR = Q[L]*(EQ[L]**QBR)

            # GRID VALUES
            for K in range(xp.NO):  # 35
                if xp.E[K] < E1*(1.-1e-5) or xp.E[K] >= E2*(1.+1e-5):
                    continue

                WT = 1./xp.F[11, K]
                WT = WT*WT
                if xp.E[K] > EQ[L]:
                    # right of energy grid point
                    if INT == 1:
                        ADD = AR*xp.E[K] + BR
                        AD = xp.S[K] + Q[L] - ADD
                    elif INT == 2:
                        ADD = QAR / (xp.E[K]**QBR)
                        AD = xp.S[K] * Q[L] / ADD

                elif xp.E[K] < EQ[L]:
                    # left o energy grid point
                    if INT == 1:
                        ADD = AL * xp.E[K] + BL
                        AD = xp.S[K] + Q[L] - ADD
                    elif INT == 2:
                        ADD = QAL / (xp.E[K]**QBL)
                        AD = xp.S[K] * Q[L] / ADD
                else:
                    # same energy as grid point
                    AD = xp.S[K]

                # check if difference is within requested limit of ULI*sigma
                if ULI != 0:
                    T1X = 100.*(AD-Q[L])/Q[L]
                    T2X = T1X*T1X
                    TEST = np.sqrt(WT*T2X)
                    if TEST >= ULI:
                        F33 = xp.F[2, K] * xp.F[2, K]
                        F44 = 1./WT - F33
                        FNEW = T2X / (ULI*ULI)
                        F33N = FNEW - F44
                        xp.F[11, K] = np.sqrt(FNEW)
                        xp.F[2, K] = np.sqrt(F33N)
                        WT = 1./FNEW

                        format511 = "(20X,' VALUE OUTSIDE ',F5.2,' SIGMA BY ',F10.2)"
                        fort_write(file_IO2, format511, [ULI, TEST])

                AV = AV + AD*WT
                WTS = WTS + WT

                # statistical uncertainty reduces if more than one value contributes,
                # average for all other uncertainties
                for M in range(11):  # 38
                    if xp.NETG[M] != 9:
                        xp.F[M, MAXF-1] = xp.F[M, MAXF-1] + xp.F[M, K]
                    elif xp.F[M, K] != 0.0:
                        xp.F[M, MAXF-1] = xp.F[M, MAXF-1] + (1./xp.F[M, K])**2

                NKOT = NKOT + 1

            AKOT = NKOT
            if AV != 0.0:
                # GRID VALUE AND OUT
                EEE = EQ[L]
                QQQ = AV / WTS
                DIF = QQQ / Q[L]
                for N in range(11):  # 39
                    if xp.NETG[N] != 9:
                        xp.F[N, MAXF-1] = xp.F[N, MAXF-1] / AKOT
                    elif xp.F[N, MAXF-1] > 0.0:
                        xp.F[N, MAXF-1] = 1. / np.sqrt(xp.F[N, MAXF-1])
                    else:
                        xp.F[N, MAXF-1] = 0.

                # OUTPUT
                fort_write(gma_file_handle, format200, [EEE, QQQ, xp.F[0:12, MAXF-1]])
                fort_write(file_IO2, format290, [EEE, QQQ, xp.F[0:12, MAXF-1], DIF])

        # end of data set
        fort_write(gma_file_handle, format200, [0, 0, xp.F[0:12, MAXF-1]])
        fort_write(file_IO2, format200, [0, 0, xp.F[0:12, MAXF-1]])

        if xp.NCO == 0:
            continue

        format6114 = '(1X,10F8.5)'
        format6115 = '(10F8.5)'
        for KL in range(xp.NCO):  # 6113
            fort_write(file_IO2, format6114, [xp.ECOR[KL, :(KL+1)]])
            fort_write(gma_file_handle, format6115, [xp.ECOR[KL, :(KL+1)]])

    # DATA FILE COMPLETE
    format256 = '(4HEND*,1X,2I5)'
    fort_write(gma_file_handle, format250, [NQMM, NQMM])
    fort_write(file_IO2, format250, [NQMM, NQMM])
    fort_write(gma_file_handle, format256, [NQMM, NQMM])
    fort_write(file_IO2, format256, [NQMM, NQMM])

    prior_file_handle.close()
    file_IO2.close()
    expdata_file_handle.close()
    gma_file_handle.close()


atexit.register(must_be_called.check_called)
reduce_data()

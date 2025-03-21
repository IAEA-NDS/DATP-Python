import numpy as np
import math
from .helpers import (
    fort_read,
    fort_write,
    unflatten,
    Bunch,
)
from .constants import (
    MAX_NUM_REACTIONS,
    MAX_NUM_POINTS,
    END_DATA_BLOCK_INDICATION_STRING,
    END_DATASET_INDICATION_STRING,
    MAXF,
)


def read_apriori(prior_file_handle):
    prior_energy_mesh = np.zeros((MAX_NUM_REACTIONS, MAX_NUM_POINTS), dtype=float)
    prior_cross_section = np.zeros((MAX_NUM_REACTIONS, MAX_NUM_POINTS), dtype=float)
    prior_label = np.empty((MAX_NUM_REACTIONS,), dtype=object)
    prior_number_points = np.zeros((35,), dtype=int)

    format99 = '(A16)'  # original: (8A2)
    format100r = '(2E10.4)'
    num_reactions = 0
    for L in range(MAX_NUM_REACTIONS):
        cur_label = fort_read(prior_file_handle, format99)[0]
        # exit loop if no more prior reactions to read
        if cur_label.strip() == '':
            break

        prior_label[L] = cur_label
        num_points = 0
        for K in range(MAX_NUM_POINTS):
            EQ9, TQ9 = fort_read(prior_file_handle, format100r)
            if EQ9 == 0:
                break
            prior_energy_mesh[L, K] = EQ9
            prior_cross_section[L, K] = TQ9
            num_points += 1

        prior_number_points[L] = num_points
        num_reactions += 1

    # shrink arrays to real data size
    return {
        'number_points': prior_number_points[:num_reactions],
        'energy_mesh': prior_energy_mesh[:num_reactions,:],
        'cross_section': prior_cross_section[:num_reactions,:],
        'label': prior_label[:num_reactions],
    }


def read_dataset(expdata_file_handle):

    # variables with local scope
    NAU: str; NREF: str; NQT: str
    NCOM: str; NXQT: str; NXAU: str

    # this declaration is not present in Fortran code
    # but assumed to be implicitly done
    SES = 0.

    # data set identification
    # original string
    # format100 = '(2I4,12A2,14A2,10A2)'
    format100 = '(2I4, A24, A28, A20)'
    NR = 0
    while NR == 0:
        NR, NY, NQT, NAU, NREF = fort_read(expdata_file_handle, format100)

    if NR == 9999:
        dataset = {'dataset_id': NR}
        datablock_complete = True
        return dataset, datablock_complete

    format103 = '(4I2,I3,I5,5I3)'
    NQ, NT, NCO, NCS, NCCO, NO, NID = unflatten(
            fort_read(expdata_file_handle, format103), [6, [5]])

    # COMMENTS
    # original: (40A2)
    format106 = '(A80)'
    NCOM = []
    for i in range(NCCO):
        NCOM.append(fort_read(expdata_file_handle, format106))

    # NORMALIZATION UNCERTAINTIES
    ENF = None
    NENF = None
    SES = 0.

    if (not (NT == 2 or NT == 4)) and (not (NT == 8 or NT == 9)):
        format107 = '(10F5.1, 10I3)'
        ENF, NENF = unflatten(fort_read(expdata_file_handle, format107), [[10], [10]])
        for K in range(10):
            SES = SES + ENF[K]*ENF[K]

    # ENERGY DEPENDENT UNCERTAINTY CORRELATIONS PARAMETERS AND TAGS
    format110 = '(3F5.2)'
    EPA = np.empty((3,11), dtype=float)
    for i in range(11):
        EPA[:,i] = fort_read(expdata_file_handle, format110)

    for k in range(11):
        absum = EPA[0,k] + EPA[1,k]
        if absum > 1.0:
            EPA[1,k] = 1.0 - EPA[0,k]

    format111 = '(11I3)'
    NETG = fort_read(expdata_file_handle, format111)

    # DATA
    E = np.empty((NO,), dtype=float)
    S = np.empty((NO,), dtype=float)
    F = np.zeros((12, MAXF), dtype=float)
    format114 = '(2E10.4,12F5.1)'
    for K in range(NO):
        E[K], S[K], F[:,K] = unflatten(
                fort_read(expdata_file_handle, format114), [2, [12]])
        SSS = 0.
        for M in range(2, 11):
            SSS = SSS + F[M, K]*F[M, K]

        F[11, K] = np.sqrt(SES+SSS)

    # CORRELATIONS WITH PRECEDING DATA SETS

    # line in fortran code not required here
    # if no cross-correlations present,
    # respective arrays will be empty
    #if NCS == 0: goto .lbl29
    NCST = np.zeros((NCS,), dtype=int)
    NEC = np.zeros((2,10,NCS), dtype=int)
    FCFC = np.zeros((10,NCS), dtype=float)
    for K in range(NCS):
        format116 = '(I5,20I2)'
        tmp = fort_read(expdata_file_handle, format116)
        # ISSUE: there are not always 20 I2 numbers
        #        in the GMDATA file but sometimes less
        tmp = [x for x in tmp if x is not None]
        tmp2 = np.zeros((20,), dtype=int)
        tmp2[:(len(tmp)-1)] = tmp[1:]
        NCST[K] = tmp[0]
        NEC[0, :, K] = tmp2[:10]
        NEC[1, :, K] = tmp2[10:]

        format452 = '(10F5.1)'
        tmp = fort_read(expdata_file_handle, format452, none_as=0.)
        if np.any([math.isnan(x) for x in tmp]):
            raise ValueError

        FCFC[:, K] = tmp
        for ji in range(10):
            if FCFC[ji, K] > 1.0:
                FCFC[ji, K] = 1.0
            if FCFC[ji, K] < -1.0:
                FCFC[ji, K] = -1.0

    # CORRELATION MATRIX INPUT
    ECOR = np.zeros((NCO, NCO), dtype=float)
    format117 = '(10F8.5)'
    for L in range(NCO):
        num_el_read = 0
        num_el_desired = L + 1
        res = []
        while num_el_read < num_el_desired:
            tmp = fort_read(expdata_file_handle, format117)
            tmp = [x for x in tmp if x is not None]
            res += tmp
            num_el_read += len(tmp)
        ECOR[L, :(L+1)] = res

    # special marker for thermal constants
    NID[3] = 1 if NR >= 910 and NR <= 934 else 0

    # read datablock/dataset termination indicator
    format118 = '(A2)'
    end_indicator = fort_read(expdata_file_handle, format118)[0]
    if end_indicator not in (END_DATASET_INDICATION_STRING, END_DATA_BLOCK_INDICATION_STRING):
        raise ValueError('Expected End-of-Datablock or End-of-Dataset indicator')
    datablock_complete = (end_indicator == END_DATA_BLOCK_INDICATION_STRING)

    dataset = {
        'dataset_id': NR,
        'year': NY,
        'quantity_name': NQT,
        'author': NAU,
        'pubref': NREF,
        'NQ': NQ,
        'quantity_type': NT,
        'cormat_dim': NCO,
        'num_comments': NCCO,
        'comments': NCOM,
        'NO': NO,
        'reaction_ids': NID,
        'ENF': ENF,
        'NENF': NENF,
        'SES': SES,
        'EPA': EPA,
        'NETG': NETG,
        'energies': E,
        'measured_values': S,
        'uncertainties': F,
        'NCST': NCST,
        'NCST_size': NCS, 
        'NEC': NEC,
        'FCFC': FCFC,
        'cormat': ECOR
    }

    return dataset, datablock_complete


def read_datablocks(expdata_file_handle):
    datablocks = []
    datasets = []
    while True:
        dataset, datablock_complete = read_dataset(expdata_file_handle)
        dataset = Bunch(dataset)
        if dataset.dataset_id == 9999:
            break

        format3733 = "(' read data set  ',i7)"
        fort_write(None, format3733, [dataset.dataset_id])
        datasets.append(dataset)

        if datablock_complete:
            datablocks.append(datasets)
            datasets = []

    if len(datasets) > 0:
        raise ValueError(
            'Encountered incomplete datablock at end of file. '
            'Termination suffix {END_DATA_BLOCK_INDICATION_STRING} missing'
        )
    return datablocks

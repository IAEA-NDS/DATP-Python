import numpy as np
import math
from .helpers import (
    fort_read,
    fort_write,
    unflatten,
)
from .constants import (
    MAX_NUM_REACTIONS,
    MAX_NUM_POINTS,
    END_DATA_BLOCK_INDICATION_STRING,
    ELIMINATION_BLOCK_INDICATION_STRING,
    DOWNWEIGHT_BLOCK_INDICATION_STRING,
    FISSION_SPECTRUM_BLOCK_INDICATION_STRING,
    MAXF,
    SHOULD_TEST_OUTPUT,
)


def copy_gma_controls(prior_file_handle, file_IO2, gma_file_handle):
    for K in range(10):
        format260 = '(A2,A2,A1,8I5)'
        format483 = "(' reading  ', a2)"
        KCO1, KCO2, KCO3, MC1, MC2, MC3, MC4, MC5, MC6, MC7, MC8 = \
            fort_read(prior_file_handle, format260, none_as=0.)
        fort_write(None, format483, [KCO1])

        # exit loop if nothing more to read
        if KCO1.strip() == '':
            break

        fort_write(gma_file_handle,  format260,
                   [KCO1, KCO2, KCO3, MC1, MC2, MC3, MC4, MC5, MC6, MC7, MC8])
        fort_write(file_IO2,  format260,
                   [KCO1, KCO2, KCO3, MC1, MC2, MC3, MC4, MC5, MC6, MC7, MC8])
        if KCO1 == DOWNWEIGHT_BLOCK_INDICATION_STRING and MC2 == 10:
            # data set numbers selected for downweighting
            while True:
                format402 = '16I5'
                MSEP = np.empty((16,), dtype=int)
                MSEP[:] = fort_read(prior_file_handle, format402)
                fort_write(file_IO2, format402, [MSEP])
                fort_write(gma_file_handle, format402, [MSEP])
                if MSEP[0] == 0:
                    break

        elif KCO1 == FISSION_SPECTRUM_BLOCK_INDICATION_STRING and MC1 != 0:
            # fission spectrum
            while True:
                format404 = '(2E13.5)'
                AE, BS = fort_read(prior_file_handle, format404, none_as=0.)
                fort_write(file_IO2, format404, [AE, BS])
                fort_write(gma_file_handle, format404, [AE, BS])
                if AE == 0.0:
                    break

        elif KCO1 == ELIMINATION_BLOCK_INDICATION_STRING:
            format408 = '(16i5)'
            format468 = "('Data Sets to be Excluded')"
            NEXL = fort_read(prior_file_handle, format408)
            fort_write(gma_file_handle, format408, NEXL)
            fort_write(file_IO2, format468, [None])
            fort_write(file_IO2, format408, NEXL)


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
    prior_number_points = prior_number_points[:num_reactions]
    prior_energy_mesh = prior_energy_mesh[:num_reactions,:]
    prior_cross_section = prior_cross_section[:num_reactions,:]
    prior_label = prior_label[:num_reactions]

    return prior_number_points, prior_label, prior_energy_mesh, prior_cross_section


def transfer_apriori_to_output_file(

    file_IO2, gma_file_handle, prior_number_points, prior_label, prior_energy_mesh, prior_cross_section
):
    num_reactions = prior_number_points.shape[0]
    ITOT = 0
    for K in range(num_reactions):
        ITOT = ITOT + prior_number_points[K]
    ITOT = ITOT - 2*num_reactions

    format264 = '(5HAPRI ,2I5)'
    fort_write(gma_file_handle, format264, [ITOT, num_reactions])
    fort_write(file_IO2, format264, [ITOT, num_reactions])

    format3731 = "('cross section',i5,' number',i5,A16)"
    format101 = '(2X,I4,2X,2E12.4)'
    format109 = '(8x,2e10.4)'
    format100 = '(2E14.6)' if SHOULD_TEST_OUTPUT else '(2E10.4)'
    format99 = '(A16)'  # original: (8A2)
    for L in range(num_reactions):
        NOR = prior_number_points[L] - 1
        NOR2 = NOR - 1
        fort_write(gma_file_handle, format99, [prior_label[L]])
        fort_write(file_IO2, format99, [prior_label[L]])
        fort_write(None, format3731, [L+1, NOR2, prior_label[L]])

        for K in range(1, NOR):
            fort_write(file_IO2, format101, [K, prior_energy_mesh[L, K], prior_cross_section[L, K]])
            fort_write(gma_file_handle, format100, [prior_energy_mesh[L, K], prior_cross_section[L, K]])
        fort_write(file_IO2, format109, [0, 0])
        fort_write(gma_file_handle, format100, [0, 0])


def read_dataset(expdata_file_handle, NZ: int, IBZ: int):

    # variables with local scope
    NAU: str; NREF: str; NQT: str
    NCOM: str; NQQ: str; NXQT: str
    NXAU: str; ICC: str; NES: str; NEB: str

    ICC = 'C '
    NES = 'ES'
    NEB = END_DATA_BLOCK_INDICATION_STRING

    # this declaration is not present in Fortran code
    # but assumed to be implicitly done
    SES = 0.

    if NZ == 5:
        for K in range(1,1000):
            NXY[K] = 0

    NBQZ = 1
    if NBQZ == 1: NQQ = NES
    if IBZ == 2: NQQ = NEB

    # data set identification
    # original string
    # format100 = '(2I4,12A2,14A2,10A2)'
    format100 = '(2I4, A24, A28, A20)'
    NR = 0
    while NR == 0:
        NR, NY, NQT, NAU, NREF = fort_read(expdata_file_handle, format100)

    if NR == 9999:
        return {'NR': NR}
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

    format118 = '(A2)'
    NQQ = fort_read(expdata_file_handle, format118)
    assert len(NQQ) == 1
    NQQ = NQQ[0]

    return({'NR': NR, 'NY': NY, 'NQT': NQT, 'NAU': NAU, 'NREF': NREF,
            'NQ': NQ, 'NT': NT, 'NCO': NCO, 'NCS': NCS, 'NCCO': NCCO, 'NO': NO, 'NID': NID,
            'NCOM': NCOM, 'ENF': ENF, 'NENF': NENF, 'SES': SES, 'EPA': EPA,
            'NETG': NETG, 'E': E, 'S': S, 'F': F,
            'NCST': NCST, 'NEC': NEC, 'NQQ': NQQ, 'FCFC': FCFC, 'ECOR': ECOR})

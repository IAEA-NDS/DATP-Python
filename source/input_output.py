import numpy as np
from helpers import (
    fort_read,
    fort_write,
)
from constants import (
    NQM,
    NOM,
    NHEL,
    NHMO,
    NHFI,
    MTY,
    SHOULD_TEST_OUTPUT,
)


def copy_gma_controls(file_IO1, file_IO2, file_IO4):
    for K in range(10):
        format260 = '(A2,A2,A1,8I5)'
        format483 = "(' reading  ', a2)"
        KCO1, KCO2, KCO3, MC1, MC2, MC3, MC4, MC5, MC6, MC7, MC8 = \
            fort_read(file_IO1, format260, none_as=0.)
        fort_write(None, format483, [KCO1])
        if KCO1 == MTY:
            break
        fort_write(file_IO4,  format260,
                   [KCO1, KCO2, KCO3, MC1, MC2, MC3, MC4, MC5, MC6, MC7, MC8])
        fort_write(file_IO2,  format260,
                   [KCO1, KCO2, KCO3, MC1, MC2, MC3, MC4, MC5, MC6, MC7, MC8])
        if KCO1 == NHMO and MC2 == 10:
            # data set numbers selected for downweighting
            while True:
                format402 = '16I5'
                MSEP = np.empty((16,), dtype=int)
                MSEP[:] = fort_read(file_IO1, format402)
                fort_write(file_IO2, format402, [MSEP])
                fort_write(file_IO4, format402, [MSEP])
                if MSEP[0] == 0:
                    break

        elif KCO1 == NHFI and MC1 != 0:
            # fission spectrum
            while True:
                format404 = '(2E13.5)'
                AE, BS = fort_read(file_IO1, format404, none_as=0.)
                fort_write(file_IO2, format404, [AE, BS])
                fort_write(file_IO4, format404, [AE, BS])
                if AE == 0.0:
                    break

        elif KCO1 == NHEL:
            format408 = '(16i5)'
            format468 = "('Data Sets to be Excluded')"
            NEXL = fort_read(file_IO1, format408)
            fort_write(file_IO4, format408, NEXL)
            fort_write(file_IO2, format468, [None])
            fort_write(file_IO2, format408, NEXL)


def read_apriori(file_IO1):
    ER = np.zeros((NQM, NOM), dtype=float)
    T = np.zeros((NQM, NOM), dtype=float)
    format99 = '(A16)'  # original: (8A2)
    format100r = '(2E10.4)'
    LAB = np.empty((NQM,), dtype=object)
    NOD = np.zeros((35,), dtype=int)
    for L in range(NQM):
        LAB[L] = fort_read(file_IO1, format99)

        for K in range(1, NOM+1):
            EQ9, TQ9 = fort_read(file_IO1, format100r)
            if EQ9 == 0:
                NOD[L] = K-1
                break
            ER[L, K-1] = EQ9
            T[L, K-1] = TQ9

        if K == NOM:
            NOD[L] = NOM

    return NOD, LAB, ER, T


def transfer_apriori_to_output_file(file_IO2, file_IO4, NOD, LAB, ER, T):
    ITOT = 0
    for K in range(NQM):
        ITOT = ITOT + NOD[K]
    ITOT = ITOT - 2*NQM

    format264 = '(5HAPRI ,2I5)'
    fort_write(file_IO4, format264, [ITOT, NQM])
    fort_write(file_IO2, format264, [ITOT, NQM])

    format3731 = "('cross section',i5,' number',i5,A16)"
    format101 = '(2X,I4,2X,2E12.4)'
    format109 = '(8x,2e10.4)'
    format100 = '(2E14.6)' if SHOULD_TEST_OUTPUT else '(2E10.4)'
    format99 = '(A16)'  # original: (8A2)
    for L in range(NQM):
        NOR = NOD[L] - 1
        NOR2 = NOR - 1
        fort_write(file_IO4, format99, [LAB[L]])
        fort_write(file_IO2, format99, [LAB[L]])
        fort_write(None, format3731, [L+1, NOR2, LAB[L]])

        for K in range(1, NOR):
            fort_write(file_IO2, format101, [K, ER[L, K], T[L, K]])
            fort_write(file_IO4, format100, [ER[L, K], T[L, K]])
        fort_write(file_IO2, format109, [0, 0])
        fort_write(file_IO4, format100, [0, 0])

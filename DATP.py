# helpful for tentative fortran to python conversion
from goto import with_goto
from fortranformat import FortranRecordReader
from fortranformat import FortranRecordWriter
# other python packages
import os
import numpy as np
import time
import math


def debugout(debugstr):
    print('DEBUG: ' + str(debugstr))

# classes goto and label are soley
# defined to suppress linting messages
# goto and label statements are handled
# by the goto package


class goto:
    dummy = 0


class label:
    dummy = 0


class Bunch(object):
    def __init__(self, adict):
        self.__dict__.update(adict)


def fort_read(fobj, formatstr, none_as=None, debug=False):

    frr = FortranRecordReader(formatstr)
    if not isinstance(fobj, str):
        fname = fobj.name
        inpline = fobj.readline()
    else:
        fname = 'console'
        inpline = fobj

    res = frr.read(inpline)
    if none_as is not None:
        res = [none_as if x is None else x for x in res]

    if debug:
        print('--- reading ---')
        print('file: ' + fname)
        print('fmt: ' + formatstr)
        print('str: ' + inpline)

    return res


def fort_write(fobj, formatstr, values, debug=False):
    vals = list(flatten(values))
    vals = [v for v in vals if v is not None]
    if debug:
        print('--- writing ---')
        try:
            print('file: ' + fobj.name)
        except AttributeError:
            print('file: console')
        print('fmt: ' + formatstr)
        print('values: ')
        print(vals)
    frw = FortranRecordWriter(formatstr)
    line = frw.write(vals)
    if fobj is None:
        print(line)
    else:
        fobj.write(line + '\n')


def flatten(container):
    for i in container:
        if isinstance(i, (list, tuple, np.ndarray)):
            for j in flatten(i):
                yield j
        else:
            yield i


def unflatten(flat_list, skeleton):
    it = iter(flat_list)

    def rec(skel):
        result = []
        for cursk in skel:
            if isinstance(cursk, int):
                for i in range(cursk):
                    result.append(next(it))
            else:
                result.append(rec(cursk))
        return result

    result = rec(skeleton)
    # check if all elements traversed
    try:
        next(it)
        raise IndexError
    except StopIteration:
        return result

# global vars


MAXF = 900
pyMAXF = MAXF - 1  # helper variable for Python indexing
NOM = 200
NQM = 33  # number of apriori tables in DAT.INP (excluding FIS*)
W = 0.

# LIMIT OF DIF/UNC FOR OUTLIERS?
ULI = 3.0


@with_goto
def reduce_data():

    NQND = 'EB'
    NHEL = 'EL'
    NQST = 'ST'
    NHMO = 'MO'
    NHFI = 'FI'
    MTY = '  '
    NQQA = NQST

    basedir = '.'
    # OPEN(14,FILE='DAT.INP')
    file_IO1 = open(os.path.join(basedir, 'DAT.INP'), 'r')
    # OPEN(15,FILE='DAT.LST')
    file_IO2 = open(os.path.join(basedir, 'DAT.LST'), 'w')
    # OPEN(12,FILE='GMDATA.CRD')
    file_IO3 = open(os.path.join(basedir, 'GMDATA.CRD'), 'r')
    file_ID3 = file_IO3
    # OPEN(13,FILE='DAT.RES')
    file_IO4 = open(os.path.join(basedir, 'DAT.RES'), 'w')

    format250 = '(4HEDBL,1X,2I5)'

    # COPY OF CONTROLS FOR GMA
    for K in range(10):
        format260 = '(A2,A2,A1,8I5)'
        format483 = "(' reading  ', a2)"
        KCO1, KCO2, KCO3, MC1, MC2, MC3, MC4, MC5, MC6, MC7, MC8 = \
            fort_read(file_IO1, format260, none_as=0.)
        fort_write(None, format483, [KCO1])
        if KCO1 == MTY:
            break  # goto .lbl84
        fort_write(file_IO4,  format260,
                   [KCO1, KCO2, KCO3, MC1, MC2, MC3, MC4, MC5, MC6, MC7, MC8])
        fort_write(file_IO2,  format260,
                   [KCO1, KCO2, KCO3, MC1, MC2, MC3, MC4, MC5, MC6, MC7, MC8])
        if KCO1 == NHMO and MC2 == 10:
            goto .lbl400
        goto .lbl401

        # data set numbers selected for downweighting
        label .lbl400
        format402 = '16I5'
        MSEP = np.empty((16,), dtype=int)
        MSEP[:] = fort_read(file_IO1, format402)
        fort_write(file_IO2, format402, [MSEP])
        fort_write(file_IO4, format402, [MSEP])
        if MSEP[0] == 0:
            goto .lbl407
        goto .lbl400

        label .lbl401
        if KCO1 == NHFI and MC1 != 0:
            goto .lbl403
        goto .lbl405

        # fission spectrum
        label .lbl403
        format404 = '(2E13.5)'
        AE, BS = fort_read(file_IO1, format404, none_as=0.)
        fort_write(file_IO2, format404, [AE, BS])
        fort_write(file_IO4, format404, [AE, BS])
        if AE == 0.0:
            goto .lbl407
        goto .lbl403

        label .lbl405
        if KCO1 == NHEL:
            goto .lbl406
        goto .lbl407

        # data set numbers for exclusion
        label .lbl406
        format408 = '(16i5)'
        format468 = "('Data Sets to be Excluded')"
        JNEX = MC1
        NEXL = fort_read(file_IO1, format408)
        fort_write(file_IO4, format408, NEXL)
        fort_write(file_IO2, format468, [None])
        fort_write(file_IO2, format408, NEXL)

        label .lbl407
        label .lbl80
    label .lbl84  # terminal label of FOR loop

    # CLEAR
    Q = np.zeros((NOM,), dtype=float)
    ER = np.zeros((NQM, NOM), dtype=float)
    T = np.zeros((NQM, NOM), dtype=float)

    # READ APRIORI
    format99 = '(A16)'  # original: (8A2)
    format100 = '(2E10.4)'
    LAB = np.empty((NQM,), dtype=object)
    NOD = np.zeros((35,), dtype=int)
    for L in range(NQM):
        LAB[L] = fort_read(file_IO1, format99)

        for K in range(1, NOM+1):
            EQ9, TQ9 = fort_read(file_IO1, format100)
            if EQ9 == 0:
                goto .lbl3
            ER[L, K-1] = EQ9
            T[L, K-1] = TQ9

        label .lbl2
        NOD[L] = NOM
        goto .lbl4

        label .lbl3
        
        NOD[L] = K-1
        label .lbl4
    label .lbl1   # terminal label of FOR loop

    # TRANSFER APRIORI TO OUTPUT FILE
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
    for L in range(NQM):
        NOR = NOD[L] - 1
        NOR2 = NOR - 1
        fort_write(file_IO4, format99, [LAB[L]])
        fort_write(file_IO2, format99, [LAB[L]])
        fort_write(None, format3731, [L+1, NOR2, LAB[L]])

        for K in range(1, NOR):
            fort_write(file_IO2, format101, [K, ER[L, K], T[L, K]])
            fort_write(file_IO4, format100, [ER[L, K], T[L, K]])
        label .lbl8  # terminal label of FOR loop
        fort_write(file_IO2, format109, [W, W])
        fort_write(file_IO4, format100, [W, W])
    label .lbl7  # terminal label of FOR loop

    # START OF REDUCTION AND TRANSFER
    label .lbl500

    # variable not used in DATP fortran code
    NX = None
    # Bunch allows to access the dictionary elements
    # returned by DATRCL using the syntax expdata.varname
    xp = Bunch(DATRCL(file_ID3, 1, NX, 1))
    if xp.NR == 9999:
        goto .lbl160
    format3733 = "(' read data set  ',i7)"
    fort_write(None, format3733, [xp.NR])

    #  Interpolation type (this is very system specific -
    # does not apply for other simultaneous evaluations
    INT = 1
    if xp.NT == 1 or xp.NT == 2:
        goto .lbl92
    if xp.NT == 5 or xp.NT == 8:
        goto .lbl92
    goto .lbl94

    label .lbl92
    if xp.NID[0] == 2 or xp.NID[0] == 5:
        goto .lbl94
    if xp.NID[0] == 10:
        goto .lbl94
    if xp.NID[0] > 10:
        goto .lbl94
    INT = 2

    label .lbl94

    # CONSTRUCT APRIORI
    EQ = np.empty((200,), dtype=float)

    for K in range(NOM):
        Q[K] = 0.

    # NOTE: computed goto of fortran replaced
    #       by if-else statements
    if xp.NT == 1:
        goto .lbl11
    if xp.NT == 2:
        goto .lbl11
    if xp.NT == 3:
        goto .lbl12
    if xp.NT == 4:
        goto .lbl12
    if xp.NT == 5:
        goto .lbl13
    if xp.NT == 6:
        goto .lbl14
    if xp.NT == 7:
        goto .lbl15
    if xp.NT == 8:
        goto .lbl13
    if xp.NT == 9:
        goto .lbl15
    assert xp.NT >= 1 and xp.NT <= 9

    # CS + CS SHAPE
    label .lbl11
    M1 = xp.NID[0]
    pyM1 = M1 - 1
    NON = NOD[pyM1]
    NON1 = NON-1
    # NOTE: M1-1 due to different start index
    #       in Python (0) and Fortran (1)
    pyM1 = M1 - 1
    pyNON = NON - 1
    pyNON1 = NON1 - 1

    E11 = (ER[pyM1, 0] + ER[pyM1, 1]) / 2.
    E22 = (ER[pyM1, pyNON] + ER[pyM1, pyNON1]) / 2.
    for K in range(NON):
        EQ[K] = ER[pyM1, K]
        Q[K] = T[pyM1, K]
    label .lbl20  # terminal label of loop
    mxm = NON
    mxm1 = mxm - 1
    format4611 = "(i6,'  cr. sec. apriori for interp. ')"
    fort_write(None, format4611, [mxm])
    goto .lbl30

    # RATIO + RATIO SHAPE
    label .lbl12
    M1 = xp.NID[0]
    M2 = xp.NID[1]
    pyM1 = M1 - 1
    pyM2 = M2 - 1
    NO1 = NOD[pyM1]
    NON = NOD[pyM2]
    mxm = 0

    for K in range(NO1):
        # find matching energies
        etst1 = ER[pyM1, K] * 0.9999
        etst2 = ER[pyM1, K] * 1.0001
        for L in range(NON):
            if ER[pyM2, L] > etst1 and ER[pyM2, L] < etst2:
                goto .lbl56
        label .lbl55  # terminal label of loop
        goto .lbl58

        label .lbl56
        mxm = mxm + 1
        pymxm = mxm - 1
        EQ[pymxm] = ER[pyM1, K]
        if T[pyM2, L] == 0.:
            goto .lbl4616

        goto .lbl4617

        label .lbl4616
        format4618 = ":(' apriori constr. ',4i5,f10.7,'*****************')"
        fortK = K + 1
        fortL = L + 1
        fort_write(None, format4618, [M1, fortK, M2, fortL, EQ[pymxm]])
        exit()

        label .lbl4617
        Q[pymxm] = T[pyM1, K] / T[pyM2, L]

        label .lbl58

    label .lbl21  # terminal label of loop

    pymxm = mxm - 1
    mxm1 = mxm - 1
    pymxm1 = mxm1 - 1
    E11 = (EQ[0] + EQ[1]) / 2.
    E22 = (EQ[pymxm] + EQ[pymxm1]) / 2.

    format4612 = "(i6,'  ratio apriori for interp. ')"
    fort_write(None, format4612, [mxm])
    goto .lbl30

    # SUM AND SHAPE OF SUM
    label .lbl13
    M1 = xp.NID[0]
    pyM1 = M1 - 1
    M2 = xp.NID[1]
    pyM2 = M2 - 1
    M3 = xp.NID[2]
    pyM3 = M3 - 1
    NO1 = NOD[pyM1]
    NON = NOD[pyM2]
    if M3 == 0:
        goto .lbl22
    NO3 = NOD[pyM3]

    label .lbl22
    mxm = 0
    for K in range(NO1):  # 23
        etst1 = ER[pyM1, K] * 0.9999
        etst2 = ER[pyM1, K] * 1.0001

        for L in range(NON):  # 16
            if ER[pyM2, L] > etst1 and ER[pyM2, L] < etst2:
                goto .lbl17
        label .lbl16  # terminal label of loop

        goto .lbl26

        label .lbl17
        if M3 == 0:
            goto .lbl18

        for J in range(NO3):  # 19
            if ER[pyM3, J] > etst1 and ER[pyM3, J] < etst2:
                goto .lbl18
        label .lbl19  # terminal label of loop

        goto .lbl26

        label .lbl18
        mxm = mxm + 1
        pymxm = mxm - 1
        EQ[pymxm] = ER[pyM1, K]
        Q[pymxm] = T[pyM1, K] + T[pyM2, L]
        if M3 == 0:
            goto .lbl26
        Q[pymxm] = Q[pymxm] + T[pyM3, J]

        label .lbl26

    label .lbl23  # terminal label of loop

    mxm1 = mxm - 1
    pymxm1 = mxm1 - 1
    E11 = (EQ[0] + EQ[1]) / 2.
    E22 = (EQ[pymxm] + EQ[pymxm1]) / 2.

    format4613 = "(i6,'  sum apriori for interp. ')"
    fort_write(None, format4613, [mxm])
    goto .lbl30

    #  RATIO OF CS VS. SUM + SHAPE
    label .lbl15
    M1 = xp.NID[0]
    pyM1 = M1 - 1
    M2 = xp.NID[1]
    pyM2 = M2 - 1
    M3 = xp.NID[2]
    pyM3 = M3 - 1
    NO1 = NOD[pyM1]
    NON = NOD[pyM2]
    NO3 = NOD[pyM3]
    mxm = 0
    pymxm = mxm - 1

    for K in range(NO1):  # 27
        etst1 = ER[pyM1, K] * 0.9999
        etst2 = ER[pyM1, K] * 1.0001
        for L in range(NON):  # 76
            if ER[pyM2, L] > etst1 and ER[pyM2, L] < etst2:
                goto .lbl77
        label .lbl76  # end loop

        goto .lbl71

        label .lbl77
        for J in range(NO3):  # 79
            if ER[pyM3, J] > etst1 and ER[pyM3, J] < etst2:
                goto .lbl78
        label .lbl79

        goto .lbl71

        label .lbl78
        mxm = mxm + 1
        pymxm = mxm - 1
        EQ[pymxm] = ER[pyM1, K]
        Q[pymxm] = T[pyM1, K] / (T[pyM2, L] + T[pyM3, J])

        label .lbl71
    label .lbl27  # end of loop

    mxm1 = mxm - 1
    pymxm1 = mxm1 - 1
    E11 = (EQ[0] + EQ[1]) / 2.
    E22 = (EQ[pymxm] + EQ[pymxm1]) / 2.

    format4614 = "(i6,'  cr. sec. vs sum apriori for interp. ')"
    fort_write(None, format4614, [mxm])

    # REDUCTION
    label .lbl30

    # FIND USEFUL DATA RANGE
    if xp.E[0] > E22:
        goto .lbl61
    pyNO = xp.NO - 1
    if xp.E[pyNO] < E11:
        goto .lbl61

    label .lbl14
    NQMM = 0
    if NQQA == NQST:
        goto .lbl62
    if NQQA == NQND:
        goto .lbl63

    goto .lbl64

    label .lbl63
    fort_write(file_IO4, format250, [NQMM, NQMM])
    fort_write(file_IO2, format250, [NQMM, NQMM])

    label .lbl62
    format251 = '(4HBLCK,1X,2I5)'
    fort_write(file_IO4, format251, [NQMM, NQMM])
    fort_write(file_IO2, format251, [NQMM, NQMM])
    goto .lbl64

    # out of range
    label .lbl61
    if NQQA != NQND and xp.NQQ == NQND:
        NQQA = xp.NQQ
    goto .lbl500

    label .lbl64
    NQQA = xp.NQQ

    # find nr of CS involved
    for MN in range(1, 5):  # 66
        pyMN = MN - 1
        if xp.NID[pyMN] == 0:
            goto .lbl67
    label .lbl66  # end of loop

    label .lbl67
    NNN = MN - 1
    xp.NID[3] = 0

    # special marker for thermal constants
    if xp.NR >= 910 and xp.NR <= 934:
        xp.NID[3] = 1

    format253 = '(5HDATA ,9I5)'
    fort_write(file_IO4, format253, [xp.NR, xp.NT, xp.NCO, NNN, xp.NID[0:4]])
    fort_write(file_IO2, format253, [xp.NR, xp.NT, xp.NCO, NNN, xp.NID[0:4]])

    format254 = '(3I5,A28,8X,A20)'
    fort_write(file_IO4, format254, [xp.NY, xp.NQ, xp.NCS, xp.NAU, xp.NREF])
    fort_write(file_IO2, format254, [xp.NY, xp.NQ, xp.NCS, xp.NAU, xp.NREF])

    if xp.NT == 2 or xp.NT == 4:
        goto .lbl81
    if xp.NT == 8 or xp.NT == 9:
        goto .lbl81

    # normalization uncertainties
    format261 = '(10F5.1,10I3)'
    fort_write(file_IO4, format261, [xp.ENF[0:10], xp.NENF[0:10]])
    fort_write(file_IO2, format261, [xp.ENF[0:10], xp.NENF[0:10]])

    label .lbl81

    # energy dep. unc. parameters
    format262 = '(3F5.2,I3)'
    # NOTE: this loop is implicit in the write statement
    #       of the fortran code
    for K in range(11):
        fort_write(file_IO2, format262, [xp.EPA[0:3, K], xp.NETG[K]])
        fort_write(file_IO4, format262, [xp.EPA[0:3, K], xp.NETG[K]])
    if xp.NCS == 0:
        goto .lbl82

    # cross correlations
    format263 = '(I5,20I3)'
    format293 = '(10F5.1)'
    for K in range(xp.NCS):  # 83
        # NOTE: during flattening in fort_write first index should
        #       change fastest
        fort_write(file_IO4, format263, [xp.NCST[K], xp.NEC[:, :, K]])
        fort_write(file_IO2, format263, [xp.NCST[K], xp.NEC[:, :, K]])
        fort_write(file_IO4, format293, [xp.FCFC[0:10, K]])
        fort_write(file_IO2, format293, [xp.FCFC[0:10, K]])
    label .lbl83

    label .lbl82

    if xp.NT == 6:
        goto .lbl1003

    # GET GRID VALUES  - try at all apriori energies to find data

    # NOTE: format200 and format290 will be used
    #       only much later
    format200 = '(2E10.4,12F5.1)'
    format290 = '(2E10.4,12F5.1,F7.3)'

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
            xp.F[N, pyMAXF] = 0.
        label .lbl133  # end loop

        if E1 > .03:
            INT = 1

        # INTERPOLATION CONST.
        if INT == 1:
            goto .lbl88
        if INT == 2:
            goto .lbl89

        # LIN LIN
        label .lbl88
        AL = (Q[L-1]-Q[L])/(EQ[L-1]-EQ[L])
        BL = Q[L]-AL*EQ[L]
        AR = (Q[L]-Q[L+1])/(EQ[L]-EQ[L+1])
        BR = Q[L]-AR*EQ[L]
        goto .lbl28

        # LOG LOG
        label .lbl89

        QBL = (np.log(Q[L-1])-np.log(Q[L]))/(np.log(EQ[L])-np.log(EQ[L-1]))
        QAL = Q[L]*(EQ[L]**QBL)
        QBR = (np.log(Q[L])-np.log(Q[L+1]))/(np.log(EQ[L+1])-np.log(EQ[L]))
        QAR = Q[L]*(EQ[L]**QBR)

        label .lbl28

        # GRID VALUES
        for K in range(xp.NO):  # 35
            if xp.E[K] < E1*(1.-1e-5) or xp.E[K] >= E2*(1.+1e-5):
                goto .lbl91

            WT = 1./xp.F[11, K]
            WT = WT*WT
            if xp.E[K] > EQ[L]:
                goto .lbl32
            if xp.E[K] == EQ[L]:
                goto .lbl33

            # left o energy grid point
            if INT == 1:
                goto .lbl36
            if INT == 2:
                goto .lbl37

            label .lbl36
            ADD = AL * xp.E[K] + BL
            AD = xp.S[K] + Q[L] - ADD
            goto .lbl34

            label .lbl37
            ADD = QAL / (xp.E[K]**QBL)
            AD = xp.S[K] * Q[L] / ADD
            goto .lbl34

            # right of energy grid point
            label .lbl32
            if INT == 1:
                goto .lbl51
            if INT == 2:
                goto .lbl52

            label .lbl51
            ADD = AR*xp.E[K] + BR
            AD = xp.S[K] + Q[L] - ADD
            goto .lbl34

            label .lbl52
            ADD = QAR / (xp.E[K]**QBR)
            AD = xp.S[K] * Q[L] / ADD
            goto .lbl34

            # same energy as grid point
            label .lbl33
            AD = xp.S[K]
            label .lbl34

            # check if difference is within requested limit of ULI*sigma
            if ULI == 0:
                goto .lbl510
            T1X = 100.*(AD-Q[L])/Q[L]
            T2X = T1X*T1X
            TEST = np.sqrt(WT*T2X)
            if TEST < ULI:
                goto .lbl510
            F33 = xp.F[2, K] * xp.F[2, K]
            F44 = 1./WT - F33
            FNEW = T2X / (ULI*ULI)
            F33N = FNEW - F44
            xp.F[11, K] = np.sqrt(FNEW)
            xp.F[2, K] = np.sqrt(F33N)
            WT = 1./FNEW

            format511 = "(20X,' VALUE OUTSIDE ',F5.2,' SIGMA BY ',F10.2)"
            fort_write(file_IO2, format511, [ULI, TEST])

            label .lbl510
            AV = AV + AD*WT
            WTS = WTS + WT

            # statistical uncertainty reduces if more than one value contributes,
            # average for all other uncertainties
            for M in range(11):  # 38
                if xp.NETG[M] == 9:
                    goto .lbl53

                xp.F[M, pyMAXF] = xp.F[M, pyMAXF] + xp.F[M, K]
                goto .lbl54

                label .lbl53
                if xp.F[M, K] == 0.0:
                    goto .lbl54

                xp.F[M, pyMAXF] = xp.F[M, pyMAXF] + (1./xp.F[M, K])**2

                label .lbl54
            label .lbl38  # end of loop

            NKOT = NKOT + 1

            label .lbl91

        label .lbl35  # end of loop

        AKOT = NKOT
        if AV == 0.0:
            goto .lbl42

        # GRID VALUE AND OUT
        EEE = EQ[L]
        QQQ = AV / WTS
        DIF = QQQ / Q[L]
        for N in range(11):  # 39
            if xp.NETG[N] == 9:
                goto .lbl44
            xp.F[N, pyMAXF] = xp.F[N, pyMAXF] / AKOT
            goto .lbl46

            label .lbl44
            if xp.F[N, pyMAXF] <= 0.0:
                goto .lbl47
            xp.F[N, pyMAXF] = 1. / np.sqrt(xp.F[N, pyMAXF])
            goto .lbl46

            label .lbl47
            xp.F[N, pyMAXF] = 0.

            label .lbl46
        label .lbl39  # end of loop

        # OUTPUT
        fort_write(file_IO4, format200, [EEE, QQQ, xp.F[0:12, pyMAXF]])
        fort_write(file_IO2, format290, [EEE, QQQ, xp.F[0:12, pyMAXF], DIF])

        label .lbl42
    label .lbl40  # end of loop

    # end of data set
    fort_write(file_IO4, format200, [W, W, xp.F[0:12, pyMAXF]])
    fort_write(file_IO2, format200, [W, W, xp.F[0:12, pyMAXF]])

    if xp.NCO == 0:
        goto .lbl500

    format6114 = '(1X,10F8.5)'
    format6115 = '(10F8.5)'
    for KL in range(xp.NCO):  # 6113
        fort_write(file_IO2, format6114, [xp.ECOR[KL, :(KL+1)]])
        fort_write(file_IO4, format6115, [xp.ECOR[KL, :(KL+1)]])
    label .lbl6113  # end of loop

    goto .lbl500

    # fission spectrum average data set
    label .lbl1003
    fort_write(file_IO2, format200, [xp.E[0], xp.S[0], xp.F[0:12, 0]])
    fort_write(file_IO4, format200, [xp.E[0], xp.S[0], xp.F[0:12, 0]])
    fort_write(file_IO2, format200, [W, W, xp.F[0:12, pyMAXF]])
    fort_write(file_IO4, format200, [W, W, xp.F[0:12, pyMAXF]])

    goto .lbl500

    # DATA FILE COMPLETE
    label .lbl160
    format256 = '(4HEND*,1X,2I5)'
    fort_write(file_IO4, format250, [NQMM, NQMM])
    fort_write(file_IO2, format250, [NQMM, NQMM])
    fort_write(file_IO4, format256, [NQMM, NQMM])
    fort_write(file_IO2, format256, [NQMM, NQMM])

    file_IO1.close()
    file_IO2.close()
    file_IO3.close()
    file_IO4.close()


@with_goto
def DATRCL(file_ID3, NZ: int, NAB: int, IBZ: int):

    # variables with local scope
    NAU: str; NREF: str; NQT: str
    NCOM: str; NQQ: str; NXQT: str
    NXAU: str; ICC: str; NES: str; NEB: str

    ICC = 'C '
    NES = 'ES'
    NEB = 'EB'

    # this declaration is not present in Fortran code
    # but assumed to be implicitly done
    SES = 0.

    if NZ == 5:
        goto .lbl5030
    goto .lbl5031

    label .lbl5030
    for K in range(1,1000):
        NXY[K] = 0

    label .lbl5031
    NBQZ = 1
    if NBQZ == 1: NQQ = NES
    if IBZ == 2: NQQ = NEB

    # data set identification
    label .lbl10
    # original string
    #format100 = '(2I4,12A2,14A2,10A2)'

    format100 = '(2I4, A24, A28, A20)'
    NR, NY, NQT, NAU, NREF = fort_read(file_ID3, format100)
    if NR == 0:
        goto .lbl10
    if NR == 9999:
        return {'NR': NR}
    format103 = '(4I2,I3,I5,5I3)'
    NQ, NT, NCO, NCS, NCCO, NO, NID = unflatten(
            fort_read(file_ID3, format103), [6, [5]])

    # COMMENTS
    # original: (40A2)
    format106 = '(A80)'
    NCOM = []
    for i in range(NCCO):
        NCOM.append(fort_read(file_ID3, format106))

    # NORMALIZATION UNCERTAINTIES
    ENF = None
    NENF = None
    SES = 0.
    if NT == 2 or NT == 4:
        goto .lbl20
    if NT == 8 or NT == 9:
        goto .lbl20
    format107 = '(10F5.1, 10I3)'
    ENF, NENF = unflatten(fort_read(file_ID3, format107), [[10], [10]])

    for K in range(10):
        SES = SES + ENF[K]*ENF[K]

    # ENERGY DEPENDENT UNCERTAINTY CORRELATIONS PARAMETERS AND TAGS
    label .lbl20
    format110 = '(3F5.2)'
    EPA = np.empty((3,11), dtype=float)
    for i in range(11):
        EPA[:,i] = fort_read(file_ID3, format110)

    for k in range(11):
        absum = EPA[0,k] + EPA[1,k]
        if absum > 1.0:
            EPA[1,k] = 1.0 - EPA[0,k]

    format111 = '(11I3)'
    NETG = fort_read(file_ID3, format111)

    # DATA
    E = np.empty((NO,), dtype=float)
    S = np.empty((NO,), dtype=float)
    F = np.zeros((12, MAXF), dtype=float)
    format114 = '(2E10.4,12F5.1)'
    for K in range(NO):
        E[K], S[K], F[:,K] = unflatten(
                fort_read(file_ID3, format114), [2, [12]])
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
        tmp = fort_read(file_ID3, format116)
        # ISSUE: there are not always 20 I2 numbers
        #        in the GMDATA file but sometimes less
        tmp = [x for x in tmp if x is not None]
        tmp2 = np.zeros((20,), dtype=int)
        tmp2[:(len(tmp)-1)] = tmp[1:]
        NCST[K] = tmp[0]
        NEC[0, :, K] = tmp2[:10]
        NEC[1, :, K] = tmp2[10:]

        format452 = '(10F5.1)'
        tmp = fort_read(file_ID3, format452, none_as=0.)
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
            tmp = fort_read(file_ID3, format117)
            tmp = [x for x in tmp if x is not None]
            res += tmp
            num_el_read += len(tmp)
        ECOR[L, :(L+1)] = res

    format118 = '(A2)'
    NQQ = fort_read(file_ID3, format118)
    assert len(NQQ) == 1
    NQQ = NQQ[0]


    return({'NR': NR, 'NY': NY, 'NQT': NQT, 'NAU': NAU, 'NREF': NREF,
            'NQ': NQ, 'NT': NT, 'NCO': NCO, 'NCS': NCS, 'NCCO': NCCO, 'NO': NO, 'NID': NID,
            'NCOM': NCOM, 'ENF': ENF, 'NENF': NENF, 'SES': SES, 'EPA': EPA,
            'NETG': NETG, 'E': E, 'S': S, 'F': F,
            'NCST': NCST, 'NEC': NEC, 'NQQ': NQQ, 'FCFC': FCFC, 'ECOR': ECOR})




reduce_data()





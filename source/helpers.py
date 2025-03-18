import numpy as np
# helpful for tentative fortran to python conversion
from fortranformat import FortranRecordReader
from fortranformat import FortranRecordWriter


def debugout(debugstr):
    print('DEBUG: ' + str(debugstr))


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


def find_indices_with_tol(a, v, atol, rtol):
    a = np.array(a)
    v = np.array(v)
    # bring into order
    ordidcs = np.argsort(a)
    a = a[ordidcs]
    # match elements
    idcs = np.minimum(np.searchsorted(a, v), len(a)-1)
    s = np.isclose(a[idcs], v, atol=atol, rtol=rtol)
    rem_idcs = np.maximum(idcs[~s]-1, 0)
    idcs[~s] = rem_idcs
    s[~s] = np.isclose(a[rem_idcs], v[~s], atol=atol, rtol=rtol)
    idcs[~s] = -1  # indication for not found
    # map to original order
    found_sel = idcs != -1
    idcs[found_sel] = ordidcs[idcs[found_sel]]
    return idcs

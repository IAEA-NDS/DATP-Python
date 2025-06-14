import functools
import inspect
from copy import deepcopy
import json
import numpy as np
import argparse
from collections.abc import (
    MutableMapping,
    MutableSequence,
)


# DEBUG: This function is only used during debugging/modernization
#        to ensure that it is indeed called and bad changes in the
#        function will impact test results.
def must_be_called(func):

    def get_current_function():
        caller_frame = inspect.stack()[1]
        caller_function_name = caller_frame.function
        caller_function = caller_frame.frame.f_globals[caller_function_name]
        return caller_function

    def check_called():
        for func in this_decorator.registered_funcs:
            if not func.called:
                raise ValueError(
                    f'function {func.__name__} was not called'
                )
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        wrapper.called = True
        return func(*args, **kwargs)
    wrapper.called = False

    this_decorator = get_current_function()
    if not hasattr(this_decorator, "registered_funcs"):
        this_decorator.registered_funcs = set()
    this_decorator.registered_funcs.add(wrapper)
    this_decorator.check_called = check_called
    return wrapper


def diff_objs(d1, d2, path):
    if isinstance(d1, np.int64):
        d1 = int(d1)
    if isinstance(d2, np.int64):
        d2 = int(d2)
    if type(d1) != type(d2):
        print(f'type mismatch at {path} ({type(d1)} vs {type(d2)}')
    elif isinstance(d1, str) and isinstance(d2, str):
        d1 = d1.strip()
        d2 = d2.strip()
        if d1 != d2:
            print(f'value mismatch at {path} ({d1} vs {d2})')
    elif isinstance(d1, float) and isinstance(d2, float):
        if any(s in path for s in ('CO', 'ENFF', 'FCFC')):
            if not np.isclose(d1, d2, atol=0.051, rtol=0.0):
                print(f'value mismatch at {path} ({d1} vs {d2})')
        else:
            if not np.isclose(d1, d2, rtol=1e-3):
                print(f'value mismatch at {path} ({d1} vs {d2})')
    elif d1 != d2:
        print(f'value mismatch at {path} ({d1} vs {d2})')


def fancy_numeric_list_diff(d1, d2, path):
    work_d1 = d1.copy()
    work_d2 = d2.copy()
    i = 0
    while len(work_d1) > 0:
        v1 = work_d1.pop(0)
        if len(work_d2) == 0:
            break
        min_j = np.argmin([abs(v1-w) for w in work_d2])
        v2 = work_d2[min_j]
        if not np.isclose(v1, v2, rtol=1e-1):
            p1 = path + '/' + str(i)
            print(f'only in first list {p1} = {v1}')
        elif v1 != v2:
            j = d2.index(v2)
            p1 = path + '/' + str(i)
            p2 = path + '/' + str(j)
            print(f'inexact match {p1}={v1} and {p2}={v2}')
            work_d2.pop(min_j)
        else:  # perfect equality
            work_d2.pop(min_j)
        i += 1

    for v1 in work_d1:
        i = d1.index(v1)
        p1 = path + '/' + str(i)
        print(f'only in first list {p1} = {v1}')

    for v2 in work_d2:
        j = d2.index(v2)
        p2 = path + '/' + str(j)
        print(f'only in second list {p2} = {v2}')


def diff_lists(d1, d2, path):
    if len(d1) != len(d2):
        print(f'list length mismatch at {path} ({len(d1)} vs {len(d2)})')

    # sophisticated comparison for lists with floats
    numeric_types = (int, float, np.float64, np.int64)
    float_types = (float, np.float64)
    if (all(isinstance(x, numeric_types) for x in d1)
            and all(isinstance(x, numeric_types) for x in d2)
            and any(isinstance(x, float_types) for x in d1)
            and any(isinstance(x, float_types) for x in d2)):
        fancy_numeric_list_diff(d1, d2, path)
        return

    if len(d1) != len(d2):
        # comparison of generic lists with unequal length not supported
        return

    # regular comparison for the rest
    for i in range(len(d1)):
        mydiff(d1[i], d2[i], path + '/' + str(i))


def diff_dicts(d1, d2, path):
    keys_only_d1 = set(d1.keys()).difference(d2.keys())
    keys_only_d2 = set(d2.keys()).difference(d1.keys())
    for k in keys_only_d1:
        if k in ('comments', 'computed'):
            continue
        print(f'{path}/{k} only in first dictionary')
    for k in keys_only_d2:
        if k in ('comments', 'computed'):
            continue
        print(f'{path}/{k} only in second dictionary')
    common_keys = set(d1.keys()).intersection(d2.keys())
    for k in common_keys:
        curpath = path + '/' + str(k)
        mydiff(d1[k], d2[k], curpath)


def mydiff(d1, d2, path=''):
    if isinstance(d1, np.ndarray):
        d1 = d1.tolist()
    if isinstance(d2, np.ndarray):
        d2 = d2.tolist()
    if isinstance(d1, MutableMapping) and isinstance(d2, MutableMapping):
        diff_dicts(d1, d2, path)
    elif isinstance(d1, MutableSequence) and isinstance(d2, MutableSequence):
        diff_lists(d1, d2, path)
    else:
        diff_objs(d1, d2, path)


def extract_datasets(datablocks):
    datasets = {}
    for d in datablocks:
        datasets.update(
            {ds['NS']: ds for ds in d['datasets']}
        )
    return datasets


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='subcommand')
    subparsers.required = True

    parser_compare = subparsers.add_parser('compare')
    parser_compare.add_argument(
        '--extract-datasets', action='store_true', help='convert datablock list to dict'
    )
    parser_compare.add_argument("files", nargs=2, help="files for comparison")

    args = parser.parse_args()

    if args.subcommand == 'compare':
        file1 = args.files[0]
        file2 = args.files[1]
        with open(file1, 'r') as f:
            d1 = json.load(f)
        with open(file2, 'r') as f:
            d2 = json.load(f)
        if args.extract_datasets:
            d1 = extract_datasets(d1['datablocks'])
            d2 = extract_datasets(d2['datablocks'])
        mydiff(d1, d2)

import functools
import inspect
import numpy as np


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
            print(f'value mismatch at {path} ({d1} vs {d2}')
    elif isinstance(d1, float) and isinstance(d2, float):
        if any(s in path for s in ('CO', 'ENFF', 'FCFC')):
            if not np.isclose(d1, d2, atol=0.051, rtol=0.0):
                print(f'value mismatch at {path} ({d1} vs {d2}')
        else:
            if not np.isclose(d1, d2, rtol=1e-3):
                print(f'value mismatch at {path} ({d1} vs {d2}')
    elif d1 != d2:
        print(f'value mismatch at {path} ({d1} vs {d2}')


def diff_lists(d1, d2, path):
    if len(d1) != len(d2):
        print(f'list length mismatch at {path} ({len(d1)} vs {len(d2)}')
        return
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
        curpath = path + '/' + k
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

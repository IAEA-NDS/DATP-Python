import functools
import inspect


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

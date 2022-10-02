from functools import wraps
from time import time


def timing(f):
    '''Dev use only. Benchmarks time to complete any function decorated with `@timing`'''
    @wraps(f)
    def wrap(*args, **kw):
        ts = time()
        result = f(*args, **kw)
        te = time()
        print(f'**BENCHMARK** Function: {f.__name__} took: {te-ts} sec')
        return result
    return wrap



import unittest
import sys
import pdb
import functools
import traceback

#From https://stackoverflow.com/questions/4398967/python-unit-testing-automatically-running-the-debugger-when-a-test-fails
def debug_on(*exceptions):
    if not exceptions:
        exceptions = (AssertionError, )
    def decorator(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            try:
                return f(*args, **kwargs)
            except exceptions:
                info = sys.exc_info()
                traceback.print_exception(*info) 
                pdb.post_mortem(info[2])
        return wrapper
    return decorator


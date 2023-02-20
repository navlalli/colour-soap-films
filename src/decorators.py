""" Decorators """
import time

def timeit(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        output = func(*args, **kwargs)
        end_time = time.time()
        print(f"{func.__name__}() --> {end_time - start_time:.2f} s execution time")
        return output
    return wrapper


""" Decorators """
import time

def timeit(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        output = func(*args, **kwargs)
        end_time = time.time()
        print(f"Function {func.__name__}() took {end_time - start_time:.2f} s to execute")
        return output
    return wrapper


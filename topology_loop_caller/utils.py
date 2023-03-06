import time
from functools import wraps
from typing import Union, List
import os
from loguru import logger

def timeit(func):
    """
    Decorator to measure the execution time of functions.
    Can be disabled during mixer object initialization with
    disable_logs=False
    """

    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        logger.info(f"Function {func.__name__} took {total_time:.4f} seconds.")
        return result

    return timeit_wrapper


def list_full_paths(directory: str) -> list:
    """
    Function to list absolute paths of file in the given directory.
    """
    return [os.path.join(directory, file) for file in os.listdir(directory)]


def filter_list_of_paths(paths_list: List[str], selected_extentions = List[str]) -> Union[list, None]:
    """
    Function to filter file paths by their formats. The passed values are returned.
    """
    return [x for x in paths_list if any([x.endswith(ext) for ext in selected_extentions])]


def get_current_path() -> str:
    """
    Returned a path from which the scripts are executed.
    """
    print(os.getcwd())
    return os.getcwd()
RESULTS_FOLDER = os.path.join(get_current_path(), 'output')


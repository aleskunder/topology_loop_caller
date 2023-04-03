import time
from functools import wraps
from typing import Union, List, Dict, Any
import os
import glob
from loguru import logger
import argparse


def timeit(func):
    """
    Decorator to measure the execution time of functions.
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


def list_full_paths(files: list) -> list:
    """
    Function to list absolute paths of files in the given list.
    """
    result = []
    for file in files:
        if os.path.isdir(file):
            result.extend([os.path.join(file, f) for f in os.listdir(file)])
        else:
            result.append(os.path.abspath(file))
    return result

def listdir_nohidden(path):
    return glob.glob(os.path.join(path, '*'))


def filter_list_of_paths(
    paths_list: List[str], selected_extentions=List[str]
) -> Union[list, None]:
    """
    Function to filter file paths by their formats. The passed values are returned.
    """
    return [
        x for x in paths_list if any([x.endswith(ext) for ext in selected_extentions])
    ]


def get_current_path() -> str:
    """
    Returned a path from which the scripts are executed.
    """
    return os.getcwd()


def get_relative_path(filename: str) -> str:
    """
    Get the path relative to this script
    :param filename: str, a filename of the script
    :return: str, the path to a script
    """
    return os.path.join(os.path.dirname(__file__), filename)


def get_args_from_groups(parser, group_nums):
    """
    Get a list of command-line arguments from the specified groups in the parser object
    :param parser: argparse.ArgumentParser, the parser object
    :param group_nums: list of integers, the group numbers to extract arguments from
    :return: list of str, the command-line arguments
    """
    arg_list = []
    for group_num in group_nums:
        group = parser._action_groups[group_num]
        arg_list += [action.dest for action in group._group_actions]

    return arg_list


def filter_out_flag_args(args: argparse.Namespace) -> Dict[str, Any]:
    """
    Filter argparse.Namespace to remove arguments with only `action_true` and `action_false` actions.
    :param args: argparse.Namespace
    :return: Dict[str, Any]
    """
    return {
        k: v
        for k, v in vars(args).items()
        if not (
            isinstance(v, argparse._StoreTrueAction)
            or isinstance(v, argparse._StoreFalseAction)
        )
    }


def filter_present_flags(args, group_nums):
    # Get the action flags that are present in the specified groups
    action_flags = []
    for group_num in group_nums:
        for action in args._action_groups[group_num]._group_actions:
            if isinstance(action, argparse._StoreTrueAction) or isinstance(
                action, argparse._StoreFalseAction
            ):
                if getattr(args, action.dest) is True:
                    action_flags.append(f"--{action.dest}")
    return action_flags


def str_int_float_args(args: dict) -> dict:
    """
    Convert integer and float arguments to strings in a dictionary of parsed arguments
    :param args: dict, the parsed arguments from argparse
    :return: dict, the parsed arguments with integer and float arguments converted to strings
    """
    converted_args = {}
    for key, value in args.items():
        if isinstance(value, int) or isinstance(value, float):
            converted_args[key] = str(value)
        else:
            converted_args[key] = value
    return converted_args


RESULTS_FOLDER = os.path.join(get_current_path(), "output")

import argparse
import subprocess
import os
import json
from typing import Union, Dict, Any
from loguru import logger
from utils import (
    get_relative_path,
    get_args_from_groups,
    filter_out_flag_args,
    list_full_paths,
    listdir_nohidden,
    filter_present_flags,
    str_int_float_args,
    RESULTS_FOLDER
)


def pipeline_arg_parsing() -> argparse.Namespace:
    """
    Step 0: Define command-line arguments for the pipeline script
    :return: argparse.Namespace object
    """

    parser = argparse.ArgumentParser(
        description="""
    Command-line tool to execute a full pipeline: \n 
    1) load Cooler files & transform them into distance matrices by selected method,\n
    2) calculate persistent homologies,\n
    3) (TO DO) filter homologies, generate features.
    """
    )
    general_group = parser.add_argument_group("Arguments used in all steps")
    general_group.add_argument(
        "-i",
        "--input-files",
        dest="input_files",
        type=str,
        nargs="+",
        help="""Required: input files or a regular expression matching multiple input files.
               Example: -i /path/to/test.cool /path/to/other/test2.cool
               or: -i ./*.cool""",
        required=True,
        metavar="<files>",
    )
    general_group.add_argument(
        "-o",
        "--output-dir",
        dest="output_dir",
        type=str,
        help="Optional: output folder. If nothing provided output will be written into ./output/",
        default=RESULTS_FOLDER,
        required=False,
        metavar="<dir>",
    )
    step1_group = parser.add_argument_group("Arguments used in step 1 only")
    step1_group.add_argument(
        "-d",
        "--distance-function",
        dest="distance_function",
        type=str,
        help="Optional: method of contacts-to-distances transition. Default: pearson",
        default="pearson",
        choices=["pearson", "log"],
        required=False,
    )
    step1_group.add_argument(
        "--do-not-save-preserved-bins",
        dest="do_not_save_preserved_bins",
        help="A flag: do not preserve indices of no nan bins, if present. The bin numbers are required for subsequent bins-to-coordinates transformations.",
        action='store_true',
        required=False,
    )

    step1_group.add_argument(
        "--saved-file-prefix",
        dest="saved_file_prefix",
        type=str,
        metavar="No/prefix",
        default="No",
        help="Optional: a prefix of saved files. If No, each input filename if splitted by '_' and the first part is taken.",
        required=False,
    )
    step1_group.add_argument(
        "--no-pearson-sqrt",
        dest="no_pearson_sqrt",
        help="For pearson distance function: if not present, a square root of (1 - corr. matrix) is taken.",
        action='store_true',
        required=False,
    )
    step1_group.add_argument(
        "--log-zero-replacement-strategy",
        dest="log_zero_replacement_strategy",
        type=str,
        help="For logarithm distance function: zero replacement is required for logarithmization. See https://www.sciencedirect.com/science/article/pii/S0169743921000162",
        choices=["martin_fernandez", "half_min", "unif"],
        default="martin_fernandez",
        required=False,
    )
    step1_group.add_argument(
        "--log-base",
        dest="log_base",
        type=float,
        help="For logarithm distance function: a base of logarithm, default is 10.",
        default=10.0,
        metavar="int/float",
        required=False,
    )
    step1_group.add_argument(
        "--mcool-resolution",
        dest="mcool_resolution",
        type=int,
        help="Optional: a resolution for .mcool files, default is 1000.",
        default=1000,
        metavar="int",
        required=False,
    )
    step1_group.add_argument(
        "--no-balance-matrix",
        dest="no_balance_matrix",
        help="A flag: whether to balance the interaction frequency matrix before distance transformation",
        action='store_true',
        required=False,
    )
    step1_group.add_argument(
        "--fetch-fragment",
        dest="fetch_fragment",
        type=str,
        help="Optional: if necessary, a fragment (chr/chr fragment) is subsetted for each file.",
        default="No",
        metavar="No/str",
        required=False,
    )
    step2_group = parser.add_argument_group("Arguments used in step 2 only")
    step2_group.add_argument(
        "--maxdim",
        dest="maxdim",
        type=int,
        default=2,
        help="Compute persistent homology in dimensions mindim, ..., k.",
    )
    step2_group.add_argument(
        "--mindim",
        dest="mindim",
        type=int,
        default=1,
        help="Compute persistent homology in dimensions j, ..., maxdim.",
    )
    step2_group.add_argument(
        "--minrad",
        dest="minrad",
        type=float,
        default=-float("inf"),
        help="Compute homology from time t onward.",
    )
    step2_group.add_argument(
        "--maxrad",
        dest="maxrad",
        type=float,
        default=float("inf"),
        help="Stop computing homology after time t.",
    )
    step2_group.add_argument(
        "--numrad",
        dest="numrad",
        type=float,
        default=float("inf"),
        help="Divide the interval from minrad to maxrad into N equally spaced steps, and compute the homology of each step. If the value of numrad is set to Inf, then homology will computed at every time point.",
    )
    step2_group.add_argument(
        "--model",
        dest="model",
        type=str,
        default="vr",
        choices=["pc", "vr", "complex"],
        help="Used Eirene model, 'pc' (point cloud), 'vr' (vietoris-rips), or 'complex'.",
    )
    return parser.parse_args()


def run_first_step(args: argparse.Namespace) -> None:
    """
    Step 1: loading Cooler files, transforming to distance matrices, saving.
    Use subprocess.run to call the Python script with its command-line arguments
    :param args: argparse.Namespace, the output of parser.parse_args()
    :return: None
    """
    # Get the path of the calculate_persistent_homologies.jl file relative to this script
    matrix_transform_path = get_relative_path("matrix_operations.py")
    flag_arg_names = [
        "do_not_save_preserved_bins",
        "no_pearson_sqrt",
        "no_balance_matrix"
    ]
    cmd_arg_names = [
        "input_files",
        "output_dir",
        "distance_function",
        "saved_file_prefix",
        "log_zero_replacement_strategy",
        "log_base",
        "mcool_resolution",
        "fetch_fragment"
    ]
    args_dict = vars(args)
    flag_args = [f"--{flag_arg_name.replace('_', '-')}" for flag_arg_name in flag_arg_names if args_dict[flag_arg_name]]
    cmd_args_dict = str_int_float_args({k:args_dict[k] for k in cmd_arg_names})
    cmd_args = [f"--{k.replace('_', '-')}={' '.join(v) if type(v) == list else v}" for k, v in cmd_args_dict.items()]
    subprocess.run(["python", matrix_transform_path, *cmd_args, *flag_args])


def run_second_step(args: argparse.Namespace) -> None:
    """
    Step 2: Run a Julia script with Topological Data Analysis.
    Use subprocess.run to call the Python script with its command-line arguments
    :param args: argparse.Namespace, the output of parser.parse_args()
    :return: None
    """
    # Define input dir and output dir from basedir:
    matrices_path = os.path.join(args.output_dir, "distance_matrices/")
    matrices_paths = list_full_paths(listdir_nohidden(matrices_path))
    results_path = os.path.join(
        args.output_dir, "persistent_homology_results"
    )

    # Create the subdirectories if they don't exist
    os.makedirs(results_path, exist_ok=True)

    # Get the path of the calculate_persistent_homologies.jl file relative to this script
    calculate_homologies_path = get_relative_path("calculate_persistent_homologies.jl")

    cmd_arg_names = [
        "maxdim",
        "mindim",
        "maxrad",
        "minrad",
        "numrad",
        "model"
    ]
    args_dict = vars(args)
    cmd_args_dict = str_int_float_args({k:args_dict[k] for k in cmd_arg_names})
    cmd_args = [f"--{k.replace('_', '-')}={v}" for k, v in cmd_args_dict.items()]
    subprocess.run(["julia", calculate_homologies_path, "--input-matrices",
            *matrices_paths,
            "--output-path",
            results_path, *cmd_args])


def main():
    logger.add("pipeline.log", rotation="10 MB", compression="zip")
    args = pipeline_arg_parsing()
    logger.info(f"Parsed arguments: {json.dumps(vars(args), indent=4)}",)
    logger.info("Pipeline started")
    run_first_step(args)
    logger.info("First step completed")
    run_second_step(args)
    logger.info("Second step completed")
# Place for step3.py
    logger.success("Pipeline completed")


if __name__ == "__main__":
    main()

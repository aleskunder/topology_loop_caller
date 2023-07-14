import cooler
from loguru import logger
import numpy as np
import pandas as pd
import os
import glob
import re
import yaml
from typing import Union, Tuple, List
import argparse
from topology_loop_caller.utils import (
    timeit,
    list_full_paths,
    filter_list_of_paths,
    convert_to_python_dict,
    RESULTS_FOLDER,
)
from topology_loop_caller.distance_transform import (
    pearson_distance,
    negative_log_transformation,
)


@timeit



def main() -> None:
    parser = argparse.ArgumentParser(
        description="""Loads persistent homology files, processes features, adds new ones and saves results (.tsv files) to a designated path.""",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input-homologies",
        dest="input_homologies",
        type=str,
        nargs="+",
        help="""Required: input homologies files or a regular expression matching multiple input files.
               Example: -i /path/to/test.tsv /path/to/other/test2.tsv
               or: -i ./*.tsv""",
        required=True,
        metavar="<files>",
    )
    parser.add_argument(
        "-o",
        "--output-path",
        dest="output_path",
        type=str,
        help="Optional: output folder. If nothing provided output will be written into ./output/persistent_homology_postprocessed",
        default=RESULTS_FOLDER,
        required=False,
        metavar="<dir>",
    )
    parser.add_argument(
        "--include-trans-chrom-contacts",
        dest="include_trans_chrom_contacts",
        help="A flag: whether to compute the distance matrix for the whole genome. If not present, distance matrices are saved for all chromosomes separately",
        action="store_true",
        required=False,
    )
    parser.add_argument(
        "--fetch-fragment",
        dest="fetch_fragment",
        type=str,
        help="Optional: if necessary, a fragment (chr/chr fragment) is subsetted for each file.",
        default="No",
        metavar="No/str",
        required=False,
    )

    args = parser.parse_args()

    logger.info(args)
    # Check if any value is an empty string, and if so, replace it with None
    for k, v in vars(args).items():
        if v == "No":
            setattr(args, k, None)

    input_files_str = " ".join(args.input_files)

    input_files = []
    for input_file in args.input_files:
        input_files += glob.glob(input_file)

    if len(input_files) == 0:
        logger.error("The given files do not exist or are not .tsv files.")
        return

    # Optional arguments
    saved_file_base_folder = args.output_dir
    cis_chrom_only = not args.include_trans_chrom_contacts
    fetch_fragment = args.fetch_fragment
    input_file_paths = list_full_paths(input_files)
    logger.success(f"{input_file_paths} are parsed.")
    input_file_paths = filter_list_of_paths(input_file_paths, [".tsv"])
    if input_file_paths:
        logger.success(f"{input_file_paths} are kept after filtering.")
    else:
        logger.error("The given folder does not contain .tsv files.")

    # Create a list of dictionaries to store the arguments for each file
    file_args = []
    for i, name in enumerate(input_file_paths, 1):
        # Determine the replica name from the filename
        if not saved_file_prefix:
            replica_name = name.split("/")[-1].split(".")[0].split("_")[0]
        elif len(input_file_paths) > 1:
            replica_name = f"{saved_file_prefix}_{str(i)}"
        else:
            replica_name = saved_file_prefix

    for i, name in enumerate(input_file_paths):

        logger.info(
            f"Start postprocessing with the following arguments: {file_args[i]}"
        )
        transform_and_save_matrix(bal, **file_args[i])


if __name__ == "__main__":
    main()

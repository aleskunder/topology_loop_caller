import cooler
from loguru import logger
import numpy as np
import os
from typing import Union
import argparse
from topology_loop_caller.utils import (
    timeit,
    list_full_paths,
    filter_list_of_paths,
    RESULTS_FOLDER,
)
from topology_loop_caller.distance_transform import (
    pearson_distance,
    negative_log_transformation,
)


@timeit
def load_cooler(
    file_path: str,
    resolution: int = 1000,
    balance_matrix: bool = True,
    fetch_fragment: str = None,
    **kwargs,
) -> np.array:
    """
    Function to load Cooler object and return interaction ferquency matrix
    :param file_path: Union[str, Path], a path to a cooler file
    :param resolution: int, a resolution for .mcool files
    :param balance_matrix: bool, whether to balance the interaction ferquency matrix
    :param fetch_fragment: if necessary, a fragment (chr/chr fragment) is subsetted
    :returns: 2D NumPy array
    """
    parsed_file_format = parse_format(file_path)
    if parsed_file_format == "cool":
        c = cooler.Cooler(file_path, **kwargs)
        logger.success(f"Loaded given .cool file {file_path.split('/')[-1]}.")
    else:
        c = cooler.Cooler(f"{file_path}::resolutions/{resolution}", **kwargs)
        logger.success(
            f"Loaded given .mcool file {file_path.split('/')[-1]} with resolution {resolution}."
        )
    if fetch_fragment:
        result_matrix = c.matrix(balance=balance_matrix).fetch(fetch_fragment)
    else:
        result_matrix = c.matrix(balance=balance_matrix)[:, :]

    return result_matrix


def parse_format(file_path: str) -> str:
    """
    Function to check the right format (.cool, .mcool) of a given Cooler file.
    :param file_path: Union[str, Path], a path to a cooler file
    :returns: NumPy matrix
    """
    parsed_format = file_path.split(".")[-1]
    if parsed_format not in ["cool", "mcool"]:
        logger.error(
            f"Check the given file path. .cool or .mcool format is required, {parsed_format} is given"
        )
    return parsed_format


@timeit
def transform_and_save_matrix(
    balanced_matrix: np.array,
    saved_file_prefix: str = "output",
    saved_file_base_folder: str = RESULTS_FOLDER,
    save_preserved_bins: bool = True,
    distance_function: str = "pearson",
    **kwargs,
) -> None:
    """
    Function transforms balanced matrix to a distance matrix by a given method
    and saves .npy and (optionally) indices of preserved bins.
    :param balanced_matrix: np.array, an interaction frequency matrix
    :param saved_file_prefix: str, a prefix of saved files
    :param saved_file_base_folder: Union[str, Path], a path to a results folder
    :param save_preserved_bins: bool, whether to preserve indices of no nan bins. The bin numbers are required for subsequent bins-to-coordinates transformations.
    :param distance_function: 'pearson' or 'log', a method of contacts - distances transition.
    :returns: none
    """
    assert distance_function in [
        "pearson",
        "log",
    ], "Wrong distance_function argument: 'pearson' or 'log' are the options"
    not_nan_bin_idx = np.logical_not(np.isnan(balanced_matrix).all(axis=1))

    # Indices saving:
    if save_preserved_bins:
        # Create (sub)folders, if not exist:
        preserved_bins_folder = os.path.join(saved_file_base_folder, "bins_indices")
        os.makedirs(preserved_bins_folder, exist_ok=True)

        deleted_bins = np.arange(balanced_matrix.shape[0])[
            ~not_nan_bin_idx
        ]  # Indices of deleted bins

        saved_bins = np.arange(balanced_matrix.shape[0])[
            not_nan_bin_idx
        ]  # Indices of saved bins

        np.save(
            os.path.join(preserved_bins_folder, f"{saved_file_prefix}_saved.npy"),
            saved_bins,
        )
        logger.success(f"{saved_file_prefix}_saved.npy is saved")

        np.save(
            os.path.join(preserved_bins_folder, f"{saved_file_prefix}_deleted.npy"),
            deleted_bins,
        )
        logger.success(f"{saved_file_prefix}_deleted.npy is saved")

    # Create (sub)folders, if not exist:
    dm_path = os.path.join(saved_file_base_folder, "distance_matrices")
    os.makedirs(dm_path, exist_ok=True)

    if distance_function == "pearson":
        np.save(
            os.path.join(dm_path, f"{saved_file_prefix}.npy"),
            pearson_distance(balanced_matrix, **kwargs),
        )
    else:
        np.save(
            os.path.join(dm_path, f"{saved_file_prefix}.npy"),
            negative_log_transformation(balanced_matrix, **kwargs),
        )
    logger.success(
        f"{saved_file_prefix}.npy is transformed with {distance_function} method and saved."
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="""Loads cooler files, converts to distance matrices with a selected method and saves results (.npy matrices optionally with NaN bin numbers) to a designated path.""",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input-dir",
        dest="input_dir",
        type=str,
        help="Path to a folder with cooler files.",
        required=True,
        metavar="<dir>",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        dest="saved_file_base_folder",
        type=str,
        help="Output folder. If nothing provided output will be written into ./output/",
        default=RESULTS_FOLDER,
        required=False,
        metavar="<dir>",
    )
    parser.add_argument(
        "-d",
        "--distance-function",
        dest="distance_function",
        type=str,
        help="Method of contacts - distances transition.",
        default="pearson",
        choices=["pearson", "log"],
        required=False,
    )
    parser.add_argument(
        "--save-preserved-bins",
        dest="save_preserved_bins",
        type=bool,
        metavar="bool",
        help="Whether to preserve indices of no nan bins. The bin numbers are required for subsequent bins-to-coordinates transformations.",
        default=True,
        required=False,
    )
    parser.add_argument(
        "--saved-file-prefix",
        dest="saved_file_prefix",
        type=Union[str, None],
        metavar="None/prefix",
        default=None,
        help="Prefix of saved files. If None, each input filename if splitted by '_' and the first part is taken.",
        required=False,
    )
    parser.add_argument(
        "--pearson-sqrt",
        dest="pearson_sqrt",
        type=bool,
        help="For pearson distance function: whether to take a square root of 1 - corr. matrix.",
        default=True,
        metavar="bool",
        required=False,
    )
    parser.add_argument(
        "--log-zero-replacement-strategy",
        dest="log_zero_replacement_strategy",
        type=str,
        help="For logarithm distance function: zero replacement is required for logarithmization. See https://www.sciencedirect.com/science/article/pii/S0169743921000162",
        choices=["martin_fernandez", "half_min", "unif"],
        default="martin_fernandez",
        required=False,
    )
    parser.add_argument(
        "--log-base",
        dest="log_base",
        type=Union[float, int],
        help="For logarithm distance function: a base of logarithm, default is 10.",
        default=10,
        metavar="int/float",
        required=False,
    )
    parser.add_argument(
        "--mcool-resolution",
        dest="mcool_resolution",
        type=int,
        help="A resolution for .mcool files, default is 1000.",
        default=1000,
        metavar="int",
        required=False,
    )
    parser.add_argument(
        "--to-balance-matrix",
        dest="to_balance_matrix",
        type=bool,
        help="bool, whether to balance the interaction ferquency matrix before distance transformation",
        metavar="bool",
        default=True,
        required=False,
    )
    parser.add_argument(
        "--fetch-fragment",
        dest="fetch_fragment",
        type=Union[None, str],
        help="if necessary, a fragment (chr/chr fragment) is subsetted for each file.",
        default=None,
        metavar="None/str",
        required=False,
    )

    args = parser.parse_args()

    # Required arguments
    input_dir = args.input_dir

    # Optional arguments
    saved_file_base_folder = args.saved_file_base_folder
    distance_function = args.distance_function
    save_preserved_bins = args.save_preserved_bins
    saved_file_prefix = args.saved_file_prefix
    distance_function_kwargs = {
        "pearson": {"sqrt": args.pearson_sqrt},
        "log": {
            "zero_replacement_strategy": args.log_zero_replacement_strategy,
            "log_base": args.log_base,
        },
    }
    resolution = args.mcool_resolution
    to_balance_matrix = args.to_balance_matrix
    fetch_fragment = args.fetch_fragment
    input_file_paths = list_full_paths(input_dir)
    logger.success(f"{input_file_paths} are parsed.")
    input_file_paths = filter_list_of_paths(input_file_paths, [".cool", ".mcool"])
    if input_file_paths:
        logger.success(f"{input_file_paths} are kept after filtering.")
    else:
        logger.error("The given folder does not contain .cool or .mcool files.")

    # Create a list of dictionaries to store the arguments for each file
    file_args = []
    for i, name in enumerate(input_file_paths, 1):
        # Determine the replica name from the filename
        if not saved_file_prefix:
            replica_name = name.split("/")[-1].split(".")[0].split("_")[0]
        else:
            replica_name = f"{saved_file_prefix}_{str(i)}"

        # Set the arguments for transform_and_save_matrix
        transform_args = {
            "saved_file_prefix": replica_name,
            "saved_file_base_folder": saved_file_base_folder,
            "save_preserved_bins": save_preserved_bins,
        }

        # Determine which distance function to use and set the appropriate arguments
        transform_args.update(
            {
                "distance_function": distance_function,
                **distance_function_kwargs[distance_function],
            }
        )

        # Append the arguments for this file to the list of file_args
        file_args.append(transform_args)

    # Use a list comprehension to call transform_and_save_matrix for each file with the appropriate arguments
    for i, name in enumerate(input_file_paths):
        # Load the cooler file
        bal = load_cooler(
            name,
            resolution=resolution,
            balance_matrix=to_balance_matrix,
            fetch_fragment=fetch_fragment,
        )

        logger.info(
            f"Start transforming and saving with the following arguments: {file_args[i]}"
        )
        transform_and_save_matrix(bal, **file_args[i])


if __name__ == "__main__":
    main()

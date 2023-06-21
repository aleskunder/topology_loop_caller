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
def load_cooler(
    file_path: str,
    resolution: int = 1000,
    balance_matrix: bool = True,
    fetch_fragment: str = None,
    cooler_info_folder: str = None,
    cis_contacts_only: bool = False,
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
    chrom_names = c.chromnames
    info = c.info
    matrices = {}
    if fetch_fragment:
        result_matrix = c.matrix(balance=balance_matrix).fetch(fetch_fragment)
        matrices['fragment'] = result_matrix
        logger.success(f"Fetched the given fragment {fetch_fragment}.")
    elif cis_contacts_only:
        for chrom_name in chrom_names:
            matrices[chrom_name] = c.matrix(balance=balance_matrix).fetch(chrom_name)
    else:
        result_matrix = c.matrix(balance=balance_matrix)[:, :]
        matrices['genome'] = result_matrix

    return matrices, info

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
    input_matrices: tuple,
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
    

    matrices, cooler_info = input_matrices

    # Indices saving:
    if save_preserved_bins:
        # Create (sub)folders, if not exist:
        preserved_bins_folder = os.path.join(saved_file_base_folder, "bins_indices",)
        os.makedirs(preserved_bins_folder, exist_ok=True)
        folders = [
            name for name in os.listdir(preserved_bins_folder)
            if os.path.isdir(os.path.join(preserved_bins_folder, name))
        ]

        if len(folders) == 0:
            next_folder_num = 1
        else:
            folder_nums = [int(re.findall(r'\d+', name)[0]) for name in folders]
            next_folder_num = max(folder_nums) + 1

        chrom_bins_folder = os.path.join(preserved_bins_folder, f"file_{next_folder_num}")
        os.makedirs(chrom_bins_folder, exist_ok=True)
        for i, (chrom_name, matrix) in enumerate(matrices.items()):
            not_nan_bin_idx = np.logical_not(np.isnan(matrix).all(axis=1))
            deleted_bins = np.arange(matrix.shape[0])[
                ~not_nan_bin_idx
            ]  # Indices of deleted bins

            saved_bins = np.arange(matrix.shape[0])[
                not_nan_bin_idx
            ]  # Indices of saved bins

            np.save(
                os.path.join(chrom_bins_folder, f"{saved_file_prefix}_{chrom_name}_saved.npy"),
                saved_bins,
            )

            np.save(
                os.path.join(chrom_bins_folder, f"{saved_file_prefix}_{chrom_name}_deleted.npy"),
                deleted_bins,
            )

    # Create (sub)folders, if not exist:
    dm_path = os.path.join(saved_file_base_folder, "distance_matrices")
    os.makedirs(dm_path, exist_ok=True)
    folders = [
        name for name in os.listdir(dm_path)
        if os.path.isdir(os.path.join(dm_path, name))
    ]

    if len(folders) == 0:
        next_folder_num = 1
    else:
        folder_nums = [int(re.findall(r'\d+', name)[0]) for name in folders]
        next_folder_num = max(folder_nums) + 1

    chrom_dm_path = os.path.join(dm_path, f"file_{next_folder_num}")
    os.makedirs(chrom_dm_path, exist_ok=True)

    for i, (chrom_name, matrix) in enumerate(matrices.items()):
        if np.all(np.isnan(matrix)):
            continue
        if distance_function == "pearson":
            np.save(
                os.path.join(chrom_dm_path, f"{saved_file_prefix}_{chrom_name}.npy"),
                pearson_distance(matrix, **kwargs),
            )
        else:
            np.save(
                os.path.join(chrom_dm_path, f"{saved_file_prefix}_{chrom_name}.npy"),
                negative_log_transformation(matrix, **kwargs),
            )

    temp_folder = os.path.join(saved_file_base_folder, "temp")
    os.makedirs(temp_folder, exist_ok=True)

    cooler_info_path = os.path.join(temp_folder, f"{saved_file_prefix}_cooler_info.yaml")
    with open(cooler_info_path, "w") as file:
        yaml.dump(convert_to_python_dict(cooler_info), file)
    logger.success(
        f"{saved_file_prefix}.npy is transformed with {distance_function} method and saved."
    )

def bin_to_coor(
    bin_num: int, organism_clr: cooler.Cooler, scale_factor: int = 2000
) -> Tuple:
    """
    Function transforms bin number to genomic coordinates
    :param bin_num: The bin number to be converted
    :param organism_clr: A cooler.Cooler object representing the organism's Hi-C data
    :param scale_factor: The scale factor to be used in computing the genomic coordinates
    :returns: A tuple representing the genomic coordinates, with the format (chromosome, start, end)
    """
    chromsizes = organism_clr.chromsizes  # get chromosome sizes from the cooler object
    k = 0  # initialize a variable to keep track of the cumulative offset in bin numbers for each chromosome
    for i, size in enumerate(
        chromsizes
    ):  # iterate over the chromosomes and their sizes
        if (
            bin_num - k
        ) * scale_factor > size:  # if the bin number is past the end of the current chromosome
            k += np.ceil(
                size / scale_factor
            )  # update the cumulative offset by the number of bins in the current chromosome
        else:  # if the bin number is within the current chromosome
            start = (
                bin_num - k
            ) * scale_factor  # compute the start coordinate of the bin
            end = min(
                start + scale_factor, size
            )  # compute the end coordinate of the bin, ensuring that it doesn't extend past the end of the chromosome
            return (
                organism_clr.chromnames[i],
                int(start),
                int(end),
            )  # return the chromosome name and computed start and end coordinates
    return "error"  # if the bin number is past the end of the last chromosome, return an error message


def coor_to_bin(chr_coor: Tuple, organism_clr: cooler.Cooler, scale_factor: int = 2000) -> Union[int, str]:
    """
    Function transforms genomic coordinates to bin number
    :param chr_coor: tuple of chromosome name, start, end coordinates
    :param organism_clr: cooler object
    :param scale_factor: scaling factor, default 2000
    :return: bin number or "error" if input is invalid
    """
    chromsizes = organism_clr.chromsizes
    chromnames = organism_clr.chromnames
    
    # Check if chromosome name is valid
    chrom_name = chr_coor[0]
    if chrom_name not in chromnames:
        return "error: invalid chromosome name"
    chrom_index = chromnames.index(chrom_name)
    chrom_size = chromsizes[chrom_index]
    
    # Check if coordinates are valid
    start_coord = chr_coor[1]
    end_coord = chr_coor[2]
    if end_coord > chrom_size or end_coord <= start_coord or end_coord - start_coord > scale_factor:
        return "error: invalid coordinates"
    
    # Calculate bin number
    start_bin = np.floor(start_coord / scale_factor)
    end_bin = np.floor(end_coord / scale_factor)
    bin_num = end_bin - start_bin
    bin_start = start_bin + organism_clr.offset(chrom_name)
    bin_end = bin_start + bin_num
    
    return int(bin_start - 1 + bin_end - bin_start)


def find_matching_loop_indices(source_df: pd.DataFrame, target_df: pd.DataFrame, 
                               window_size: int = 0, scaling_factor: int = 2000) -> List[int]:
    '''
    Find indices of loops in source_df that have matching loop starts in target_df.
    
    :params source_df: pd.DataFrame with columns chrom1 & 2, start1 & 2, replica
    :params target_df: pd.DataFrame with same columns as source_df
    :param window_size: int, window for loop, loops are intersected
                        iff abs(start1 for source_df - start1 for target_df) <= window_size & 
                        same for start2 with corresponding chromosomes
    :param scaling_factor: int, scaling factor to convert window size to base pairs
    :returns: a list of indices in source_df with intersecting loop starts' intersection
    '''
    matching_indices = []
    for source_ind, source_row in source_df.iterrows():
        target_subset = target_df[(target_df['replica'] == source_row['replica']) &
                                  (target_df['chrom1'] == source_row['chrom1']) &
                                  (target_df['chrom2'] == source_row['chrom2']) &
                                  (abs(source_row['start1'] - target_df['start1']) <= window_size * scaling_factor) &
                                  (abs(source_row['start2'] - target_df['start2']) <= window_size * scaling_factor)]
        if not target_subset.empty:
            matching_indices.append(source_ind)
    return matching_indices


def main() -> None:
    parser = argparse.ArgumentParser(
        description="""Loads cooler files, converts to distance matrices with a selected method and saves results (.npy matrices optionally with NaN bin numbers) to a designated path.""",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
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
    parser.add_argument(
        "-o",
        "--output-dir",
        dest="output_dir",
        type=str,
        help="Optional: output folder. If nothing provided output will be written into ./output/",
        default=RESULTS_FOLDER,
        required=False,
        metavar="<dir>",
    )
    parser.add_argument(
        "-d",
        "--distance-function",
        dest="distance_function",
        type=str,
        help="Optional: method of contacts-to-distances transition. Default: Pearson",
        default="pearson",
        choices=["pearson", "log"],
        required=False,
    )
    parser.add_argument(
        "--do-not-save-preserved-bins",
        dest="do_not_save_preserved_bins",
        help="A flag: do not preserve indices of no nan bins, if present. The bin numbers are required for subsequent bins-to-coordinates transformations.",
        action="store_true",
        required=False,
    )

    parser.add_argument(
        "--saved-file-prefix",
        dest="saved_file_prefix",
        type=str,
        metavar="No/prefix",
        default="No",
        help="Optional: a prefix of saved files. If No, each input filename if splitted by '_' and the first part is taken.",
        required=False,
    )
    parser.add_argument(
        "--no-pearson-sqrt",
        dest="no_pearson_sqrt",
        help="For pearson distance function: if not present, a square root of (1 - corr. matrix) is taken.",
        action="store_true",
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
        type=float,
        help="For logarithm distance function: a base of logarithm, default is 10.",
        default=10.0,
        metavar="int/float",
        required=False,
    )
    parser.add_argument(
        "--mcool-resolution",
        dest="mcool_resolution",
        type=int,
        help="Optional: a resolution for .mcool files, default is 1000.",
        default=1000,
        metavar="int",
        required=False,
    )
    parser.add_argument(
        "--no-balance-matrix",
        dest="no_balance_matrix",
        help="A flag: whether to balance the interaction frequency matrix before distance transformation",
        action="store_true",
        required=False,
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
        logger.error("The given files do not exist or are not .cool or .mcool files.")
        return

    # Optional arguments
    saved_file_base_folder = args.output_dir
    distance_function = args.distance_function
    save_preserved_bins = not args.do_not_save_preserved_bins
    saved_file_prefix = args.saved_file_prefix
    distance_function_kwargs = {
        "pearson": {"sqrt": not args.no_pearson_sqrt},
        "log": {
            "zero_replacement_strategy": args.log_zero_replacement_strategy,
            "log_base": float(args.log_base),
        },
    }
    resolution = args.mcool_resolution
    to_balance_matrix = not args.no_balance_matrix
    cis_chrom_only = not args.include_trans_chrom_contacts
    fetch_fragment = args.fetch_fragment
    input_file_paths = list_full_paths(input_files)
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
        elif len(input_file_paths) > 1:
            replica_name = f"{saved_file_prefix}_{str(i)}"
        else:
            replica_name = saved_file_prefix

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
            cis_contacts_only=cis_chrom_only,
            fetch_fragment=fetch_fragment,
        )

        logger.info(
            f"Start transforming and saving with the following arguments: {file_args[i]}"
        )
        transform_and_save_matrix(bal, **file_args[i])


if __name__ == "__main__":
    main()

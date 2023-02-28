import cooler
from loguru import logger
import numpy as np
from topology_loop_caller.utils import timeit
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
    if parsed_file_format == ".cool":
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
    if parsed_format not in [".cool", ".mcool"]:
        logger.error(
            f"Check the given file path. .cool or .mcool format is required, {parsed_format} is given"
        )
    return parsed_format


@timeit
def transform_and_save_matrix(
    balanced_matrix: np.array,
    saved_file_prefix: str = "output",
    saved_file_base_folder: str = "../results/NoNAN_DM_new/",
    save_preserved_idx: bool = True,
    preserved_idx_folder: str = "../results/DM_indices/",
    distance_function: str = "pearson",
    **kwargs,
) -> None:
    """
    Function transforms balanced matrix to a distance matrix by a given method
    and saves .npy and (optionally) indices of preserved bins.
    :param balanced_matrix: np.array, an interaction frequency matrix
    :param saved_file_prefix: str, a prefix of saved files
    :param saved_file_base_folder: Union[str, Path], a path to a results folder
    :param save_preserved_idx: bool, whether to preserve indices of no nan bins
    :param preserved_idx_folder: Union[str, Path], a path to a bins folder
    :param distance_function: 'pearson' or 'log', a method of contacts - distances transition.
    :returns: none
    """
    assert distance_function in [
        "pearson",
        "log",
    ], "Wrong distance_function argument: 'pearson' or 'log' are the options"
    not_nan_bin_idx = np.logical_not(np.isnan(balanced_matrix).all(axis=1))
    if save_preserved_idx:
        deleted_bins = np.arange(balanced_matrix.shape[0])[
            ~not_nan_bin_idx
        ]  # Indices of deleted bins
        saved_bins = np.arange(balanced_matrix.shape[0])[
            not_nan_bin_idx
        ]  # Indices of saved bins
        np.save(f"{preserved_idx_folder}{saved_file_prefix}_saved.npy", saved_bins)
        logger.success(f"{saved_file_prefix}_saved.npy is saved")
        np.save(f"{preserved_idx_folder}{saved_file_prefix}_deleted.npy", deleted_bins)
        logger.success(f"{saved_file_prefix}_deleted.npy is saved")
    if distance_function == "pearson":
        np.save(
            f"{saved_file_base_folder}{saved_file_prefix}.npy",
            pearson_distance(balanced_matrix, **kwargs),
        )
    else:
        np.save(
            f"{saved_file_base_folder}{saved_file_prefix}.npy",
            negative_log_transformation(balanced_matrix, **kwargs),
        )
    logger.success(
        f"{saved_file_prefix}.npy is transformed with {distance_function} method and saved."
    )

import cooler
from loguru import logger
import numpy as np
import pandas as pd
import os
import glob
from typing import Union, Tuple, List, Dict
import argparse
from topology_loop_caller.utils import (
    timeit,
    list_full_paths,
    filter_list_of_paths,
    group_files_by_prefix,
    RESULTS_FOLDER,
)


def dist_to_cont_bin(saved_indices: np.ndarray, bin_num: int) -> Union[int, float]:
    """
    Converts the bin number in a filtered matrix to the corresponding bin number in the original matrix.
    
    Args:
    - saved_indices: a numpy array storing the information of preserved Hi-C bins
    - bin_num: the bin number in the filtered matrix
    
    Returns:
    - the corresponding bin number in the original matrix
    - np.nan if the bin number is not found in the saved indices
    
    """
    if len(saved_indices) > bin_num:
        return saved_indices[bin_num]
    else: 
        print(f"Bin number {bin_num} is not found in the saved indices.")
        return np.nan 


def cont_to_dist_bin(saved_indices: np.ndarray, contig_bin_num: int) -> Union[int, str]:
    """
    Converts a contiguous bin number to its corresponding index in the saved_indices array.
    :param saved_indices: numpy array of indices of preserved Hi-C bins
    :param contig_bin_num: the bin number in the contiguous matrix
    :return: the corresponding index in saved_indices, or an error message
    """
    try:
        dist_bin_num = int(np.where(saved_indices == contig_bin_num)[0])
        return dist_bin_num
    except:
        return "Invalid bin number or saved indices array."


def bottleneck_death(dist_matrix, loops_df):
    d = []
    for ind in loops_df.index:
        vert = loops_df.loc[ind, 'vertices']
        temp_dist_m = dist_matrix[vert, :][:, vert]
        d.append(np.amax(temp_dist_m, axis=1).min())
    return d


def maxrad(dist_matrix, loops_df):
    maxrad = []
    for ind in loops_df.index:
        vert = loops_df.loc[ind, 'vertices']
        temp_dist_m = dist_matrix[vert, :][:, vert]
        maxrad.append(temp_dist_m.max())
    return maxrad


def loop_perimeter(dist_matrix, loops_df):
    per = []
    for ind in loops_df.index:
        vert = loops_df.loc[ind, 'vertices']
        temp_dist_m = dist_matrix[vert, :][:, vert]
        per.append(np.sum(np.diagonal(temp_dist_m, offset=1))+temp_dist_m[0, -1])
    return per


def internal_edges_sum(dist_matrix, loops_df):
    diam = []
    for ind in loops_df.index:
        vert = loops_df.loc[ind, 'vertices']
        temp_dist_m = dist_matrix[vert, :][:, vert]
        diam.append(temp_dist_m.sum() / 2 - loops_df.loc[ind, 'perimeter'])
    return diam

def main() -> None:
    parser = argparse.ArgumentParser(
        description="""Loads homologies dataframes, postprocesses them: adding new features, filter old ones.""",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input-files",
        dest="input_files",
        type=str,
        nargs="+",
        help="""Required: input files or a regular expression matching multiple input files.
               Example: -i /path/to/res.csv /path/to/other/res.csv
               or: -i ./*.csv""",
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
        "--mcool-resolution",
        dest="mcool_resolution",
        type=int,
        help="Optional: a resolution for .mcool files, default is 1000.",
        default=1000,
        metavar="int",
        required=False,
    )
    parser.add_argument(
        "--include-trans-chrom-contacts",
        dest="include_trans_chrom_contacts",
        help="A flag: the distance matrix was computed for the whole genome. If not present, distance matrices are treated for all chromosomes separately",
        action="store_true",
        required=False,
    )
    parser.add_argument(
        "--fetch-fragment",
        dest="fetch_fragment",
        type=str,
        help="Optional: if present, a fragment (chr/chr fragment) is treated for each file.",
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
        logger.error("The given files do not exist or are not .csv files.")
        return

    # Optional arguments
    saved_file_base_folder = args.output_dir
    ony_cis_chrom_contacts = not args.include_trans_chrom_contacts
    resolution = args.mcool_resolution
    fetch_fragment = args.fetch_fragment
    input_file_paths = list_full_paths(input_files)
    logger.success(f"{input_file_paths} are parsed.")
    input_file_paths = filter_list_of_paths(input_file_paths, [".csv"])
    if input_file_paths:
        logger.success(f"{input_file_paths} are kept after filtering.")
    else:
        logger.error("The given folder does not contain .csv files.")

    if ony_cis_chrom_contacts:
        filegroups = group_files_by_prefix(input_file_paths)
    else:
        filegroups = {os.path.basename(file_path):file_path for file_path in input_file_paths}
    # Use a list comprehension to call transform_and_save_matrix for each file with the appropriate arguments
    for prefix_name, full_paths in filegroups.items():
        for full_path in full_paths:
            temp_df = pd.read_csv(full_path)
            saved_indices_path,  all_distances_path = get_new_paths(full_path)
            saved_indices = np.load(saved_indices_path)
            all_distances = np.load(all_distances_path)
            temp_df = temp_df.loc[temp_df.Dim == 1]
            temp_df.Vertices = temp_df.Vertices.apply(lambda x: literal_eval(x))
            temp_df.Vertices = temp_df.Vertices.apply(lambda x: [k-1 for k in x])

            temp_df['start_bin'] = temp_df.Vertices.apply(lambda x: x[0])
            temp_df['end_bin'] = temp_df.Vertices.apply(lambda x: x[-1])
            temp_df['start1'] = np.nan
            temp_df['end1'] = np.nan
            temp_df['chrom1'] = ''
            temp_df['start2'] = np.nan
            temp_df['end2'] = np.nan
            temp_df['chrom2'] = ''
            temp_df.loc[:, 'start'] = temp_df.loc[:, 'start_bin'].apply(lambda x: (chr_name, 10000*dist_to_cont_bin(saved_indices, x), 10000*dist_to_cont_bin(saved_indices, x)+10000))
            temp_df.loc[:, 'end'] = temp_df.loc[:, 'end_bin'].apply(lambda x: (chr_name, 10000*dist_to_cont_bin(saved_indices, x), 10000*dist_to_cont_bin(saved_indices, x)+10000))
            temp_df['chrom1'], temp_df['start1'], temp_df['end1']  = zip(*temp_df['start'].apply(lambda x: ( x[0], x[1], x[2])))
            temp_df['chrom2'], temp_df['start2'], temp_df['end2']  = zip(*temp_df['end'].apply(lambda x: (x[0], x[1], x[2])))
            temp_df.loc[:, 'Death'] = bottleneck_death(all_distances, temp_df)
            temp_df['Lifetime'] = temp_df['Death'] - temp_df['Birth']
            temp_df['max_rad'] = 0
            temp_df.loc[:, 'max_rad'] = maxrad(all_distances, temp_df)
            temp_df['lifetime_cech'] = 0
            temp_df['lifetime_cech'] = temp_df['max_rad'] - temp_df['Birth']
            temp_df['perimeter'] = 0
            temp_df.loc[:, 'perimeter'] = loop_perimeter(all_distances, temp_df)
            temp_df['perimeter_Normalized'] = 0
            temp_df['perimeter_Normalized'] = temp_df['perimeter'] / temp_df['Numvert']
            temp_df['internal_edges'] = 0
            temp_df['internal_edges'] = internal_edges_sum(all_distances, temp_df)
            temp_df['Mean_Diameter'] = 0
            temp_df['Mean_Diameter'] = 2 * temp_df['internal_edges'] / (temp_df['Numvert'] - 3) / temp_df['Numvert']
            temp_df['Perimeter_log'] = np.log10(temp_df['perimeter'])
            temp_df['Internal_Edges_log'] = np.log10(temp_df['internal_edges'])
            temp_df['Numvert_log'] = np.log10(temp_df['Numvert'])
            temp_df['Range_log'] = np.log10(temp_df['Range'])
            results.append(temp_df)
        results = pd.concat(results)
        results

if __name__ == "__main__":
    main()

import argparse
import subprocess
import os
from topology_loop_caller.utils import RESULTS_FOLDER


def pipeline_arg_parsing() -> argparse.Namespace:
    """
    Step 0: Define command-line arguments for the pipeline script
    :return: argparse.Namespace object
    """
    parser = argparse.ArgumentParser(
        description="""
    Command-line tool to execute full pipeline:
    1) load Cooler files & transform them into distance matrices by selected method,
    2) calculate persistent homologies,
    3) (TO DO) filter homologies, generate features.
    """
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
    parser.add_argument("--maxdim", dest="maxdim",type=int, default=2,
                        help="Compute persistent homology in dimensions 0, ..., k.")
    parser.add_argument("--minrad", dest="minrad", type=float, default=-float('inf'),
                        help="Compute homology from time t onward.")
    parser.add_argument("--maxrad", dest="maxrad", type=float, default=float('inf'),
                        help="Stop computing homology after time t.")
    parser.add_argument("--numrad", dest="numrad", type=float, default=float('inf'),
                        help="Divide the interval from minrad to maxrad into N equally spaced steps, and compute the homology of each step. If the value of numrad is set to Inf, then homology will computed at every time point.")
    parser.add_argument("--model", dest="model", type=str, default="vr", choices=["pc", "vr", "complex"],
                        help="Used Eirene model, 'pc' (point cloud), 'vr' (vietoris-rips), or 'complex'.")
    parser.add_argument("--zero-order-homologies-skip", dest="zero_order_homologies_skip", type=bool, default=True,
                        help="Whether to skip zero order homologies.")
    return parser.parse_args()
    return args


def run_first_step(args: argparse.Namespace) -> None:
    """
    Step 1: loading Cooler files, transforming to distance matrices, saving.
    Use subprocess.run to call the Python script with its command-line arguments
    :param args: argparse.Namespace, the output of parser.parse_args()
    :return: None
    """
    subprocess.run(
        [
            "python",
            "matrix_transform.py",
            "--input-dir",
            args.input_dir,
            "--output-dir",
            args.saved_file_base_folder,
            "--distance-function",
            args.distance_function,
            "--save-preserved-bins",
            str(args.save_preserved_bins),
            "--saved-file-prefix",
            str(args.saved_file_prefix),
            "--pearson-sqrt",
            str(args.pearson_sqrt),
            "--log-zero-replacement-strategy",
            args.log_zero_replacement_strategy,
            "--log-base",
            str(args.log_base),
            "--mcool-resolution",
            str(args.mcool_resolution),
            "--to-balance-matrix",
            str(args.to_balance_matrix),
            "--fetch-fragment",
            str(args.fetch_fragment),
        ]
    )


def run_second_step(args: argparse.Namespace) -> None:
    """
    Step 2: Run a Julia script with Topological Data Analysis.
    Use subprocess.run to call the Python script with its command-line arguments
    :param args: argparse.Namespace, the output of parser.parse_args()
    :return: None
    """

    # Define input dir and output dir from basedir:
    matrices_path = os.path.join(args.saved_file_base_folder, 'distance_matrices')
    results_path = os.path.join(args.saved_file_base_folder, 'persistent_homology_results')

    # Create the subdirectories if they don't exist
    os.makedirs(results_path, exist_ok=True)

    subprocess.run(
        [
            "julia",
            "calculate_persistent_homologies.jl",
            "--matrices-path",
            matrices_path,
            "--results-path",
            results_path,
            "--maxdim",
            str(args.maxdim),
            "--minrad",
            str(args.minrad),
            "--maxrad",
            str(args.maxrad),
            "--numrad",
            str(args.numrad),
            "--model",
            args.model,
            "--zero-order-homologies-skip",
            str(args.zero_order_homologies_skip),
        ],
    )


def main():
    args = pipeline_arg_parsing()
    run_first_step(args)
    run_second_step(args)


if __name__ == "__main__":
    main()

import argparse
import subprocess
from topology_loop_caller.utils import RESULTS_FOLDER


def pipeline_arg_parsing():
    # Define command-line arguments for the pipeline script
    parser = argparse.ArgumentParser(
        description="""
    Command-line tool to execute full pipeline:
    1) load Cooler files & transform them into distance matrices by selected method,
    3) calculate persistent homologies,
    4) (TO DO) filter homologies, generate features.
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
    args = parser.parse_args()
    return args


def run_first_step(args):
    # Step 1: loading Cooler files, transforming to distance matrices, saving.
    # Use subprocess.run to call the Python script with its command-line arguments
    subprocess.run(
        [
            "python",
            "step1.py",
            "--input-folder",
            args.input_folder,
            "--string-arg",
            args.string_arg,
            "--float-arg",
            str(args.float_arg),
        ]
    )


def run_second_step(args):
    # Step 2: Run a Julia script
    # Use subprocess.run to call the Julia script with its command-line arguments
    subprocess.run(
        [
            "julia",
            "step2.jl",
            "--input-folder",
            args.input_folder,
            "--string-arg",
            args.string_arg,
            "--float-arg",
            str(args.float_arg),
        ],
        stdout=open(args.output_file, "w"),
    )


def main():
    args = pipeline_arg_parsing()
    run_first_step(args)
    run_second_step(args)


if __name__ == "__main__":
    main()

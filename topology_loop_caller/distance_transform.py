import numpy as np
from topology_loop_caller.utils import timeit


@timeit
def negative_log_transformation(
    M: np.array,
    zero_replacement_strategy: str = "martin_fernandez",
    log_base: float = 10,
):
    """
    Function transforms contact map to (pseudo)distance matrix by log transformation:
    Distances satisfy the first and second metric space axioms and are normalized to [0;1] interval.
    :param M: 2-D np.array with quadratic shape, the given contact map;
    :param zero_replacement_strategy: "martin_fernandez", "half_min" or "unif": zero replacement is
        required for logarithmization. See https://www.sciencedirect.com/science/article/pii/S0169743921000162
    :param log_base: float, a base of logarithm, default is 10.
    :returns: np.array with the same shape as M.
    """
    # Replace all NaNs with zeros
    result = np.nan_to_num(M)
    minimal_val = np.unique(result)[1]
    if zero_replacement_strategy == "half_min":
        # Replace all zeros with the minimal nonzero value * 0.65;
        result = np.where(result == 0.0, minimal_val * 0.65, result)
    elif zero_replacement_strategy == "half_min":
        # Replace all zeros with the minimal nonzero value / 2;
        result = np.where(result == 0.0, minimal_val / 2, result)
    else:
        # Replace all zeros with uniform values from [0.1*DL;DL]
        result = np.where(
            result == 0.0,
            np.random.uniform(0.1 * minimal_val, minimal_val, result.shape),
            result,
        )
    # Set all the elements <= 1;
    result = result / result.max() if result.max() > 1.0 else result
    # Log-transformation: all values are set <= 0;
    result = np.log(result) / np.log(log_base)
    # Now all values are set >= 0;
    result = result * (-1)
    # Now all values are in the [0; 1] interval;
    result = result / result.max()
    # Zero main diagonal for subsequent topological analysis.
    np.fill_diagonal(result, 0.0)
    # A step to get surely symmetric matrix:
    result = (result + result.transpose()) / 2

    return result


@timeit
def pearson_distance(
    balanced_matrix: np.array,
    sqrt: bool = True,
):
    """
    Function transforms contact map to (pseudo)distance matrix by Pearson corr calculation:
    Distances satisfy the first and second metric space axioms and are normalized to [0;1] interval.
    :param balanced_matrix: 2-D np.array with quadratic shape;
    :param sqrt: bool, whether to take a square root of 1 - corr. matrix.
    :returns: np.array with no NaN shape similarity np.array.
    """
    nan_idx = np.isnan(balanced_matrix).all(axis=1)
    no_nan_balanced_matrix = balanced_matrix[~nan_idx][:, ~nan_idx]
    pearson_corr_matrix = np.corrcoef(no_nan_balanced_matrix)
    pearson_corr_matrix = np.nan_to_num(pearson_corr_matrix)
    pearson_corr_matrix = 1 - pearson_corr_matrix
    if sqrt:
        M1_e = np.sqrt(pearson_corr_matrix)
    # A step to get surely symmetric matrix:
    pearson_corr_matrix = (pearson_corr_matrix + pearson_corr_matrix.transpose()) / 2
    return pearson_corr_matrix

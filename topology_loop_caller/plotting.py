import numpy as np
import matplotlib.pyplot as plt


def plot_hic_matrix_fragment(
    matrix: np.array,
    subset=None,
    log_scaled: bool = True,
    lop1p_add: float = 1,
    ax=None,
    figsize: tuple = (10, 10),
    cmap: str = "RdBu",
    interpolation: str = "none",
    plot_title: str = None,
    title_fontsize: float = 12,
):
    """
    The function visualizes Hi-C contact map fragment.
    :param matrix: numpy 2-D array
    :param subset: the region of interest, tuple(x_start, y_start, n_bins). If None, the whole matrix is plotted
    :param log_scaled: bool, whether to perform log-normalization
    :param lop1p_add: float, the value to add before log-normalization
    :param ax: matplotlib.axes
    :param figsize: figsize if ax=None
    :param cmap: plot color map
    :param interpolation: interpolation method
    :param plot_title: str, plot title text
    :param title_fontsize: plot title fontsize
    :returns: ax
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    if subset:
        matrix_subset = matrix[
            subset[0] : subset[0] + subset[2], subset[1] : subset[1] + subset[2]
        ]
    else:
        matrix_subset = matrix[:, :]
    if log_scaled:
        ax.imshow(
            np.log(matrix_subset + lop1p_add), cmap=cmap, interpolation=interpolation
        )
    else:
        ax.imshow(matrix_subset, cmap=cmap, interpolation=interpolation)
    if plot_title:
        ax.set_title(plot_title, fontsize=title_fontsize)
    return ax


def genomic_dist_interaction_freq_plot(
    matrix: np.array,
    # chromsizes, TODO: add only cis contacts
    ax=None,
    x_log_scaled: bool = True,
    y_log_scaled: bool = True,
    figsize: tuple = (10, 10),
    color="b",
    lw: float = 2,
    plot_title: str = None,
    title_fontsize: float = 12,
):
    """
    The function visualizes a relation between interaction frequencies and genomic distance (cischrom).
    :param matrix: numpy 2-D array, normally balanced matrix
    :param chromsizes: chromsizes in cooler.chromsizes format in n_bins
    :param x_log_scaled: bool, whether to perform log-normalization of x-axis
    :param y_log_scaled: bool, whether to perform log-normalization of y-axis
    :param ax: matplotlib.axes
    :param figsize: figsize if ax=None
    :param plot_title: str, plot title text
    :param title_fontsize: plot title fontsize
    :returns: ax
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    # Deleting NaN rows in case of balanced matrix
    ind = np.isnan(matrix).all(axis=1)
    no_nan_matrix = matrix[~ind][:, ~ind]
    mean_frequences = [
        np.mean(no_nan_matrix.diagonal(i)) for i in range(no_nan_matrix.shape[0])
    ]
    ax.plot(mean_frequences, color=color, lw=lw)
    if x_log_scaled:
        ax.set_xscale("log")
    if y_log_scaled:
        ax.set_yscale("log")

    if plot_title:
        ax.set_title(plot_title, fontsize=title_fontsize)
    return ax

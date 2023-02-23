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

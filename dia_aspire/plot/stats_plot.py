from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_target_decoy_dist(
    data: pd.DataFrame,
    metric: str,
    decoy_col: str = "decoy",
    target_label: str = "Target",
    decoy_label: str = "Decoy",
    bins: int = 50,
    kde: bool = True,
    stat: str = "density",
    element: str = "step",
    alpha: float = 0.4,
    ax: Optional[plt.Axes] = None,
    **kwargs,
):
    """
    Plot histogram and KDE for target and decoy distributions.

    This function creates a comparison plot showing the distribution
    of a metric for target and decoy PSMs/peptides.

    Parameters
    ----------
    data : pandas.DataFrame
        The input data containing both target and decoy entries.
    metric : str
        The column name of the metric to plot.
    decoy_col : str, optional
        The column name indicating decoy status (0 for target, 1 for decoy).
        Default is "decoy".
    target_label : str, optional
        Label for target entries in the legend. Default is "Target".
    decoy_label : str, optional
        Label for decoy entries in the legend. Default is "Decoy".
    bins : int, optional
        Number of bins for histogram. Default is 50.
    kde : bool, optional
        Whether to overlay KDE curve. Default is True.
    stat : str, optional
        Aggregate statistic to compute in each bin. Default is "density".
        Options: "count", "frequency", "probability", "density", "percent".
    element : str, optional
        Visual representation of the histogram. Default is "step".
        Options: "bars", "step", "poly".
    alpha : float, optional
        Transparency level for histogram. Default is 0.4.
    ax : matplotlib.pyplot.Axes, optional
        The matplotlib Axes on which to plot. If `None`, the current
        Axes instance is used.
    **kwargs : dict, optional
        Additional arguments passed to :py:func:`seaborn.histplot`.

    Returns
    -------
    matplotlib.pyplot.Axes
        An :py:class:`matplotlib.axes.Axes` with the target-decoy
        distribution plot.

    Examples
    --------
    >>> import pandas as pd
    >>> import matplotlib.pyplot as plt
    >>> from dia_aspire_rescore.plot import plot_target_decoy_dist
    >>>
    >>> # Create sample data
    >>> data = pd.DataFrame({
    ...     'score': [0.8, 0.9, 0.3, 0.4, 0.85, 0.35],
    ...     'decoy': [0, 0, 1, 1, 0, 1]
    ... })
    >>>
    >>> # Single plot
    >>> fig, ax = plt.subplots()
    >>> plot_target_decoy_dist(data, 'score', ax=ax)
    >>> plt.show()
    >>>
    >>> # Multiple subplots (2x2 grid)
    >>> metrics = ['score1', 'score2', 'score3', 'score4']
    >>> fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    >>> for idx, metric in enumerate(metrics):
    ...     ax = axes[idx // 2, idx % 2]
    ...     plot_target_decoy_dist(data, metric, ax=ax)
    >>> plt.tight_layout()
    >>> plt.show()
    """
    import seaborn as sns

    if ax is None:
        ax = plt.gca()

    default_colors = sns.color_palette()
    palette = {target_label: default_colors[0], decoy_label: default_colors[1]}

    # Prepare data with Type column
    plot_data = data.copy()
    plot_data["Type"] = plot_data[decoy_col].map({0: target_label, 1: decoy_label})

    # Plot histogram with KDE
    sns.histplot(
        data=plot_data,
        x=metric,
        hue="Type",
        hue_order=[target_label, decoy_label],
        palette=palette,
        ax=ax,
        bins=bins,
        stat=stat,
        common_norm=False,  # Normalize each group separately
        kde=kde,
        element=element,
        alpha=alpha,
        **kwargs,
    )

    ax.set_xlabel(metric)
    ax.set_ylabel(stat.capitalize())

    return ax


def plot_qvalues(
    qvalues: np.ndarray,
    threshold: float = 0.1,
    ax: Optional[plt.Axes] = None,
    **kwargs,
):
    """
    Plot the cumulative number of discoveries over range of q-values.

    This implementation is adapted from the Mokapot package.
    https://github.com/wfondrie/mokapot

    Parameters
    ----------
    qvalues : numpy.ndarray
        The q-values to plot.
    threshold : float, optional
        Indicates the maximum q-value to plot.
    ax : matplotlib.pyplot.Axes, optional
        The matplotlib Axes on which to plot. If `None` the current
        Axes instance is used.
    **kwargs : dict, optional
        Arguments passed to :py:func:`matplotlib.axes.Axes.plot`.

    Returns
    -------
    matplotlib.pyplot.Axes
        An :py:class:`matplotlib.axes.Axes` with the cumulative
        number of accepted target PSMs or peptides.
    """

    if ax is None:
        ax = plt.gca()

    # Calculate cumulative targets at each q-value
    qvals = pd.Series(qvalues, name="qvalue")
    qvals = qvals.sort_values(ascending=True).to_frame()
    qvals["target"] = 1
    qvals["num"] = qvals["target"].cumsum()
    qvals = qvals.groupby(["qvalue"]).max().reset_index()
    qvals = qvals[["qvalue", "num"]]

    zero = pd.DataFrame({"qvalue": qvals["qvalue"][0], "num": 0}, index=[-1])
    qvals = pd.concat([zero, qvals], sort=True).reset_index(drop=True)

    xmargin = threshold * 0.05
    ymax = qvals.num[qvals["qvalue"] <= (threshold + xmargin)].max()
    ymargin = ymax * 0.05

    # Set margins
    curr_ylims = ax.get_ylim()
    if curr_ylims[1] < ymax + ymargin:
        ax.set_ylim(0 - ymargin, ymax + ymargin)

    ax.set_xlim(0 - xmargin, threshold + xmargin)
    ax.set_xlabel("q-value")
    ax.set_ylabel("Discoveries")

    ax.step(qvals["qvalue"].values, qvals.num.values, where="post", **kwargs)

    return ax

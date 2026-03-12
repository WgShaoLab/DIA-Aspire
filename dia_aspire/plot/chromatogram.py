from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from alphabase.constants.atom import MASS_ISOTOPE
from alphabase.peptide.fragment import create_fragment_mz_dataframe, get_charged_frag_types

from dia_aspire.extraction import extract_xic

# TODO: color palette?
# TODO: ms2 top k fragments?


def plot_xic_ms1(
    psm: dict | pd.Series,
    spectrum_df: pd.DataFrame,
    peak_df: pd.DataFrame,
    n_isotopes: int = 2,
    rt_extension: float = 10.0,
    show_apex: bool = True,
    show_boundary: bool = True,
    show_legend: bool = True,
    ppm_tolerance: float = 20.0,
    ax: plt.Axes | None = None,
    **kwargs,
) -> plt.Axes:
    """
    Plot MS1 XIC for precursor isotope peaks.

    This function extracts and visualizes MS1 XIC traces for multiple
    isotope peaks (M, M+1, M+2, etc.) of a given PSM.

    Parameters
    ----------
    psm : dict | pd.Series
        Single PSM containing:
        - `precursor_mz`: Precursor m/z value
        - `charge`: Charge state
        - `rt`: Retention time in minutes
        - `rt_start`: Peak start RT in minutes
        - `rt_stop`: Peak stop RT in minutes
    spectrum_df : pd.DataFrame
        AlphaRaw spectrum DataFrame with RT in minutes.
    peak_df : pd.DataFrame
        AlphaRaw peak DataFrame.
    n_isotopes : int, optional
        Number of isotope peaks to extract. For example:
        - `n_isotopes=2`: Extract M and M+1
        - `n_isotopes=3`: Extract M, M+1, M+2
        Default is 2.
    rt_extension : float, optional
        RT window extension in seconds. The extraction window will be
        `[rt_start*60 - rt_extension, rt_stop*60 + rt_extension]` seconds.
        Default is 10.0 seconds.
    show_apex : bool, optional
        Whether to mark the apex (maximum intensity point) for each isotope.
        Default is True.
    show_boundary : bool, optional
        Whether to show vertical lines at `rt_start` and `rt_stop` boundaries.
        Default is True.
    show_legend : bool, optional
        Whether to display legend with isotope labels (M, M+1, M+2, etc.).
        Default is True.
    ppm_tolerance : float, optional
        m/z tolerance in ppm for peak matching. Default is 20.0.
    ax : matplotlib.pyplot.Axes, optional
        The matplotlib Axes on which to plot. If `None`, the current
        Axes instance is used.
    **kwargs : dict, optional
        Additional arguments passed to `matplotlib.axes.Axes.plot`.

    Returns
    -------
    matplotlib.pyplot.Axes
        The matplotlib Axes with the MS1 XIC plot. X-axis is in seconds.

    Examples
    --------
    >>> import pandas as pd
    >>> import matplotlib.pyplot as plt
    >>> from dia_aspire_rescore.plot import plot_xic_ms1
    >>>
    >>> # Single PSM MS1 XIC
    >>> psm = psm_df.iloc[0]
    >>> ax = plot_xic_ms1(psm, spectrum_df, peak_df, n_isotopes=3, rt_extension=10.0)
    >>> ax.set_xlabel('Retention Time (s)')
    >>> plt.show()
    >>>
    >>> # In subplot
    >>> fig, ax = plt.subplots()
    >>> plot_xic_ms1(psm, spectrum_df, peak_df, n_isotopes=2, ax=ax)
    >>> plt.tight_layout()
    >>> plt.show()
    """
    if ax is None:
        ax = plt.gca()

    precursor_mz = float(psm["precursor_mz"])
    charge = int(psm["charge"])
    rt_min = float(psm["rt"])
    rt_start_min = float(psm["rt_start"])
    rt_stop_min = float(psm["rt_stop"])

    isotope_mzs = np.array([precursor_mz + i * MASS_ISOTOPE / charge for i in range(n_isotopes)])

    rt_start_sec = rt_start_min * 60 - rt_extension
    rt_stop_sec = rt_stop_min * 60 + rt_extension

    rt_start_extract = rt_start_sec / 60
    rt_stop_extract = rt_stop_sec / 60

    # TODO: color palette?
    colors = plt.get_cmap("tab10")(np.linspace(0, 1, n_isotopes))
    rt_sec = rt_min * 60  # in seconds

    for i, isotope_mz in enumerate(isotope_mzs):
        rt_values_min, intensities = extract_xic(
            spectrum_df=spectrum_df,
            peak_df=peak_df,
            query_mzs=np.array([isotope_mz]),
            rt_start=rt_start_extract,
            rt_stop=rt_stop_extract,
            precursor_mz=None,  # not needed for MS1
            ppm_tolerance=ppm_tolerance,
            ms_level=1,
        )

        if len(rt_values_min) == 0:
            continue

        rt_values_sec = rt_values_min * 60
        intensity = intensities[0]

        label = f"[M+{i}]^{charge} ({isotope_mzs[i]:.2f} m/z)" if i > 0 else f"[M]^{charge} ({isotope_mzs[0]:.2f} m/z)"
        ax.plot(
            rt_values_sec,
            intensity,
            color=colors[i],
            label=label,
            linewidth=kwargs.get("linewidth", 1.5),
            **{k: v for k, v in kwargs.items() if k != "linewidth"},
        )

    if show_apex:
        ax.axvline(
            rt_sec,
            color="gray",
            linestyle="--",
            linewidth=2,
            alpha=0.7,
        )

    if show_boundary:
        rt_start_sec_boundary = rt_start_min * 60
        rt_stop_sec_boundary = rt_stop_min * 60
        ax.axvline(
            rt_start_sec_boundary,
            color="red",
            linestyle="--",
            linewidth=2,
            alpha=0.7,
        )
        ax.axvline(
            rt_stop_sec_boundary,
            color="red",
            linestyle="--",
            linewidth=2,
            alpha=0.7,
        )

    ax.set_xlabel("Retention Time (s)", fontsize=12)
    ax.set_ylabel("Intensity", fontsize=12)

    if show_legend:
        ax.legend(fontsize=10, loc="best")

    ax.grid(True, alpha=0.3)

    return ax


def plot_xic_ms2(
    psm: dict | pd.Series,
    spectrum_df: pd.DataFrame,
    peak_df: pd.DataFrame,
    rt_extension: float = 10.0,
    show_apex: bool = True,
    show_boundary: bool = True,
    show_legend: bool = True,
    charged_frag_types: list[str] | None = None,
    ppm_tolerance: float = 20.0,
    ax: plt.Axes | None = None,
    **kwargs,
) -> plt.Axes:
    """
    Plot MS2 XIC (Extracted Ion Chromatogram) for fragment ions.

    This function extracts and visualizes MS2 XIC traces for all fragment
    ions of a given PSM. Each fragment type is plotted as a separate line.

    Parameters
    ----------
    psm : dict | pd.Series
        Single PSM containing:
        - `sequence`: Peptide sequence
        - `mod`: Modification string (can be empty)
        - `mod_sites`: Modification sites string (can be empty)
        - `charge`: Charge state
        - `precursor_mz`: Precursor m/z value
        - `rt`: Retention time in minutes
        - `rt_start`: Peak start RT in minutes
        - `rt_stop`: Peak stop RT in minutes
    spectrum_df : pd.DataFrame
        AlphaRaw spectrum DataFrame with RT in minutes.
    peak_df : pd.DataFrame
        AlphaRaw peak DataFrame.
    rt_extension : float, optional
        RT window extension in seconds. The extraction window will be
        `[rt_start*60 - rt_extension, rt_stop*60 + rt_extension]` seconds.
        Default is 10.0 seconds.
    show_apex : bool, optional
        Whether to mark the apex (maximum intensity point) for each fragment.
        Default is True.
    show_boundary : bool, optional
        Whether to show vertical lines at `rt_start` and `rt_stop` boundaries.
        Default is True.
    show_legend : bool, optional
        Whether to display legend with fragment type labels (e.g., b_z1, y_z1).
        Default is True.
    charged_frag_types : list[str] | None, optional
        Fragment types to generate. If None, defaults to
        `get_charged_frag_types(["b", "y"], 2)` (b/z1, b/z2, y/z1, y/z2).
    ppm_tolerance : float, optional
        m/z tolerance in ppm for peak matching. Default is 20.0.
    ax : matplotlib.pyplot.Axes, optional
        The matplotlib Axes on which to plot. If `None`, the current
        Axes instance is used.
    **kwargs : dict, optional
        Additional arguments passed to `matplotlib.axes.Axes.plot`.

    Returns
    -------
    matplotlib.pyplot.Axes
        The matplotlib Axes with the MS2 XIC plot. X-axis is in seconds.

    Examples
    --------
    >>> import pandas as pd
    >>> import matplotlib.pyplot as plt
    >>> from dia_aspire_rescore.plot import plot_xic_ms2
    >>>
    >>> # Single PSM MS2 XIC
    >>> psm = psm_df.iloc[0]
    >>> ax = plot_xic_ms2(psm, spectrum_df, peak_df, rt_extension=10.0, show_apex=True)
    >>> ax.set_xlabel('Retention Time (s)')
    >>> plt.show()
    >>>
    >>> # In subplot
    >>> fig, ax = plt.subplots()
    >>> plot_xic_ms2(psm, spectrum_df, peak_df, ax=ax)
    >>> plt.tight_layout()
    >>> plt.show()
    """
    if ax is None:
        ax = plt.gca()

    precursor_mz = float(psm["precursor_mz"])
    rt_min = float(psm["rt"])
    rt_start_min = float(psm["rt_start"])
    rt_stop_min = float(psm["rt_stop"])

    psm_for_frag = pd.DataFrame([psm]).copy()

    # ensure required columns exist and are string type for create_fragment_mz_dataframe
    if "mods" not in psm_for_frag.columns:
        psm_for_frag["mods"] = ""
    else:
        psm_for_frag["mods"] = psm_for_frag["mods"].fillna("")
    psm_for_frag["mods"] = psm_for_frag["mods"].astype(str)

    if "mod_sites" not in psm_for_frag.columns:
        psm_for_frag["mod_sites"] = ""
    else:
        psm_for_frag["mod_sites"] = psm_for_frag["mod_sites"].fillna("")
    psm_for_frag["mod_sites"] = psm_for_frag["mod_sites"].astype(str)

    if charged_frag_types is None:
        charged_frag_types = get_charged_frag_types(["b", "y"], 2)

    fragment_mz_df = create_fragment_mz_dataframe(psm_for_frag, charged_frag_types)

    frag_start = int(psm_for_frag.iloc[0]["frag_start_idx"])
    frag_stop = int(psm_for_frag.iloc[0]["frag_stop_idx"])

    frag_mzs = fragment_mz_df.iloc[frag_start:frag_stop]

    all_frag_mzs = frag_mzs.values.flatten()
    query_mzs = all_frag_mzs[all_frag_mzs > 0]

    if len(query_mzs) == 0:
        ax.text(
            0.5,
            0.5,
            "No fragments found",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=14,
        )
        return ax

    rt_start_sec = rt_start_min * 60 - rt_extension
    rt_stop_sec = rt_stop_min * 60 + rt_extension

    rt_start_extract = rt_start_sec / 60
    rt_stop_extract = rt_stop_sec / 60

    rt_values_min, intensities = extract_xic(
        spectrum_df=spectrum_df,
        peak_df=peak_df,
        query_mzs=query_mzs,
        rt_start=rt_start_extract,
        rt_stop=rt_stop_extract,
        precursor_mz=precursor_mz,
        ppm_tolerance=ppm_tolerance,
        ms_level=2,
    )

    if len(rt_values_min) == 0:
        ax.text(
            0.5,
            0.5,
            "No XIC data found",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=14,
        )
        return ax

    rt_values_sec = rt_values_min * 60
    colors = plt.get_cmap("tab20")(np.linspace(0, 1, len(query_mzs)))
    rt_sec = rt_min * 60

    _plot_fragment_traces(ax, frag_mzs, query_mzs, intensities, rt_values_sec, colors, **kwargs)

    if show_apex:
        ax.axvline(
            rt_sec,
            color="gray",
            linestyle="--",
            linewidth=2,
            alpha=0.7,
        )

    if show_boundary:
        rt_start_sec_boundary = rt_start_min * 60
        rt_stop_sec_boundary = rt_stop_min * 60
        ax.axvline(
            rt_start_sec_boundary,
            color="red",
            linestyle="--",
            linewidth=2,
            alpha=0.7,
        )
        ax.axvline(
            rt_stop_sec_boundary,
            color="red",
            linestyle="--",
            linewidth=2,
            alpha=0.7,
        )

    ax.set_xlabel("Retention Time (s)", fontsize=12)
    ax.set_ylabel("Intensity", fontsize=12)

    if show_legend:
        ax.legend(fontsize=9, loc="best", ncol=2)

    ax.grid(True, alpha=0.3)

    return ax


def _parse_fragment_label(col_name: str, row_idx: int, frag_mz: float) -> str:
    """Parse fragment column name to generate label."""
    if "_z" in col_name:
        frag_type, charge_str = col_name.split("_z")
        charge_num = int(charge_str)
    else:
        frag_type = col_name[0]
        charge_num = 1

    frag_num = row_idx + 1
    return f"{frag_type}{frag_num}^{charge_num} ({frag_mz:.2f} m/z)"


def _plot_fragment_traces(
    ax: plt.Axes,
    frag_mzs: pd.DataFrame,
    query_mzs: np.ndarray,
    intensities: np.ndarray,
    rt_values_sec: np.ndarray,
    colors: np.ndarray,
    **kwargs,
) -> None:
    """Plot XIC traces for all fragments."""
    frag_cols = frag_mzs.columns.tolist()
    frag_idx = 0

    for col_name in frag_cols:
        col_frag_mzs = frag_mzs[col_name].values
        for row_idx, frag_mz in enumerate(col_frag_mzs):
            if frag_mz <= 0:
                continue

            frag_idx_in_query = np.where(np.abs(query_mzs - frag_mz) < 1e-6)[0]
            if len(frag_idx_in_query) == 0:
                continue
            frag_idx_in_query = frag_idx_in_query[0]

            intensity = intensities[frag_idx_in_query]
            label = _parse_fragment_label(col_name, row_idx, frag_mz)

            ax.plot(
                rt_values_sec,
                intensity,
                color=colors[frag_idx],
                label=label,
                linewidth=kwargs.get("linewidth", 1),
                **{k: v for k, v in kwargs.items() if k != "linewidth"},
            )

            frag_idx += 1

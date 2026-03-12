"""Core numba-accelerated functions for XIC extraction."""

from __future__ import annotations

import numba as nb
import numpy as np
import pandas as pd
from alpharaw.match.match_utils import match_closest_peaks, match_highest_peaks

from dia_aspire.constants.spectrum import SpectrumDfCols


@nb.njit()
def get_spec_idxes_in_rt_window(
    spec_rts: np.ndarray,
    rt_start: float,
    rt_stop: float,
) -> np.ndarray:
    """
    Get spectrum indices within the RT window.

    Parameters
    ----------
    spec_rts : np.ndarray
        RT values of spectra (must be sorted).
    rt_start : float
        Start of RT window (in the same unit as spec_rts).
    rt_stop : float
        End of RT window (in the same unit as spec_rts).

    Returns
    -------
    np.ndarray
        Array of spectrum indices within the RT window.
    """
    start_idx = np.searchsorted(spec_rts, rt_start, side="left")
    stop_idx = np.searchsorted(spec_rts, rt_stop, side="left")
    return np.arange(start_idx, stop_idx, dtype=np.int64)


@nb.njit()
def filter_spec_idxes_by_isolation(
    spec_idxes: np.ndarray,
    isolation_lower_mzs: np.ndarray,
    isolation_upper_mzs: np.ndarray,
    precursor_mz: float,
) -> np.ndarray:
    """
    Filter spectrum indices by isolation window containing the precursor m/z.

    Parameters
    ----------
    spec_idxes : np.ndarray
        Spectrum indices to filter.
    isolation_lower_mzs : np.ndarray
        Lower bounds of isolation windows for all spectra.
    isolation_upper_mzs : np.ndarray
        Upper bounds of isolation windows for all spectra.
    precursor_mz : float
        Precursor m/z to check.

    Returns
    -------
    np.ndarray
        Filtered spectrum indices where isolation window contains precursor_mz.
    """
    # TODO: prealloc len(spec_idxes) and return out[:j] instead?
    count = 0
    for idx in spec_idxes:
        if isolation_lower_mzs[idx] <= precursor_mz <= isolation_upper_mzs[idx]:
            count += 1
    result = np.empty(count, dtype=np.int64)

    j = 0
    for idx in spec_idxes:
        if isolation_lower_mzs[idx] <= precursor_mz <= isolation_upper_mzs[idx]:
            result[j] = idx
            j += 1

    return result


@nb.njit()
def filter_spec_idxes_by_ms_level(
    spec_idxes: np.ndarray,
    ms_levels: np.ndarray,
    target_ms_level: int,
) -> np.ndarray:
    """
    Filter spectrum indices by MS level.

    Parameters
    ----------
    spec_idxes : np.ndarray
        Spectrum indices to filter.
    ms_levels : np.ndarray
        MS levels for all spectra.
    target_ms_level : int
        Target MS level (1 for MS1, 2 for MS2).

    Returns
    -------
    np.ndarray
        Filtered spectrum indices with the target MS level.
    """
    # TODO: prealloc len(spec_idxes) and return out[:j] instead?
    count = 0
    for idx in spec_idxes:
        if ms_levels[idx] == target_ms_level:
            count += 1
    result = np.empty(count, dtype=np.int64)

    j = 0
    for idx in spec_idxes:
        if ms_levels[idx] == target_ms_level:
            result[j] = idx
            j += 1

    return result


@nb.njit()
def extract_xic_for_mzs(
    spec_idxes: np.ndarray,
    query_mzs: np.ndarray,
    query_mz_tols: np.ndarray,
    peak_mzs: np.ndarray,
    peak_intensities: np.ndarray,
    peak_start_idxes: np.ndarray,
    peak_stop_idxes: np.ndarray,
    match_closest: bool = True,
) -> np.ndarray:
    """
    Extract XIC intensity matrix for multiple m/z queries.

    Parameters
    ----------
    spec_idxes : np.ndarray
        Spectrum indices to extract from.
    query_mzs : np.ndarray
        Query m/z values, shape (n_ions,).
    query_mz_tols : np.ndarray
        Absolute m/z tolerances, shape (n_ions,).
    peak_mzs : np.ndarray
        All peak m/z values from peak_df.
    peak_intensities : np.ndarray
        All peak intensity values from peak_df.
    peak_start_idxes : np.ndarray
        Peak start indices for each spectrum.
    peak_stop_idxes : np.ndarray
        Peak stop indices for each spectrum.
    match_closest : bool, default=True
        Whether to match the closest peak instead of the highest intensity peak.

    Returns
    -------
    np.ndarray
        Intensity matrix, shape (n_ions, n_specs).
    """
    n_ions = len(query_mzs)
    n_specs = len(spec_idxes)

    intensities = np.zeros((n_ions, n_specs), dtype=np.float32)

    for i in range(n_specs):
        spec_idx = spec_idxes[i]
        peak_start = peak_start_idxes[spec_idx]
        peak_stop = peak_stop_idxes[spec_idx]

        spec_peak_mzs = peak_mzs[peak_start:peak_stop]
        spec_peak_intensities = peak_intensities[peak_start:peak_stop]

        if len(spec_peak_mzs) == 0:
            continue

        if match_closest:
            matched_indices = match_closest_peaks(
                spec_peak_mzs,
                spec_peak_intensities,
                query_mzs,
                query_mz_tols,
            )
        else:
            matched_indices = match_highest_peaks(
                spec_peak_mzs,
                spec_peak_intensities,
                query_mzs,
                query_mz_tols,
            )

        for j in range(n_ions):
            matched_idx = matched_indices[j]
            if matched_idx >= 0:
                intensities[j, i] = spec_peak_intensities[matched_idx]

    return intensities


def extract_xic(
    spectrum_df: pd.DataFrame,
    peak_df: pd.DataFrame,
    query_mzs: np.ndarray,
    rt_start: float,
    rt_stop: float,
    precursor_mz: float | None = None,
    ppm_tolerance: float = 20.0,
    ms_level: int = 2,
    match_closest: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """
    High-level API for XIC extraction.

    Parameters
    ----------
    spectrum_df : pd.DataFrame
        AlphaRaw spectrum_df containing RT, isolation windows, peak indices.
    peak_df : pd.DataFrame
        AlphaRaw peak_df containing m/z and intensity values.
    query_mzs : np.ndarray
        Query m/z values (can be from alphabase for any ion type).
    rt_start : float
        Start of RT window (in the same unit as spectrum_df.rt).
    rt_stop : float
        End of RT window (in the same unit as spectrum_df.rt).
    precursor_mz : float, optional
        Precursor m/z for filtering by DIA isolation window (required for MS2).
    ppm_tolerance : float, default=20.0
        m/z tolerance in ppm.
    ms_level : int, default=2
        1 for MS1 (precursor), 2 for MS2 (fragment).
    match_closest : bool, default=True
        Whether to match the closest peak instead of the highest intensity peak.

    Returns
    -------
    rt_values : np.ndarray
        RT values for each extracted point.
    intensities : np.ndarray
        Intensity matrix, shape (n_ions, n_specs).
    """
    query_mzs = np.asarray(query_mzs, dtype=np.float64)
    query_mz_tols = query_mzs * ppm_tolerance * 1e-6

    spec_rts = spectrum_df[SpectrumDfCols.RT].values
    ms_levels = spectrum_df[SpectrumDfCols.MS_LEVEL].values
    peak_start_idxes = spectrum_df[SpectrumDfCols.PEAK_START_IDX].values
    peak_stop_idxes = spectrum_df[SpectrumDfCols.PEAK_STOP_IDX].values

    peak_mzs = peak_df["mz"].values
    peak_intensities = peak_df["intensity"].values

    spec_idxes = get_spec_idxes_in_rt_window(spec_rts, rt_start, rt_stop)

    if len(spec_idxes) == 0:
        return np.array([], dtype=np.float64), np.zeros((len(query_mzs), 0), dtype=np.float32)

    spec_idxes = filter_spec_idxes_by_ms_level(spec_idxes, ms_levels, ms_level)

    if len(spec_idxes) == 0:
        return np.array([], dtype=np.float64), np.zeros((len(query_mzs), 0), dtype=np.float32)

    if ms_level == 2:
        if precursor_mz is None:
            raise ValueError("precursor_mz is required for MS2 XIC extraction")

        isolation_lower_mzs = spectrum_df[SpectrumDfCols.ISOLATION_LOWER_MZ].values
        isolation_upper_mzs = spectrum_df[SpectrumDfCols.ISOLATION_UPPER_MZ].values

        spec_idxes = filter_spec_idxes_by_isolation(
            spec_idxes,
            isolation_lower_mzs,
            isolation_upper_mzs,
            precursor_mz,
        )

    if len(spec_idxes) == 0:
        return np.array([], dtype=np.float64), np.zeros((len(query_mzs), 0), dtype=np.float32)

    intensities = extract_xic_for_mzs(
        spec_idxes,
        query_mzs,
        query_mz_tols,
        peak_mzs,
        peak_intensities,
        peak_start_idxes,
        peak_stop_idxes,
        match_closest,
    )

    rt_values = spec_rts[spec_idxes]

    return rt_values, intensities

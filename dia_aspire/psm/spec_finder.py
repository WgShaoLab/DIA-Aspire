import numba
import numpy as np
import pandas as pd

from dia_aspire.constants.spectrum import SpectrumDfCols


@numba.njit
def find_dia_spec_idxes_same_window(
    spec_rt_values: np.ndarray,
    query_rt_values: np.ndarray,
    max_spec_per_query: int,
) -> np.ndarray:
    """
    For given array of query RT values, find spectrum indices
    from the subset of spectra within the same normal DIA m/z window.
    This function is numba accelerated.

    Adapted from alpharaw: https://github.com/MannLabs/alpharaw

    Parameters
    ----------
    spec_rt_values : np.ndarray
        RT values of given DIA spectra.
    query_rt_values : np.ndarray
        Query RT values.
    max_spec_per_query : int
        Return maximal spectrum indices (scan windows) for the given query.
        Must be less than or equal to the number of spectra.

    Returns
    -------
    ndarray[int32]
        Result spectrum indices with shape (query num, max_spec_per_query).
    """
    rt_idxes = np.searchsorted(spec_rt_values, query_rt_values)

    spec_idxes = np.full((len(rt_idxes), max_spec_per_query), -1, dtype=np.int32)
    n = max_spec_per_query // 2

    # TODO: closest RT
    for iquery in range(len(rt_idxes)):
        if rt_idxes[iquery] < n:
            spec_idxes[iquery, :] = np.arange(0, max_spec_per_query)
        elif rt_idxes[iquery] + n >= len(spec_rt_values):  # fix the boundary issue
            spec_idxes[iquery, :] = np.arange(len(spec_rt_values) - max_spec_per_query, len(spec_rt_values))
        else:
            spec_idxes[iquery, :] = np.arange(rt_idxes[iquery] - n, rt_idxes[iquery] - n + max_spec_per_query)
    return spec_idxes


def find_DIA_spec_idxes_by_rt(
    spectrum_df: pd.DataFrame,
    query_rts: np.ndarray,
    query_precursor_mzs: np.ndarray,
    rt_tolerance: float = -1.0,
) -> np.ndarray:
    """
    Find MS2 spectrum indices (int32) from the `spectrum_df`,
    given the `query_rts` and `query_precursor_mzs`,
    ensuring the isolation window contains the precursor m/z.

    Parameters
    ----------
    spectrum_df : pd.DataFrame
        DataFrame containing spectrum information with columns:
        SpectrumDfCols.RT, SpectrumDfCols.ISOLATION_LOWER_MZ, SpectrumDfCols.ISOLATION_UPPER_MZ.
    query_rts : np.ndarray
        Query RT values.
    query_precursor_mzs : np.ndarray
        Query precursor m/z values.
    rt_tolerance : float, optional
        RT tolerance. If negative (default -1.0), no RT check is performed.

    Returns
    -------
    np.ndarray
        Array of spectrum indices. -1 indicates no valid spectrum found.
    """
    spec_rts = spectrum_df[SpectrumDfCols.RT].to_numpy()
    spec_isolation_lower_mzs = spectrum_df[SpectrumDfCols.ISOLATION_LOWER_MZ].to_numpy()
    spec_isolation_upper_mzs = spectrum_df[SpectrumDfCols.ISOLATION_UPPER_MZ].to_numpy()

    return find_batch_DIA_spec_idxes_by_rt(
        spec_rts,
        spec_isolation_lower_mzs,
        spec_isolation_upper_mzs,
        query_rts,
        query_precursor_mzs,
        rt_tolerance,
    )


@numba.njit
def find_single_DIA_spec_idx_by_rt(
    spec_rts: np.ndarray,
    spec_isolation_lower_mzs: np.ndarray,
    spec_isolation_upper_mzs: np.ndarray,
    query_rt: float,
    query_precursor_mz: float,
    rt_tolerance: float = -1.0,
) -> np.int32:
    """
    Find the closest MS2 spectrum index for DIA data by RT,
    ensuring the isolation window contains the precursor m/z.
    This function is numba accelerated.

    Parameters
    ----------
    spec_rts : np.ndarray
        RT values of the spectra (must be sorted).
    spec_isolation_lower_mzs : np.ndarray
        Left m/z values of the isolation windows.
    spec_isolation_upper_mzs : np.ndarray
        Right m/z values of the isolation windows.
    query_rt : float
        Query RT value.
    query_precursor_mz : float
        Query precursor m/z value.
    rt_tolerance : float, optional
        RT tolerance. If negative (default -1.0), no RT check is performed.

    Returns
    -------
    np.int32
        The spectrum index. Returns -1 if no valid spectrum found.
    """
    rt_idx = int(np.searchsorted(spec_rts, query_rt))

    # Determine the two closest spectrum indices
    # TODO: check if is MS1
    # TODO: in ms1, both upper and lower isolation m/z are -1
    if rt_idx == 0:
        candidates = [0, 1] if len(spec_rts) > 1 else [0]
    elif rt_idx >= len(spec_rts):
        candidates = [len(spec_rts) - 2, len(spec_rts) - 1] if len(spec_rts) > 1 else [len(spec_rts) - 1]
    else:
        left_idx = rt_idx - 1
        right_idx = rt_idx
        left_dist = abs(spec_rts[left_idx] - query_rt)
        right_dist = abs(spec_rts[right_idx] - query_rt)
        candidates = [left_idx, right_idx] if left_dist <= right_dist else [right_idx, left_idx]

    for idx in candidates:
        if idx < 0 or idx >= len(spec_rts):
            continue
        # Check if isolation window contains precursor m/z
        if spec_isolation_lower_mzs[idx] <= query_precursor_mz <= spec_isolation_upper_mzs[idx]:
            # Check RT tolerance (only if rt_tolerance >= 0)
            if rt_tolerance >= 0 and abs(spec_rts[idx] - query_rt) > rt_tolerance:
                continue

            return np.int32(idx)

    return np.int32(-1)


@numba.njit
def find_batch_DIA_spec_idxes_by_rt(
    spec_rts: np.ndarray,
    spec_isolation_lower_mzs: np.ndarray,
    spec_isolation_upper_mzs: np.ndarray,
    query_rts: np.ndarray,
    query_precursor_mzs: np.ndarray,
    rt_tolerance: float = -1.0,
) -> np.ndarray:
    """
    Batch version.

    Parameters
    ----------
    spec_rts : np.ndarray
        RT values of the spectra (must be sorted).
    spec_isolation_lower_mzs : np.ndarray
        Left m/z values of the isolation windows.
    spec_isolation_upper_mzs : np.ndarray
        Right m/z values of the isolation windows.
    query_rts : np.ndarray
        Query RT values.
    query_precursor_mzs : np.ndarray
        Query precursor m/z values.
    rt_tolerance : float, optional
        RT tolerance. If negative (default -1.0), no RT check is performed.

    Returns
    -------
    np.ndarray
        Array of spectrum indices. -1 indicates no valid spectrum found.
    """
    n_queries = query_rts.shape[0]
    result_indices = np.empty(n_queries, dtype=np.int32)

    for i in range(n_queries):
        result_indices[i] = find_single_DIA_spec_idx_by_rt(
            spec_rts,
            spec_isolation_lower_mzs,
            spec_isolation_upper_mzs,
            query_rts[i],
            query_precursor_mzs[i],
            rt_tolerance,
        )

    return result_indices

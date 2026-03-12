import logging
from typing import ClassVar, Optional

import numba
import numpy as np
import pandas as pd
from alphabase.peptide.fragment import create_fragment_mz_dataframe
from alpharaw.match.psm_match import load_ms_data
from tqdm import tqdm

from dia_aspire_rescore.config.ms2_match import MS2MatchConfig
from dia_aspire_rescore.constants.spectrum import PsmDfColsExt
from dia_aspire_rescore.extraction.utils import extract_xic
from dia_aspire_rescore.features.base import BaseFeatureGenerator

logger = logging.getLogger(__name__)


class XICFeatureGenerator(BaseFeatureGenerator):
    """
    Generate xic features.
    """

    DEFAULT_FRAG_TYPES: ClassVar[list[str]] = ["b_z1", "b_z2", "y_z1", "y_z2"]
    DEFAULT_PPM_TOLERANCE: ClassVar[float] = 20.0

    def __init__(
        self,
        ms_files: dict[str, str],
        ms_file_type: str = "mzml",
        ppm_tolerance: float = DEFAULT_PPM_TOLERANCE,
        ms2_match_config: Optional[MS2MatchConfig] = None,  # TODO:deprecated
    ):
        """
        Initialize XICFeatureGenerator.

        Parameters
        ----------
        ms_files : dict[str, str]
            Dict mapping raw_name to file path.
        ms_file_type : str, optional
            MS file type, by default "mzml".
        ppm_tolerance : float, optional
            m/z tolerance in ppm for XIC extraction, by default 20.0.
        ms2_match_config : Optional[MS2MatchConfig], optional
            MS2 matching config
        """
        self.ms_files = ms_files
        self.ms_file_type = ms_file_type
        self.ppm_tolerance = ppm_tolerance
        self.frag_types = self.DEFAULT_FRAG_TYPES
        self.ms2_match_config = ms2_match_config or MS2MatchConfig()

    @property
    def feature_names(self) -> list[str]:
        return [
            # cosine correlation features (unfiltered)
            "cos_corr_top3_median",
            "cos_corr_top3_mean",
            "cos_corr_top6_median",
            "cos_corr_top6_mean",
            "cos_corr_top12_median",
            "cos_corr_top12_mean",
            # Pearson correlation features (unfiltered)
            "pearson_corr_top3_median",
            "pearson_corr_top3_mean",
            "pearson_corr_top6_median",
            "pearson_corr_top6_mean",
            "pearson_corr_top12_median",
            "pearson_corr_top12_mean",
            # cosine correlation features (noise filtered)
            "cos_corr_filtered_top3_median",
            "cos_corr_filtered_top3_mean",
            "cos_corr_filtered_top6_median",
            "cos_corr_filtered_top6_mean",
            "cos_corr_filtered_top12_median",
            "cos_corr_filtered_top12_mean",
            # pearson correlation features (noise filtered)
            "pearson_corr_filtered_top3_median",
            "pearson_corr_filtered_top3_mean",
            "pearson_corr_filtered_top6_median",
            "pearson_corr_filtered_top6_mean",
            "pearson_corr_filtered_top12_median",
            "pearson_corr_filtered_top12_mean",
            # fragment quality metrics
            "n_filtered_fragments",
            "fragment_filter_ratio",
            # peak shape quality features
            "peak_symmetry",
            "peak_fwhm",
            "peak_jaggedness",
            "peak_modality",
            "apex_rt_ratio",
            "peak_base_width",
            "peak_tailing_factor",
            "total_xic_intensity",
            "max_xic_intensity",
        ]

    def generate(
        self,
        psm_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Generate XIC correlation features.

        Parameters
        ----------
        psm_df : pd.DataFrame
            PSM DataFrame with columns: sequence, charge, mods, mod_sites,
            rt_start, rt_stop, precursor_mz, raw_name.

        Returns
        -------
        pd.DataFrame
            psm_df with added XIC feature columns.
        """
        fragment_mz_df = create_fragment_mz_dataframe(psm_df, self.frag_types)

        for feat_name in self.feature_names:
            psm_df[feat_name] = 0.0

        for raw_name, file_path in self.ms_files.items():
            logger.info(f"Processing XIC features for {raw_name}")

            ms_data = load_ms_data(file_path, self.ms_file_type)
            spectrum_df = ms_data.spectrum_df
            peak_df = ms_data.peak_df

            mask = psm_df.raw_name == raw_name
            psm_indices = psm_df[mask].index.tolist()

            if len(psm_indices) == 0:
                continue

            feature_array = np.zeros((len(psm_indices), len(self.feature_names)), dtype=np.float32)

            for i, idx in enumerate(tqdm(psm_indices, desc=f"Processing {raw_name}")):
                psm = psm_df.loc[idx]
                features = self._calculate_xic_features(
                    psm,
                    fragment_mz_df,
                    spectrum_df,
                    peak_df,
                )
                feature_array[i, :] = features

            psm_df.loc[psm_indices, self.feature_names] = feature_array

        return psm_df

    def _calculate_xic_features(
        self,
        psm: pd.Series,
        fragment_mz_df: pd.DataFrame,
        spectrum_df: pd.DataFrame,
        peak_df: pd.DataFrame,
    ) -> np.ndarray:
        """
        Calculate XIC co-elution features for a single PSM.

        Returns
        -------
        np.ndarray
            Array of 35 XIC features (dtype=float32):
            - 12 unfiltered MS2 correlation features
            - 12 filtered MS2 correlation features
            - 2 fragment quality metrics
            - 9 peak shape features
        """
        frag_start = int(psm[PsmDfColsExt.FRAG_START_IDX])
        frag_stop = int(psm[PsmDfColsExt.FRAG_STOP_IDX])
        frag_mzs = fragment_mz_df.iloc[frag_start:frag_stop]

        query_mzs = frag_mzs.values.flatten()
        query_mzs = query_mzs[query_mzs > 0]

        n_valid_fragments = len(query_mzs)

        if len(query_mzs) < 2:
            return np.zeros(len(self.feature_names), dtype=np.float32)

        try:
            rt_values, intensity_matrix = extract_xic(
                spectrum_df=spectrum_df,
                peak_df=peak_df,
                query_mzs=query_mzs,
                rt_start=psm[PsmDfColsExt.RT_START],
                rt_stop=psm[PsmDfColsExt.RT_STOP],
                precursor_mz=psm[PsmDfColsExt.PRECURSOR_MZ],
                ppm_tolerance=self.ppm_tolerance,
                ms_level=2,
                match_closest=True,
            )
        except Exception:
            return np.zeros(len(self.feature_names), dtype=np.float32)

        if len(rt_values) < 3 or intensity_matrix.shape[1] < 3:
            return np.zeros(len(self.feature_names), dtype=np.float32)

        cos_corr_matrix = _calc_cos_corr_matrix(intensity_matrix)
        pearson_corr_matrix = _calc_pearson_corr_matrix(intensity_matrix)

        cos_features = _extract_top_k_corr_features(cos_corr_matrix)
        pearson_features = _extract_top_k_corr_features(pearson_corr_matrix)

        # noise filtering based on intensity threshold
        filtered_intensity_matrix, valid_indices = _filter_fragments_by_intensity(
            intensity_matrix, min_intensity_ratio=0.1, top_n=12
        )

        n_filtered_fragments = len(valid_indices)
        fragment_filter_ratio = n_filtered_fragments / n_valid_fragments if n_valid_fragments > 0 else 0.0

        if n_filtered_fragments >= 2:
            cos_corr_matrix_filtered = _calc_cos_corr_matrix(filtered_intensity_matrix)
            pearson_corr_matrix_filtered = _calc_pearson_corr_matrix(filtered_intensity_matrix)

            cos_features_filtered = _extract_top_k_corr_features(cos_corr_matrix_filtered)
            pearson_features_filtered = _extract_top_k_corr_features(pearson_corr_matrix_filtered)
        else:
            cos_features_filtered = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            pearson_features_filtered = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        # peak shape features
        peak_shape_features = _calculate_peak_shape_features(
            filtered_intensity_matrix if n_filtered_fragments >= 2 else intensity_matrix
        )

        all_features = (
            cos_features
            + pearson_features
            + cos_features_filtered
            + pearson_features_filtered
            + (float(n_filtered_fragments), fragment_filter_ratio)
            + peak_shape_features
        )

        return np.array(all_features, dtype=np.float32)


@numba.njit
def _calc_cos_corr_matrix(intensity_matrix: np.ndarray) -> np.ndarray:
    """
    Calculate all cosine correlations between fragments.

    Parameters
    ----------
    intensity_matrix : np.ndarray
        Intensity matrix, shape (n_ions, n_specs).

    Returns
    -------
    np.ndarray
        Cosine correlation matrix, shape (n_ions, n_ions).
    """
    n_ions = intensity_matrix.shape[0]
    corr_matrix = np.zeros((n_ions, n_ions), dtype=np.float32)

    for i in range(n_ions):
        corr_matrix[i, i] = 1.0  # diagonal is always 1
        for j in range(i + 1, n_ions):
            corr = _calc_cos_corr(intensity_matrix[i, :], intensity_matrix[j, :])
            corr_matrix[i, j] = corr
            corr_matrix[j, i] = corr

    return corr_matrix


@numba.njit
def _calc_cos_corr(
    array_1: np.ndarray,
    array_2: np.ndarray,
) -> float:
    """
    Calculate cosine correlation between two arrays.
    """
    norm_1 = np.linalg.norm(array_1)
    norm_2 = np.linalg.norm(array_2)
    if norm_1 == 0 or norm_2 == 0:
        return 0.0
    return np.dot(array_1, array_2) / (norm_1 * norm_2)


@numba.njit
def _calc_pearson_corr_matrix(intensity_matrix: np.ndarray) -> np.ndarray:
    """
    Calculate all Pearson correlations between fragments.

    Parameters
    ----------
    intensity_matrix : np.ndarray
        Intensity matrix, shape (n_ions, n_specs).

    Returns
    -------
    np.ndarray
        Pearson correlation matrix, shape (n_ions, n_ions).
    """
    n_ions = intensity_matrix.shape[0]
    corr_matrix = np.zeros((n_ions, n_ions), dtype=np.float32)

    for i in range(n_ions):
        corr_matrix[i, i] = 1.0  # diagonal is always 1
        for j in range(i + 1, n_ions):
            corr = _calc_pearson_corr(intensity_matrix[i, :], intensity_matrix[j, :])
            corr_matrix[i, j] = corr
            corr_matrix[j, i] = corr

    return corr_matrix


@numba.njit
def _calc_pearson_corr(
    array_1: np.ndarray,
    array_2: np.ndarray,
) -> float:
    """
    Calculate Pearson correlation between two arrays.
    """
    if np.std(array_1) == 0 or np.std(array_2) == 0:
        return 0.0

    mean_1 = np.mean(array_1)
    mean_2 = np.mean(array_2)
    centered_1 = array_1 - mean_1
    centered_2 = array_2 - mean_2
    numerator = np.sum(centered_1 * centered_2)
    denominator = np.sqrt(np.sum(centered_1**2) * np.sum(centered_2**2))

    if denominator == 0:
        return 0.0

    return numerator / denominator


def _filter_fragments_by_intensity(
    intensity_matrix: np.ndarray,
    min_intensity_ratio: float = 0.1,
    top_n: int = 12,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Filter fragments by intensity to remove noise.

    Strategy:
    1. Calculate max intensity for each fragment across all scans
    2. Filter fragments with max_intensity > min_intensity_ratio * global_max
    3. Keep only top_n fragments by max intensity

    Parameters
    ----------
    intensity_matrix : np.ndarray
        Intensity matrix, shape (n_ions, n_specs).
    min_intensity_ratio : float, default=0.1
        Minimum intensity ratio relative to global maximum.
    top_n : int, default=12
        Maximum number of fragments to keep.

    Returns
    -------
    filtered_matrix : np.ndarray
        Filtered intensity matrix, shape (n_filtered, n_specs).
    valid_indices : np.ndarray
        Indices of fragments that passed filtering.
    """
    n_ions, n_specs = intensity_matrix.shape

    if n_ions == 0:
        return np.zeros((0, n_specs), dtype=np.float32), np.array([], dtype=np.int64)

    max_intensities = np.zeros(n_ions, dtype=np.float32)
    max_intensities = np.max(intensity_matrix, axis=1)

    global_max = np.max(max_intensities)

    if global_max == 0:
        return np.zeros((0, n_specs), dtype=np.float32), np.array([], dtype=np.int64)

    intensity_threshold = global_max * min_intensity_ratio
    valid_mask = max_intensities >= intensity_threshold

    valid_indices = np.where(valid_mask)[0]  # (n_valid_fragments,)

    if len(valid_indices) == 0:
        return np.zeros((0, n_specs), dtype=np.float32), np.array([], dtype=np.int64)

    valid_max_intensities = max_intensities[valid_indices]
    sorted_idx = np.argsort(valid_max_intensities)[::-1]
    n_to_keep = min(top_n, len(sorted_idx))
    top_indices = sorted_idx[:n_to_keep]  # (n_to_keep)
    final_indices = valid_indices[top_indices]
    filtered_matrix = intensity_matrix[final_indices, :]

    return filtered_matrix, final_indices


@numba.njit
def _extract_top_k_corr_features(
    corr_matrix: np.ndarray,
) -> tuple[float, float, float, float, float, float]:
    """
    Extract top-k correlation features from correlation matrix.

    Takes the upper triangle (excluding diagonal) of the correlation matrix,
    selects the top 3, 6, and 12 values, and calculates median and mean.

    Parameters
    ----------
    corr_matrix : np.ndarray
        Correlation matrix, shape (n_ions, n_ions).

    Returns
    -------
    tuple[float, float, float, float, float, float]
        (top3_median, top3_mean, top6_median, top6_mean,
         top12_median, top12_mean)
    """
    n_ions = corr_matrix.shape[0]

    # excluding diagonal
    upper_tri_values = []
    for i in range(n_ions):
        for j in range(i + 1, n_ions):
            upper_tri_values.append(corr_matrix[i, j])

    if len(upper_tri_values) == 0:
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    upper_tri_array = np.array(upper_tri_values, dtype=np.float32)

    sorted_corr = np.sort(upper_tri_array)[::-1]

    top3 = sorted_corr[: min(3, len(sorted_corr))]
    top3_median = np.median(top3)
    top3_mean = np.mean(top3)

    top6 = sorted_corr[: min(6, len(sorted_corr))]
    top6_median = np.median(top6)
    top6_mean = np.mean(top6)

    top12 = sorted_corr[: min(12, len(sorted_corr))]
    top12_median = np.median(top12)
    top12_mean = np.mean(top12)

    return (
        float(top3_median),
        float(top3_mean),
        float(top6_median),
        float(top6_mean),
        float(top12_median),
        float(top12_mean),
    )


def _calculate_peak_shape_features(
    intensity_matrix: np.ndarray,
) -> tuple[float, float, float, float, float, float, float, float, float]:
    """
    Calculate peak shape quality features from intensity matrix.

    Parameters
    ----------
    intensity_matrix : np.ndarray
        Intensity matrix, shape (n_ions, n_specs).

    Returns
    -------
    tuple of 9 floats
        (peak_symmetry, peak_fwhm, peak_jaggedness, peak_modality, apex_rt_ratio,
         peak_base_width, peak_tailing_factor, total_xic_intensity, max_xic_intensity)
    """
    if intensity_matrix.shape[0] == 0 or intensity_matrix.shape[1] == 0:
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    summed_xic = np.sum(intensity_matrix, axis=0)

    if len(summed_xic) == 0 or np.max(summed_xic) == 0:
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    peak_symmetry = float(_calc_peak_symmetry(summed_xic))
    peak_fwhm = float(_calc_peak_fwhm(summed_xic))
    peak_jaggedness = float(_calc_peak_jaggedness(summed_xic))
    peak_modality = float(_calc_peak_modality(summed_xic))
    apex_rt_ratio = float(_calc_apex_rt_ratio(summed_xic))
    peak_base_width = float(_calc_peak_base_width(summed_xic))
    peak_tailing_factor = float(_calc_peak_tailing_factor(summed_xic))
    total_xic_intensity = float(np.log(np.sum(intensity_matrix) + 1))
    max_xic_intensity = float(np.log(np.max(intensity_matrix) + 1))

    return (
        peak_symmetry,
        peak_fwhm,
        peak_jaggedness,
        peak_modality,
        apex_rt_ratio,
        peak_base_width,
        peak_tailing_factor,
        total_xic_intensity,
        max_xic_intensity,
    )


@numba.njit
def _calc_peak_symmetry(xic: np.ndarray) -> float:
    """Calculate peak symmetry as ratio of left/right half areas."""
    if len(xic) < 3:
        return 0.0

    apex_idx = np.argmax(xic)
    if apex_idx == 0 or apex_idx == len(xic) - 1:
        return 0.0

    left_area = np.sum(xic[:apex_idx])
    right_area = np.sum(xic[apex_idx + 1 :])

    if left_area == 0 and right_area == 0:
        return 1.0

    return min(left_area, right_area) / max(left_area, right_area) if max(left_area, right_area) > 0 else 0.0


@numba.njit
def _calc_peak_fwhm(xic: np.ndarray) -> float:
    """Calculate full width at half maximum."""
    if len(xic) < 3:
        return 0.0

    max_intensity = np.max(xic)
    if max_intensity == 0:
        return 0.0

    half_max = max_intensity / 2.0

    above_half_max = xic >= half_max
    n_above = np.sum(above_half_max)

    return float(n_above)


@numba.njit
def _calc_peak_jaggedness(xic: np.ndarray) -> float:
    """Calculate jaggedness as smoothness measure (1 - sign_changes / n_points)."""
    if len(xic) < 3:
        return 0.0

    diffs = np.diff(xic)

    sign_changes = 0
    for i in range(len(diffs) - 1):
        if diffs[i] * diffs[i + 1] < 0:  # indicate sign change
            sign_changes += 1

    max_changes = len(xic) - 2
    if max_changes > 0:
        return 1.0 - (sign_changes / max_changes)
    return 1.0


@numba.njit
def _calc_peak_modality(xic: np.ndarray) -> float:
    """Calculate number of local maxima."""
    if len(xic) < 3:
        return 1.0

    local_maxima = 0
    for i in range(1, len(xic) - 1):
        if xic[i] > xic[i - 1] and xic[i] > xic[i + 1]:
            local_maxima += 1

    return float(max(1, local_maxima))


@numba.njit
def _calc_apex_rt_ratio(xic: np.ndarray) -> float:
    """Calculate position of apex relative to peak boundaries."""
    if len(xic) < 2:
        return 0.5

    apex_idx = np.argmax(xic)
    return float(apex_idx) / float(len(xic) - 1)


@numba.njit
def _calc_peak_base_width(xic: np.ndarray) -> float:
    """Calculate width at 10% of max intensity."""
    if len(xic) < 3:
        return 0.0

    max_intensity = np.max(xic)
    if max_intensity == 0:
        return 0.0

    threshold = max_intensity * 0.1

    above_threshold = xic >= threshold
    n_above = np.sum(above_threshold)

    return float(n_above)


@numba.njit
def _calc_peak_tailing_factor(xic: np.ndarray) -> float:
    """Calculate tailing factor as right_width / left_width at 10% height."""
    if len(xic) < 3:
        return 1.0

    apex_idx = np.argmax(xic)
    max_intensity = xic[apex_idx]

    if max_intensity == 0:
        return 1.0

    threshold = max_intensity * 0.1

    left_width = 0
    for i in range(apex_idx, -1, -1):
        if xic[i] >= threshold:
            left_width += 1
        else:
            break

    right_width = 0
    for i in range(apex_idx, len(xic)):
        if xic[i] >= threshold:
            right_width += 1
        else:
            break

    if left_width == 0:
        return 1.0

    return float(right_width) / float(left_width)

from __future__ import annotations

import logging

import numba
import numpy as np
import pandas as pd
from alphabase.peptide.fragment import get_charged_frag_types
from alpharaw.dia.normal_dia import NormalDIAGrouper
from alpharaw.match.match_utils import (
    match_closest_peaks,
    match_highest_peaks,
)
from alpharaw.match.psm_match import PepSpecMatch, load_ms_data
from alpharaw.ms_data_base import (
    PEAK_INTENSITY_DTYPE,
    PEAK_MZ_DTYPE,
)
from alpharaw.utils.ms_path_utils import parse_ms_files_to_dict

from dia_aspire.constants.spectrum import PsmDfColsExt
from dia_aspire.psm.spec_finder import find_dia_spec_idxes_same_window

logger = logging.getLogger(__name__)


class DIAPeptideSpectrumMatcher(PepSpecMatch):
    """
    Peptide-spectrum matcher for DIA data.
    Using retention time to find spectra in the same DIA window.

    Adapted from alpharaw: https://github.com/MannLabs/alpharaw
    """

    def __init__(
        self,
        charged_frag_types: list | None = None,
        match_closest: bool = True,
        use_ppm: bool = True,
        tol_value: float = 20.0,
        n_neighbors: int = 0,
    ):
        self.charged_frag_types = (
            get_charged_frag_types(["b", "y"], 2) if charged_frag_types is None else charged_frag_types
        )
        self.match_closest = match_closest
        self.use_ppm = use_ppm
        self.tolerance = tol_value
        self.max_spec_per_query = 1 + 2 * n_neighbors

    def _add_missing_columns_to_psm_df(self, psm_df: pd.DataFrame, raw_data=None):
        # DIA results do not have spec_idx/scan_num in psm_df, nothing to merge
        return psm_df

    def _prepare_matching_dfs(self):
        fragment_mz_df = self.get_fragment_mz_df()
        fragment_mz_df = pd.concat([fragment_mz_df] * self.max_spec_per_query, ignore_index=True)
        if self.use_ppm:
            self.all_frag_mz_tols = fragment_mz_df.values * self.tolerance * 1e-6
        else:
            self.all_frag_mz_tols = np.full_like(fragment_mz_df.values, self.tolerance)

        psm_df_list = []
        len_frags = len(fragment_mz_df) // self.max_spec_per_query
        for i in range(self.max_spec_per_query):
            psm_df = self.psm_df.copy()
            psm_df[PsmDfColsExt.FRAG_START_IDX] = psm_df.frag_start_idx + i * len_frags
            psm_df[PsmDfColsExt.FRAG_STOP_IDX] = psm_df.frag_stop_idx + i * len_frags
            psm_df_list.append(psm_df)
        self.psm_df = pd.concat(psm_df_list, ignore_index=True)

        matched_intensity_df = pd.DataFrame(
            np.zeros_like(fragment_mz_df.values, dtype=PEAK_INTENSITY_DTYPE),
            columns=fragment_mz_df.columns,
        )

        matched_mz_err_df = pd.DataFrame(
            np.zeros_like(fragment_mz_df.values, dtype=PEAK_MZ_DTYPE),
            columns=fragment_mz_df.columns,
        )

        return (fragment_mz_df, matched_intensity_df, matched_mz_err_df)

    def _match_ms2_one_raw_numba(self, raw_name: str, psm_df_one_raw: pd.DataFrame) -> pd.DataFrame:
        """
        Internal method to extract peak information with numba as backend.

        Parameters
        ----------
        raw_name : str
            The raw name of the raw file. `psm_df_one_raw` dataframe should also
            contain the same raw name in `raw_name` column.
        psm_df_one_raw : pd.DataFrame
            The dataframe for PSMs.

        Returns
        -------
        pd.DataFrame
            `psm_df_one_raw`
        """
        psm_df_one_raw = psm_df_one_raw.reset_index(drop=True)
        psm_df_one_raw[PsmDfColsExt.SPEC_IDX] = -1

        if raw_name in self._ms_file_dict:
            raw_data = load_ms_data(
                self._ms_file_dict[raw_name],
                self._ms_file_type,
                process_count=self.ms_loader_thread_num,
            )
            if raw_data is None:
                return

            psm_origin_len = len(psm_df_one_raw) // self.max_spec_per_query

            grouper = NormalDIAGrouper(raw_data)

            psm_groups = grouper.assign_dia_groups(psm_df_one_raw.precursor_mz.values[:psm_origin_len])

            all_spec_idxes = np.full(len(psm_df_one_raw), -1, dtype=np.int32)

            for dia_group, group_df in grouper.dia_group_dfs:
                logger.debug(f"Processing DIA group: {dia_group}")
                psm_idxes = psm_groups[dia_group]
                logger.debug(f"PSM idxes in this DIA group: {psm_idxes}")
                if len(psm_idxes) == 0:
                    continue
                psm_idxes = np.array(psm_idxes, dtype=np.int32)

                # spec_idxes shape: (len(psm_idxes), max_spec_per_query)
                spec_idxes = find_dia_spec_idxes_same_window(
                    group_df.rt.values,
                    psm_df_one_raw.rt.values[psm_idxes],
                    max_spec_per_query=self.max_spec_per_query,
                )

                logger.debug(f"RTs in psm_df: {psm_df_one_raw.rt.values[psm_idxes]}")
                logger.debug(f"# Raw specs in this DIA group: {len(group_df.rt.values)}")
                logger.debug(f"Relative spec_idxes return by find_dia_spec_idxes_same_window: {spec_idxes}")

                # extract the spec_idx in group_df by the results from find_dia_spec_idxes_same_window
                absolute_spec_idxes = group_df.spec_idx.values[spec_idxes]
                logger.debug(f"True spec_idxes in raw file: {absolute_spec_idxes}")

                spec_rts = raw_data.spectrum_df.rt.values[absolute_spec_idxes]
                logger.debug(f"Corresponding RTs in raw file: {spec_rts}")
                logger.debug(f"Length of psm_one_raw: {len(psm_df_one_raw)}")

                # TODO: assign spec_idx back to psm_df_one_raw
                for i in range(spec_idxes.shape[-1]):
                    logger.debug(
                        f"Assigning spec_idxes for neighbor {i}, spec_idxes: {spec_idxes[:, i]}, absolute_spec_idxes: {absolute_spec_idxes[:, i]}"
                    )
                    all_spec_idxes[psm_idxes + psm_origin_len * i] = spec_idxes[:, i]
                    # assign spec_idx to psm_df_one_raw
                    psm_df_one_raw.loc[psm_idxes + psm_origin_len * i, PsmDfColsExt.SPEC_IDX] = absolute_spec_idxes[
                        :, i
                    ]

                # Collect indices of all PSMs in current group (including all neighbors)
                current_group_psm_indices_list = []
                for psm_idx in psm_idxes:
                    for neighbor_idx in range(self.max_spec_per_query):
                        idx_in_all = psm_idx + psm_origin_len * neighbor_idx
                        current_group_psm_indices_list.append(idx_in_all)

                current_group_psm_indices = np.array(current_group_psm_indices_list, dtype=np.int32)

                # Extract data for current group only
                current_spec_idxes = all_spec_idxes[current_group_psm_indices]
                current_frag_start = psm_df_one_raw.frag_start_idx.values[current_group_psm_indices]
                current_frag_stop = psm_df_one_raw.frag_stop_idx.values[current_group_psm_indices]

                logger.debug(f"Current group PSM indices: {current_group_psm_indices}")
                logger.debug(f"Current group spec_idxes: {current_spec_idxes}")

                # Only match PSMs in current group with current group's peaks
                match_one_raw_with_numba(
                    current_spec_idxes,
                    current_frag_start,
                    current_frag_stop,
                    self.fragment_mz_df.values,
                    self.all_frag_mz_tols,
                    raw_data.peak_df.mz.values,
                    raw_data.peak_df.intensity.values,
                    group_df.peak_start_idx.values,
                    group_df.peak_stop_idx.values,
                    self.matched_intensity_df.values,
                    self.matched_mz_err_df.values,
                    self.match_closest,
                )
        else:
            print(f"`{raw_name}` is not found in ms_file_dict.")
            return
        return psm_df_one_raw

    def match_ms2_multi_raw(
        self,
        psm_df: pd.DataFrame,
        ms_files: dict | list,
        ms_file_type: str = "alpharaw_hdf",
        process_num: int = 8,
    ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Match peaks for the given `psm_df` against the corresponding MS spectrum files.

        Parameters
        ----------
        psm_df : pd.DataFrame
            Peptide-spectrum matches in alphabase dataframe format.
        ms_files : dict | list
            The absolute or relative paths of MS files.
            If `dict`, the format should be `{'raw_name1': 'path/to/file1', ...}`.
            If `list`, it should be a list of file paths.
        ms_file_type : str, optional
            MS file type that is already registered in
            :obj:`alpharaw.ms_data_base.ms_reader_provider`.
            By default "alpharaw_hdf".
        process_num : int, optional
            Match peaks by using multiprocessing, by default 8

        Returns
        -------
        tuple
            psm_df : pd.DataFrame
                PSM DataFrame with added columns. PSMs are grouped by 'raw_name'.
                Fragment indices via `frag_start_idx`, `frag_stop_idx`.

            fragment_mz_df : pd.DataFrame
                Fragment m/z in alphabase wide format.

            matched_intensity_df : pd.DataFrame
                Matched fragment intensities. Rows indexed by `frag_start_idx:frag_stop_idx`.

            matched_mz_err_df : pd.DataFrame
                Matched mass errors (ppm or Da). np.inf if a fragment is not matched.

        Notes
        -----
        - Output PSMs are grouped by `raw_name` for efficient processing.
        - If `n_neighbors > 0`, each PSM is replicated `(1 + 2*n_neighbors)` times.
        - Fragment DataFrames align with psm_df via `frag_start_idx`, `frag_stop_idx`.
        """
        if isinstance(ms_files, list):
            ms_files = parse_ms_files_to_dict(ms_files)
        psm_df = psm_df[psm_df.raw_name.isin(ms_files)].reset_index(drop=True)
        super().match_ms2_multi_raw(psm_df, ms_files, ms_file_type, process_num)

        return (
            self.psm_df,
            self.fragment_mz_df,
            self.matched_intensity_df,
            self.matched_mz_err_df,
        )

    def match_ms2_one_raw(
        self,
        psm_df_one_raw: pd.DataFrame,
        ms_file: str,
        ms_file_type: str = "alpharaw_hdf",
    ) -> pd.DataFrame:
        """
        Match psm_df_one_raw against ms2_file

        Parameters
        ----------
        psm_df_one_raw : pd.DataFrame
            PSM DataFrame that contains only one raw file.
        ms_file : str
            The path to the MS2 file.
        ms_file_type : str, optional
            The type of the MS2 file.

        Returns
        -------
        tuple
            psm_df_one_raw : pd.DataFrame
                PSM DataFrame that contains only one raw file.
            fragment_mz_df : pd.DataFrame
                Fragment m/z in alphabase wide format.
            matched_intensity_df : pd.DataFrame
                Matched fragment intensities. Rows indexed by `frag_start_idx:frag_stop_idx`.
            matched_mz_err_df : pd.DataFrame
                Matched mass errors (ppm or Da). np.inf if a fragment is not matched.
        """
        if len(psm_df_one_raw.raw_name.unique()) > 1:
            raise ValueError("psm_df_one_raw should contain only one raw file.")

        ms_file_dict = {psm_df_one_raw.raw_name.values[0]: ms_file}
        return self.match_ms2_multi_raw(psm_df_one_raw, ms_file_dict, ms_file_type)


@numba.njit
def match_one_raw_with_numba(
    spec_idxes: np.ndarray,
    frag_start_idxes: np.ndarray,
    frag_stop_idxes: np.ndarray,
    all_frag_mzs: np.ndarray,
    all_frag_mz_tols: np.ndarray,
    all_spec_mzs: np.ndarray,
    all_spec_intensities: np.ndarray,
    peak_start_idxes: np.ndarray,
    peak_stop_idxes: np.ndarray,
    matched_intensities: np.ndarray,
    matched_mz_errs: np.ndarray,
    match_closest: bool = True,
) -> None:
    """
    Internel function to match fragment mz values to spectrum mz values.
    Matched_mz_errs[i] = np.inf if no peaks are matched.

    Results will saved in place of matched_intensities
    and matched_mz_errs.
    """
    for spec_idx, frag_start, frag_stop in zip(spec_idxes, frag_start_idxes, frag_stop_idxes):
        if spec_idx == -1:
            continue
        peak_start = peak_start_idxes[spec_idx]
        peak_stop = peak_stop_idxes[spec_idx]
        if peak_stop == peak_start:
            continue
        spec_mzs = all_spec_mzs[peak_start:peak_stop]
        spec_intens = all_spec_intensities[peak_start:peak_stop]

        frag_mzs = all_frag_mzs[frag_start:frag_stop, :].copy()
        frag_mz_tols = all_frag_mz_tols[frag_start:frag_stop, :].copy()

        if match_closest:
            matched_idxes = match_closest_peaks(spec_mzs, spec_intens, frag_mzs, frag_mz_tols).reshape(-1)
        else:
            matched_idxes = match_highest_peaks(spec_mzs, spec_intens, frag_mzs, frag_mz_tols).reshape(-1)

        matched_intens = spec_intens[matched_idxes]
        matched_intens[matched_idxes == -1] = 0

        matched_mass_errs = np.abs(spec_mzs[matched_idxes.reshape(-1)] - frag_mzs.reshape(-1))
        matched_mass_errs[matched_idxes == -1] = np.inf

        matched_intensities[frag_start:frag_stop, :] = matched_intens.reshape(frag_mzs.shape)

        matched_mz_errs[frag_start:frag_stop, :] = matched_mass_errs.reshape(frag_mzs.shape)

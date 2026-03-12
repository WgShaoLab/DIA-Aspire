from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from alphabase.constants.atom import MASS_ISOTOPE
from alpharaw.match.psm_match import load_ms_data
from tqdm import tqdm

from dia_aspire.extraction.utils import extract_xic

logger = logging.getLogger(__name__)


class XICExtractor:
    """
    Extract XIC traces for MS1 isotopes and MS2 transitions from DIA raw files based on the spectral library.

    Parameters
    ----------
    precursor_df : pd.DataFrame
        Spectral library precursor DataFrame from SpecLibReader.
    transition_df : pd.DataFrame
        Spectral library transition DataFrame from SpecLibReader.
    n_isotopes : int, default=3
        Number of MS1 isotopes to extract (M, M+1, M+2, ...).
    rt_extension : float, default=30.0
        RT window extension in seconds (applied to both sides).
    ppm_tolerance : float, default=15.0
        m/z tolerance in ppm for XIC extraction.

    Examples
    --------
    >>> from dia_aspire_rescore.io import read_speclib, read_diann2
    >>> precursor_df, transition_df = read_speclib("library.tsv")
    >>> psm_df = read_diann2("report.parquet")
    >>> extractor = XICExtractor(precursor_df, transition_df)
    >>> extractor.extract_all(psm_df, ms_files, output_dir="./output")
    """

    def __init__(
        self,
        precursor_df: pd.DataFrame,
        transition_df: pd.DataFrame,
        n_isotopes: int = 3,
        rt_extension: float = 30.0,
        ppm_tolerance: float = 15.0,
    ):
        self.precursor_df = precursor_df
        self.transition_df = transition_df
        self.n_isotopes = n_isotopes
        self.rt_extension = rt_extension
        self.ppm_tolerance = ppm_tolerance

    def match_psm_to_speclib(self, psm_df: pd.DataFrame) -> pd.DataFrame:
        """
        Match PSMs to spectral library to get transition information.

        Parameters
        ----------
        psm_df : pd.DataFrame
            PSM DataFrame from DIA-NN report.

        Returns
        -------
        pd.DataFrame
            Merged DataFrame with PSM info + speclib info (transition indices, precursor_id).
        """
        merge_keys = ["sequence", "mods", "mod_sites", "charge"]
        speclib_cols = [*merge_keys, "precursor_id", "transition_start_idx", "transition_stop_idx"]
        psm_cols = [*merge_keys, "rt", "rt_start", "rt_stop", "precursor_mz", "raw_name"]

        matched_df = psm_df[psm_cols].merge(self.precursor_df[speclib_cols], on=merge_keys, how="inner")

        logger.info(f"Matched {len(matched_df)}/{len(psm_df)} PSMs to spectral library")
        return matched_df

    def extract_one_raw(
        self,
        matched_df: pd.DataFrame,
        ms_file_path: str,
        ms_file_type: str = "mzml",
    ) -> pd.DataFrame:
        """
        Extract XIC traces for all matched precursors from one raw file.

        Parameters
        ----------
        matched_df : pd.DataFrame
            Matched DataFrame (PSM + speclib info) for this raw file.
        ms_file_path : str
            Path to the MS raw file.
        ms_file_type : str, default="mzml"
            MS file type ("mzml" or "alpharaw_hdf").

        Returns
        -------
        pd.DataFrame
            XIC DataFrame with columns: precursor_id, level, mz, type, charge, rt, intensity.
        """
        if len(matched_df) == 0:
            return pd.DataFrame(columns=["precursor_id", "level", "mz", "type", "charge", "rt", "intensity"])

        ms_data = load_ms_data(ms_file_path, ms_file_type)
        if ms_data is None:
            logger.error(f"Failed to load MS data from {ms_file_path}")
            return pd.DataFrame(columns=["pr", "feature", "info", "rt", "value"])

        spectrum_df = ms_data.spectrum_df
        peak_df = ms_data.peak_df

        xic_records: list[dict] = []

        for _, row in tqdm(matched_df.iterrows(), total=len(matched_df), desc="Extracting XICs"):
            precursor_id = row["precursor_id"]
            precursor_mz = float(row["precursor_mz"])
            charge = int(row["charge"])
            rt_start = float(row["rt_start"])
            rt_stop = float(row["rt_stop"])
            rt_extension_min = self.rt_extension / 60.0  # to minutes
            rt_start_extract = rt_start - rt_extension_min
            rt_stop_extract = rt_stop + rt_extension_min

            ms1_records = self._extract_ms1_xic(
                precursor_id=precursor_id,
                precursor_mz=precursor_mz,
                charge=charge,
                rt_start=rt_start_extract,
                rt_stop=rt_stop_extract,
                spectrum_df=spectrum_df,
                peak_df=peak_df,
            )
            xic_records.extend(ms1_records)

            # transitions in speclib
            transition_start = int(row["transition_start_idx"])
            transition_stop = int(row["transition_stop_idx"])
            ms2_records = self._extract_ms2_xic(
                precursor_id=precursor_id,
                precursor_mz=precursor_mz,
                transition_start=transition_start,
                transition_stop=transition_stop,
                rt_start=rt_start_extract,
                rt_stop=rt_stop_extract,
                spectrum_df=spectrum_df,
                peak_df=peak_df,
            )
            xic_records.extend(ms2_records)

        if not xic_records:
            return pd.DataFrame(columns=["precursor_id", "level", "mz", "type", "charge", "rt", "intensity"])

        xic_df = pd.DataFrame(xic_records)
        xic_df["precursor_id"] = xic_df["precursor_id"].astype("object")
        xic_df["level"] = xic_df["level"].astype("object")
        xic_df["mz"] = xic_df["mz"].astype(np.float32)
        xic_df["type"] = xic_df["type"].astype("object")
        xic_df["charge"] = xic_df["charge"].astype(np.int8)
        xic_df["rt"] = xic_df["rt"].astype(np.float32)
        xic_df["intensity"] = xic_df["intensity"].astype(np.float32)
        return xic_df

    def _extract_ms1_xic(
        self,
        precursor_id: str,
        precursor_mz: float,
        charge: int,
        rt_start: float,
        rt_stop: float,
        spectrum_df: pd.DataFrame,
        peak_df: pd.DataFrame,
    ) -> list[dict]:
        records: list[dict] = []
        isotope_mzs = np.array([precursor_mz + i * MASS_ISOTOPE / charge for i in range(self.n_isotopes)])
        # TODO: extract all isotopes at once
        for isotope_idx, isotope_mz in enumerate(isotope_mzs):
            rt_values, intensities = extract_xic(
                spectrum_df=spectrum_df,
                peak_df=peak_df,
                query_mzs=np.array([isotope_mz]),
                rt_start=rt_start,
                rt_stop=rt_stop,
                precursor_mz=None,
                ppm_tolerance=self.ppm_tolerance,
                ms_level=1,
            )
            # intensities (1, timepoints)
            if len(rt_values) == 0:
                continue
            intensity_trace = intensities[0]  # (timepoints,)
            isotope_type = "M" if isotope_idx == 0 else f"M+{isotope_idx}"
            for rt, intensity in zip(rt_values, intensity_trace):
                records.append({
                    "precursor_id": precursor_id,
                    "level": "ms1",
                    "mz": float(isotope_mz),
                    "type": isotope_type,
                    "charge": charge,
                    "rt": rt,
                    "intensity": intensity,
                })

        return records

    def _extract_ms2_xic(
        self,
        precursor_id: str,
        precursor_mz: float,
        transition_start: int,
        transition_stop: int,
        rt_start: float,
        rt_stop: float,
        spectrum_df: pd.DataFrame,
        peak_df: pd.DataFrame,
    ) -> list[dict]:
        records: list[dict] = []
        if transition_start >= transition_stop:
            return records
        transitions = self.transition_df.iloc[transition_start:transition_stop]

        for _, trans_row in transitions.iterrows():
            trans_mz = float(trans_row["mz"])
            trans_type = str(trans_row["type"])
            trans_charge = int(trans_row["charge"])

            rt_values, intensities = extract_xic(
                spectrum_df=spectrum_df,
                peak_df=peak_df,
                query_mzs=np.array([trans_mz]),
                rt_start=rt_start,
                rt_stop=rt_stop,
                precursor_mz=precursor_mz,
                ppm_tolerance=self.ppm_tolerance,
                ms_level=2,
            )

            if len(rt_values) == 0:
                continue

            intensity_trace = intensities[0]
            for rt, intensity in zip(rt_values, intensity_trace):
                records.append({
                    "precursor_id": precursor_id,
                    "level": "ms2",
                    "mz": trans_mz,
                    "type": trans_type,
                    "charge": trans_charge,
                    "rt": rt,
                    "intensity": intensity,
                })
        return records

    def extract_all(
        self,
        psm_df: pd.DataFrame,
        ms_files: dict[str, str],
        ms_file_type: str = "mzml",
        output_dir: str = "./xic",
        output_format: str = "parquet",
    ) -> dict[str, Path]:
        """
        Extract XIC traces for all raw files.

        Parameters
        ----------
        psm_df : pd.DataFrame
            PSM DataFrame from DIA-NN report.
        ms_files : dict[str, str]
            Dict mapping raw_name to file path.
        ms_file_type : str, default="mzml"
            MS file type.
        output_dir : str, default="./xic"
            Output directory for output xic files.
        output_format : str, default="parquet"
            Output format: "parquet" or "csv".

        Returns
        -------
        dict[str, Path]
            Dict mapping raw_name to output file path.
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        matched_df = self.match_psm_to_speclib(psm_df)
        if len(matched_df) == 0:
            logger.warning("No PSMs matched to spectral library, nothing to extract")
            return {}

        output_files = {}
        for raw_name, ms_file_path in ms_files.items():
            logger.info(f"Processing raw file: {raw_name}")
            mask = matched_df["raw_name"] == raw_name
            matched_one_raw = matched_df[mask].reset_index(drop=True)
            if len(matched_one_raw) == 0:
                logger.warning(f"No matched PSMs for {raw_name}, skipping")
                continue
            logger.info(f"Extracting XICs for {len(matched_one_raw)} precursors from {raw_name}")
            xic_df = self.extract_one_raw(
                matched_df=matched_one_raw,
                ms_file_path=ms_file_path,
                ms_file_type=ms_file_type,
            )
            if output_format.lower() == "csv":
                output_file = output_path / f"{raw_name}.xic.csv"
                xic_df.to_csv(output_file, index=False)
            elif output_format.lower() == "parquet":
                output_file = output_path / f"{raw_name}.xic.parquet"
                xic_df.to_parquet(output_file, index=False)
            else:
                raise ValueError(f"Unsupported output format: {output_format}. Use 'parquet' or 'csv'.")
            logger.info(f"Saved {len(xic_df)} XIC records to {output_file}")
            output_files[raw_name] = output_file
        return output_files

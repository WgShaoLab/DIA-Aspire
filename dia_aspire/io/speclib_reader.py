# dia_aspire_rescore/io/speclib_reader.py

from pathlib import Path
from typing import Optional, Dict, List
import numpy as np
import pandas as pd
from alphabase.psm_reader.maxquant_reader import ModifiedSequenceReader
from alphabase.psm_reader.psm_reader import psm_reader_provider, psm_reader_yaml
from alphabase.psm_reader.utils import get_column_mapping_for_df
from alphabase.yaml_utils import load_yaml

PACKAGE_DIR = Path(__file__).parent.parent
CONFIG_DIR = PACKAGE_DIR / "constants"
SPECLIB_CONFIG_PATH = CONFIG_DIR / "speclib_reader.yaml"


def _load_speclib_config():
    if not SPECLIB_CONFIG_PATH.exists():
        raise FileNotFoundError(f"Configuration file not found: {SPECLIB_CONFIG_PATH}")

    custom_config = load_yaml(SPECLIB_CONFIG_PATH)
    speclib_config = custom_config.get("spec_lib")

    if not speclib_config:
        raise ValueError("'spec_lib' section not found in config file")

    psm_reader_yaml["spec_lib"] = speclib_config
    return speclib_config


_speclib_config = _load_speclib_config()


class SpecLibReader(ModifiedSequenceReader):
    """Reader for spectral library TSV files.

    This reader parses spectral library files where each row represents a transition.
    It produces two DataFrames:
    - precursor_df: One row per precursor with transition_start_idx and transition_stop_idx
    - transition_df: All transitions with mz, intensity, charge, and type columns

    Examples
    --------
    >>> reader = SpecLibReader()
    >>> reader.import_file("library.tsv")
    >>> precursor_df = reader.precursor_df
    >>> transition_df = reader.transition_df
    """

    _reader_type = "spec_lib"
    _add_unimod_to_mod_mapping = True

    def __init__(
        self,
        *,
        column_mapping: Optional[Dict] = None,
        modification_mapping: Optional[Dict] = None,
        mod_seq_columns: Optional[List[str]] = None,
        rt_unit: Optional[str] = None,
        **kwargs,
    ):
        """Initialize the SpecLibReader.

        Parameters
        ----------
        column_mapping : dict, optional
            Column mapping for precursor columns. If None, read from config.
        modification_mapping : dict, optional
            Modification mapping. If None, use default UniMod mappings.
        mod_seq_columns : list[str], optional
            Columns containing modified sequences. If None, read from config.
        rt_unit : str, optional
            RT unit. If None, read from config.
        """
        super().__init__(
            column_mapping=column_mapping,
            modification_mapping=modification_mapping,
            mod_seq_columns=mod_seq_columns,
            fdr=1.0,  # no fdr
            keep_decoy=True,  # no decoy
            rt_unit=rt_unit,
            **kwargs,
        )

        self._transition_df = pd.DataFrame()
        self._transition_column_mapping = psm_reader_yaml[self._reader_type].get("transition_column_mapping", {})

    @property
    def precursor_df(self) -> pd.DataFrame:
        return self._psm_df

    @property
    def transition_df(self) -> pd.DataFrame:
        """
        Transition DataFrame.

        Columns:
            - mz: m/z
            - intensity: Library intensity
            - charge: Fragment charge
            - type: Fragment type, e.g., "y7", "b3", "y5-17", "m3:8"
        """
        return self._transition_df

    def import_file(self, file_path: str) -> pd.DataFrame:
        """
        Import a spectral library TSV file.

        Parameters
        ----------
        file_path : str
            Path to the spectral library TSV file.

        Returns
        -------
        tuple[pd.DataFrame, pd.DataFrame]
            Tuple of (precursor_df, transition_df).
        """
        origin_df = self._load_file(file_path)

        if len(origin_df) == 0:
            self._psm_df = pd.DataFrame()
            self._transition_df = pd.DataFrame()
            return self._psm_df

        precursor_id_col = self.column_mapping.get("precursor_id", "ions")
        if precursor_id_col not in origin_df.columns:
            raise ValueError(
                f"Precursor ID column '{precursor_id_col}' not found in file. "
                f"Available columns: {list(origin_df.columns)}"
            )

        self.mod_seq_column = self._get_actual_column(self._mod_seq_columns, origin_df)

        origin_df = origin_df.sort_values(precursor_id_col).reset_index(drop=True)

        self._build_transition_df(origin_df)
        self._build_precursor_df(origin_df, precursor_id_col)

        if self.mod_seq_column is not None:
            precursor_origin_df = origin_df.drop_duplicates(subset=[precursor_id_col]).reset_index(drop=True)
            self._load_modifications(precursor_origin_df)
            self._translate_modifications()

        return self.precursor_df, self.transition_df

    def _build_transition_df(self, origin_df: pd.DataFrame) -> None:
        """
        Build the transition DataFrame from the origin DataFrame.

        Parameters
        ----------
        origin_df : pd.DataFrame
            Raw spectral library DataFrame sorted by precursor.
        """
        rename_dict = {v: k for k, v in self._transition_column_mapping.items()}
        self._transition_df = origin_df.rename(columns=rename_dict)[list(self._transition_column_mapping.keys())]
        self._transition_df["mz"] = self._transition_df["mz"].astype(np.float32)
        self._transition_df["intensity"] = self._transition_df["intensity"].astype(np.float32)
        self._transition_df["charge"] = self._transition_df["charge"].astype(np.int8)

    def _build_precursor_df(self, origin_df: pd.DataFrame, precursor_id_col: str) -> None:
        """
        Build the precursor DataFrame with transition indices.

        Parameters
        ----------
        origin_df : pd.DataFrame
            Raw spectral library DataFrame sorted by precursor.
        precursor_id_col : str
            Column name for precursor identification.
        """
        grouped = origin_df.groupby(precursor_id_col)
        first_idx = grouped.head(1).index
        group_sizes = grouped.size()

        column_mapping_for_df = get_column_mapping_for_df(self.column_mapping, origin_df)
        precursor_origin = origin_df.loc[first_idx, list(column_mapping_for_df.values())]

        self._psm_df = precursor_origin.rename(columns={v: k for k, v in column_mapping_for_df.items()}).reset_index(
            drop=True
        )

        transition_stop_idx = np.cumsum(group_sizes.values)
        transition_start_idx = np.concatenate([[0], transition_stop_idx[:-1]])

        self._psm_df["transition_start_idx"] = transition_start_idx
        self._psm_df["transition_stop_idx"] = transition_stop_idx


psm_reader_provider.register_reader("spec_lib", SpecLibReader)


def read_speclib(file_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read a spectral library TSV file.

    Parameters
    ----------
    file_path : str
        Path to the spectral library TSV file.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        Tuple of (precursor_df, transition_df).
    """
    speclib_config = psm_reader_yaml["spec_lib"]

    reader = SpecLibReader(
        column_mapping=speclib_config.get("column_mapping"),
        modification_mapping=speclib_config.get("modification_mapping"),
    )
    return reader.import_file(file_path)

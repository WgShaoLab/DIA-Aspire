# dia_aspire_rescore/io/diann2_reader.py

from pathlib import Path

import numpy as np
import pandas as pd
from alphabase.psm_reader.keys import PsmDfCols
from alphabase.psm_reader.maxquant_reader import ModifiedSequenceReader
from alphabase.psm_reader.psm_reader import psm_reader_provider, psm_reader_yaml
from alphabase.yaml_utils import load_yaml
from .utils import correct_fdr

PACKAGE_DIR = Path(__file__).parent.parent
CONFIG_DIR = PACKAGE_DIR / "constants"
DIANN2_CONFIG_PATH = CONFIG_DIR / "diann2_reader.yaml"


def _load_diann2_config():
    """Load DIANN2 config and inject into alphabase's global config."""
    if not DIANN2_CONFIG_PATH.exists():
        raise FileNotFoundError(f"Configuration file not found: {DIANN2_CONFIG_PATH}")

    custom_config = load_yaml(DIANN2_CONFIG_PATH)
    diann2_config = custom_config.get("diann2")

    if not diann2_config:
        raise ValueError("'diann2' section not found in config file")

    psm_reader_yaml["diann2"] = diann2_config
    return diann2_config


_diann2_config = _load_diann2_config()


class Diann2Reader(ModifiedSequenceReader):
    """Reader for DIANN 2.x parquet output."""

    _reader_type = "diann2"

    def _translate_decoy(self) -> None:
        """DIANN 2.x has 'Decoy' column directly."""
        self._psm_df[PsmDfCols.DECOY] = (self._psm_df[PsmDfCols.DECOY] == 1).astype(np.int8)


psm_reader_provider.register_reader("diann2", Diann2Reader)


def read_diann2(
    file_path: str,
    fdr: float = 1,
    fdrv: float = 0.01,
    keep_decoy: bool = True,
    MBR: bool = False,
) -> pd.DataFrame:
    """
    Read DIANN 2.x parquet report file.

    Parameters
    ----------
    file_path : str
        Path to the DIANN 2.x parquet file.
    fdr : float, optional
        FDR threshold for filtering PSMs. Default is 1.0 (no filtering).
    keep_decoy : bool, optional
        Whether to keep decoy PSMs. Default is False.

    Returns
    -------
    pd.DataFrame
        Processed PSM DataFrame.

    """
    diann2_config = psm_reader_yaml["diann2"]

    reader = psm_reader_provider.get_reader(
        reader_type="diann2",
        column_mapping=diann2_config.get("column_mapping"),
        modification_mapping=diann2_config.get("modification_mapping"),
        fdr=fdr,
        keep_decoy=keep_decoy,
    )
    psm_df = reader.import_file(file_path)
    psm_dfs = correct_fdr(psm_df,fdrv,MBR)
    return psm_dfs

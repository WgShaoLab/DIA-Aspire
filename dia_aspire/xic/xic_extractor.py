"""Command-line interface for DIA-Aspire-Rescore."""

import logging
import sys
import warnings
#from dia_aspire_rescore.config import FineTuneConfig, IOConfig  # noqa: E402
from alpharaw import register_all_readers
from dia_aspire.extraction import XICExtractor
from dia_aspire.io import find_ms_files, read_diann2, read_speclib

logger = logging.getLogger(__name__)

def setup_warnings() -> None:
    """Configure warning filters before importing heavy dependencies."""
    warnings.filterwarnings(
        "ignore",
        message="mask_modloss is deprecated",
        category=UserWarning,
        module="peptdeep.model.ms2",
    )

    warnings.filterwarnings(
        "ignore",
        message="Dotnet-based dependencies could not be loaded",
        category=UserWarning,
        module="alpharaw",  # only involve mzml files in this project
    )

setup_warnings()

def setup_logging(verbose: bool) -> None:
    """
    Configure logging for the application.

    Parameters
    ----------
    verbose : bool
        If True, set log level to DEBUG; otherwise INFO.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
        force=True,
    )

    logging.getLogger("alphabase").setLevel(logging.ERROR)
    logging.getLogger("peptdeep").setLevel(logging.ERROR)
    logging.getLogger("alpharaw").setLevel(logging.ERROR)
    # Re-apply warning filters after logging setup (in case logging reset them)
    setup_warnings()

def extract_xic(
    report: str,
    speclib: str,
    ms_file_dir: str,
    ms_file_type: str,
    output_dir: str = './output',
    output_format: str = 'parquet',
    n_isotopes: int = 1,
    rt_extension: float = 50,
    ppm_tolerance: float = 20,
    verbose: bool = False,
):
    """
    Extract XIC traces for MS1 isotopes and MS2 transitions.
    Outputs one file per raw file with columns: precursor_id, level, mz, type, charge, rt, intensity.
    """
    setup_logging(verbose)
    register_all_readers()

    ms_files = find_ms_files(ms_file_dir, ms_file_type)

    #logger.info(f"Loading spectral library from {speclib}")
    precursor_df, transition_df = read_speclib(speclib)
    #logger.info(f"Loaded {len(precursor_df)} precursors, {len(transition_df)} transitions")

    #logger.info(f"Loading DIA-NN report from {report}")
    psm_df = read_diann2(report)
    #logger.info(f"Loaded {len(psm_df)} PSMs")

    extractor = XICExtractor(
        precursor_df=precursor_df,
        transition_df=transition_df,
        n_isotopes=n_isotopes,
        rt_extension=rt_extension,
        ppm_tolerance=ppm_tolerance,
    )
    extractor.extract_all(
        psm_df=psm_df,
        ms_files=ms_files,
        ms_file_type=ms_file_type,
        output_dir=output_dir,
        output_format=output_format,
    )
    logger.info("XIC extraction completed.")

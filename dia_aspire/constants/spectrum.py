"""Column name constants for DIA data."""

from alphabase.psm_reader.keys import PsmDfCols


class SpectrumDfCols:
    """Column names for DIA spectrum DataFrame."""

    SPEC_IDX = PsmDfCols.SPEC_IDX
    RT = PsmDfCols.RT
    PRECURSOR_MZ = PsmDfCols.PRECURSOR_MZ
    ISOLATION_LOWER_MZ = "isolation_lower_mz"
    ISOLATION_UPPER_MZ = "isolation_upper_mz"
    MS_LEVEL = "ms_level"
    PEAK_START_IDX = "peak_start_idx"
    PEAK_STOP_IDX = "peak_stop_idx"
    FRAG_START_IDX = "frag_start_idx"
    FRAG_STOP_IDX = "frag_stop_idx"

class PsmDfColsExt(PsmDfCols):
    """Column names for extended PSM dataframe."""

    SPEC_START_IDX = "spec_start_idx"
    SPEC_STOP_IDX = "spec_stop_idx"
    FRAG_START_IDX = "frag_start_idx"
    FRAG_STOP_IDX = "frag_stop_idx"

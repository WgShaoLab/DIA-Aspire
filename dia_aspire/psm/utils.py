import pandas as pd
from alphabase.peptide.precursor import refine_precursor_df

from dia_aspire.constants.spectrum import PsmDfColsExt


def refine_matcher_results(
    psm_df: pd.DataFrame,
    fragment_mz_df: pd.DataFrame,
    matched_intensity_df: pd.DataFrame,
    matched_mz_err_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Refine the matcher results by sorting psm_df by nAA and resetting the fragment indices.

    Parameters
    ----------
    psm_df : pd.DataFrame
        The PSM dataframe to refine.
    fragment_mz_df : pd.DataFrame
        The fragment mz dataframe to refine.
    matched_intensity_df : pd.DataFrame
        The matched intensity dataframe to refine.
    matched_mz_err_df : pd.DataFrame
        The matched mz error dataframe to refine.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]
        The refined PSM dataframe, fragment mz dataframe, matched intensity dataframe, and matched mz error dataframe.
    """
    psm_df = refine_precursor_df(psm_df)
    reordered_psm_df = reset_frag_idx(psm_df)
    reordered_fragment_mz_df = order_matched_matrix(
        reordered_psm_df,
        psm_df,
        fragment_mz_df,
    )
    reordered_matched_intensity_df = order_matched_matrix(
        reordered_psm_df,
        psm_df,
        matched_intensity_df,
    )
    reordered_matched_mz_err_df = order_matched_matrix(
        reordered_psm_df,
        psm_df,
        matched_mz_err_df,
    )
    return reordered_psm_df, reordered_fragment_mz_df, reordered_matched_intensity_df, reordered_matched_mz_err_df


def prepare_finetuning_input(
    psm_df: pd.DataFrame,
    intensity_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Prepare the input for PeptDeep finetuning.
    The psm_df should be sorted by nAA and reset the fragment indices.

    Parameters
    ----------
    psm_df : pd.DataFrame
        The PSM dataframe to prepare.
    intensity_df : pd.DataFrame
        The intensity dataframe to prepare.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        The reordered PSM dataframe and intensity dataframe.
    """
    psm_df = refine_precursor_df(psm_df)
    reordered_psm_df = reset_frag_idx(psm_df)
    reordered_intensity_df = order_matched_matrix(
        reordered_psm_df,
        psm_df,
        intensity_df,
    )
    return reordered_psm_df, reordered_intensity_df


def reset_frag_idx(
    psm_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Reset the frag_start_idx and frag_stop_idx of the dataframe so both columns will be monotonically increasing.

    Parameters
    ----------
    psm_df : pd.DataFrame
        The dataframe to reset the indices.

    Returns
    -------
    pd.DataFrame
        The dataframe with the reset indices.

    Notes
    -----
    This function is an internal helper function.
    """
    psm_df_new = psm_df.copy()
    number_of_frags = psm_df_new.frag_stop_idx - psm_df_new.frag_start_idx
    accumulated_frags = number_of_frags.cumsum()

    new_frag_start_idx = accumulated_frags - number_of_frags
    new_frag_stop_idx = accumulated_frags

    psm_df_new[PsmDfColsExt.FRAG_START_IDX] = new_frag_start_idx
    psm_df_new[PsmDfColsExt.FRAG_STOP_IDX] = new_frag_stop_idx
    return psm_df_new


def order_matched_matrix(
    refined_psm_df: pd.DataFrame,
    original_psm_df: pd.DataFrame,
    original_matched_matrix: pd.DataFrame,
) -> pd.DataFrame:
    """
    Rearrange the fragment intensities to match the order used by the start and stop indices in refined_psm_df.
    The goal of this is to reorder the fragment intensities using a newer psm_df that only has different start and stop indices.

    Parameters
    ----------
    refined_psm_df : pd.DataFrame
        The dataframe with the new frag_start_idx and frag_stop_idx to respect.
    original_psm_df : pd.DataFrame
        The dataframe with the old frag_start_idx and frag_stop_idx.
    original_matched_matrix : pd.DataFrame
        'fragment_mz_df', 'matched_intensity_df', 'matched_mz_err_df' from matcher.

    Returns
    -------
    pd.DataFrame
        The reordered fragment intensity dataframe.
    """
    reordered = original_matched_matrix.copy()
    for i in range(len(original_psm_df)):
        new_start_idx = refined_psm_df.iloc[i][PsmDfColsExt.FRAG_START_IDX]
        new_end_idx = refined_psm_df.iloc[i][PsmDfColsExt.FRAG_STOP_IDX]

        old_start_idx = original_psm_df.iloc[i][PsmDfColsExt.FRAG_START_IDX]
        old_end_idx = original_psm_df.iloc[i][PsmDfColsExt.FRAG_STOP_IDX]

        reordered.iloc[new_start_idx:new_end_idx, :] = original_matched_matrix.iloc[old_start_idx:old_end_idx, :]
    return reordered

from typing import ClassVar, Optional

import numpy as np
import pandas as pd
from peptdeep.model.ms2 import calc_ms2_similarity
from peptdeep.pretrained_models import ModelManager

from dia_aspire.config import MS2MatchConfig
from dia_aspire.features.base import BaseFeatureGenerator
from dia_aspire.psm.matcher import DIAPeptideSpectrumMatcher
from dia_aspire.psm.utils import refine_matcher_results


class MS2FeatureGenerator(BaseFeatureGenerator):
    """
    Generate MS2 spectrum similarity and fragment matching features.

    Features include:
    - Similarity metrics (cos, sa, spc, pcc) for all fragments, b-ions, and y-ions
    - Weighted fragment scores
    - Fragment matching statistics

    TODO: include spectral entropy
    TODO: https://oktoberfest.readthedocs.io/en/stable/svm_features.html
    """

    DEFAULT_FRAG_TYPES: ClassVar[list[str]] = ["b_z1", "b_z2", "y_z1", "y_z2"]
    DEFAULT_SPC_TOP_K: ClassVar[int] = 0  # use all fragments for SPC calculation

    def __init__(
        self,
        model_mgr: ModelManager,
        ms_files: dict[str, str],
        ms_file_type: str = "mzml",
        ms2_match_config: Optional[MS2MatchConfig] = None,
    ):
        """
        Initialize MS2FeatureGenerator.

        Parameters
        ----------
        model_mgr : ModelManager
            ModelManager from peptdeep for intensity prediction.
        ms_files : dict[str, str]
            Dict mapping raw_name to file path.
        ms_file_type : str, optional
            MS file type, by default "mzml".
        ms2_match_config : MS2MatchConfig, optional
            MS2 matching config. If None, uses defaults.
        """
        self.model_mgr = model_mgr
        self.ms_files = ms_files
        self.ms_file_type = ms_file_type
        self.ms2_match_config = ms2_match_config or MS2MatchConfig()

        self.frag_types = self.DEFAULT_FRAG_TYPES
        self.spc_top_k = self.DEFAULT_SPC_TOP_K
        self.used_frag_types = [
            frag_type for frag_type in self.frag_types if frag_type in self.model_mgr.ms2_model.charged_frag_types
        ]

    @property
    def feature_names(self) -> list[str]:
        return [
            # Overall similarity metrics
            "cos",
            "sa",
            "spc",
            "pcc",
            # B-ion metrics
            "cos_bion",
            "sa_bion",
            "spc_bion",
            "pcc_bion",
            # Y-ion metrics
            "cos_yion",
            "sa_yion",
            "spc_yion",
            "pcc_yion",
            # Weighted scores
            "merr_weighted_frag_score",
            "pred_weighted_frag_score",
            "merr_weighted_bion_score",
            "pred_weighted_bion_score",
            "merr_weighted_yion_score",
            "pred_weighted_yion_score",
            # Fragment statistics - all fragments
            "matched_frag_num",
            "matched_frag_ratio",
            "both_matched_pred_frag_num",
            "both_matched_pred_frag_to_matched",
            "both_matched_pred_frag_to_pred",
            "matched_not_pred_frag_num",
            "matched_not_pred_frag_ratio",
            "pred_not_matched_frag_num",
            "pred_not_matched_frag_ratio",
            "matched_frag_rel_to_pred",
            "pred_frag_rel_to_matched",
            # Fragment statistics - b ions
            "matched_bion_num",
            "matched_bion_ratio",
            "both_matched_pred_bion_num",
            "both_matched_pred_bion_to_matched",
            "both_matched_pred_bion_to_pred",
            "matched_not_pred_bion_num",
            "matched_not_pred_bion_ratio",
            "pred_not_matched_bion_num",
            "pred_not_matched_bion_ratio",
            "matched_bion_rel_to_pred",
            "pred_bion_rel_to_matched",
            # Fragment statistics - y ions
            "matched_yion_num",
            "matched_yion_ratio",
            "both_matched_pred_yion_num",
            "both_matched_pred_yion_to_matched",
            "both_matched_pred_yion_to_pred",
            "matched_not_pred_yion_num",
            "matched_not_pred_yion_ratio",
            "pred_not_matched_yion_num",
            "pred_not_matched_yion_ratio",
            "matched_yion_rel_to_pred",
            "pred_yion_rel_to_matched",
        ]

    def generate(self, psm_df: pd.DataFrame) -> pd.DataFrame:
        """
        Generate MS2 features. Handles MS2 matching internally.

        Parameters
        ----------
        psm_df : pd.DataFrame
            PSM DataFrame.

        Returns
        -------
        pd.DataFrame
            psm_df with added feature columns.
        """
        matcher = DIAPeptideSpectrumMatcher(
            match_closest=self.ms2_match_config.match_closest,
            use_ppm=self.ms2_match_config.use_ppm,
            tol_value=self.ms2_match_config.tolerance,
            n_neighbors=0,
        )

        psm_df, _, matched_intensity_df, matched_mz_err_df = matcher.match_ms2_multi_raw(
            psm_df, self.ms_files, self.ms_file_type
        )

        psm_df, _, matched_intensity_df, matched_mz_err_df = refine_matcher_results(
            psm_df, _, matched_intensity_df, matched_mz_err_df
        )

        predict_intensity_df = self.model_mgr.predict_ms2(psm_df)
        predict_intensity_df = predict_intensity_df[self.used_frag_types]

        psm_df = self._calculate_similarity_all_frags(
            psm_df, predict_intensity_df, matched_intensity_df, matched_mz_err_df
        )
        psm_df = self._calculate_similarity_bions(psm_df, predict_intensity_df, matched_intensity_df, matched_mz_err_df)
        psm_df = self._calculate_similarity_yions(psm_df, predict_intensity_df, matched_intensity_df, matched_mz_err_df)

        return psm_df

    def _calculate_similarity_all_frags(
        self,
        psm_df: pd.DataFrame,
        predict_intensity_df: pd.DataFrame,
        matched_intensity_df: pd.DataFrame,
        matched_mz_err_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Calculate similarity metrics for all fragments."""
        used_frag_types = self.used_frag_types

        psm_df, _ = calc_ms2_similarity(
            psm_df,
            predict_intensity_df,
            matched_intensity_df,
            charged_frag_types=used_frag_types,
            metrics=["COS", "SA", "SPC", "PCC"],
            spc_top_k=self.spc_top_k,
        )

        psm_df.rename(
            columns={
                "COS": "cos",
                "SA": "sa",
                "SPC": "spc",
                "PCC": "pcc",
            },
            inplace=True,
        )

        psm_df = self._get_psm_scores(
            psm_df,
            predict_intensity_df=predict_intensity_df[used_frag_types],
            matched_intensity_df=matched_intensity_df[used_frag_types],
            matched_mass_err_df=matched_mz_err_df[used_frag_types],
        )

        psm_df.rename(
            columns={
                "merr_weighted_score": "merr_weighted_frag_score",
                "pred_weighted_score": "pred_weighted_frag_score",
            },
            inplace=True,
        )

        has_matched_intens = matched_intensity_df[used_frag_types].values > 0
        has_predicted_intens = predict_intensity_df[used_frag_types].values > 0.001
        has_both_matched_predicted = has_matched_intens & has_predicted_intens

        (
            psm_df["matched_frag_num"],
            psm_df["matched_frag_ratio"],
            psm_df["both_matched_pred_frag_num"],
            psm_df["both_matched_pred_frag_to_matched"],
            psm_df["both_matched_pred_frag_to_pred"],
            psm_df["matched_not_pred_frag_num"],
            psm_df["matched_not_pred_frag_ratio"],
            psm_df["pred_not_matched_frag_num"],
            psm_df["pred_not_matched_frag_ratio"],
            psm_df["matched_frag_rel_to_pred"],
            psm_df["pred_frag_rel_to_matched"],
        ) = zip(
            *psm_df[["frag_start_idx", "frag_stop_idx"]].apply(
                self._get_frag_features,
                axis=1,
                matched_inten_values=matched_intensity_df[used_frag_types].values,
                predicted_inten_values=predict_intensity_df[used_frag_types].values,
                has_matched_intens=has_matched_intens,
                has_predicted_intens=has_predicted_intens,
                has_both_matched_predicted=has_both_matched_predicted,
            )
        )

        return psm_df

    def _calculate_similarity_bions(
        self,
        psm_df: pd.DataFrame,
        predict_intensity_df: pd.DataFrame,
        matched_intensity_df: pd.DataFrame,
        matched_mz_err_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Calculate similarity metrics for b-ions."""
        b_frag_types = [_t for _t in self.used_frag_types if _t.startswith("b")]

        if len(b_frag_types) > 0:
            # Calculate similarity metrics
            psm_df, _ = calc_ms2_similarity(
                psm_df,
                predict_intensity_df,
                matched_intensity_df,
                charged_frag_types=b_frag_types,
                metrics=["COS", "SA", "SPC", "PCC"],
            )

            psm_df.rename(
                columns={
                    "COS": "cos_bion",
                    "SA": "sa_bion",
                    "SPC": "spc_bion",
                    "PCC": "pcc_bion",
                },
                inplace=True,
            )

            psm_df = self._get_psm_scores(
                psm_df,
                predict_intensity_df=predict_intensity_df[b_frag_types],
                matched_intensity_df=matched_intensity_df[b_frag_types],
                matched_mass_err_df=matched_mz_err_df[b_frag_types],
            )

            psm_df.rename(
                columns={
                    "merr_weighted_score": "merr_weighted_bion_score",
                    "pred_weighted_score": "pred_weighted_bion_score",
                },
                inplace=True,
            )

            has_matched_intens = matched_intensity_df[b_frag_types].values > 0
            has_predicted_intens = predict_intensity_df[b_frag_types].values > 0
            has_both_matched_predicted = has_matched_intens & has_predicted_intens

            (
                psm_df["matched_bion_num"],
                psm_df["matched_bion_ratio"],
                psm_df["both_matched_pred_bion_num"],
                psm_df["both_matched_pred_bion_to_matched"],
                psm_df["both_matched_pred_bion_to_pred"],
                psm_df["matched_not_pred_bion_num"],
                psm_df["matched_not_pred_bion_ratio"],
                psm_df["pred_not_matched_bion_num"],
                psm_df["pred_not_matched_bion_ratio"],
                psm_df["matched_bion_rel_to_pred"],
                psm_df["pred_bion_rel_to_matched"],
            ) = zip(
                *psm_df[["frag_start_idx", "frag_stop_idx"]].apply(
                    self._get_frag_features,
                    axis=1,
                    matched_inten_values=matched_intensity_df[b_frag_types].values,
                    predicted_inten_values=predict_intensity_df[b_frag_types].values,
                    has_matched_intens=has_matched_intens,
                    has_predicted_intens=has_predicted_intens,
                    has_both_matched_predicted=has_both_matched_predicted,
                )
            )
        else:
            # Set all b-ion features to 0
            psm_df[
                [
                    "cos_bion",
                    "sa_bion",
                    "spc_bion",
                    "pcc_bion",
                    "merr_weighted_bion_score",
                    "pred_weighted_bion_score",
                    "matched_bion_num",
                    "matched_bion_ratio",
                    "both_matched_pred_bion_num",
                    "both_matched_pred_bion_to_matched",
                    "both_matched_pred_bion_to_pred",
                    "matched_not_pred_bion_num",
                    "matched_not_pred_bion_ratio",
                    "pred_not_matched_bion_num",
                    "pred_not_matched_bion_ratio",
                    "matched_bion_rel_to_pred",
                    "pred_bion_rel_to_matched",
                ]
            ] = 0

        return psm_df

    def _calculate_similarity_yions(
        self,
        psm_df: pd.DataFrame,
        predict_intensity_df: pd.DataFrame,
        matched_intensity_df: pd.DataFrame,
        matched_mz_err_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Calculate similarity metrics for y-ions."""
        y_frag_types = [_t for _t in self.used_frag_types if _t.startswith("y")]

        if len(y_frag_types) > 0:
            psm_df, _ = calc_ms2_similarity(
                psm_df,
                predict_intensity_df,
                matched_intensity_df,
                charged_frag_types=y_frag_types,
                metrics=["COS", "SA", "SPC", "PCC"],
            )

            psm_df.rename(
                columns={
                    "COS": "cos_yion",
                    "SA": "sa_yion",
                    "SPC": "spc_yion",
                    "PCC": "pcc_yion",
                },
                inplace=True,
            )

            psm_df = self._get_psm_scores(
                psm_df,
                predict_intensity_df=predict_intensity_df[y_frag_types],
                matched_intensity_df=matched_intensity_df[y_frag_types],
                matched_mass_err_df=matched_mz_err_df[y_frag_types],
            )

            psm_df.rename(
                columns={
                    "merr_weighted_score": "merr_weighted_yion_score",
                    "pred_weighted_score": "pred_weighted_yion_score",
                },
                inplace=True,
            )

            has_matched_intens = matched_intensity_df[y_frag_types].values > 0
            has_predicted_intens = predict_intensity_df[y_frag_types].values > 0
            has_both_matched_predicted = has_matched_intens & has_predicted_intens

            (
                psm_df["matched_yion_num"],
                psm_df["matched_yion_ratio"],
                psm_df["both_matched_pred_yion_num"],
                psm_df["both_matched_pred_yion_to_matched"],
                psm_df["both_matched_pred_yion_to_pred"],
                psm_df["matched_not_pred_yion_num"],
                psm_df["matched_not_pred_yion_ratio"],
                psm_df["pred_not_matched_yion_num"],
                psm_df["pred_not_matched_yion_ratio"],
                psm_df["matched_yion_rel_to_pred"],
                psm_df["pred_yion_rel_to_matched"],
            ) = zip(
                *psm_df[["frag_start_idx", "frag_stop_idx"]].apply(
                    self._get_frag_features,
                    axis=1,
                    matched_inten_values=matched_intensity_df[y_frag_types].values,
                    predicted_inten_values=predict_intensity_df[y_frag_types].values,
                    has_matched_intens=has_matched_intens,
                    has_predicted_intens=has_predicted_intens,
                    has_both_matched_predicted=has_both_matched_predicted,
                )
            )
        else:
            # Set all y-ion features to 0
            psm_df[
                [
                    "cos_yion",
                    "sa_yion",
                    "spc_yion",
                    "pcc_yion",
                    "merr_weighted_yion_score",
                    "pred_weighted_yion_score",
                    "matched_yion_num",
                    "matched_yion_ratio",
                    "both_matched_pred_yion_num",
                    "both_matched_pred_yion_to_matched",
                    "both_matched_pred_yion_to_pred",
                    "matched_not_pred_yion_num",
                    "matched_not_pred_yion_ratio",
                    "pred_not_matched_yion_num",
                    "pred_not_matched_yion_ratio",
                    "matched_yion_rel_to_pred",
                    "pred_yion_rel_to_matched",
                ]
            ] = 0

        return psm_df

    def _get_psm_scores(
        self,
        psm_df: pd.DataFrame,
        predict_intensity_df: pd.DataFrame,
        matched_intensity_df: pd.DataFrame,
        matched_mass_err_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Calculate weighted scores from predictions and matches.

        This calculates two types of scores:
        1. merr_weighted_score: Score weighted by mass error
        2. pred_weighted_score: Score weighted by both mass error and predictions
        """
        matched_norm_intensity_df = pd.DataFrame(
            np.log(matched_intensity_df.values + 1),
            columns=matched_intensity_df.columns.values,
        )

        matched_merr_weight_df = matched_mass_err_df.mask(matched_mass_err_df > 1000000, 0).abs()
        max_merr = matched_merr_weight_df.values.max()
        if max_merr > 0:
            matched_merr_weight_df /= max_merr
        matched_merr_weight_df = 1 - matched_merr_weight_df.pow(4)

        # intensity weighted by mass error
        peak_score_df = matched_norm_intensity_df * matched_merr_weight_df
        # intensity weighted by both mass error and predictions
        pred_weighted_score_df = peak_score_df * predict_intensity_df

        def _get_one_score(
            frag_start_end,
            peak_score_values,
            pred_weighted_score_values,
        ):
            frag_start, frag_end = frag_start_end
            frag_ratio = (peak_score_values[frag_start:frag_end] > 0).mean() ** 0.5
            return (
                peak_score_values[frag_start:frag_end].sum() * frag_ratio,
                pred_weighted_score_values[frag_start:frag_end].sum() * frag_ratio,
            )

        (
            psm_df["merr_weighted_score"],
            psm_df["pred_weighted_score"],
        ) = zip(
            *psm_df[["frag_start_idx", "frag_stop_idx"]].apply(
                _get_one_score,
                axis=1,
                peak_score_values=peak_score_df.values,
                pred_weighted_score_values=pred_weighted_score_df.values,
            )
        )

        return psm_df

    @staticmethod
    def _get_frag_features(
        frag_start_end,
        matched_inten_values,
        predicted_inten_values,
        has_matched_intens,
        has_predicted_intens,
        has_both_matched_predicted,
    ):
        """
        Calculate fragment matching statistics.

        Returns various metrics about matched and predicted fragments including:
        - Number and ratio of matched fragments
        - Overlap between matched and predicted fragments
        - Relative intensities
        """
        frag_start, frag_end = frag_start_end

        matched_frag_num = has_matched_intens[frag_start:frag_end].sum(dtype=np.float32)
        pred_frag_num = has_predicted_intens[frag_start:frag_end].sum(dtype=np.float32)
        matched_frag_ratio = matched_frag_num / (matched_inten_values.shape[1] * (frag_end - frag_start))

        both_matched_pred_frag_num = has_both_matched_predicted[frag_start:frag_end].sum(dtype=np.float32)

        matched_not_pred_frag_num = (
            has_matched_intens[frag_start:frag_end] & ~has_both_matched_predicted[frag_start:frag_end]
        ).sum(dtype=np.float32)

        pred_not_matched_frag_num = (
            has_predicted_intens[frag_start:frag_end] & ~has_both_matched_predicted[frag_start:frag_end]
        ).sum(dtype=np.float32)

        if matched_frag_num > 0:
            both_matched_pred_frag_to_matched = both_matched_pred_frag_num / matched_frag_num
            matched_not_pred_frag_ratio = matched_not_pred_frag_num / matched_frag_num
        else:
            both_matched_pred_frag_to_matched = 0
            matched_not_pred_frag_ratio = 0

        if pred_frag_num > 0:
            both_matched_pred_frag_to_pred = both_matched_pred_frag_num / pred_frag_num
            pred_not_matched_frag_ratio = pred_not_matched_frag_num / pred_frag_num
        else:
            both_matched_pred_frag_to_pred = 0
            pred_not_matched_frag_ratio = 0

        matched_frag_rel_to_pred = matched_inten_values[frag_start:frag_end][
            has_predicted_intens[frag_start:frag_end]
        ].sum()
        if matched_frag_rel_to_pred > 0:
            matched_frag_rel_to_pred /= matched_inten_values[frag_start:frag_end].sum()

        pred_frag_rel_to_matched = predicted_inten_values[frag_start:frag_end][
            has_matched_intens[frag_start:frag_end]
        ].sum()
        if pred_frag_rel_to_matched > 0:
            pred_frag_rel_to_matched /= predicted_inten_values[frag_start:frag_end].sum()

        return (
            matched_frag_num,
            matched_frag_ratio,
            both_matched_pred_frag_num,
            both_matched_pred_frag_to_matched,
            both_matched_pred_frag_to_pred,
            matched_not_pred_frag_num,
            matched_not_pred_frag_ratio,
            pred_not_matched_frag_num,
            pred_not_matched_frag_ratio,
            matched_frag_rel_to_pred,
            pred_frag_rel_to_matched,
        )

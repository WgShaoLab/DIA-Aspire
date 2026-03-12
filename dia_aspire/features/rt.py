import numpy as np
import pandas as pd
from peptdeep.pretrained_models import ModelManager

from dia_aspire.features.base import BaseFeatureGenerator


class RTFeatureGenerator(BaseFeatureGenerator):
    """
    Generate retention time prediction and delta features.

    Features:
    - rt_pred: Predicted retention time
    - rt_delta: Difference between predicted and observed (calibrated)
    - rt_delta_abs: Absolute difference between predicted and observed
    - rt_ratio: Ratio of predicted and observed
    """

    def __init__(self, model_mgr: ModelManager):
        """
        Initialize RTFeatureGenerator.

        Parameters
        ----------
        model_mgr : ModelManager
            ModelManager from peptdeep for RT prediction
        """
        self.model_mgr = model_mgr

    @property
    def feature_names(self) -> list[str]:
        """Get the names of the features."""
        return ["rt_pred", "delta_rt", "abs_rt_delta", "rt_ratio"]

    def generate(self, psm_df: pd.DataFrame) -> pd.DataFrame:
        """
        Generate RT features and add them to psm_df.

        Parameters
        ----------
        psm_df : pd.DataFrame

        Returns
        -------
        pd.DataFrame
            psm_df with RT feature columns added
        """
        if "rt_norm" in psm_df.columns:
            psm_df = self.model_mgr.predict_rt(psm_df)
            psm_df["delta_rt"] = psm_df.rt_pred - psm_df.rt_norm
            psm_df["abs_rt_delta"] = psm_df["delta_rt"].abs()
            psm_df["rt_ratio"] = np.minimum(psm_df.rt_pred, psm_df.rt_norm) / np.maximum(psm_df.rt_pred, psm_df.rt_norm)

        else:
            psm_df["rt_pred"] = 0
            psm_df["delta_rt"] = 0
            psm_df["abs_rt_delta"] = 0
            psm_df["rt_ratio"] = 0

        return psm_df

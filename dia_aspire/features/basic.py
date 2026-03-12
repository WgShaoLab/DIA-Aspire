import numpy as np
import pandas as pd

from dia_aspire.features.base import BaseFeatureGenerator


class BasicFeatureGenerator(BaseFeatureGenerator):
    """
    Generate basic peptide-level features independent of predictions.

    Features:
    - charge_1 through charge_7: Charge one-hot encoding (charge_gt_6 for charges > 6)
    - mod_num: Count of modifications (excluding Carbamidomethyl@C)
    - pep_length: Length of the peptide sequence
    - precursor_mz: Precursor m/z of the peptide
    """

    def __init__(self):
        pass

    @property
    def feature_names(self) -> list[str]:
        """Get the names of the features."""
        return [
            "charge_1",
            "charge_2",
            "charge_3",
            "charge_4",
            "charge_5",
            "charge_6",
            "charge_gt_6",
            "mod_num",
            "pep_length",
            "precursor_mz",  # already in psm_df
        ]

    def generate(self, psm_df: pd.DataFrame) -> pd.DataFrame:
        """
        Generate basic features and add them to psm_df.

        Parameters
        ----------
        psm_df : pd.DataFrame
            PSM DataFrame containing 'charge', 'mods', and 'precursor_mz' columns

        Returns
        -------
        pd.DataFrame
            psm_df with basic feature columns added
        """

        def _charge_one_hot(ch):
            x = [0] * 7
            if ch > 6:
                x[-1] = 1
            else:
                x[ch - 1] = 1
            return tuple(x)

        (
            psm_df["charge_1"],
            psm_df["charge_2"],
            psm_df["charge_3"],
            psm_df["charge_4"],
            psm_df["charge_5"],
            psm_df["charge_6"],
            psm_df["charge_gt_6"],
        ) = zip(*psm_df.charge.astype(np.int8).apply(_charge_one_hot))

        def _mod_count(mods):
            if not mods:
                return 0
            mod_count = 0
            for mod in mods.split(";"):
                if mod != "Carbamidomethyl@C":
                    mod_count += 1
            return mod_count

        psm_df["mod_num"] = psm_df.mods.apply(_mod_count)
        psm_df["pep_length"] = psm_df.sequence.apply(len)

        return psm_df

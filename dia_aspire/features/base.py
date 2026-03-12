from abc import ABC, abstractmethod

import pandas as pd


class BaseFeatureGenerator(ABC):
    """
    Interface for feature generators.
    """

    @property
    @abstractmethod
    def feature_names(self):
        """
        Get the names of the features.
        """
        pass

    @abstractmethod
    def generate(self, psm_df: pd.DataFrame, *args, **kwargs) -> pd.DataFrame:
        """
        Generate and add features to psm_df.
        """
        pass

from dataclasses import asdict, dataclass
from typing import Any

import yaml  # type: ignore[import-untyped]


@dataclass
class ConfigBase:
    """
    Base class for all configuration classes.
    """

    def to_dict(self) -> dict[str, Any]:
        """
        Convert config to dictionary.

        Returns
        -------
        dict
            Dictionary representation of the config
        """
        return asdict(self)

    def to_yaml(self, path: str) -> None:
        """
        Save configuration to YAML file.

        Parameters
        ----------
        path : str
            Path to save YAML file
        """
        with open(path, "w") as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False)

    @classmethod
    def from_dict(cls, config_dict: dict[str, Any]):
        """
        Create config from dictionary.

        Parameters
        ----------
        config_dict : dict
            Dictionary containing configuration values

        Returns
        -------
        ConfigBase
            Configuration instance
        """
        return cls(**config_dict)

    @classmethod
    def from_yaml(cls, path: str):
        """
        Load configuration from YAML file.

        Parameters
        ----------
        path : str
            Path to YAML file

        Returns
        -------
        ConfigBase
            Configuration instance
        """
        with open(path) as f:
            config_dict = yaml.safe_load(f)
        return cls.from_dict(config_dict)

"""Configuration for MS2 spectrum matching."""

from dataclasses import dataclass

from dia_aspire.config.base import ConfigBase


@dataclass
class MS2MatchConfig(ConfigBase):
    """MS2 matching configuration."""

    # TODO: add option to calibrate mz error
    match_closest: bool = True
    use_ppm: bool = True
    tolerance: float = 20.0

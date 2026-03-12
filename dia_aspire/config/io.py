"""Configuration for input/output paths and formats."""

from dataclasses import dataclass

from dia_aspire.config.base import ConfigBase


@dataclass
class IOConfig(ConfigBase):
    """Input/output configuration."""

    report_file: str  # DIA-NN parquet file path
    ms_file_dir: str  # Directory containing MS files
    ms_file_type: str = "mzml"
    output_dir: str = "./output"

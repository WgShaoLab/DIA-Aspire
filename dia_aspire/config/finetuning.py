"""Configuration for fine-tuning AlphaPeptDeep models."""

from dataclasses import dataclass
from typing import Optional

from dia_aspire.config import ConfigBase


@dataclass
class FineTuneConfig(ConfigBase):
    """Fine-tuning configuration for MS2 and RT models."""

    # FDR threshold for selecting training PSMs
    fdr_threshold: float = 0.01

    # Model settings
    instrument: str = "QE"
    nce: float = 27
    mask_modloss: bool = True
    device: Optional[str] = None

    # MS2 training parameters
    psm_num_to_train_ms2: int = 10000
    psm_num_to_test_ms2: int = 0
    epoch_to_train_ms2: int = 20
    warmup_epoch_to_train_ms2: int = 10
    batch_size_to_train_ms2: int = 512
    lr_to_train_ms2: float = 0.0001
    psm_num_per_mod_to_train_ms2: int = 50
    top_n_mods_to_train: int = 10

    # RT training parameters
    psm_num_to_train_rt_ccs: int = 10000
    epoch_to_train_rt_ccs: int = 40
    warmup_epoch_to_train_rt_ccs: int = 10
    batch_size_to_train_rt_ccs: int = 1024
    lr_to_train_rt_ccs: float = 0.0001
    psm_num_per_mod_to_train_rt_ccs: int = 50

    # General settings
    train_verbose: bool = False
    thread_num: int = 36
    # model_path: Optional[str] = None # TODO: support external models

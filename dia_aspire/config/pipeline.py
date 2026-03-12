"""Pipeline for DIA rescoring workflow."""

import logging
from pathlib import Path
from typing import Optional

import pandas as pd
from alphabase.peptide.precursor import refine_precursor_df
from alpharaw import register_all_readers

from dia_aspire.config.finetuning import FineTuneConfig
from dia_aspire.config.io import IOConfig
from dia_aspire.config.ms2_match import MS2MatchConfig
from dia_aspire.features import (
    BasicFeatureGenerator,
    MS2FeatureGenerator,
    RTFeatureGenerator,
    XICFeatureGenerator,
)
from dia_aspire.finetuning import FineTuner
from dia_aspire.io import find_ms_files, read_diann2
from dia_aspire.rescore.perform_rescore import PerformRescore

logger = logging.getLogger(__name__)


class Pipeline:
    """
    DIA rescoring pipeline.

    Workflow:
    1. load_psm()
    2. finetune() (optional)
    3. generate_features()
    4. rescore() (placeholder)
    5. write_report()
    """

    # TODO: unify and clean up the configs
    # ms2_match_config is not used in the generators

    def __init__(
        self,
        io_config: IOConfig,
        ms2_match_config: Optional[MS2MatchConfig] = None,
        finetune_config: Optional[FineTuneConfig] = None,
        feature_generators: Optional[list[str]] = None,
        model: str='DNN',
    ):
        self.io_config = io_config
        self.ms2_match_config = ms2_match_config or MS2MatchConfig()
        self.finetune_config = finetune_config or FineTuneConfig()
        self.feature_generators = feature_generators or ["basic", "xic", "ms2", "rt"]

        self.psm_df: Optional[pd.DataFrame] = None
        self.finetuner: Optional[FineTuner] = None
        self.ms_files: dict[str, str] = {}
        self.model = model

        self._setup()

    def _setup(self):
        register_all_readers()
        Path(self.io_config.output_dir).mkdir(parents=True, exist_ok=True)
        self.ms_files = find_ms_files(self.io_config.ms_file_dir, self.io_config.ms_file_type)
        if not self.ms_files:
            raise ValueError(f"No {self.io_config.ms_file_type} files found in {self.io_config.ms_file_dir}")
        logger.info(f"Found {len(self.ms_files)} MS files")

    def run(self, skip_finetuning: bool = False) -> pd.DataFrame:
        """Full pipeline."""
        self.load_psm()
        if not skip_finetuning:
            self.finetune()
        self.generate_features()
        self.rescore()
        return self.psm_df

    def run_feature_generation(self, skip_finetuning: bool = False) -> pd.DataFrame:
        """Run pipeline up to feature generation (no rescoring)."""
        self.load_psm()
        if not skip_finetuning:
            self.finetune()
        self.generate_features()
        return self.psm_df

    def load_psm(self):
        """Load PSM data from report file."""
        self.psm_df = read_diann2(self.io_config.report_file)
        self.psm_df = refine_precursor_df(self.psm_df)
        logger.info(f"Loaded {len(self.psm_df)} PSMs")

    def finetune(self):
        """Finetune peptdeep models."""
        self.finetuner = FineTuner(self.finetune_config)
        self.finetuner.load_pretrained("generic")
        self.finetuner.train(
            self.psm_df,
            self.ms_files,
            ms_file_type=self.io_config.ms_file_type,
            ms2_match_config=self.ms2_match_config,
        )
        logger.info("Finetuning completed")

    def generate_features(self):
        """Generate rescoring features."""
        if self.finetuner is None:
            self.finetuner = FineTuner(self.finetune_config)
            self.finetuner.load_pretrained("generic")

        for gen_name in self.feature_generators:
            generator = self._create_generator(gen_name)
            if generator is None:
                logger.warning(f"Unknown generator: {gen_name}")
                continue

            self.psm_df = generator.generate(self.psm_df)
            logger.info(f"Generating {gen_name} features...")
            #logger.info(f"Added {len(generator.feature_names)} features")
        logger.info(f"psm with features generated!")

        output_path = Path(self.io_config.output_dir) / "psm_with_features_xic.csv"
        self.psm_df.to_csv(output_path, index=False)
        logger.info(f"Saved to {output_path}")

    def _create_generator(self, name: str):
        if name == "basic":
            return BasicFeatureGenerator()

        if name == "xic":
            return XICFeatureGenerator(
                ms_files=self.ms_files,
                ms_file_type=self.io_config.ms_file_type,
                ppm_tolerance=self.ms2_match_config.tolerance,
                ms2_match_config=self.ms2_match_config,
                enable_smoothing=False,
            )

        # ms2 and rt generators require finetuner
        if self.finetuner is None:
            raise RuntimeError(f"Finetuner is required for '{name}' generator but was not initialized")

        if name == "ms2":
            return MS2FeatureGenerator(
                model_mgr=self.finetuner.model_manager,
                ms_files=self.ms_files,
                ms_file_type=self.io_config.ms_file_type,
                ms2_match_config=self.ms2_match_config,
            )
        elif name == "rt":
            return RTFeatureGenerator(model_mgr=self.finetuner.model_manager)
        return None

    def rescore(self):
        self.psm_df = PerformRescore(self.psm_df, Path(self.io_config.output_dir), self.model)
        outp = Path(self.io_config.output_dir) / "rescore_pr_matrix.csv"
        logger.info(f"Matrix result saved to {outp}")
        return True

    # def write_report(self):
    #     output_path = Path(self.io_config.output_dir) / "psm_with_features_xic.csv"
    #     self.psm_df.to_csv(output_path, index=False)
    #     logger.info(f"Saved to {output_path}")

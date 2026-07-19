#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MLHuaiNanzi temperature prediction analyzer for MethArCT

Predicts microbial growth temperature range (T_min, T_opt, T_max) from
genome-wide amino acid composition (AAC + dipeptide composition) features
using the embedded MLHuaiNanzi PLS regression models.
"""

import os
import pandas as pd
from pathlib import Path
from typing import Dict, Optional, Union
from ..utils.config import Config
from ..utils.logger import get_logger
from ..utils.file_utils import FileUtils


class HuainanziAnalyzer:
    """MLHuaiNanzi growth temperature range predictor"""

    def __init__(self, config: Optional[Config] = None):
        self.config = config or Config()
        self.logger = get_logger("huainanzi_analyzer")

        # Results directory
        self.results_dir = self.config.get('output.base_dir', 'results')
        FileUtils.ensure_dir(self.results_dir)

        self._predict_mod = None  # lazy-loaded module reference
        self.tool_available = self._check_availability()

    # ------------------------------------------------------------------
    # Availability check
    # ------------------------------------------------------------------
    def _check_availability(self) -> bool:
        """Check whether the embedded MLHuaiNanzi Python module can be imported."""
        try:
            from .mlhuainanzi.predict import predict_temperatures  # noqa: F401
            self.logger.info("MLHuaiNanzi module available (embedded)")
            return True
        except ImportError as e:
            self.logger.warning(f"MLHuaiNanzi module not available: {e}")
            return False

    def _ensure_loaded(self):
        """Lazy-load MLHuaiNanzi prediction functions."""
        if self._predict_mod is None:
            from .mlhuainanzi import predict as _mod
            self._predict_mod = _mod

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def predict_temperature(self,
                            input_file: Union[str, Path],
                            output_prefix: Optional[str] = None) -> Dict:
        """
        Predict growth temperature range for a genome FASTA file.

        Args:
            input_file: Path to protein FASTA file (.faa / .fasta)
            output_prefix: Prefix for output files

        Returns:
            Analysis results dictionary with keys:
                status, input_file, analysis_timestamp, tool,
                prediction: {T_min, T_opt, T_max} in °C,
                summary: {T_min, T_opt, T_max, temperature_category},
                output_files: {tsv}
        """
        input_file = Path(input_file)

        if not input_file.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")

        if not FileUtils.validate_fasta(input_file):
            raise ValueError(f"Invalid FASTA file: {input_file}")

        if output_prefix is None:
            output_prefix = input_file.stem

        self.logger.info(f"Starting MLHuaiNanzi temperature prediction for {input_file.name}")
        self._ensure_loaded()

        return self._run_prediction(input_file, output_prefix)

    # ------------------------------------------------------------------
    # Internal prediction wrapper
    # ------------------------------------------------------------------
    def _run_prediction(self, input_file: Path, output_prefix: str) -> Dict:
        """Execute MLHuaiNanzi prediction and return structured results."""
        try:
            mod = self._predict_mod

            results = mod.predict_temperatures(str(input_file))
            t_min = results['T_min']
            t_opt = results['T_opt']
            t_max = results['T_max']

            # Determine temperature category
            temperature_category = self._categorize_temperature(t_opt)

            self.logger.info(
                f"MLHuaiNanzi prediction: T_min={t_min:.1f}, T_opt={t_opt:.1f}, "
                f"T_max={t_max:.1f} °C ({temperature_category})"
            )

            # Save TSV output
            FileUtils.ensure_dir(self.results_dir)
            tsv_path = os.path.join(self.results_dir, f"{output_prefix}_Huainanzi_Summary.tsv")
            with open(tsv_path, 'w', encoding='utf-8') as f:
                f.write("Genome\tT_min (°C)\tT_opt (°C)\tT_max (°C)\tCategory\n")
                f.write(f"{input_file.name}\t{t_min:.2f}\t{t_opt:.2f}\t{t_max:.2f}\t{temperature_category}\n")

            result = {
                'status': 'success',
                'input_file': str(input_file),
                'analysis_timestamp': pd.Timestamp.now().isoformat(),
                'tool': 'MLHuaiNanzi',
                'prediction': {
                    'T_min': round(t_min, 2),
                    'T_opt': round(t_opt, 2),
                    'T_max': round(t_max, 2),
                },
                'output_files': {
                    'tsv': tsv_path,
                },
                'summary': {
                    'T_min': round(t_min, 2),
                    'T_opt': round(t_opt, 2),
                    'T_max': round(t_max, 2),
                    'temperature_category': temperature_category,
                },
            }
            return result

        except Exception as e:
            self.logger.error(f"MLHuaiNanzi prediction failed: {e}")
            return {
                'status': 'failed',
                'error': str(e),
            }

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _categorize_temperature(t_opt: float) -> str:
        """Classify organism by optimal growth temperature."""
        if t_opt < 20:
            return 'Psychrophilic'
        elif t_opt < 40:
            return 'Mesophilic'
        elif t_opt < 60:
            return 'Thermophilic'
        else:
            return 'Hyperthermophilic'

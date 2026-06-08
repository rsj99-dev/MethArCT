#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pH preference analyzer for MethArCT

Predicts microbial growth pH preference (optimum, maximum, minimum) from
genome-wide amino acid composition features using pre-trained GenomeSpot
Lasso regression models.
"""

import os
import json
import pandas as pd
from pathlib import Path
from typing import Dict, Optional, Union

from ..utils.config import Config
from ..utils.logger import get_logger
from ..utils.file_utils import FileUtils


class PHAnalyzer:
    """
    pH preference analyzer based on GenomeSpot pre-trained models.

    Predicts pH optimum, pH maximum, and pH minimum from whole-genome
    protein sequences using amino acid composition features.
    """

    def __init__(self, config: Optional[Config] = None):
        self.config = config or Config()
        self.logger = get_logger("ph_analyzer")

        # Results directory
        self.results_dir = self.config.get('output.base_dir', 'results')
        FileUtils.ensure_dir(self.results_dir)

        self._ph_predictor = None  # lazy-loaded module reference
        self.tool_available = self._check_ph_availability()

    # ------------------------------------------------------------------
    # Availability check
    # ------------------------------------------------------------------
    def _check_ph_availability(self) -> bool:
        """Check whether the embedded pH predictor module can be imported."""
        try:
            from .ph_predictor.predictor import PHPredictor  # noqa: F401
            from .ph_predictor.feature_extractor import extract_features  # noqa: F401
            self.logger.info("pH predictor module available (embedded)")
            return True
        except ImportError as e:
            self.logger.warning(f"pH predictor module not available: {e}")
            return False

    def _ensure_loaded(self):
        """Lazy-load pH predictor functions."""
        if self._ph_predictor is None:
            from .ph_predictor import predictor as _mod
            self._ph_predictor = _mod

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def predict_ph(self,
                   input_file: Union[str, Path],
                   output_prefix: Optional[str] = None) -> Dict:
        """
        Predict pH preference for a genome FASTA file.

        Args:
            input_file: Path to protein FASTA file (.faa / .fasta)
            output_prefix: Prefix for output files

        Returns:
            Analysis results dictionary
        """
        input_file = Path(input_file)

        if not input_file.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")

        if not FileUtils.validate_fasta(input_file):
            raise ValueError(f"Invalid FASTA file: {input_file}")

        if output_prefix is None:
            output_prefix = input_file.stem

        self.logger.info(f"Starting pH prediction for {input_file.name}")
        self._ensure_loaded()

        # Run prediction
        prediction_result = self._run_prediction(input_file, output_prefix)
        return prediction_result

    # ------------------------------------------------------------------
    # Internal prediction wrapper
    # ------------------------------------------------------------------
    def _run_prediction(self, input_file: Path, output_prefix: str) -> Dict:
        """
        Execute pH prediction and return structured results.
        """
        try:
            from .ph_predictor import DEFAULT_MODEL_DIR
            from .ph_predictor.feature_extractor import extract_features
            from .ph_predictor.predictor import PHPredictor, format_predictions_tsv

            # 1. Extract features
            features = extract_features(str(input_file))
            if features is None:
                return {
                    'status': 'failed',
                    'error': 'Failed to extract features from FASTA file',
                }

            # 2. Load model and predict
            predictor = PHPredictor(model_dir=DEFAULT_MODEL_DIR)
            predictions = predictor.predict(features)

            # 3. Build structured results
            ph_values = {}
            for target, result in predictions.items():
                ph_values[target] = {
                    'value': float(result['value']) if result.get('value') is not None else None,
                    'error': float(result['error']) if result.get('error') is not None else None,
                    'is_novel': bool(result['is_novel']) if result.get('is_novel') is not None else False,
                    'warning': str(result['warning']) if result.get('warning') is not None else None,
                    'units': result.get('units', 'pH'),
                }

            self.logger.info(
                f"pH prediction complete - optimum: {ph_values.get('ph_optimum', {}).get('value')}, "
                f"max: {ph_values.get('ph_max', {}).get('value')}, "
                f"min: {ph_values.get('ph_min', {}).get('value')}"
            )

            # 4. Save outputs
            output_dir = self.results_dir
            FileUtils.ensure_dir(output_dir)

            # TSV output
            tsv_path = os.path.join(output_dir, f"{output_prefix}_pH_Summary.tsv")
            tsv_content = format_predictions_tsv(predictions)
            with open(tsv_path, 'w', encoding='utf-8') as f:
                f.write(tsv_content + "\n")

            # JSON output (full features + predictions)
            json_path = os.path.join(output_dir, f"{output_prefix}_pH_Details.json")
            with open(json_path, 'w', encoding='utf-8') as f:
                json.dump({
                    'input_file': str(input_file),
                    'predictions': ph_values,
                }, f, indent=2, ensure_ascii=False)

            # 5. Build result dict
            is_novel = any(r.get('is_novel', False) for r in ph_values.values())

            result = {
                'status': 'success',
                'input_file': str(input_file),
                'analysis_timestamp': pd.Timestamp.now().isoformat(),
                'tool': 'pH_predictor (GenomeSpot)',
                'prediction': ph_values,
                'is_novel': is_novel,
                'output_files': {
                    'tsv': tsv_path,
                    'json': json_path,
                },
                'summary': {
                    'ph_optimum': ph_values.get('ph_optimum', {}).get('value'),
                    'ph_max': ph_values.get('ph_max', {}).get('value'),
                    'ph_min': ph_values.get('ph_min', {}).get('value'),
                    'is_novel': is_novel,
                },
            }
            return result

        except Exception as e:
            self.logger.error(f"pH prediction failed: {e}")
            return {
                'status': 'failed',
                'error': str(e),
            }

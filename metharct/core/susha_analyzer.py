#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SuSha salinity adaptation analyzer for MethArCT

Predicts microbial salinity adaptation from genome-wide amino acid composition
features using the SuSha ensemble learning model.
"""

import os
import pandas as pd
from pathlib import Path
from typing import Dict, Optional, Union
from ..utils.config import Config
from ..utils.logger import get_logger
from ..utils.file_utils import FileUtils


# Salinity label mapping (mirrors susha.config.LABEL_MAP)
SALINITY_LABEL_MAP = {
    0: "Salt-sensitive",
    1: "Halotolerant",
    2: "Slight halophilic",
    3: "Moderate halophilic",
    4: "Extreme halophilic",
}


class SuShaAnalyzer:
    """SuSha salinity adaptation analyzer"""

    def __init__(self, config: Optional[Config] = None):
        self.config = config or Config()
        self.logger = get_logger("susha_analyzer")

        # Results directory
        self.results_dir = self.config.get('output.base_dir', 'results')
        FileUtils.ensure_dir(self.results_dir)

        self._susha_predict = None  # lazy-loaded module reference
        self.tool_available = self._check_susha_availability()

    # ------------------------------------------------------------------
    # Availability check
    # ------------------------------------------------------------------
    def _check_susha_availability(self) -> bool:
        """Check whether the embedded SuSha Python module can be imported."""
        try:
            from .susha.predict import run_prediction, process_fasta  # noqa: F401
            self.logger.info("SuSha module available (embedded)")
            return True
        except ImportError as e:
            self.logger.warning(f"SuSha module not available: {e}")
            return False

    def _ensure_loaded(self):
        """Lazy-load SuSha prediction functions."""
        if self._susha_predict is None:
            from .susha import predict as _mod
            self._susha_predict = _mod

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def predict_salinity(self,
                         input_file: Union[str, Path],
                         output_prefix: Optional[str] = None) -> Dict:
        """
        Predict salinity adaptation for a genome FASTA file.

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

        self.logger.info(f"Starting SuSha salinity prediction for {input_file.name}")
        self._ensure_loaded()

        # Run prediction
        prediction_result = self._run_prediction(input_file, output_prefix)
        return prediction_result

    # ------------------------------------------------------------------
    # Internal prediction wrapper
    # ------------------------------------------------------------------
    def _run_prediction(self, input_file: Path, output_prefix: str) -> Dict:
        """
        Execute SuSha prediction and return structured results.
        """
        try:
            mod = self._susha_predict

            # 1. Extract features
            features = mod.process_fasta(str(input_file))
            if features is None:
                return {
                    'status': 'failed',
                    'error': 'Failed to extract features from FASTA file',
                }

            X_df = pd.DataFrame([features])[mod.FEATURE_COLS]

            # 2. Load model and predict
            import joblib
            model = joblib.load(mod.DEFAULT_MODEL_PATH)

            pred_idx = int(model.predict(X_df)[0])
            pred_label = SALINITY_LABEL_MAP.get(pred_idx, str(pred_idx))

            probs = model.predict_proba(X_df)[0]
            max_prob = float(max(probs))

            top3_indices = list(map(int, sorted(range(len(probs)), key=lambda i: probs[i], reverse=True)[:3]))
            top3 = [
                {'rank': rank, 'label': SALINITY_LABEL_MAP[int(idx)], 'probability': float(probs[idx])}
                for rank, idx in enumerate(top3_indices, 1)
            ]

            self.logger.info(f"SuSha prediction: {pred_label} (Confidence: {max_prob:.2%})")

            # 3. SHAP interpretation (optional, non-blocking)
            shap_contributions = {}
            try:
                shap_contributions = mod.interpret_model(model, X_df, pred_idx)
            except Exception as e:
                self.logger.warning(f"SHAP interpretation skipped: {e}")

            # 4. Save outputs
            output_dir = self.results_dir
            FileUtils.ensure_dir(output_dir)

            tsv_path = os.path.join(output_dir, f"{output_prefix}_SuSha_Summary.tsv")
            with open(tsv_path, 'w', encoding='utf-8') as f:
                f.write("Genome\tPredicted_Salinity\tConfidence\n")
                f.write(f"{input_file.name}\t{pred_label}\t{max_prob:.4f}\n")

            # Optionally save Excel (non-blocking if openpyxl missing)
            excel_path = os.path.join(output_dir, f"{output_prefix}_SuSha_Result.xlsx")
            try:
                mod.run_prediction(input_file, os.path.join(output_dir, output_prefix))
            except Exception as e:
                self.logger.debug(f"SuSha Excel output skipped: {e}")

            # 5. Build result dict
            result = {
                'status': 'success',
                'input_file': str(input_file),
                'analysis_timestamp': pd.Timestamp.now().isoformat(),
                'tool': 'SuSha',
                'prediction': {
                    'salinity_label': pred_label,
                    'salinity_index': pred_idx,
                    'confidence': max_prob,
                },
                'top3_predictions': top3,
                'output_files': {
                    'tsv': tsv_path,
                },
                'summary': {
                    'predicted_salinity': pred_label,
                    'confidence': max_prob,
                },
            }
            return result

        except Exception as e:
            self.logger.error(f"SuSha prediction failed: {e}")
            return {
                'status': 'failed',
                'error': str(e),
            }

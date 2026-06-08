"""
pH predictor module.

Loads GenomeSpot pre-trained models and performs pH prediction on extracted features.

Supported prediction targets:
  - ph_optimum: optimal growth pH
  - ph_max: maximum tolerated pH
  - ph_min: minimum tolerated pH

Each prediction result contains:
  - value: predicted value
  - error: error estimate (RMSE)
  - is_novel: novelty flag (True means input features differ significantly from training set)
  - warning: warning flag (set when predicted value exceeds valid range and is clipped)
  - units: unit (pH)
"""

import logging
import os
from typing import Dict, Optional, Tuple

import joblib
import numpy as np

from .feature_extractor import PH_FEATURE_NAMES, features_to_array

logger = logging.getLogger(__name__)

# Valid range for pH predictions (consistent with GenomeSpot)
PH_BOUNDS = (0.5, 14.0)

# Supported prediction targets
PH_TARGETS = ["ph_optimum", "ph_max", "ph_min"]


class PHPredictor:
    """
    pH predictor based on GenomeSpot pre-trained models.

    Usage:
        predictor = PHPredictor(model_dir="GenomeSPOT_pH")
        results = predictor.predict(features_dict)

    Args:
        model_dir: Path to model directory containing the following files:
            - ph_optimum.joblib, ph_max.joblib, ph_min.joblib
            - novelty_ph.joblib
            - error_ph_optimum.joblib, error_ph_max.joblib, error_ph_min.joblib
    """

    def __init__(self, model_dir: str):
        self.model_dir = model_dir
        self._validate_model_dir()
        self._load_models()

    def _validate_model_dir(self):
        """Validate that model directory contains all required files."""
        required_files = [
            "ph_optimum.joblib",
            "ph_max.joblib",
            "ph_min.joblib",
            "novelty_ph.joblib",
            "error_ph_optimum.joblib",
            "error_ph_max.joblib",
            "error_ph_min.joblib",
        ]
        missing = [f for f in required_files if not os.path.exists(os.path.join(self.model_dir, f))]
        if missing:
            raise FileNotFoundError(
                f"Model directory {self.model_dir} is missing the following files:\n" + "\n".join(f"  - {f}" for f in missing)
            )

    def _load_models(self):
        """Load all model files into memory."""
        logger.info("Loading models from %s...", self.model_dir)

        # Load prediction models (sklearn Pipeline)
        self.models = {}
        for target in PH_TARGETS:
            model_path = os.path.join(self.model_dir, f"{target}.joblib")
            self.models[target] = joblib.load(model_path)
            logger.info("  Loaded: %s", target)

        # Load error models (numpy arrays)
        self.error_models = {}
        for target in PH_TARGETS:
            error_path = os.path.join(self.model_dir, f"error_{target}.joblib")
            self.error_models[target] = joblib.load(error_path)
            logger.info("  Loaded: error_%s", target)

        # Load novelty detection model (OneClassSVM)
        novelty_path = os.path.join(self.model_dir, "novelty_ph.joblib")
        self.novelty_model = joblib.load(novelty_path)
        logger.info("  Loaded: novelty_ph")

    def predict_error(self, y_pred: float, error_arr: np.ndarray) -> float:
        """
        Look up the error estimate corresponding to the predicted value.

        The error model is an N x 2 array:
          Column 0: reference predicted value
          Column 1: corresponding RMSE
        Finds the nearest reference point to y_pred and returns its RMSE.

        Args:
            y_pred: current predicted value
            error_arr: error reference array (N, 2)

        Returns:
            Estimated RMSE
        """
        closest_idx = np.argmin(np.abs(error_arr[:, 0] - y_pred))
        _, ref_err = error_arr[closest_idx]
        return float(ref_err)

    def predict_novelty(self, X: np.ndarray) -> bool:
        """
        Novelty detection: determine if input features differ significantly from training set.
    
        OneClassSVM returns 1 for 'not novel' (within training set range),
        and -1 for 'novel' (anomalous sample outside training set).
    
        Args:
            X: feature array of shape (1, 60)
    
        Returns:
            True if novel (prediction may be less reliable), False if normal
        """
        result = self.novelty_model.predict(X)
        return result[0] != 1

    def check_prediction_range(self, y_pred: float) -> Tuple[float, Optional[str]]:
        """
        Check if predicted value is within valid range; clip and warn if exceeded.

        Args:
            y_pred: raw predicted value

        Returns:
            (clipped prediction value, warning string or None)
        """
        min_val, max_val = PH_BOUNDS
        if y_pred < min_val:
            return min_val, "min_exceeded"
        elif y_pred > max_val:
            return max_val, "max_exceeded"
        else:
            return y_pred, None

    def predict(self, features: Dict[str, float]) -> Dict[str, dict]:
        """
        Perform pH prediction on extracted features.

        Args:
            features: feature dictionary (60 features, keys are feature names)

        Returns:
            Prediction result dictionary with structure:
            {
                "ph_optimum": {"value": float, "error": float, "is_novel": bool, "warning": str|None, "units": "pH"},
                "ph_max":     {"value": float, "error": float, "is_novel": bool, "warning": str|None, "units": "pH"},
                "ph_min":     {"value": float, "error": float, "is_novel": bool, "warning": str|None, "units": "pH"},
            }
        """
        # Convert feature dictionary to model input array
        X = features_to_array(features, PH_FEATURE_NAMES)

        # Check feature completeness
        if np.any(np.isnan(X[0])):
            missing_feats = [
                name for name, val in zip(PH_FEATURE_NAMES, X[0])
                if np.isnan(val)
            ]
            logger.warning("Missing features (NaN values): %s", missing_feats[:5])
            warning = "genome missing features"
            return {
                target: {"value": None, "error": None, "is_novel": None, "warning": warning, "units": "pH"}
                for target in PH_TARGETS
            }

        # Novelty detection (shared across all targets)
        is_novel = self.predict_novelty(X)
        if is_novel:
            logger.warning("Input genome features differ significantly from training set; predictions may be less reliable.")

        # Predict for each target
        predictions = {}
        for target in PH_TARGETS:
            model = self.models[target]
            error_model = self.error_models[target]

            # Predict
            y_pred = float(model.predict(X)[0])

            # Range check
            y_pred, warning = self.check_prediction_range(y_pred)

            # Error estimate
            error = self.predict_error(y_pred, error_model)

            predictions[target] = {
                "value": y_pred,
                "error": error,
                "is_novel": is_novel,
                "warning": warning,
                "units": "pH",
            }

            logger.info(
                "%s = %.4f (error: %.4f, novel: %s, warning: %s)",
                target, y_pred, error, is_novel, warning,
            )

        return predictions

    def predict_and_format(self, features: Dict[str, float]) -> str:
        """
        Predict and format results as a TSV string (ready for file writing).

        Returns:
            TSV-formatted prediction result string
        """
        predictions = self.predict(features)
        return format_predictions_tsv(predictions)


def format_predictions_tsv(predictions: Dict[str, dict]) -> str:
    """
    Format prediction results as a TSV string.

    Args:
        predictions: return value of predict()

    Returns:
        TSV-formatted string with header and data rows
    """
    cols = ["target", "value", "error", "is_novel", "warning", "units"]
    lines = ["\t".join(cols)]
    for target in sorted(predictions):
        row = predictions[target]
        line = [
            target,
            str(row["value"]),
            str(row["error"]),
            str(row["is_novel"]),
            str(row["warning"]),
            str(row["units"]),
        ]
        lines.append("\t".join(line))
    return "\n".join(lines)

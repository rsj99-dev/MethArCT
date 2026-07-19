"""
Prediction engine: loads pre-trained PLS models and predicts growth temperatures.
"""

import os
import warnings
import numpy as np
import joblib

from .features import read_fasta, extract_features

# Suppress sklearn cross-version unpickle warnings (models are version-compatible)
warnings.filterwarnings('ignore', message='Trying to unpickle estimator')

_MODELS_DIR = os.path.join(os.path.dirname(__file__), 'models')
TARGETS = ['T_min', 'T_opt', 'T_max']


def predict_temperatures(faa_path):
    """Predict T_min, T_opt, T_max for a given .faa proteome file.

    Returns dict with keys: genome_name, valid_residues, T_min, T_opt, T_max.
    """
    seq = read_fasta(faa_path)
    if len(seq) == 0:
        raise ValueError(f"No valid protein residues found in: {faa_path}")

    features = extract_features(seq)
    X = np.array(features, dtype=float).reshape(1, -1)

    results = {
        'genome_name': os.path.basename(faa_path).rsplit('.', 1)[0],
        'valid_residues': len(seq),
    }

    for target in TARGETS:
        pkl_path = os.path.join(_MODELS_DIR, f'model_{target}.pkl')
        pipe = joblib.load(pkl_path)
        results[target] = float(pipe.predict(X).ravel()[0])

    return results

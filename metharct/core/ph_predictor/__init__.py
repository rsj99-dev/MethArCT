"""
pH_predictor: genome-wide pH preference prediction from protein sequences.

Based on GenomeSpot (https://github.com/cultivarium/genomespot),
using pre-trained Lasso regression models to predict microbial pH preference
from amino acid composition features.

This module is embedded within MethArCT as a core analysis component.
"""

# Default model directory (bundled alongside this package)
import os

DEFAULT_MODEL_DIR = os.path.join(os.path.dirname(__file__), "models")

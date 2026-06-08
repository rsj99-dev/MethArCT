"""
Signal peptide prediction module.

Uses a Hidden Markov Model (HMM) to predict signal peptides in bacterial and archaeal proteins.
Adapted from GenomeSpot's SignalPeptideHMM class.

References:
- Eddy, S. What is a hidden Markov model?. Nat Biotechnol 22, 1315–1316 (2004).
- Bagos et al., Combined prediction of Tat and Sec signal peptides with HMMs,
  Bioinformatics (2010)
- Nielsen, H., & Krogh, A. Prediction of signal peptides and signal anchors
  by a hidden Markov model. In ISMB (1998).
"""

import os
from typing import Tuple

import joblib
import numpy as np

# Signal peptide HMM model path (distributed with this package)
TRAINED_MODEL = os.path.join(os.path.dirname(__file__), "hmm", "hmm_signal_peptide.joblib")

# Standard amino acid symbols (consistent with HMM training)
SYMBOLS = [
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
]

# HMM hidden states
# B: start region, C/C1/C2/C3: around cleavage site, H: hydrophobic core, M: transmembrane, N/N1/N2/N3: N-terminal region
STATES = ["B", "C", "C1", "C2", "C3", "H", "M", "N", "N1", "N2", "N3"]

# Prediction parameters (consistent with GenomeSpot)
N_TERMINUS_LENGTH = 50          # Only analyze first 50 amino acids at N-terminus
SIGNAL_PEPTIDE_END_STATE = "C1" # C1 state marks signal peptide cleavage site
THRESHOLD_LOG_PROB = -134.0     # Log probability threshold


class SignalPeptideHMM:
    """
    Predict signal peptides in bacterial/archaeal proteins using a pre-trained HMM.

    Prediction logic:
    1. Take the first 50 amino acids at the N-terminus of the protein sequence
    2. Decode with HMM to get the most likely hidden state sequence
    3. If the sequence contains a C1 state and log probability > threshold, signal peptide is present
    4. The index of the C1 state is the signal peptide cleavage site

    Args:
        model_file: path to pre-trained HMM model (defaults to built-in model in this package)
    """

    def __init__(self, model_file: str = TRAINED_MODEL):
        if not os.path.exists(model_file):
            raise FileNotFoundError(
                f"Signal peptide HMM model not found: {model_file}\n"
                "Please ensure the hmm/hmm_signal_peptide.joblib file exists."
            )
        self.model = joblib.load(model_file)
        self.symbols = SYMBOLS
        self.states = STATES
        self.threshold_log_prob = THRESHOLD_LOG_PROB
        self.signal_end_state = SIGNAL_PEPTIDE_END_STATE
        self.nterminus_length = N_TERMINUS_LENGTH
        self.symbol_to_idx = dict(zip(self.symbols, range(len(self.symbols))))
        self.state_to_index = dict(zip(self.states, range(len(self.states))))
        self.idx_to_state = dict(zip(range(len(self.states)), self.states))

    def _format_protein_sequence(self, protein_sequence: str) -> np.ndarray:
        """
        Convert amino acid sequence to HMM input format (integer index array).
        Only the first 50 N-terminal amino acids are used; non-standard amino acids are replaced with glycine (G).
        """
        default_symbol = self.symbol_to_idx.get("G")
        protein_nterminus = protein_sequence[0: self.nterminus_length]
        arr_sequence = np.array(
            [self.symbol_to_idx.get(aa, default_symbol) for aa in protein_nterminus]
        ).reshape(-1, 1)
        return arr_sequence

    def _predict_hidden_states(self, formatted_sequence: np.ndarray):
        """Decode with HMM and return predicted hidden state sequence and log probability."""
        log_prob, pred_state_indices = self.model.decode(formatted_sequence)
        pred_states = [self.idx_to_state[idx] for idx in pred_state_indices]
        return pred_states, log_prob

    def predict_signal_peptide(self, protein_sequence: str) -> Tuple[bool, int]:
        """
        Predict whether a protein sequence contains a signal peptide.

        Args:
            protein_sequence: full amino acid sequence string

        Returns:
            (is_exported, signal_end_index):
                is_exported: True if a signal peptide is detected (protein may be secreted/localized extracellularly)
                signal_end_index: signal peptide cleavage site index (-1 if no signal peptide)
        """
        input_sequence = self._format_protein_sequence(protein_sequence)

        if len(input_sequence) < self.nterminus_length:
            # Sequence too short for reliable prediction
            is_exported = False
            signal_end_index = -1
        else:
            pred_states, log_prob = self._predict_hidden_states(input_sequence)
            has_cut_site = self.signal_end_state in pred_states
            is_exported = (log_prob > self.threshold_log_prob) and has_cut_site

            if is_exported:
                signal_end_index = pred_states.index(self.signal_end_state)
            else:
                signal_end_index = -1

        return is_exported, signal_end_index

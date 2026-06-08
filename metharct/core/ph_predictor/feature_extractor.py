"""
Feature extraction module.

Extracts the 60-dimensional feature vector required for pH prediction from whole-genome protein sequences.

Feature sources (consistent with GenomeSpot):
  - 20 extracellular soluble protein amino acid frequencies (extracellular_soluble_aa_X)
  - 20 intracellular soluble protein amino acid frequencies (intracellular_soluble_aa_X)
  - 20 membrane protein amino acid frequencies (membrane_aa_X)

Protein localization rules:
  - Membrane: GRAVY > mean GRAVY + 0.5
  - Extracellular soluble: has signal peptide AND not membrane
  - Intracellular soluble: no signal peptide AND not membrane

Feature value calculation: length-weighted mean of amino acid frequencies for each protein class.
"""

import logging
from collections import Counter
from typing import Dict, List, Tuple

import numpy as np

from .fasta_parser import read_fasta
from .signal_peptide import SignalPeptideHMM

logger = logging.getLogger(__name__)

# Standard amino acids (20 types, excluding non-standard characters)
STANDARD_AMINO_ACIDS = sorted({
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
})

# Kyte & Doolittle (1982) hydropathy scale
HYDROPHOBICITY = {
    "A":  1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C":  2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I":  4.5,
    "L":  3.8, "K": -3.9, "M":  1.9, "F":  2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V":  4.2,
}

# Membrane protein threshold: GRAVY difference above mean protein GRAVY
DIFF_HYDROPHOBICITY_MEMBRANE = 0.5


def format_sequence(sequence: str) -> str:
    """
    Clean amino acid sequence: remove non-standard characters and convert to uppercase.
    """
    return "".join(
        aa for aa in sequence.strip().upper()
        if aa in set(STANDARD_AMINO_ACIDS)
    )


def compute_aa_1mer_frequencies(sequence: str, start_pos: int = 1) -> Dict[str, float]:
    """
    Compute single amino acid frequencies, skipping the N-terminal Met (position 0).

    Args:
        sequence: cleaned amino acid sequence
        start_pos: start position (default 1, skipping Met)

    Returns:
        {amino_acid: frequency} for 20 standard amino acids
    """
    trimmed = sequence[start_pos:]
    if len(trimmed) < 1:
        return {aa: 0.0 for aa in STANDARD_AMINO_ACIDS}

    counts = Counter(trimmed)
    total = len(trimmed)
    return {aa: counts.get(aa, 0) / total for aa in STANDARD_AMINO_ACIDS}


def compute_gravy(sequence: str, start_pos: int = 1) -> float:
    """
    Compute Grand Average of Hydropathy (GRAVY), i.e. mean protein hydropathy.
    """
    trimmed = sequence[start_pos:]
    if len(trimmed) < 1:
        return np.nan
    return float(np.mean([HYDROPHOBICITY[aa] for aa in trimmed if aa in HYDROPHOBICITY]))


def analyze_proteins(
    sequences: Dict[str, str],
    signal_peptide_model: SignalPeptideHMM,
) -> Dict[str, dict]:
    """
    Analyze each protein: signal peptide prediction, sequence cleaning, feature computation.

    Returns:
        {protein_id: {aa_A, aa_C, ..., gravy, length, is_exported, signal_end_index}}
    """
    protein_data = {}

    for protein_id, raw_seq in sequences.items():
        seq = format_sequence(raw_seq)
        if len(seq) < 2:
            continue

        # Signal peptide prediction
        is_exported, signal_end_index = signal_peptide_model.predict_signal_peptide(seq)

        # Determine start position after cleavage
        start_pos = signal_end_index + 1 if is_exported and signal_end_index >= 0 else 1
        # Ensure start_pos does not exceed sequence length
        if start_pos >= len(seq):
            start_pos = 1

        trimmed_seq = seq[start_pos:]
        length = len(trimmed_seq)
        if length < 1:
            continue

        # Amino acid frequencies (skip region before start_pos)
        aa_freqs = compute_aa_1mer_frequencies(seq, start_pos=start_pos)
        gravy = compute_gravy(seq, start_pos=start_pos)

        metrics = {f"aa_{aa}": freq for aa, freq in aa_freqs.items()}
        metrics["gravy"] = gravy
        metrics["length"] = length
        metrics["is_exported"] = is_exported
        protein_data[protein_id] = metrics

    return protein_data


def assign_localization(protein_data: Dict[str, dict]) -> Dict[str, str]:
    """
    Classify proteins by localization based on signal peptide and hydrophobicity.

    Returns:
        {protein_id: 'membrane' | 'extra_soluble' | 'intra_soluble'}
    """
    # Compute mean GRAVY across all proteins
    gravy_values = [v["gravy"] for v in protein_data.values() if not np.isnan(v["gravy"])]
    if not gravy_values:
        mean_gravy = 0.0
    else:
        mean_gravy = float(np.mean(gravy_values))

    localization = {}
    for protein_id, metrics in protein_data.items():
        gravy = metrics.get("gravy", 0.0)
        if np.isnan(gravy):
            gravy = mean_gravy

        if (gravy - mean_gravy) >= DIFF_HYDROPHOBICITY_MEMBRANE:
            localization[protein_id] = "membrane"
        elif metrics.get("is_exported", False):
            localization[protein_id] = "extra_soluble"
        else:
            localization[protein_id] = "intra_soluble"

    return localization


def length_weighted_average(values: List[float], lengths: List[int]) -> float:
    """Compute length-weighted average."""
    total_length = sum(lengths)
    if total_length == 0:
        return 0.0
    return float(sum(l * v for l, v in zip(lengths, values)) / total_length)


def compute_class_features(
    protein_data: Dict[str, dict],
    protein_ids: List[str],
) -> Dict[str, float]:
    """
    Compute length-weighted amino acid frequency features for a group of proteins.

    Args:
        protein_data: metrics dictionary for all proteins
        protein_ids: list of protein IDs in this class

    Returns:
        {aa_A: freq, aa_C: freq, ..., aa_Y: freq} (20 features)
    """
    if not protein_ids:
        return {f"aa_{aa}": np.nan for aa in STANDARD_AMINO_ACIDS}

    values_by_aa = {aa: [] for aa in STANDARD_AMINO_ACIDS}
    lengths = []

    for pid in protein_ids:
        metrics = protein_data[pid]
        length = metrics.get("length", 0)
        if length <= 0:
            continue
        lengths.append(length)
        for aa in STANDARD_AMINO_ACIDS:
            values_by_aa[aa].append(metrics.get(f"aa_{aa}", 0.0))

    if not lengths:
        return {f"aa_{aa}": np.nan for aa in STANDARD_AMINO_ACIDS}

    result = {}
    for aa in STANDARD_AMINO_ACIDS:
        result[f"aa_{aa}"] = length_weighted_average(values_by_aa[aa], lengths)

    return result


def extract_features(faa_path: str) -> Dict[str, float]:
    """
    Extract 60-dimensional features required for pH prediction from a protein FASTA file.

    Features are ordered as follows (consistent with GenomeSpot models):
      extracellular_soluble_aa_A ... extracellular_soluble_aa_Y  (20)
      intracellular_soluble_aa_A ... intracellular_soluble_aa_Y  (20)
      membrane_aa_A ... membrane_aa_Y                            (20)

    Args:
        faa_path: path to protein FASTA file (supports .gz)

    Returns:
        Feature dictionary with feature names as keys and values
    """
    logger.info("Reading protein sequences: %s", faa_path)
    sequences = read_fasta(faa_path)
    logger.info("Read %d protein sequences in total", len(sequences))

    if not sequences:
        raise ValueError(f"No protein sequences read from file {faa_path}.")

    # Initialize signal peptide model
    logger.info("Loading signal peptide HMM model...")
    sp_model = SignalPeptideHMM()

    # Analyze each protein
    logger.info("Analyzing protein sequences (signal peptide prediction, amino acid frequency computation)...")
    protein_data = analyze_proteins(sequences, sp_model)
    logger.info("Valid protein count: %d", len(protein_data))

    if not protein_data:
        raise ValueError("No valid proteins remaining after processing. Please check the input file.")

    # Protein localization classification
    logger.info("Performing protein localization classification...")
    localization = assign_localization(protein_data)

    # Group by class
    extra_ids = [pid for pid, loc in localization.items() if loc == "extra_soluble"]
    intra_ids = [pid for pid, loc in localization.items() if loc == "intra_soluble"]
    memb_ids  = [pid for pid, loc in localization.items() if loc == "membrane"]

    logger.info(
        "Localization results: extracellular_soluble=%d, intracellular_soluble=%d, membrane=%d",
        len(extra_ids), len(intra_ids), len(memb_ids),
    )

    # Compute features per class
    extra_feats = compute_class_features(protein_data, extra_ids)
    intra_feats = compute_class_features(protein_data, intra_ids)
    memb_feats  = compute_class_features(protein_data, memb_ids)

    # Assemble full feature dictionary (key order consistent with GenomeSpot: sorted by localization, then aa)
    features = {}
    for aa in STANDARD_AMINO_ACIDS:
        features[f"extracellular_soluble_aa_{aa}"] = extra_feats[f"aa_{aa}"]
    for aa in STANDARD_AMINO_ACIDS:
        features[f"intracellular_soluble_aa_{aa}"] = intra_feats[f"aa_{aa}"]
    for aa in STANDARD_AMINO_ACIDS:
        features[f"membrane_aa_{aa}"] = memb_feats[f"aa_{aa}"]

    return features


# 60 feature names required by pH models (in model input order)
PH_FEATURE_NAMES = [
    f"{loc}_aa_{aa}"
    for loc in ["extracellular_soluble", "intracellular_soluble", "membrane"]
    for aa in STANDARD_AMINO_ACIDS
]


def features_to_array(features: Dict[str, float], feature_names: List[str] = None) -> np.ndarray:
    """
    Convert feature dictionary to numpy array for model input.

    Args:
        features: feature dictionary
        feature_names: list of feature names (defaults to PH_FEATURE_NAMES)

    Returns:
        numpy array of shape (1, 60)
    """
    if feature_names is None:
        feature_names = PH_FEATURE_NAMES
    X = np.array([features.get(name, np.nan) for name in feature_names]).reshape(1, -1)
    return X

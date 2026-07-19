"""
Amino acid composition (AAC) and dipeptide composition (DPC) feature extraction.

Extracts a 420-dimensional feature vector from a protein sequence:
  - 20 AAC features: normalized frequency of each standard amino acid
  - 400 DPC features: normalized frequency of each dipeptide (AA, AC, ..., YY)
"""

import re
from collections import Counter

# 20 standard amino acids
AMINO_ACIDS = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# All 400 possible dipeptides
DIPEPTIDES = [a + b for a in AMINO_ACIDS for b in AMINO_ACIDS]

AA_SET = set(AMINO_ACIDS)


def read_fasta(filepath):
    """Read a FASTA file, return concatenated sequence of standard AAs only."""
    sequences = []
    current_seq = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_seq:
                    sequences.append(''.join(current_seq))
                    current_seq = []
            else:
                current_seq.append(line.upper())
        if current_seq:
            sequences.append(''.join(current_seq))

    full_seq = ''.join(sequences)
    full_seq = re.sub(f'[^{"".join(AMINO_ACIDS)}]', '', full_seq)
    return full_seq


def extract_aac(sequence):
    """Extract 20-dimensional amino acid composition."""
    total = len(sequence)
    if total == 0:
        return [0.0] * 20
    counter = Counter(sequence)
    return [counter.get(aa, 0) / total for aa in AMINO_ACIDS]


def extract_dpc(sequence):
    """Extract 400-dimensional dipeptide composition."""
    total_dipeptides = len(sequence) - 1
    if total_dipeptides <= 0:
        return [0.0] * 400

    dp_counter = Counter()
    for i in range(len(sequence) - 1):
        dp = sequence[i:i + 2]
        if dp[0] in AA_SET and dp[1] in AA_SET:
            dp_counter[dp] += 1

    return [dp_counter.get(dp, 0) / total_dipeptides for dp in DIPEPTIDES]


def extract_features(sequence):
    """Extract full 420-dimensional vector (20 AAC + 400 DPC)."""
    return extract_aac(sequence) + extract_dpc(sequence)

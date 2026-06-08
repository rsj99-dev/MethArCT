"""
FASTA file parser module.

Supports reading plain text and gzip-compressed FASTA files.
"""

import gzip
import io
from typing import IO, Iterator, Tuple


def iterate_fasta(fasta_file: IO) -> Iterator[Tuple[str, str]]:
    """
    Iterate over a FASTA file, yielding (header, sequence) tuples one at a time.

    Adapted from a classic Biostars implementation:
    https://www.biostars.org/p/710/

    Args:
        fasta_file: an open FASTA file object

    Yields:
        (header, sequence): header is the description line without '>', sequence is the full amino acid sequence
    """
    from itertools import groupby

    faiter = (x[1] for x in groupby(fasta_file, lambda line: line[0] == ">"))
    for header in faiter:
        header_str = header.__next__()[1:].strip()
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (header_str, seq)


def read_fasta(filepath: str) -> dict:
    """
    Read a FASTA file and return a {protein_id: sequence} dictionary.

    Automatically detects file extension and supports .gz compressed format.
    protein_id is taken as the first whitespace-delimited token in the header.

    Args:
        filepath: path to FASTA file (supports .faa / .faa.gz)

    Returns:
        dict: {protein_id: amino_acid_sequence}
    """
    sequences = {}

    if filepath.endswith(".gz"):
        fh = io.TextIOWrapper(io.BufferedReader(gzip.open(filepath, "r")))
    else:
        fh = open(filepath, "r", encoding="utf-8")

    try:
        fh.seek(0)
        for header, sequence in iterate_fasta(fh):
            protein_id = header.split(" ")[0]
            sequences[protein_id] = sequence
    finally:
        fh.close()

    return sequences

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sequence processing utilities for MethArCT
"""

import re
from typing import Dict, List, Tuple, Optional, Union
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
try:
    # New version Biopython (>=1.80) uses gc_fraction
    from Bio.SeqUtils import gc_fraction, molecular_weight
    def GC(sequence):
        """Wrapper for compatibility with old GC function"""
        return 100 * gc_fraction(sequence, ambiguous="ignore")
except ImportError:
    # Old version Biopython (<1.80) uses GC
    from Bio.SeqUtils import GC, molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis

class SequenceUtils:
    """Utility class for sequence operations"""
    
    @staticmethod
    def validate_sequence_type(sequence: str) -> str:
        """
        Determine if sequence is DNA, RNA, or protein
        
        Args:
            sequence: Input sequence string
            
        Returns:
            'dna', 'rna', or 'protein'
        """
        sequence = sequence.upper().replace('-', '').replace('N', '')
        
        if not sequence:
            return 'unknown'
        
        # Count nucleotides and amino acids
        dna_chars = set('ATCG')
        rna_chars = set('AUCG')
        protein_chars = set('ACDEFGHIKLMNPQRSTVWY')
        
        dna_count = sum(1 for c in sequence if c in dna_chars)
        rna_count = sum(1 for c in sequence if c in rna_chars)
        protein_count = sum(1 for c in sequence if c in protein_chars)
        
        total_len = len(sequence)
        
        # If >90% are valid DNA characters and no U
        if dna_count / total_len > 0.9 and 'U' not in sequence:
            return 'dna'
        # If >90% are valid RNA characters and has U
        elif rna_count / total_len > 0.9 and 'U' in sequence:
            return 'rna'
        # If >70% are valid protein characters
        elif protein_count / total_len > 0.7:
            return 'protein'
        else:
            return 'unknown'
    
    @staticmethod
    def clean_sequence(sequence: str, seq_type: str = 'auto') -> str:
        """
        Clean sequence by removing invalid characters
        
        Args:
            sequence: Input sequence
            seq_type: Sequence type ('dna', 'rna', 'protein', 'auto')
            
        Returns:
            Cleaned sequence
        """
        if seq_type == 'auto':
            seq_type = SequenceUtils.validate_sequence_type(sequence)
        
        sequence = sequence.upper().strip()
        
        if seq_type == 'dna':
            # Keep only valid DNA characters
            sequence = re.sub(r'[^ATCGN-]', '', sequence)
        elif seq_type == 'rna':
            # Keep only valid RNA characters
            sequence = re.sub(r'[^AUCGN-]', '', sequence)
        elif seq_type == 'protein':
            # Keep only valid amino acid characters
            sequence = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY*-]', '', sequence)
        
        return sequence
    
    @staticmethod
    def translate_dna(dna_sequence: str, frame: int = 1) -> str:
        """
        Translate DNA sequence to protein
        
        Args:
            dna_sequence: DNA sequence string
            frame: Reading frame (1, 2, 3, -1, -2, -3)
            
        Returns:
            Translated protein sequence
        """
        seq = Seq(dna_sequence.upper())
        
        if frame < 0:
            seq = seq.reverse_complement()
            frame = abs(frame)
        
        # Adjust for 0-based indexing
        start = frame - 1
        
        # Translate
        protein = seq[start:].translate()
        
        return str(protein)
    
    @staticmethod
    def find_orfs(dna_sequence: str, min_length: int = 100) -> List[Dict]:
        """
        Find open reading frames in DNA sequence
        
        Args:
            dna_sequence: DNA sequence string
            min_length: Minimum ORF length in nucleotides
            
        Returns:
            List of ORF dictionaries with start, end, frame, and protein sequence
        """
        orfs = []
        seq = Seq(dna_sequence.upper())
        
        # Check all 6 reading frames
        for strand in [1, -1]:
            if strand == -1:
                seq = seq.reverse_complement()
            
            for frame in range(3):
                start_pos = frame
                
                while start_pos < len(seq) - min_length:
                    # Find start codon
                    start_codon_pos = seq.find('ATG', start_pos)
                    if start_codon_pos == -1:
                        break
                    
                    # Find stop codon
                    for i in range(start_codon_pos + 3, len(seq) - 2, 3):
                        codon = seq[i:i+3]
                        if codon in ['TAA', 'TAG', 'TGA']:
                            orf_length = i + 3 - start_codon_pos
                            if orf_length >= min_length:
                                orf_seq = seq[start_codon_pos:i+3]
                                protein = orf_seq.translate()
                                
                                orfs.append({
                                    'start': start_codon_pos + 1,  # 1-based
                                    'end': i + 3,  # 1-based
                                    'length': orf_length,
                                    'frame': frame + 1 if strand == 1 else -(frame + 1),
                                    'strand': strand,
                                    'dna_sequence': str(orf_seq),
                                    'protein_sequence': str(protein)
                                })
                            break
                    
                    start_pos = start_codon_pos + 3
        
        return orfs
    
    @staticmethod
    def calculate_sequence_stats(sequence: str, seq_type: str = 'auto') -> Dict:
        """
        Calculate basic sequence statistics
        
        Args:
            sequence: Input sequence
            seq_type: Sequence type ('dna', 'rna', 'protein', 'auto')
            
        Returns:
            Dictionary with sequence statistics
        """
        if seq_type == 'auto':
            seq_type = SequenceUtils.validate_sequence_type(sequence)
        
        sequence = SequenceUtils.clean_sequence(sequence, seq_type)
        stats = {
            'length': len(sequence),
            'type': seq_type
        }
        
        if seq_type in ['dna', 'rna']:
            seq_obj = Seq(sequence)
            stats.update({
                'gc_content': round(GC(seq_obj), 2),
                'molecular_weight': round(molecular_weight(seq_obj), 2),
                'composition': {
                    'A': sequence.count('A'),
                    'T': sequence.count('T') if seq_type == 'dna' else 0,
                    'U': sequence.count('U') if seq_type == 'rna' else 0,
                    'C': sequence.count('C'),
                    'G': sequence.count('G'),
                    'N': sequence.count('N')
                }
            })
        
        elif seq_type == 'protein':
            try:
                protein_analysis = ProteinAnalysis(sequence)
                stats.update({
                    'molecular_weight': round(protein_analysis.molecular_weight(), 2),
                    'isoelectric_point': round(protein_analysis.isoelectric_point(), 2),
                    'instability_index': round(protein_analysis.instability_index(), 2),
                    'gravy': round(protein_analysis.gravy(), 2),
                    'composition': protein_analysis.get_amino_acids_percent()
                })
            except Exception:
                stats['composition'] = {}
        
        return stats
    
    @staticmethod
    def search_motifs(sequence: str, motifs: List[str]) -> Dict[str, List[int]]:
        """
        Search for sequence motifs
        
        Args:
            sequence: Input sequence
            motifs: List of motif patterns (can include IUPAC codes)
            
        Returns:
            Dictionary mapping motif to list of positions (1-based)
        """
        results = {}
        sequence = sequence.upper()
        
        # IUPAC nucleotide codes
        iupac_codes = {
            'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
            'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
            'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
        }
        
        for motif in motifs:
            motif = motif.upper()
            
            # Convert IUPAC codes to regex
            regex_pattern = motif
            for code, replacement in iupac_codes.items():
                regex_pattern = regex_pattern.replace(code, replacement)
            
            # Find all matches
            matches = []
            for match in re.finditer(regex_pattern, sequence):
                matches.append(match.start() + 1)  # 1-based position
            
            results[motif] = matches
        
        return results
    
    @staticmethod
    def reverse_complement(dna_sequence: str) -> str:
        """
        Get reverse complement of DNA sequence
        
        Args:
            dna_sequence: DNA sequence string
            
        Returns:
            Reverse complement sequence
        """
        seq = Seq(dna_sequence.upper())
        return str(seq.reverse_complement())
    
    @staticmethod
    def transcribe(dna_sequence: str) -> str:
        """
        Transcribe DNA to RNA
        
        Args:
            dna_sequence: DNA sequence string
            
        Returns:
            RNA sequence
        """
        seq = Seq(dna_sequence.upper())
        return str(seq.transcribe())
    
    @staticmethod
    def back_transcribe(rna_sequence: str) -> str:
        """
        Back-transcribe RNA to DNA
        
        Args:
            rna_sequence: RNA sequence string
            
        Returns:
            DNA sequence
        """
        seq = Seq(rna_sequence.upper())
        return str(seq.back_transcribe())
    
    @staticmethod
    def format_sequence(sequence: str, line_length: int = 80) -> str:
        """
        Format sequence with line breaks
        
        Args:
            sequence: Input sequence
            line_length: Characters per line
            
        Returns:
            Formatted sequence string
        """
        lines = []
        for i in range(0, len(sequence), line_length):
            lines.append(sequence[i:i+line_length])
        return '\n'.join(lines)
    
    @staticmethod
    def create_sequence_record(seq_id: str, sequence: str, description: str = "") -> SeqRecord:
        """
        Create BioPython SeqRecord object
        
        Args:
            seq_id: Sequence identifier
            sequence: Sequence string
            description: Sequence description
            
        Returns:
            SeqRecord object
        """
        return SeqRecord(
            Seq(sequence),
            id=seq_id,
            description=description
        )
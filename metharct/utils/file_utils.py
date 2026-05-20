#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File handling utilities for MethArCT
"""

import os
import shutil
import tempfile
from pathlib import Path
from typing import List, Optional, Union, Iterator
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class FileUtils:
    """Utility class for file operations"""
    
    @staticmethod
    def ensure_dir(path: Union[str, Path]) -> Path:
        """
        Ensure directory exists, create if not
        
        Args:
            path: Directory path
            
        Returns:
            Path object
        """
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)
        return path
    
    @staticmethod
    def validate_fasta(file_path: Union[str, Path]) -> bool:
        """
        Validate if file is a valid FASTA format
        
        Args:
            file_path: Path to FASTA file
            
        Returns:
            True if valid FASTA, False otherwise
        """
        try:
            with open(file_path, 'r') as handle:
                records = SeqIO.parse(handle, "fasta")
                # Try to read at least one record
                next(records)
                return True
        except (StopIteration, Exception):
            return False
    
    @staticmethod
    def count_sequences(file_path: Union[str, Path]) -> int:
        """
        Count number of sequences in FASTA file
        
        Args:
            file_path: Path to FASTA file
            
        Returns:
            Number of sequences
        """
        try:
            with open(file_path, 'r') as handle:
                return sum(1 for _ in SeqIO.parse(handle, "fasta"))
        except Exception:
            return 0
    
    @staticmethod
    def get_sequence_ids(file_path: Union[str, Path]) -> List[str]:
        """
        Get list of sequence IDs from FASTA file
        
        Args:
            file_path: Path to FASTA file
            
        Returns:
            List of sequence IDs
        """
        ids = []
        try:
            with open(file_path, 'r') as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    ids.append(record.id)
        except Exception:
            pass
        return ids
    
    @staticmethod
    def filter_sequences_by_length(input_file: Union[str, Path],
                                 output_file: Union[str, Path],
                                 min_length: int = 50,
                                 max_length: Optional[int] = None) -> int:
        """
        Filter sequences by length
        
        Args:
            input_file: Input FASTA file
            output_file: Output FASTA file
            min_length: Minimum sequence length
            max_length: Maximum sequence length (optional)
            
        Returns:
            Number of sequences written
        """
        count = 0
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for record in SeqIO.parse(infile, "fasta"):
                seq_len = len(record.seq)
                if seq_len >= min_length:
                    if max_length is None or seq_len <= max_length:
                        SeqIO.write(record, outfile, "fasta")
                        count += 1
        return count
    
    @staticmethod
    def merge_fasta_files(input_files: List[Union[str, Path]],
                         output_file: Union[str, Path],
                         add_source_prefix: bool = False) -> int:
        """
        Merge multiple FASTA files into one
        
        Args:
            input_files: List of input FASTA files
            output_file: Output merged FASTA file
            add_source_prefix: Whether to add source filename as prefix to sequence IDs
            
        Returns:
            Total number of sequences merged
        """
        total_count = 0
        with open(output_file, 'w') as outfile:
            for input_file in input_files:
                if not os.path.exists(input_file):
                    continue
                    
                source_name = Path(input_file).stem if add_source_prefix else None
                
                with open(input_file, 'r') as infile:
                    for record in SeqIO.parse(infile, "fasta"):
                        if add_source_prefix and source_name:
                            record.id = f"{source_name}_{record.id}"
                            record.description = f"{source_name}_{record.description}"
                        
                        SeqIO.write(record, outfile, "fasta")
                        total_count += 1
        
        return total_count
    
    @staticmethod
    def split_fasta_file(input_file: Union[str, Path],
                        output_dir: Union[str, Path],
                        sequences_per_file: int = 1000) -> List[Path]:
        """
        Split large FASTA file into smaller files
        
        Args:
            input_file: Input FASTA file
            output_dir: Output directory
            sequences_per_file: Number of sequences per output file
            
        Returns:
            List of output file paths
        """
        output_dir = Path(output_dir)
        FileUtils.ensure_dir(output_dir)
        
        input_path = Path(input_file)
        base_name = input_path.stem
        
        output_files = []
        file_count = 0
        seq_count = 0
        current_file = None
        current_handle = None
        
        try:
            with open(input_file, 'r') as infile:
                for record in SeqIO.parse(infile, "fasta"):
                    if seq_count % sequences_per_file == 0:
                        if current_handle:
                            current_handle.close()
                        
                        file_count += 1
                        current_file = output_dir / f"{base_name}_part{file_count:03d}.fasta"
                        current_handle = open(current_file, 'w')
                        output_files.append(current_file)
                    
                    SeqIO.write(record, current_handle, "fasta")
                    seq_count += 1
        
        finally:
            if current_handle:
                current_handle.close()
        
        return output_files
    
    @staticmethod
    def create_temp_file(suffix: str = ".tmp", prefix: str = "metharct_") -> str:
        """
        Create temporary file
        
        Args:
            suffix: File suffix
            prefix: File prefix
            
        Returns:
            Path to temporary file
        """
        fd, path = tempfile.mkstemp(suffix=suffix, prefix=prefix)
        os.close(fd)  # Close file descriptor, but keep the file
        return path
    
    @staticmethod
    def create_temp_dir(prefix: str = "metharct_") -> str:
        """
        Create temporary directory
        
        Args:
            prefix: Directory prefix
            
        Returns:
            Path to temporary directory
        """
        return tempfile.mkdtemp(prefix=prefix)
    
    @staticmethod
    def cleanup_temp_files(paths: List[str]):
        """
        Clean up temporary files and directories
        
        Args:
            paths: List of paths to clean up
        """
        for path in paths:
            try:
                if os.path.isfile(path):
                    os.unlink(path)
                elif os.path.isdir(path):
                    shutil.rmtree(path)
            except Exception:
                pass  # Ignore cleanup errors
    
    @staticmethod
    def get_file_size(file_path: Union[str, Path]) -> int:
        """
        Get file size in bytes
        
        Args:
            file_path: Path to file
            
        Returns:
            File size in bytes
        """
        try:
            return os.path.getsize(file_path)
        except Exception:
            return 0
    
    @staticmethod
    def format_file_size(size_bytes: int) -> str:
        """
        Format file size in human readable format
        
        Args:
            size_bytes: Size in bytes
            
        Returns:
            Formatted size string
        """
        if size_bytes == 0:
            return "0 B"
        
        size_names = ["B", "KB", "MB", "GB", "TB"]
        i = 0
        while size_bytes >= 1024 and i < len(size_names) - 1:
            size_bytes /= 1024.0
            i += 1
        
        return f"{size_bytes:.1f} {size_names[i]}"
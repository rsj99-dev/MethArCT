#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Diamond sequence alignment analyzer for MethArCT

Performs protein sequence alignment using Diamond BLAST for metabolic pathway prediction
and cultivability evaluation.
"""

import os
import shlex
import subprocess
import tempfile
import shutil
import json
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from ..utils.config import Config
from ..utils.logger import get_logger
from ..utils.file_utils import FileUtils
from ..utils.sequence_utils import SequenceUtils
from .pathway_strategies import evaluate_pathway

# ============================================================
# Module-level constants: database category classifications
# ============================================================
METHANE_DATABASES = frozenset([
    'CO2-CH4', 'Methylamine', 'Methanol', 'Methanethiol',
    'Acetate', 'Alkane', 'CO-CH4', 'Formate',
    'Methoxy', 'Fatty_acid', 'Dimethylamine', 'Trimethylamine',
    'Glycine betaine methanogenesis', 'Methylthiopropionate methanogenesis',
    'Tetramethylammonium methanogenesis', 'Methanol dismutation methanogenesis',
])

SULFUR_DATABASES = frozenset(['ASR', 'SO', 'SOX', 'S4I', 'SR', 'DSR'])

NITROGEN_DATABASES = frozenset(['ANR', 'DEN', 'DNR', 'NIT', 'AMX', 'CNIT', 'NFX1', 'NFX2'])

CULTIVATION_DATABASE = 'CULTIVATION'


class DiamondAnalyzer:
    """Diamond sequence alignment analyzer"""
    
    def __init__(self, config: Optional[Config] = None):
        self.config = config or Config()
        self.logger = get_logger("diamond_analyzer")
        
        # Tool configuration - WSL support
        self.use_wsl = self.config.get('tools.diamond.use_wsl', False)

        # Always keep the command as an argument list.  Joining "wsl diamond"
        # into a single string makes subprocess look for an executable literally
        # named "wsl diamond".
        if self.use_wsl:
            wsl_command = self.config.get('tools.diamond.wsl_path', 'wsl diamond')
            self.diamond_cmd = shlex.split(wsl_command) or ['wsl', 'diamond']
        else:
            self.diamond_cmd = [self.config.get('tools.diamond.path', 'diamond')]
        self.diamond_path = ' '.join(self.diamond_cmd)
            
        self.threads = self.config.get('tools.diamond.threads', 4)
        self.evalue = self.config.get('tools.diamond.evalue', 1e-5)
        self.max_target_seqs = self.config.get('tools.diamond.max_target_seqs', 1)
        self.identity_threshold = self.config.get('tools.diamond.identity_threshold', 30.0)
        
        # Dual threshold analysis parameters
        # HIGH_E_VALUE_THRESHOLD is now dynamic based on reference sequence length
        # Fallback E-value when average length cannot be determined:
        self._FALLBACK_HIGH_E_VALUE_THRESHOLD = 1e-100
        self.LOW_E_VALUE_THRESHOLD = 1e-5
        self.LOW_BITSCORE_THRESHOLD = 80
        self.HIGH_QUALITY_THRESHOLD = 60
        
        # Cultivability assessment thresholds
        self.CULTIVABILITY_E_VALUE_THRESHOLD = 1e-5
        self.CULTIVABILITY_BITSCORE_THRESHOLD = 50
        
        # Database paths
        self.db_base_dir = self.config.get('databases.base_dir', 'data/databases')
        self.db_dir = os.path.join(self.db_base_dir, 'diamond')
        
        # Results directory
        self.results_dir = self.config.get('output.base_dir', 'results')
        FileUtils.ensure_dir(self.results_dir)
        FileUtils.ensure_dir(self.db_dir)
        
        # Database information
        self.pathway_names = self.config.get('pathway_names', {})
        self.reference_counts = self.config.get('reference_sequence_counts', {})
        self._reference_avg_lengths = {}  # db_name -> average amino acid length
        self._reference_subject_counts = {}  # db_name -> unique subject IDs
        
        # Check tool availability
        self._check_diamond_availability()
        
        # Database preparation is deferred to analyze_sequence() for lazy initialization
        self._databases_prepared = False
        self.database_preparation_failures = []
    
    def _ensure_databases(self):
        """Ensure databases are prepared before analysis (lazy initialization)."""
        if not self._databases_prepared:
            self._prepare_databases()
            self._databases_prepared = True
    
    def _check_diamond_availability(self) -> bool:
        """
        Check if Diamond tool is available - WSL environment support
        
        Returns:
            True if Diamond is available, False otherwise
        """
        try:
            result = subprocess.run(
                self.diamond_cmd + ['version'],
                capture_output=True,
                text=True,
                timeout=30
            )
            
            if result.returncode == 0:
                version = result.stdout.strip()
                env_type = "WSL" if self.use_wsl else "Local"
                self.logger.info(f"Diamond tool available ({env_type}): {version}")
                return True
            else:
                self.logger.error(f"Diamond tool not available: {result.stderr}")
                return False
                
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            self.logger.error(f"Error checking Diamond availability: {str(e)}")
            return False

    def _host_path(self, path) -> str:
        """Convert a local path to the format expected by the DIAMOND process."""
        if self.use_wsl:
            return FileUtils.to_wsl_path(path)
        return str(path)
    
    def _prepare_databases(self):
        """
        Prepare Diamond databases from reference sequences
        """
        self.logger.info("Preparing Diamond databases...")
        failed_databases = []
        
        # Get all database paths from config
        db_paths = self.config.get_all_database_paths()
        
        for db_name, fasta_path in db_paths.items():
            # Cultivation is a directory of per-pathway FASTA files, not a
            # single searchable database. Skip any non-file entry.
            if os.path.isdir(fasta_path):
                self.logger.debug(
                    f"Skipping directory database entry: {db_name} -> {fasta_path}"
                )
                continue

            if not os.path.exists(fasta_path):
                self.logger.warning(f"Reference file not found: {fasta_path}")
                continue

            # Compute average sequence length for dynamic threshold (all databases)
            avg_len = self._compute_avg_sequence_length(fasta_path)
            if avg_len > 0:
                self._reference_avg_lengths[db_name] = avg_len
                self.logger.debug(
                    f"{db_name}: avg sequence length = {avg_len:.1f} aa"
                )

            # Record the actual number of unique reference subjects so
            # sulfur/nitrogen completeness cannot exceed 100 % while a
            # retrievable subject is still missing.
            subject_count = self._count_unique_subjects(fasta_path)
            if subject_count > 0:
                self._reference_subject_counts[db_name] = subject_count
                configured = self.reference_counts.get(db_name)
                if (db_name not in METHANE_DATABASES
                        and configured is not None
                        and int(configured) != subject_count):
                    self.logger.warning(
                        f"{db_name}: configured reference_sequence_counts="
                        f"{configured} differs from {subject_count} unique FASTA "
                        f"subjects; using {subject_count} for this run"
                    )
            
            db_file = os.path.join(self.db_dir, f"{db_name}.dmnd")
            
            # Check if database exists and is up to date.  Archives extracted
            # from zip often carry identical timestamps, so use >= instead of >.
            if (os.path.exists(db_file) and
                os.path.getmtime(db_file) >= os.path.getmtime(fasta_path)):
                self.logger.debug(f"Database {db_name} is up to date")
                continue
            
            # Create or update database
            if not self._create_diamond_database(fasta_path, db_name):
                failed_databases.append(db_name)

        self.database_preparation_failures = failed_databases
        if failed_databases:
            self.logger.error(
                "Diamond database preparation failed for: "
                + ", ".join(failed_databases)
            )

    @staticmethod
    def _count_unique_subjects(fasta_path: str) -> int:
        """Count unique subject IDs (first token of each FASTA header)."""
        subjects = set()
        try:
            with open(fasta_path, 'r', encoding='utf-8', errors='ignore') as handle:
                for line in handle:
                    if line.startswith('>'):
                        header = line[1:].strip()
                        if header:
                            subjects.add(header.split()[0])
        except Exception:
            return 0
        return len(subjects)
    
    def _compute_avg_sequence_length(self, fasta_path: str) -> float:
        """
        Compute average amino acid sequence length from a FASTA file.
        
        Args:
            fasta_path: Path to the FASTA file.
            
        Returns:
            Average sequence length, or 0 if the file cannot be parsed.
        """
        try:
            lengths = []
            with open(fasta_path, 'r') as f:
                current_len = 0
                in_seq = False
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if in_seq and current_len > 0:
                            lengths.append(current_len)
                        current_len = 0
                        in_seq = True
                    elif in_seq:
                        current_len += len(line)
                if in_seq and current_len > 0:
                    lengths.append(current_len)
            return sum(lengths) / len(lengths) if lengths else 0.0
        except Exception as e:
            self.logger.warning(f"Failed to compute avg length for {fasta_path}: {e}")
            return 0.0
    
    def _get_dynamic_high_evalue(self, db_name: str) -> float:
        """
        Get dynamic high-threshold E-value based on average reference sequence length.
        
        Mapping (average amino acid length → log10(E-value) threshold):
            >550  → 1e-140
            >500  → 1e-100
            >200  → 1e-25
            >100  → 1e-10
            >50   → 1e-5
        
        Args:
            db_name: Database name.
            
        Returns:
            Dynamic E-value threshold for the high-strictness tier.
        """
        avg_len = self._reference_avg_lengths.get(db_name, 0)
        if avg_len > 550:
            return 1e-140
        elif avg_len > 500:
            return 1e-100
        elif avg_len > 200:
            return 1e-25
        elif avg_len > 100:
            return 1e-10
        elif avg_len > 50:
            return 1e-5
        else:
            # Fallback to original fixed threshold when length is unknown
            self.logger.warning(
                f"{db_name}: avg length {avg_len:.0f} too short or unknown, "
                f"using fallback high E-value threshold"
            )
            return self._FALLBACK_HIGH_E_VALUE_THRESHOLD

    def _create_diamond_database(self, fasta_path: str, db_name: str) -> bool:
        """
        Create Diamond database from FASTA file
        
        Args:
            fasta_path: Path to input FASTA file
            db_name: Database name
            
        Returns:
            True if successful, False otherwise
        """
        try:
            db_path = os.path.join(self.db_dir, db_name)
            
            cmd = self.diamond_cmd + [
                'makedb',
                '--in', self._host_path(fasta_path),
                '--db', self._host_path(db_path),
                '--threads', str(self.threads)
            ]
            
            self.logger.info(f"Creating Diamond database: {db_name}")
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300  # 5 minutes timeout
            )
            
            if result.returncode == 0:
                self.logger.info(f"Successfully created database: {db_name}")
                return True
            else:
                self.logger.error(f"Failed to create database {db_name}: {result.stderr}")
                return False
                
        except subprocess.TimeoutExpired:
            self.logger.error(f"Timeout creating database: {db_name}")
            return False
        except Exception as e:
            self.logger.error(f"Error creating database {db_name}: {str(e)}")
            return False
    
    def analyze_sequence(self, 
                        input_file: Union[str, Path],
                        output_prefix: Optional[str] = None,
                        databases: Optional[List[str]] = None) -> Dict[str, any]:
        """
        Analyze protein sequences against reference databases
        
        Args:
            input_file: Path to input FASTA file
            output_prefix: Prefix for output files
            databases: List of databases to search against (None for all)
            
        Returns:
            Analysis results dictionary
        """
        input_file = Path(input_file)
        
        if not input_file.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")
        
        if not FileUtils.validate_fasta(input_file):
            raise ValueError(f"Invalid FASTA file: {input_file}")
        
        # Set output prefix
        if output_prefix is None:
            output_prefix = input_file.stem
        
        # Get databases to search
        if databases is None:
            databases = [
                name for name, path in self.config.get_all_database_paths().items()
                if os.path.isfile(path)
            ]
        
        self.logger.info(f"Starting Diamond analysis for {input_file.name}")
        self.logger.info(f"Searching against {len(databases)} databases")
        
        # Ensure databases are prepared
        self._ensure_databases()
        
        # Run Diamond searches
        search_results = {}
        for db_name in databases:
            result = self._run_diamond_search(input_file, db_name, output_prefix)
            if result:
                search_results[db_name] = result
        
        # Process and summarize results
        analysis_results = self._process_results(search_results, input_file, databases)
        
        # Save results
        self._save_results(analysis_results, output_prefix)
        
        return analysis_results
    
    @staticmethod
    def _add_header_to_tsv(tsv_path: str) -> None:
        """
        Prepend a header line to a Diamond outfmt 6 TSV file.
        
        The 12 standard BLAST tabular columns are:
        qseqid, sseqid, pident, length, mismatch, gapopen,
        qstart, qend, sstart, send, evalue, bitscore
        
        Args:
            tsv_path: Path to the TSV file to modify in-place.
        """
        header_line = (
            "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\t"
            "qstart\tqend\tsstart\tsend\tevalue\tbitscore\n"
        )
        try:
            with open(tsv_path, 'r') as f:
                original_content = f.read()
            # Only add header if not already present
            if not original_content.startswith('qseqid\t'):
                with open(tsv_path, 'w') as f:
                    f.write(header_line)
                    f.write(original_content)
        except Exception:
            pass  # Silently skip if file cannot be modified

    def _run_diamond_search(self, 
                           input_file: Path, 
                           db_name: str, 
                           output_prefix: str) -> Optional[Dict]:
        """
        Run Diamond search against a specific database with dual threshold analysis
        
        Args:
            input_file: Input FASTA file
            db_name: Database name
            output_prefix: Output file prefix
            
        Returns:
            Search results dictionary or None if failed
        """
        db_path = os.path.join(self.db_dir, db_name)
        
        if not os.path.exists(f"{db_path}.dmnd"):
            self.logger.warning(f"Database not found: {db_path}.dmnd")
            return None
        
        # Create output file in the same directory as output_prefix
        output_dir = os.path.dirname(output_prefix) if os.path.dirname(output_prefix) else self.results_dir
        output_file = os.path.join(
            output_dir, 
            f"{os.path.basename(output_prefix)}_{db_name}_diamond.tsv"
        )
        
        # Select E-value threshold based on database type
        if db_name == CULTIVATION_DATABASE:
            search_evalue = self.CULTIVABILITY_E_VALUE_THRESHOLD
        else:
            search_evalue = self.LOW_E_VALUE_THRESHOLD
        
        cmd = self.diamond_cmd + [
            'blastp',
            '--query', self._host_path(input_file),
            '--db', self._host_path(db_path),
            '--out', self._host_path(output_file),
            '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 
                       'mismatch', 'gapopen', 'qstart', 'qend', 
                       'sstart', 'send', 'evalue', 'bitscore',
            '--evalue', str(search_evalue),
            '--max-target-seqs', str(self.max_target_seqs),
            '--threads', str(self.threads),
            '--sensitive'
        ]
        
        try:
            self.logger.debug(f"Running Diamond search: {db_name}")
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600  # 10 minutes timeout
            )
            
            if result.returncode == 0:
                # Parse results
                if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                    # Add header line to the raw TSV file for readability
                    self._add_header_to_tsv(output_file)
                    
                    df = pd.read_csv(
                        output_file,
                        sep='\t',
                        header=0,
                        names=['qseqid', 'sseqid', 'pident', 'length',
                               'mismatch', 'gapopen', 'qstart', 'qend',
                               'sstart', 'send', 'evalue', 'bitscore']
                    )
                    
                    # Return different result structures based on database type
                    if db_name == CULTIVATION_DATABASE:
                        # Cultivability database only performs cultivability filtering
                        cultivability_hits = df[
                            (df['evalue'] <= self.CULTIVABILITY_E_VALUE_THRESHOLD) & 
                            (df['bitscore'] >= self.CULTIVABILITY_BITSCORE_THRESHOLD)
                        ]
                        
                        return {
                            'database': db_name,
                            'total_hits': len(df),
                            'cultivability_hits': cultivability_hits,
                            'output_file': output_file
                        }
                    
                    else:
                        # Metabolic pathway database performs dual threshold filtering
                        # Check if it's a methane metabolism pathway, if so use single threshold analysis
                        if db_name in METHANE_DATABASES:
                            # Methane metabolism pathway: single threshold analysis (only keep low threshold)
                            low_threshold_hits = df[
                                (df['evalue'] <= self.LOW_E_VALUE_THRESHOLD) & 
                                (df['bitscore'] >= self.LOW_BITSCORE_THRESHOLD)
                            ]
                            
                            return {
                                'database': db_name,
                                'total_hits': len(df),
                                'high_threshold_hits': pd.DataFrame(),
                                'low_threshold_hits': low_threshold_hits,
                                'output_file': output_file
                            }
                        else:
                            # Other metabolic pathways (sulfur, nitrogen) maintain dual threshold analysis
                            # High threshold E-value is dynamic based on avg reference sequence length
                            dynamic_high_evalue = self._get_dynamic_high_evalue(db_name)
                            self.logger.debug(
                                f"{db_name}: avg_len={self._reference_avg_lengths.get(db_name, 0):.0f}, "
                                f"dynamic high E-value={dynamic_high_evalue:.2e}"
                            )
                            high_threshold_hits = df[
                                df['evalue'] <= dynamic_high_evalue
                            ]
                            
                            low_threshold_hits = df[
                                (df['evalue'] <= self.LOW_E_VALUE_THRESHOLD) & 
                                (df['bitscore'] >= self.LOW_BITSCORE_THRESHOLD)
                            ]
                            
                            return {
                                'database': db_name,
                                'total_hits': len(df),
                                'high_threshold_hits': high_threshold_hits,
                                'low_threshold_hits': low_threshold_hits,
                                'output_file': output_file
                            }
                        
                else:
                    self.logger.info(f"No hits found for database: {db_name}")
                    # CULTIVATION uses a distinct field name; all other pathways share the unified structure
                    if db_name == CULTIVATION_DATABASE:
                        return {
                            'database': db_name,
                            'total_hits': 0,
                            'cultivability_hits': pd.DataFrame(),
                            'output_file': output_file
                        }
                    return {
                        'database': db_name,
                        'total_hits': 0,
                        'high_threshold_hits': pd.DataFrame(),
                        'low_threshold_hits': pd.DataFrame(),
                        'output_file': output_file
                    }
            else:
                self.logger.error(f"Diamond search failed for {db_name}: {result.stderr}")
                return None
                
        except subprocess.TimeoutExpired:
            self.logger.error(f"Diamond search timeout for database: {db_name}")
            return None
        except Exception as e:
            self.logger.error(f"Error in Diamond search for {db_name}: {str(e)}")
            return None
    
    def _get_reference_count(self, db_name: str) -> int:
        """Return the completeness denominator for a database.

        Methane pathways use the configured number of functional units because
        alternative genes are merged by ``PATHWAY_STRATEGIES``. Sulfur and
        nitrogen pathways use the number of unique reference subjects actually
        present in the FASTA file, so 100 % cannot be reported while a
        retrievable subject is still missing.
        """
        if db_name in METHANE_DATABASES:
            return int(self.reference_counts.get(db_name, 1)) or 1

        subject_count = self._reference_subject_counts.get(db_name)
        if subject_count:
            return int(subject_count)
        return int(self.reference_counts.get(db_name, 1)) or 1

    def _empty_pathway_result(self, db_name: str, reference_count: int,
                              search_status: str) -> Dict:
        """Build a zero-hit result so every searched database is represented."""
        entry = {
            'pathway_name': self.pathway_names.get(db_name, db_name),
            'reference_sequences': reference_count,
            'search_status': search_status,
        }

        if db_name == CULTIVATION_DATABASE:
            entry.update({
                'cultivability_hits': 0,
                'cultivability_percentage': 0.0,
                'cultivability_status': 'Not detected',
                'average_identity': 0,
                'best_hits': [],
            })
        elif db_name in METHANE_DATABASES:
            entry.update({
                'low_threshold_hits': 0,
                'low_completeness_percentage': 0.0,
                'average_identity_low': 0,
                'best_hits_low': [],
            })
        else:
            entry.update({
                'high_threshold_hits': 0,
                'low_threshold_hits': 0,
                'high_completeness_percentage': 0.0,
                'low_completeness_percentage': 0.0,
                'metabolic_pathway_completeness': '0.00%~0.00%',
                'average_identity_high': 0,
                'average_identity_low': 0,
                'best_hits_high': [],
                'best_hits_low': [],
            })

        return entry

    def _build_pathway_result(self, db_name: str, result: Dict,
                              reference_count: int, total_sequences: int) -> Dict:
        """Build a pathway result from a DIAMOND search result."""
        if db_name == CULTIVATION_DATABASE:
            # Calculate how many different query sequences matched
            # cultivability-related genes (deduplicated)
            cultivability_hits = result['cultivability_hits']
            cultivability_matched = len(cultivability_hits['qseqid'].drop_duplicates())
            cultivability_percentage = (
                (cultivability_matched / total_sequences) * 100
                if total_sequences else 0.0
            )

            if cultivability_percentage <= 11:
                cultivability_status = "Cultivable"
            elif cultivability_percentage <= 13:
                cultivability_status = "Co-cultivation dependent"
            else:
                cultivability_status = "Obligate symbiont"

            return {
                'pathway_name': self.pathway_names.get(db_name, db_name),
                'reference_sequences': reference_count,
                'search_status': 'searched',
                'cultivability_hits': cultivability_matched,
                'cultivability_percentage': round(cultivability_percentage, 2),
                'cultivability_status': cultivability_status,
                'average_identity': round(
                    cultivability_hits['pident'].mean(), 2
                ) if len(cultivability_hits) > 0 else 0,
                'best_hits': cultivability_hits.nlargest(5, 'bitscore').to_dict('records')
            }

        # Metabolic pathway database processing (dual threshold analysis)
        if db_name in METHANE_DATABASES:
            # Methane metabolism pathway: single threshold analysis (only keep
            # low threshold). Homologous genes (e.g., FwdF/FwdG/FwdH, Hmd/Mtd)
            # have similar sequences but represent different gene functions and
            # should be counted separately.
            low_threshold_matched = evaluate_pathway(
                db_name, result['low_threshold_hits'], reference_count, self.logger
            )
            low_threshold_matched = min(low_threshold_matched, reference_count)
            low_completeness = (low_threshold_matched / reference_count) * 100

            return {
                'pathway_name': self.pathway_names.get(db_name, db_name),
                'low_threshold_hits': low_threshold_matched,
                'reference_sequences': reference_count,
                'search_status': 'searched',
                'low_completeness_percentage': round(low_completeness, 2),
                'average_identity_low': round(
                    result['low_threshold_hits']['pident'].mean(), 2
                ) if len(result['low_threshold_hits']) > 0 else 0,
                'best_hits_low': result['low_threshold_hits'].nlargest(5, 'bitscore').to_dict('records')
            }

        # Sulfur and nitrogen pathways keep dual threshold analysis and count
        # unique reference subjects. The denominator is the number of unique
        # FASTA subjects (see _get_reference_count), never a smaller configured
        # value.
        high_threshold_matched = len(
            result['high_threshold_hits'].drop_duplicates(subset=['sseqid'])
        )
        low_threshold_matched = len(
            result['low_threshold_hits'].drop_duplicates(subset=['sseqid'])
        )
        high_threshold_matched = min(high_threshold_matched, reference_count)
        low_threshold_matched = min(low_threshold_matched, reference_count)

        high_completeness = (high_threshold_matched / reference_count) * 100
        low_completeness = (low_threshold_matched / reference_count) * 100
        metabolic_pathway_completeness = (
            f"{high_completeness:.2f}%~{low_completeness:.2f}%"
        )

        return {
            'pathway_name': self.pathway_names.get(db_name, db_name),
            'high_threshold_hits': high_threshold_matched,
            'low_threshold_hits': low_threshold_matched,
            'reference_sequences': reference_count,
            'search_status': 'searched',
            'high_completeness_percentage': round(high_completeness, 2),
            'low_completeness_percentage': round(low_completeness, 2),
            'metabolic_pathway_completeness': metabolic_pathway_completeness,
            'average_identity_high': round(
                result['high_threshold_hits']['pident'].mean(), 2
            ) if len(result['high_threshold_hits']) > 0 else 0,
            'average_identity_low': round(
                result['low_threshold_hits']['pident'].mean(), 2
            ) if len(result['low_threshold_hits']) > 0 else 0,
            'best_hits_high': result['high_threshold_hits'].nlargest(5, 'bitscore').to_dict('records'),
            'best_hits_low': result['low_threshold_hits'].nlargest(5, 'bitscore').to_dict('records')
        }

    def _process_results(self, search_results: Dict, input_file: Path,
                         databases: Optional[List[str]] = None) -> Dict:
        """
        Process and summarize Diamond search results with dual threshold analysis
        
        Args:
            search_results: Raw search results (missing keys mean the search
                did not run or failed)
            input_file: Input FASTA file
            databases: Databases requested for this analysis
            
        Returns:
            Processed analysis results
        """
        # Count input sequences
        total_sequences = FileUtils.count_sequences(input_file)
        requested = list(databases) if databases is not None else list(search_results.keys())

        pathway_results = {}
        summary_stats = {
            'total_input_sequences': total_sequences,
            'databases_requested': len(requested),
            'databases_searched': sum(1 for d in requested if d in search_results),
            'raw_hit_pathways': 0,
            'threshold_passed_pathways': 0,
            'pathways_detected': 0,
            'complete_pathways': 0,
            'total_hits': 0
        }

        for db_name in requested:
            result = search_results.get(db_name)
            reference_count = self._get_reference_count(db_name)

            if result is None or result.get('total_hits', 0) == 0:
                # Emit a zero-hit entry so downstream summaries use the full
                # configured denominator instead of silently dropping pathways.
                pathway_results[db_name] = self._empty_pathway_result(
                    db_name,
                    reference_count,
                    search_status='searched' if result is not None else 'not_searched'
                )
                continue

            summary_stats['raw_hit_pathways'] += 1
            entry = self._build_pathway_result(
                db_name, result, reference_count, total_sequences
            )
            pathway_results[db_name] = entry

            summary_stats['total_hits'] += len(result.get('low_threshold_hits', []))

            if db_name == CULTIVATION_DATABASE:
                passed = entry['cultivability_hits'] > 0
            else:
                passed = entry.get('low_threshold_hits', 0) > 0

            if passed:
                summary_stats['threshold_passed_pathways'] += 1
                summary_stats['pathways_detected'] += 1

            if entry.get('low_completeness_percentage', 0) >= 100:
                summary_stats['complete_pathways'] += 1

        return {
            'input_file': str(input_file),
            'analysis_timestamp': pd.Timestamp.now().isoformat(),
            'summary': summary_stats,
            'pathway_results': pathway_results,
            'raw_results': search_results
        }
    
    def _save_results(self, results: Dict, output_prefix: str):
        """
        Save analysis results to files with simplified output
        
        Args:
            results: Analysis results
            output_prefix: Output file prefix
        """
        # Create serializable result copy, convert DataFrame to list of dictionaries
        serializable_results = results.copy()
        
        # Process DataFrame objects in pathway_results
        for db_name, pathway_data in serializable_results['pathway_results'].items():
            # Convert DataFrame in best_hits to list of dictionaries
            if 'best_hits_high' in pathway_data and hasattr(pathway_data['best_hits_high'], 'to_dict'):
                pathway_data['best_hits_high'] = pathway_data['best_hits_high'].to_dict('records')
            if 'best_hits_low' in pathway_data and hasattr(pathway_data['best_hits_low'], 'to_dict'):
                pathway_data['best_hits_low'] = pathway_data['best_hits_low'].to_dict('records')
        
        # Process DataFrame objects in raw_results
        if 'raw_results' in serializable_results:
            for db_name, raw_data in serializable_results['raw_results'].items():
                if raw_data:
                    # Convert DataFrame to list of dictionaries
                    for key in ['high_threshold_hits', 'low_threshold_hits', 'cultivability_hits']:
                        if key in raw_data and hasattr(raw_data[key], 'to_dict'):
                            raw_data[key] = raw_data[key].to_dict('records')
        
        # Save only one summary CSV file, remove duplicate JSON and CSV file output
        output_dir = os.path.dirname(output_prefix) if os.path.dirname(output_prefix) else self.results_dir
        summary_csv_file = os.path.join(output_dir, f"{os.path.basename(output_prefix)}_metabolic_summary.csv")
        
        # Prepare summary data
        summary_data = []
        
        # Get input file information
        input_file = results.get('input_file', 'Unknown')
        analysis_timestamp = results.get('analysis_timestamp', '')
        
        # Process results for each metabolic pathway
        for db_name, pathway_data in results['pathway_results'].items():
            # Determine metabolic pathway type
            pathway_type = 'Other'
            if db_name in METHANE_DATABASES:
                pathway_type = 'Methane'
            elif db_name in SULFUR_DATABASES:
                pathway_type = 'Sulfur'
            elif db_name in NITROGEN_DATABASES:
                pathway_type = 'Nitrogen'
            elif db_name == CULTIVATION_DATABASE:
                pathway_type = 'Cultivation'
            
            # Extract corresponding data based on metabolic pathway type
            if pathway_type == 'Methane':
                # Methane metabolism pathway
                summary_data.append({
                    'Input_File': input_file,
                    'Analysis_Timestamp': analysis_timestamp,
                    'Pathway_Type': pathway_type,
                    'Database': db_name,
                    'Pathway_Name': pathway_data['pathway_name'],
                    'Reference_Sequences': pathway_data['reference_sequences'],
                    'Low_Threshold_Hits': pathway_data.get('low_threshold_hits', 0),
                    'Low_Completeness_Percentage': pathway_data.get('low_completeness_percentage', 0),
                    'Average_Identity_Low': pathway_data.get('average_identity_low', 0),
                    'Search_Status': pathway_data.get('search_status', ''),
                    'Detection_Status': 'Detected' if pathway_data.get('low_threshold_hits', 0) > 0 else 'Not Detected'
                })
            elif pathway_type in ['Sulfur', 'Nitrogen']:
                # Sulfur and nitrogen metabolism pathways
                summary_data.append({
                    'Input_File': input_file,
                    'Analysis_Timestamp': analysis_timestamp,
                    'Pathway_Type': pathway_type,
                    'Database': db_name,
                    'Pathway_Name': pathway_data['pathway_name'],
                    'Reference_Sequences': pathway_data['reference_sequences'],
                    'High_Threshold_Hits': pathway_data.get('high_threshold_hits', 0),
                    'Low_Threshold_Hits': pathway_data.get('low_threshold_hits', 0),
                    'High_Completeness_Percentage': pathway_data.get('high_completeness_percentage', 0),
                    'Low_Completeness_Percentage': pathway_data.get('low_completeness_percentage', 0),
                    'Average_Identity_High': pathway_data.get('average_identity_high', 0),
                    'Average_Identity_Low': pathway_data.get('average_identity_low', 0),
                    'Metabolic_Pathway_Completeness': pathway_data.get('metabolic_pathway_completeness', ''),
                    'Search_Status': pathway_data.get('search_status', ''),
                    'Detection_Status': 'Detected' if pathway_data.get('low_threshold_hits', 0) > 0 else 'Not Detected'
                })
            elif pathway_type == 'Cultivation':
                # Cultivation assessment
                summary_data.append({
                    'Input_File': input_file,
                    'Analysis_Timestamp': analysis_timestamp,
                    'Pathway_Type': pathway_type,
                    'Database': db_name,
                    'Pathway_Name': pathway_data['pathway_name'],
                    'Reference_Sequences': pathway_data['reference_sequences'],
                    'Cultivability_Hits': pathway_data.get('cultivability_hits', 0),
                    'Cultivability_Percentage': pathway_data.get('cultivability_percentage', 0),
                    'Cultivability_Status': pathway_data.get('cultivability_status', ''),
                    'Average_Identity': pathway_data.get('average_identity', 0),
                    'Search_Status': pathway_data.get('search_status', ''),
                    'Detection_Status': 'Detected' if pathway_data.get('cultivability_hits', 0) > 0 else 'Not Detected'
                })
        
        if summary_data:
            df_summary = pd.DataFrame(summary_data)
            df_summary.to_csv(summary_csv_file, index=False)
        
        self.logger.info(f"Results saved to {summary_csv_file}")
    
    def get_pathway_summary(self, results: Dict) -> Dict[str, Dict]:
        """
        Get pathway analysis summary
        
        Args:
            results: Analysis results
            
        Returns:
            Pathway summary dictionary
        """
        pathway_summary = {
            'methane_pathways': {},
            'sulfur_pathways': {},
            'nitrogen_pathways': {},
            'other_features': {}
        }
        
        # Categorize pathways - use consistent logic with _save_results method
        for db_name, pathway_data in results['pathway_results'].items():
            if db_name in METHANE_DATABASES:
                pathway_summary['methane_pathways'][db_name] = pathway_data
            elif db_name in SULFUR_DATABASES:
                pathway_summary['sulfur_pathways'][db_name] = pathway_data
            elif db_name in NITROGEN_DATABASES:
                pathway_summary['nitrogen_pathways'][db_name] = pathway_data
            else:
                pathway_summary['other_features'][db_name] = pathway_data
        
        return pathway_summary
    
    def cleanup_temp_files(self, output_prefix: str):
        """
        Clean up temporary files
        
        Args:
            output_prefix: Output file prefix
        """
        temp_files = []
        for file in os.listdir(self.results_dir):
            if file.startswith(output_prefix) and file.endswith('_diamond.tsv'):
                temp_files.append(os.path.join(self.results_dir, file))
        
        FileUtils.cleanup_temp_files(temp_files)
        self.logger.debug(f"Cleaned up {len(temp_files)} temporary files")
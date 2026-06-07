#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Diamond sequence alignment analyzer for MethArCT

Performs protein sequence alignment using Diamond BLAST for metabolic pathway prediction
and cultivability evaluation.
"""

import os
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

class DiamondAnalyzer:
    """Diamond sequence alignment analyzer"""
    
    def __init__(self, config: Optional[Config] = None):
        self.config = config or Config()
        self.logger = get_logger("diamond_analyzer")
        
        # Tool configuration - WSL support
        self.use_wsl = self.config.get('tools.diamond.use_wsl', False)
        
        if self.use_wsl:
            self.diamond_path = self.config.get('tools.diamond.wsl_path', 'wsl diamond')
        else:
            self.diamond_path = self.config.get('tools.diamond.path', 'diamond')
            
        self.threads = self.config.get('tools.diamond.threads', 4)
        self.evalue = self.config.get('tools.diamond.evalue', 1e-5)
        self.max_target_seqs = self.config.get('tools.diamond.max_target_seqs', 1)
        self.identity_threshold = self.config.get('tools.diamond.identity_threshold', 30.0)
        
        # Dual threshold analysis parameters
        self.HIGH_E_VALUE_THRESHOLD = 1e-100
        self.HIGH_BITSCORE_THRESHOLD = 400
        self.LOW_E_VALUE_THRESHOLD = 1e-5
        self.LOW_BITSCORE_THRESHOLD = 40
        self.HIGH_QUALITY_THRESHOLD = 60
        
        # Cultivability assessment thresholds
        self.CULTIVABILITY_E_VALUE_THRESHOLD = 1e-3
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
        
        # Check tool availability
        self._check_diamond_availability()
        
        # Prepare databases
        self._prepare_databases()
    
    def _check_diamond_availability(self) -> bool:
        """
        Check if Diamond tool is available - WSL environment support
        
        Returns:
            True if Diamond is available, False otherwise
        """
        try:
            if self.use_wsl:
                # Check in WSL environment
                result = subprocess.run(
                    ['wsl', 'diamond', 'version'],
                    capture_output=True,
                    text=True,
                    timeout=30
                )
            else:
                # Check in local environment
                result = subprocess.run(
                    [self.diamond_path, 'version'],
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
    
    def _prepare_databases(self):
        """
        Prepare Diamond databases from reference sequences
        """
        self.logger.info("Preparing Diamond databases...")
        
        # Get all database paths from config
        db_paths = self.config.get_all_database_paths()
        
        for db_name, fasta_path in db_paths.items():
            db_file = os.path.join(self.db_dir, f"{db_name}.dmnd")
            
            # Check if database exists and is newer than source file
            if (os.path.exists(db_file) and 
                os.path.exists(fasta_path) and
                os.path.getmtime(db_file) > os.path.getmtime(fasta_path)):
                self.logger.debug(f"Database {db_name} is up to date")
                continue
            
            # Create or update database
            if os.path.exists(fasta_path):
                self._create_diamond_database(fasta_path, db_name)
            else:
                self.logger.warning(f"Reference file not found: {fasta_path}")
    
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
            
            cmd = [
                self.diamond_path, 'makedb',
                '--in', fasta_path,
                '--db', db_path,
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
            databases = list(self.config.get_all_database_paths().keys())
        
        self.logger.info(f"Starting Diamond analysis for {input_file.name}")
        self.logger.info(f"Searching against {len(databases)} databases")
        
        # Run Diamond searches
        search_results = {}
        for db_name in databases:
            result = self._run_diamond_search(input_file, db_name, output_prefix)
            if result:
                search_results[db_name] = result
        
        # Process and summarize results
        analysis_results = self._process_results(search_results, input_file)
        
        # Save results
        self._save_results(analysis_results, output_prefix)
        
        return analysis_results
    
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
        db_path = os.path.join(self.db_dir, f"{db_name}.dmnd")
        
        if not os.path.exists(db_path):
            self.logger.warning(f"Database not found: {db_path}")
            return None
        
        # Create output file in the same directory as output_prefix
        output_dir = os.path.dirname(output_prefix) if os.path.dirname(output_prefix) else self.results_dir
        output_file = os.path.join(
            output_dir, 
            f"{os.path.basename(output_prefix)}_{db_name}_diamond.tsv"
        )
        
        # Select E-value threshold based on database type
        if db_name == 'CULTIVATION':
            search_evalue = self.CULTIVABILITY_E_VALUE_THRESHOLD
        else:
            search_evalue = self.LOW_E_VALUE_THRESHOLD
        
        cmd = [
            self.diamond_path, 'blastp',
            '--query', str(input_file),
            '--db', db_path,
            '--out', output_file,
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
                    df = pd.read_csv(
                        output_file,
                        sep='\t',
                        header=None,
                        names=['qseqid', 'sseqid', 'pident', 'length',
                               'mismatch', 'gapopen', 'qstart', 'qend',
                               'sstart', 'send', 'evalue', 'bitscore']
                    )
                    
                    # Return different result structures based on database type
                    if db_name == 'CULTIVATION':
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
                        if db_name in ['CO2-CH4', 'JIAAN-CH4', 'JIACHUN-CH4', 'JIALIUCHUN-CH4', 
                                     'YISUAN-CH4', 'C16-CH4', 'CO-CH4', 'JIASUAN-CH4', 
                                     'JIAYANGJI-CH4', 'ZHIFANGSUAN-CH4', '2JIAAN-CH4', '3JIAAN-CH4',
                                     'Glycine betaine methanogenesis', 'Methylthiopropionate methanogenesis', 'Tetramethylammonium methanogenesis',
                                     'Methanol dismutation methanogenesis']:
                            # Methane metabolism pathway: single threshold analysis (only keep low threshold)
                            low_threshold_hits = df[
                                (df['evalue'] <= self.LOW_E_VALUE_THRESHOLD) & 
                                (df['bitscore'] >= self.LOW_BITSCORE_THRESHOLD)
                            ]
                            
                            return {
                                'database': db_name,
                                'total_hits': len(df),
                                'low_threshold_hits': low_threshold_hits,
                                'output_file': output_file
                            }
                        else:
                            # Other metabolic pathways (sulfur, nitrogen) maintain dual threshold analysis
                            high_threshold_hits = df[
                                (df['evalue'] <= self.HIGH_E_VALUE_THRESHOLD) & 
                                (df['bitscore'] >= self.HIGH_BITSCORE_THRESHOLD)
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
                    # Return corresponding empty DataFrame structure based on database type
                    if db_name == 'CULTIVATION':
                        return {
                            'database': db_name,
                            'total_hits': 0,
                            'cultivability_hits': pd.DataFrame(),
                            'output_file': output_file
                        }
                    else:
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
    
    def _process_results(self, search_results: Dict, input_file: Path) -> Dict:
        """
        Process and summarize Diamond search results with dual threshold analysis
        
        Args:
            search_results: Raw search results
            input_file: Input FASTA file
            
        Returns:
            Processed analysis results
        """
        # Count input sequences
        total_sequences = FileUtils.count_sequences(input_file)
        
        # Initialize results
        pathway_results = {}
        summary_stats = {
            'total_input_sequences': total_sequences,
            'databases_searched': len(search_results),
            'pathways_detected': 0,
            'total_hits': 0
        }
        
        # Process each database result
        for db_name, result in search_results.items():
            if result and result['total_hits'] > 0:
                reference_count = self.reference_counts.get(db_name, 1)
                
                if db_name == 'CULTIVATION':
                    # Cultivability database processing
                    # Calculate how many different query sequences matched cultivability-related genes (deduplicated)
                    cultivability_matched = len(result['cultivability_hits']['qseqid'].drop_duplicates())
                    # Modified calculation logic: matches / total genes * 100%
                    cultivability_percentage = (cultivability_matched / total_sequences) * 100
                    
                    if cultivability_percentage <= 11:
                        cultivability_status = "Cultivable"
                    elif cultivability_percentage <= 13:
                        cultivability_status = "Co-cultivation dependent"
                    else:
                        cultivability_status = "Obligate symbiont"
                    
                    pathway_results[db_name] = {
                        'pathway_name': self.pathway_names.get(db_name, db_name),
                        'reference_sequences': reference_count,
                        'cultivability_hits': cultivability_matched,
                        'cultivability_percentage': round(cultivability_percentage, 2),
                        'cultivability_status': cultivability_status,
                        'average_identity': round(
                            result['cultivability_hits']['pident'].mean(), 2
                        ) if len(result['cultivability_hits']) > 0 else 0,
                        'best_hits': result['cultivability_hits'].nlargest(5, 'bitscore').to_dict('records')
                    }
                    
                else:
                    # Metabolic pathway database processing (dual threshold analysis)
                    # Check if it's a methane metabolism pathway, if so use single threshold analysis
                    if db_name in ['CO2-CH4', 'JIAAN-CH4', 'JIACHUN-CH4', 'JIALIUCHUN-CH4', 
                                 'YISUAN-CH4', 'C16-CH4', 'CO-CH4', 'JIASUAN-CH4', 
                                 'JIAYANGJI-CH4', 'ZHIFANGSUAN-CH4', '2JIAAN-CH4', '3JIAAN-CH4',
                                 'Glycine betaine methanogenesis', 'Methylthiopropionate methanogenesis', 'Tetramethylammonium methanogenesis',
                                 'Methanol dismutation methanogenesis']:
                        # Methane metabolism pathway: single threshold analysis (only keep low threshold)
                        # Note: Homologous genes (e.g., FwdF/FwdG/FwdH, Hmd/Mtd) have similar sequences but represent different gene functions and should be counted separately
                        low_threshold_matched = len(result['low_threshold_hits'].drop_duplicates(subset=['sseqid']))
                        
                        # For formate methanogenesis pathway, add special handling
                        if db_name == 'JIASUAN-CH4':
                            # Debug info: confirm code is being executed
                            self.logger.debug("=== JIASUAN-CH4 pathway evaluation logic executing ===")
                            
                            # Get all unique matched gene names (second column)
                            unique_hits = result['low_threshold_hits'].iloc[:, 1].unique()
                            
                            # Check key gene matching
                            fwd_fgh_present = any(gene in unique_hits for gene in ['FwdF', 'FwdG', 'FwdH'])
                            hmd_mtd_present = any(gene in unique_hits for gene in ['Hmd', 'Mtd'])
                            
                            # Check high-strict matching conditions for FdhA and FdhB genes
                            fdhA_hits = result['low_threshold_hits'][
                                (result['low_threshold_hits'].iloc[:, 1] == 'FdhA') &
                                (result['low_threshold_hits']['bitscore'] > 380) &
                                (result['low_threshold_hits']['evalue'] <= 1e-100)
                            ]
                            fdhB_hits = result['low_threshold_hits'][
                                (result['low_threshold_hits'].iloc[:, 1] == 'FdhB') &
                                (result['low_threshold_hits']['bitscore'] > 380) &
                                (result['low_threshold_hits']['evalue'] <= 1e-100)
                            ]
                            fdhA_present = len(fdhA_hits) > 0
                            fdhB_present = len(fdhB_hits) > 0
                            
                            # Debug info: output key gene detection status
                            self.logger.debug(f"JIASUAN-CH4 key gene detection status:")
                            self.logger.debug(f"  FwdF/FwdG/FwdH: {fwd_fgh_present}")
                            self.logger.debug(f"  Hmd/Mtd: {hmd_mtd_present}")
                            self.logger.debug(f"  FdhA: {fdhA_present} (hits: {len(fdhA_hits)})")
                            self.logger.debug(f"  FdhB: {fdhB_present} (hits: {len(fdhB_hits)})")
                            self.logger.debug(f"  Base gene coverage: {len(unique_hits)}/{reference_count}")
                            
                            # Optimized evaluation logic: use weighted scoring method
                            # Base gene coverage (weight 40%): number of detected genes
                            base_score = min(len(unique_hits) / reference_count * 40, 40)
                            
                            # Key functional genes (weight 60%):
                            # - FwdF/FwdG/FwdH (weight 15%)
                            # - Hmd/Mtd (weight 15%)
                            # - FdhA (weight 15%)
                            # - FdhB (weight 15%)
                            fwd_fgh_score = 15 if fwd_fgh_present else 0
                            hmd_mtd_score = 15 if hmd_mtd_present else 0
                            fdhA_score = 15 if fdhA_present else 0
                            fdhB_score = 15 if fdhB_present else 0
                            
                            # Calculate total score
                            total_score = base_score + fwd_fgh_score + hmd_mtd_score + fdhA_score + fdhB_score
                            
                            # Determine matched gene count based on score
                            # If all key functional genes are satisfied, directly mark as 100% completeness
                            if fwd_fgh_present and hmd_mtd_present and fdhA_present and fdhB_present:
                                self.logger.debug(f"JIASUAN-CH4: All key functional genes satisfied, using reference sequence count: {reference_count}")
                                low_threshold_matched = reference_count
                            elif total_score >= 80:  # Score >= 80 considered pathway complete
                                self.logger.debug(f"JIASUAN-CH4: Total score {total_score} >= 80, using reference sequence count: {reference_count}")
                                low_threshold_matched = reference_count
                            elif total_score >= 60:  # 60-79 use actual detection count
                                self.logger.debug(f"JIASUAN-CH4: Total score {total_score} between 60-79, using actual detection count: {len(unique_hits)}")
                                low_threshold_matched = len(unique_hits)
                            else:  # Below 60 use actual detection count
                                self.logger.debug(f"JIASUAN-CH4: Total score {total_score} < 60, using actual detection count: {len(unique_hits)}")
                                low_threshold_matched = len(unique_hits)
                        
                        # For CO2 reduction methanogenesis pathway, add special handling
                        elif db_name == 'CO2-CH4':
                            # Debug info: confirm code is being executed
                            self.logger.debug("=== CO2-CH4 pathway evaluation logic executing ===")
                            
                            # Get all unique matched gene names (second column)
                            unique_hits = result['low_threshold_hits'].iloc[:, 1].unique()
                            
                            # Check key gene matching
                            # CO2-CH4 pathway key genes are FwdF/FwdG/FwdH and Hmd/Mtd, no need for FdhA/FdhB
                            fwd_fgh_present = any(gene in unique_hits for gene in ['FwdF', 'FwdG', 'FwdH'])
                            hmd_mtd_present = any(gene in unique_hits for gene in ['Hmd', 'Mtd'])
                            
                            # Debug info: output key gene detection status
                            self.logger.debug(f"CO2-CH4 key gene detection status:")
                            self.logger.debug(f"  FwdF/FwdG/FwdH: {fwd_fgh_present}")
                            self.logger.debug(f"  Hmd/Mtd: {hmd_mtd_present}")
                            self.logger.debug(f"  Base gene coverage: {len(unique_hits)}/{reference_count}")
                            
                            # Optimized evaluation logic: use weighted scoring method
                            # Base gene coverage (weight 60%): number of detected genes
                            base_score = min(len(unique_hits) / reference_count * 60, 60)
                            
                            # Key functional genes (weight 40%):
                            # - FwdF/FwdG/FwdH (weight 20%)
                            # - Hmd/Mtd (weight 20%)
                            fwd_fgh_score = 20 if fwd_fgh_present else 0
                            hmd_mtd_score = 20 if hmd_mtd_present else 0
                            
                            # Calculate total score
                            total_score = base_score + fwd_fgh_score + hmd_mtd_score
                            
                            # Determine matched gene count based on score
                            # If all key functional genes are satisfied, directly mark as 100% completeness
                            if fwd_fgh_present and hmd_mtd_present:
                                self.logger.debug(f"CO2-CH4: All key functional genes satisfied, using reference sequence count: {reference_count}")
                                low_threshold_matched = reference_count
                            elif total_score >= 80:  # Score >= 80 considered pathway complete
                                self.logger.debug(f"CO2-CH4: Total score {total_score} >= 80, using reference sequence count: {reference_count}")
                                low_threshold_matched = reference_count
                            elif total_score >= 60:  # 60-79 use actual detection count
                                self.logger.debug(f"CO2-CH4: Total score {total_score} between 60-79, using actual detection count: {len(unique_hits)}")
                                low_threshold_matched = len(unique_hits)
                            else:  # Below 60 use actual detection count
                                self.logger.debug(f"CO2-CH4: Total score {total_score} < 60, using actual detection count: {len(unique_hits)}")
                                low_threshold_matched = len(unique_hits)
                        
                        # For methanethiol methanogenesis pathway, add special handling
                        elif db_name == 'JIALIUCHUN-CH4':
                            # Debug info: confirm code is being executed
                            self.logger.debug("=== JIALIUCHUN-CH4 pathway evaluation logic executing ===")
                            
                            # Get all unique matched gene names (second column)
                            unique_hits = result['low_threshold_hits'].iloc[:, 1].unique()
                            
                            # Check key gene matching
                            # Methanethiol-specific genes: MtsA1, MtsA2 (add strict matching conditions: bitscore > 200 and evalue <= 1e-100)
                            mtsA1_hits = result['low_threshold_hits'][
                                (result['low_threshold_hits'].iloc[:, 1] == 'MtsA1') &
                                (result['low_threshold_hits']['bitscore'] > 200) &
                                (result['low_threshold_hits']['evalue'] <= 1e-100)
                            ]
                            mtsA2_hits = result['low_threshold_hits'][
                                (result['low_threshold_hits'].iloc[:, 1] == 'MtsA2') &
                                (result['low_threshold_hits']['bitscore'] > 200) &
                                (result['low_threshold_hits']['evalue'] <= 1e-100)
                            ]
                            mtsA1_present = len(mtsA1_hits) > 0
                            mtsA2_present = len(mtsA2_hits) > 0
                            
                            # Methyl-coenzyme M reductase related genes: KYC55281.1, KYC55283.1, KYC55284.1, KYC55314.1
                            mcr_genes = ['KYC55281.1', 'KYC55283.1', 'KYC55284.1', 'KYC55314.1']
                            mcr_genes_present = [gene in unique_hits for gene in mcr_genes]
                            mcr_genes_count = sum(mcr_genes_present)
                            
                            # Check if any 3 of the methyl-coenzyme M reductase related genes are matched
                            mcr_3_present = mcr_genes_count >= 3
                            
                            # Debug info: output key gene detection status
                            self.logger.debug(f"JIALIUCHUN-CH4 key gene detection status:")
                            self.logger.debug(f"  MtsA1: {mtsA1_present} (hits: {len(mtsA1_hits)})")
                            self.logger.debug(f"  MtsA2: {mtsA2_present} (hits: {len(mtsA2_hits)})")
                            self.logger.debug(f"  Methyl-coenzyme M reductase gene matching: {mcr_genes_count}/4")
                            self.logger.debug(f"  Any 3 methyl-coenzyme M reductase genes matched: {mcr_3_present}")
                            self.logger.debug(f"  Base gene coverage: {len(unique_hits)}/{reference_count}")
                            
                            # Output detailed MtsA1 and MtsA2 bitscore info for debugging
                            if len(mtsA1_hits) > 0:
                                self.logger.debug(f"  MtsA1 max bitscore: {mtsA1_hits['bitscore'].max()}")
                            if len(mtsA2_hits) > 0:
                                self.logger.debug(f"  MtsA2 max bitscore: {mtsA2_hits['bitscore'].max()}")
                            
                            # Determine matched gene count based on new strategy
                            # If both MtsA1 and MtsA2 are matched, and any 3 of methyl-coenzyme M reductase related genes are matched, mark as 100% completeness
                            if mtsA1_present and mtsA2_present and mcr_3_present:
                                self.logger.debug(f"JIALIUCHUN-CH4: Key gene conditions satisfied, using reference sequence count: {reference_count}")
                                low_threshold_matched = reference_count
                            else:
                                # When conditions not met, if MtsA1 doesn't meet strict conditions, exclude MtsA1 from actual detection count
                                # Because MtsA1 is a key gene for methanethiol pathway, must meet strict conditions
                                adjusted_hits = len(unique_hits)
                                if not mtsA1_present:
                                    # If MtsA1 doesn't meet conditions, subtract 1 from detection count (MtsA1 is essential gene)
                                    adjusted_hits = max(0, len(unique_hits) - 1)
                                    self.logger.debug(f"JIALIUCHUN-CH4: MtsA1 doesn't meet strict conditions, adjusting detection count: {len(unique_hits)} -> {adjusted_hits}")
                                else:
                                    self.logger.debug(f"JIALIUCHUN-CH4: Key gene conditions not met, using actual detection count: {len(unique_hits)}")
                                low_threshold_matched = adjusted_hits
                        
                        # For glycine betaine methanogenesis pathway, add special handling
                        elif db_name == 'Glycine betaine methanogenesis':
                            # Debug info: confirm code is being executed
                            self.logger.debug("=== Glycine betaine methanogenesis pathway evaluation logic executing ===")
                            
                            # Get all unique matched gene names (second column)
                            unique_hits = result['low_threshold_hits'].iloc[:, 1].unique()
                            
                            # Check key gene matching
                            # Glycine betaine methanogenesis pathway key genes: MtgB, dimethylamine_corrinoid_protein_3, MV10360
                            mtgB_present = 'MtgB' in unique_hits
                            dimethylamine_corrinoid_present = 'dimethylamine_corrinoid_protein_3' in unique_hits
                            mv10360_present = 'MV10360' in unique_hits
                            
                            # Check if all key genes are matched
                            all_key_genes_present = mtgB_present and dimethylamine_corrinoid_present and mv10360_present
                            
                            # Debug info: output key gene detection status
                            self.logger.debug(f"Glycine betaine methanogenesis key gene detection status:")
                            self.logger.debug(f"  MtgB: {mtgB_present}")
                            self.logger.debug(f"  dimethylamine_corrinoid_protein_3: {dimethylamine_corrinoid_present}")
                            self.logger.debug(f"  MV10360: {mv10360_present}")
                            self.logger.debug(f"  All key genes matched: {all_key_genes_present}")
                            self.logger.debug(f"  Base gene coverage: {len(unique_hits)}/{reference_count}")
                            
                            # Determine matched gene count based on new strategy
                            # If all 3 key genes are matched, mark as 100% completeness
                            if all_key_genes_present:
                                self.logger.debug(f"Glycine betaine methanogenesis: All key gene conditions satisfied, using reference sequence count: {reference_count}")
                                low_threshold_matched = reference_count
                            else:
                                # When conditions not met, use actual detection count
                                self.logger.debug(f"Glycine betaine methanogenesis: Key gene conditions not met, using actual detection count: {len(unique_hits)}")
                                low_threshold_matched = len(unique_hits)
                        
                        # For methylthiopropionate methanogenesis pathway, add special handling
                        elif db_name == 'Methylthiopropionate methanogenesis':
                            # Debug info: confirm code is being executed
                            self.logger.debug("=== Methylthiopropionate methanogenesis pathway evaluation logic executing ===")
                            
                            # Get all unique matched gene names (second column)
                            unique_hits = result['low_threshold_hits'].iloc[:, 1].unique()
                            
                            # Check key gene matching, require bitscore > 100 and evalue <= 1e-5
                            # Methylthiopropionate methanogenesis pathway key genes: mtpA1, mtsA1, mtpA2, mtsA2
                            mtpA1_hits = result['low_threshold_hits'][
                                (result['low_threshold_hits'].iloc[:, 1] == 'mtpA1') &
                                (result['low_threshold_hits']['bitscore'] > 100) &
                                (result['low_threshold_hits']['evalue'] <= 1e-5)
                            ]
                            mtsA1_hits = result['low_threshold_hits'][
                                (result['low_threshold_hits'].iloc[:, 1] == 'mtsA1') &
                                (result['low_threshold_hits']['bitscore'] > 100) &
                                (result['low_threshold_hits']['evalue'] <= 1e-5)
                            ]
                            mtpA2_hits = result['low_threshold_hits'][
                                (result['low_threshold_hits'].iloc[:, 1] == 'mtpA2') &
                                (result['low_threshold_hits']['bitscore'] > 100) &
                                (result['low_threshold_hits']['evalue'] <= 1e-5)
                            ]
                            mtsA2_hits = result['low_threshold_hits'][
                                (result['low_threshold_hits'].iloc[:, 1] == 'mtsA2') &
                                (result['low_threshold_hits']['bitscore'] > 100) &
                                (result['low_threshold_hits']['evalue'] <= 1e-5)
                            ]
                            
                            mtpA1_present = len(mtpA1_hits) > 0
                            mtsA1_present = len(mtsA1_hits) > 0
                            mtpA2_present = len(mtpA2_hits) > 0
                            mtsA2_present = len(mtsA2_hits) > 0
                            
                            # Check if all key genes are matched
                            all_key_genes_present = mtpA1_present and mtsA1_present and mtpA2_present and mtsA2_present
                            
                            # Debug info: output key gene detection status
                            self.logger.debug(f"Methylthiopropionate methanogenesis key gene detection status:")
                            self.logger.debug(f"  mtpA1: {mtpA1_present} (hits: {len(mtpA1_hits)})")
                            self.logger.debug(f"  mtsA1: {mtsA1_present} (hits: {len(mtsA1_hits)})")
                            self.logger.debug(f"  mtpA2: {mtpA2_present} (hits: {len(mtpA2_hits)})")
                            self.logger.debug(f"  mtsA2: {mtsA2_present} (hits: {len(mtsA2_hits)})")
                            self.logger.debug(f"  All key genes matched: {all_key_genes_present}")
                            self.logger.debug(f"  Base gene coverage: {len(unique_hits)}/{reference_count}")
                            
                            # Determine matched gene count based on new strategy
                            # If all 4 key genes are matched (bitscore>100 and evalue<1e-5), mark as 100% completeness
                            if all_key_genes_present:
                                self.logger.debug(f"Methylthiopropionate methanogenesis: All key gene conditions satisfied, using reference sequence count: {reference_count}")
                                low_threshold_matched = reference_count
                            else:
                                # When conditions not met, use actual detection count
                                self.logger.debug(f"Methylthiopropionate methanogenesis: Key gene conditions not met, using actual detection count: {len(unique_hits)}")
                                low_threshold_matched = len(unique_hits)
                        
                        # For tetramethylammonium methanogenesis pathway, add special handling
                        elif db_name == 'Tetramethylammonium methanogenesis':
                            # Debug info: confirm code is being executed
                            self.logger.debug("=== Tetramethylammonium methanogenesis pathway evaluation logic executing ===")
                            
                            # Get all unique matched gene names (second column)
                            unique_hits = result['low_threshold_hits'].iloc[:, 1].unique()
                            
                            # Check key gene matching, require bitscore > 200 and evalue <= 1e-100
                            # Key genes for Tetramethylammonium methanogenesis pathway: MtqA/MT2, MtqB, MtqC
                            mtqA_hits = result['low_threshold_hits'][
                                (result['low_threshold_hits'].iloc[:, 1] == 'MtqA/MT2') &
                                (result['low_threshold_hits']['bitscore'] > 200) &
                                (result['low_threshold_hits']['evalue'] <= 1e-100)
                            ]
                            mtqB_hits = result['low_threshold_hits'][
                                (result['low_threshold_hits'].iloc[:, 1] == 'MtqB') &
                                (result['low_threshold_hits']['bitscore'] > 200) &
                                (result['low_threshold_hits']['evalue'] <= 1e-100)
                            ]
                            mtqC_hits = result['low_threshold_hits'][
                                (result['low_threshold_hits'].iloc[:, 1] == 'MtqC') &
                                (result['low_threshold_hits']['bitscore'] > 200) &
                                (result['low_threshold_hits']['evalue'] <= 1e-100)
                            ]
                            
                            mtqA_present = len(mtqA_hits) > 0
                            mtqB_present = len(mtqB_hits) > 0
                            mtqC_present = len(mtqC_hits) > 0
                            
                            # Check if all key genes are matched
                            all_key_genes_present = mtqA_present and mtqB_present and mtqC_present
                            
                            # Debug info: output key gene detection status
                            self.logger.debug(f"Tetramethylammonium methanogenesis key gene detection status:")
                            self.logger.debug(f"  MtqA/MT2: {mtqA_present} (hits: {len(mtqA_hits)})")
                            self.logger.debug(f"  MtqB: {mtqB_present} (hits: {len(mtqB_hits)})")
                            self.logger.debug(f"  MtqC: {mtqC_present} (hits: {len(mtqC_hits)})")
                            self.logger.debug(f"  All key genes matched: {all_key_genes_present}")
                            self.logger.debug(f"  Base gene coverage: {len(unique_hits)}/{reference_count}")
                            
                            # Determine matched gene count based on new strategy
                            # If all 3 key genes are matched (bitscore>200 and evalue<1e-100), mark as 100% completeness
                            if all_key_genes_present:
                                self.logger.debug(f"Tetramethylammonium methanogenesis: All key gene conditions satisfied, using reference sequence count: {reference_count}")
                                low_threshold_matched = reference_count
                            else:
                                # When conditions not met, use actual detection count
                                self.logger.debug(f"Tetramethylammonium methanogenesis: Key gene conditions not met, using actual detection count: {len(unique_hits)}")
                                low_threshold_matched = len(unique_hits)
                        
                        # For methanol dismutation methanogenesis pathway, add special handling
                        elif db_name == 'Methanol dismutation methanogenesis':
                            # Debug info: confirm code is being executed
                            self.logger.debug("=== Methanol dismutation methanogenesis pathway evaluation logic executing ===")
                            
                            # Get all unique matched gene names (second column)
                            unique_hits = result['low_threshold_hits'].iloc[:, 1].unique()
                            
                            # Check key gene matching
                            # Key genes for Methanol dismutation methanogenesis pathway: MvhA and elpA (one of the two genes matched is sufficient)
                            mvhA_present = 'MvhA' in unique_hits
                            elpA_present = 'elpA' in unique_hits
                            
                            # Check other gene matching (all genes except MvhA and elpA)
                            # Other key genes for Methanol dismutation methanogenesis pathway: FwdA-FwdH, Ftr, Mch, Hmd, Mtd, elpB, elpC, etc.
                            other_genes = ['FwdA', 'FwdB', 'FwdC', 'FwdD', 'FwdE', 'FwdF', 'FwdG', 'FwdH', 
                                         'Ftr', 'Mch', 'Hmd', 'Mtd', 'elpB', 'elpC']
                            other_genes_present = [gene in unique_hits for gene in other_genes]
                            other_genes_count = sum(other_genes_present)
                            
                            # Check if all other genes are matched (except MvhA and elpA)
                            all_other_genes_present = other_genes_count == len(other_genes)
                            
                            # Debug info: output key gene detection status
                            self.logger.debug(f"Methanol dismutation methanogenesis key gene detection status:")
                            self.logger.debug(f"  MvhA: {mvhA_present}")
                            self.logger.debug(f"  elpA: {elpA_present}")
                            self.logger.debug(f"  MvhA or elpA matched: {mvhA_present or elpA_present}")
                            self.logger.debug(f"  Other genes matching: {other_genes_count}/{len(other_genes)}")
                            self.logger.debug(f"  All other genes matched: {all_other_genes_present}")
                            self.logger.debug(f"  Base gene coverage: {len(unique_hits)}/{reference_count}")
                            
                            # Determine matched gene count based on new strategy
                            # If all other genes are matched, and one of MvhA/elpA is matched, mark as 100% completeness
                            if all_other_genes_present and (mvhA_present or elpA_present):
                                self.logger.debug(f"Methanol dismutation methanogenesis: All other genes matched and one of MvhA/elpA matched, using reference sequence count: {reference_count}")
                                low_threshold_matched = reference_count
                            else:
                                # When conditions not met, use actual detection count
                                self.logger.debug(f"Methanol dismutation methanogenesis: Conditions not met, using actual detection count: {len(unique_hits)}")
                                low_threshold_matched = len(unique_hits)
                            
                            # Ensure completeness does not exceed 100%
                            if low_threshold_matched > reference_count:
                                self.logger.debug(f"Methanol dismutation methanogenesis: Detected gene count ({low_threshold_matched}) exceeds reference sequence count ({reference_count}), limited to 100%")
                                low_threshold_matched = reference_count
                        
                        low_completeness = (low_threshold_matched / reference_count) * 100
                        
                        # Ensure completeness percentage does not exceed 100%
                        if low_completeness > 100:
                            self.logger.debug(f"Completeness percentage exceeds 100%, limited to 100%")
                            low_completeness = 100
                        
                        pathway_results[db_name] = {
                            'pathway_name': self.pathway_names.get(db_name, db_name),
                            'low_threshold_hits': low_threshold_matched,
                            'reference_sequences': reference_count,
                            'low_completeness_percentage': round(low_completeness, 2),
                            'average_identity_low': round(
                                result['low_threshold_hits']['pident'].mean(), 2
                            ) if len(result['low_threshold_hits']) > 0 else 0,
                            'best_hits_low': result['low_threshold_hits'].nlargest(5, 'bitscore').to_dict('records')
                        }
                    else:
                        # Other metabolic pathways (sulfur, nitrogen) maintain dual threshold analysis
                        high_threshold_matched = len(result['high_threshold_hits'].drop_duplicates(subset=['sseqid']))
                        low_threshold_matched = len(result['low_threshold_hits'].drop_duplicates(subset=['sseqid']))
                        
                        # Ensure detected gene count does not exceed reference sequence count
                        if high_threshold_matched > reference_count:
                            self.logger.debug(f"{db_name}: High threshold detected gene count ({high_threshold_matched}) exceeds reference sequence count ({reference_count}), limited to 100%")
                            high_threshold_matched = reference_count
                        if low_threshold_matched > reference_count:
                            self.logger.debug(f"{db_name}: Low threshold detected gene count ({low_threshold_matched}) exceeds reference sequence count ({reference_count}), limited to 100%")
                            low_threshold_matched = reference_count
                        
                        high_completeness = (high_threshold_matched / reference_count) * 100
                        low_completeness = (low_threshold_matched / reference_count) * 100
                        
                        # Ensure completeness percentage does not exceed 100%
                        if high_completeness > 100:
                            self.logger.debug(f"{db_name}: High threshold completeness exceeds 100%, limited to 100%")
                            high_completeness = 100
                        if low_completeness > 100:
                            self.logger.debug(f"{db_name}: Low threshold completeness exceeds 100%, limited to 100%")
                            low_completeness = 100
                        
                        # Metabolic pathway completeness format: "high_threshold%~low_threshold%"
                        metabolic_pathway_completeness = f"{high_completeness:.2f}%~{low_completeness:.2f}%"
                        
                        pathway_results[db_name] = {
                            'pathway_name': self.pathway_names.get(db_name, db_name),
                            'high_threshold_hits': high_threshold_matched,
                            'low_threshold_hits': low_threshold_matched,
                            'reference_sequences': reference_count,
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
                
                summary_stats['pathways_detected'] += 1
                if 'low_threshold_hits' in result:
                    summary_stats['total_hits'] += len(result['low_threshold_hits'])
                elif 'cultivability_hits' in result:
                    pass
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
            if db_name in ['CO2-CH4', 'JIAAN-CH4', 'JIACHUN-CH4', 'JIALIUCHUN-CH4', 
                         'YISUAN-CH4', 'C16-CH4', 'CO-CH4', 'JIASUAN-CH4', 
                         'JIAYANGJI-CH4', 'ZHIFANGSUAN-CH4', '2JIAAN-CH4', '3JIAAN-CH4',
                         'Glycine betaine methanogenesis', 'Methylthiopropionate methanogenesis', 'Tetramethylammonium methanogenesis', 'Methanol dismutation methanogenesis']:
                pathway_type = 'Methane'
            elif db_name in ['ASR', 'SO', 'SOX', 'S4I', 'SR', 'DSR']:
                pathway_type = 'Sulfur'
            elif db_name in ['ANR', 'DEN', 'DNR', 'NIT']:
                pathway_type = 'Nitrogen'
            elif db_name == 'CULTIVATION':
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
        methane_dbs = ['CO2-CH4', 'JIAAN-CH4', 'JIACHUN-CH4', 'JIALIUCHUN-CH4', 
                     'YISUAN-CH4', 'C16-CH4', 'CO-CH4', 'JIASUAN-CH4', 
                     'JIAYANGJI-CH4', 'ZHIFANGSUAN-CH4', '2JIAAN-CH4', '3JIAAN-CH4',
                     'Glycine betaine methanogenesis', 'Methylthiopropionate methanogenesis', 'Tetramethylammonium methanogenesis', 'Methanol dismutation methanogenesis']
        sulfur_dbs = ['ASR', 'SO', 'SOX', 'S4I', 'SR', 'DSR']
        nitrogen_dbs = ['ANR', 'DEN', 'DNR', 'NIT']
        
        for db_name, pathway_data in results['pathway_results'].items():
            if db_name in methane_dbs:
                pathway_summary['methane_pathways'][db_name] = pathway_data
            elif db_name in sulfur_dbs:
                pathway_summary['sulfur_pathways'][db_name] = pathway_data
            elif db_name in nitrogen_dbs:
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
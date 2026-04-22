#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Diamond sequence alignment analyzer for MethArCT

Performs protein sequence alignment using Diamond BLAST for metabolic pathway prediction,
salt tolerance assessment, and cultivability evaluation.
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

        # Cultivability evaluation thresholds
        self.CULTIVABILITY_E_VALUE_THRESHOLD = 1e-3
        self.CULTIVABILITY_BITSCORE_THRESHOLD = 50

        # Salt tolerance evaluation thresholds
        self.SALT_TOLERANCE_E_VALUE_THRESHOLD = 1e-5
        self.SALT_TOLERANCE_BITSCORE_THRESHOLD = 100
        
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
        Check if Diamond tool is available

        Returns:
            True if Diamond is available, False otherwise
        """
        try:
            if self.use_wsl:
                # WSL environment check
                result = subprocess.run(
                    ['wsl', 'diamond', 'version'],
                    capture_output=True,
                    text=True,
                    timeout=30
                )
            else:
                # Local environment check
                result = subprocess.run(
                    [self.diamond_path, 'version'],
                    capture_output=True,
                    text=True,
                    timeout=30
                )

            if result.returncode == 0:
                version = result.stdout.strip()
                env_type = "WSL" if self.use_wsl else "Local"
                self.logger.info(f"Diamond available ({env_type}): {version}")
                return True
            else:
                self.logger.error(f"Diamond not available: {result.stderr}")
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
        
        # Create temporary output file
        output_file = os.path.join(
            self.results_dir, 
            f"{output_prefix}_{db_name}_diamond.tsv"
        )
        
        # Select E-value threshold based on database type
        if db_name == 'CULTIVATION':
            search_evalue = self.CULTIVABILITY_E_VALUE_THRESHOLD
        elif db_name == 'SALT_TOLERANCE':
            search_evalue = self.SALT_TOLERANCE_E_VALUE_THRESHOLD
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
                    elif db_name == 'SALT_TOLERANCE':
                        salt_tolerance_hits = df[
                            (df['evalue'] <= self.SALT_TOLERANCE_E_VALUE_THRESHOLD) &
                            (df['bitscore'] >= self.SALT_TOLERANCE_BITSCORE_THRESHOLD)
                        ]

                        return {
                            'database': db_name,
                            'total_hits': len(df),
                            'salt_tolerance_hits': salt_tolerance_hits,
                            'output_file': output_file
                        }
                    else:
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
                    # Return empty DataFrame structure based on database type
                    if db_name == 'CULTIVATION':
                        return {
                            'database': db_name,
                            'total_hits': 0,
                            'cultivability_hits': pd.DataFrame(),
                            'output_file': output_file
                        }
                    elif db_name == 'SALT_TOLERANCE':
                        return {
                            'database': db_name,
                            'total_hits': 0,
                            'salt_tolerance_hits': pd.DataFrame(),
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
        methane_dbs = ['CO2-CH4', 'METHYLAMINE-CH4', 'METHANOL-CH4', 'METHANETHIOL-CH4', 'ACETIC_ACID-CH4',
                       'C16-CH4', 'CO-CH4', 'FORMIC_ACID-CH4', 'METHOXY-CH4', 'FATTY_ACID-CH4',
                       'DIMETHYLAMINE-CH4', 'TRIMETHYLAMINE-CH4']
        sulfur_dbs = ['ASR', 'SO', 'SOX', 'S4I', 'SR', 'DSR']
        nitrogen_dbs = ['ANR', 'DEN', 'DNR', 'NIT']

        methane_pathways_detected = 0
        sulfur_pathways_detected = 0
        nitrogen_pathways_detected = 0
        total_pathways_detected = 0
        salt_tolerance_hits = 0
        cultivability_hits = 0

        summary_stats = {
            'total_input_sequences': total_sequences,
            'databases_searched': len(search_results),
            'pathways_detected': 0,
            'total_hits': 0,
            'methane_pathways': 0,
            'sulfur_pathways': 0,
            'nitrogen_pathways': 0,
            'salt_tolerance_hits': 0,
            'cultivability_hits': 0
        }

        # Process each database result
        for db_name, result in search_results.items():
            if result and result['total_hits'] > 0:
                reference_count = self.reference_counts.get(db_name, 1)
                
                if db_name == 'CULTIVATION':
                    # Cultivability database processing
                    # Calculate unique query sequences matched to cultivation-related genes
                    cultivability_matched = len(result['cultivability_hits']['qseqid'].drop_duplicates())
                    # Modified calculation: matched items / own gene count * 100%
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
                    
                elif db_name == 'SALT_TOLERANCE':
                    # Salt tolerance database processing
                    # Calculate unique query sequences matched to salt tolerance-related genes
                    salt_tolerance_matched = len(result['salt_tolerance_hits']['qseqid'].drop_duplicates())
                    
                    if 0 <= salt_tolerance_matched <= 2:
                        salt_tolerance_level = "No Resistance"
                    elif 3 <= salt_tolerance_matched <= 5:
                        salt_tolerance_level = "Moderate"
                    elif 6 <= salt_tolerance_matched <= 8:
                        salt_tolerance_level = "Strong"
                    elif salt_tolerance_matched >= 9:
                        salt_tolerance_level = "Very Strong"
                    else:
                        salt_tolerance_level = "Unknown"
                    
                    pathway_results[db_name] = {
                        'pathway_name': self.pathway_names.get(db_name, db_name),
                        'reference_sequences': reference_count,
                        'salt_tolerance_hits': salt_tolerance_matched,
                        'salt_tolerance_level': salt_tolerance_level,
                        'average_identity': round(
                            result['salt_tolerance_hits']['pident'].mean(), 2
                        ) if len(result['salt_tolerance_hits']) > 0 else 0,
                        'best_hits': result['salt_tolerance_hits'].nlargest(5, 'bitscore').to_dict('records')
                    }
                    
                else:
                    # Metabolic pathway database processing (dual threshold analysis)
                    high_threshold_matched = len(result['high_threshold_hits'].drop_duplicates(subset=['sseqid']))
                    low_threshold_matched = len(result['low_threshold_hits'].drop_duplicates(subset=['sseqid']))

                    high_completeness = (high_threshold_matched / reference_count) * 100
                    low_completeness = (low_threshold_matched / reference_count) * 100

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

                    if db_name in methane_dbs:
                        methane_pathways_detected += 1
                    elif db_name in sulfur_dbs:
                        sulfur_pathways_detected += 1
                    elif db_name in nitrogen_dbs:
                        nitrogen_pathways_detected += 1

                summary_stats['pathways_detected'] += 1
                if db_name == 'CULTIVATION':
                    summary_stats['total_hits'] += len(result['cultivability_hits'])
                    cultivability_hits = len(result['cultivability_hits'])
                elif db_name == 'SALT_TOLERANCE':
                    summary_stats['total_hits'] += len(result['salt_tolerance_hits'])
                    salt_tolerance_hits = len(result['salt_tolerance_hits'])
                elif 'low_threshold_hits' in result:
                    summary_stats['total_hits'] += len(result['low_threshold_hits'])

        summary_stats['methane_pathways'] = methane_pathways_detected
        summary_stats['sulfur_pathways'] = sulfur_pathways_detected
        summary_stats['nitrogen_pathways'] = nitrogen_pathways_detected
        summary_stats['salt_tolerance_hits'] = salt_tolerance_hits
        summary_stats['cultivability_hits'] = cultivability_hits

        return {
            'input_file': str(input_file),
            'analysis_timestamp': pd.Timestamp.now().isoformat(),
            'summary': summary_stats,
            'pathway_results': pathway_results,
            'raw_results': search_results
        }
    
    def _save_results(self, results: Dict, output_prefix: str):
        """
        Save analysis results to files with dual threshold analysis
        
        Args:
            results: Analysis results
            output_prefix: Output file prefix
        """
        # Create serializable copy of results, convert DataFrames to list of dicts
        serializable_results = results.copy()

        # Process DataFrame objects in pathway_results
        for db_name, pathway_data in serializable_results['pathway_results'].items():
            # Convert best_hits DataFrame to list of dicts
            if 'best_hits_high' in pathway_data and hasattr(pathway_data['best_hits_high'], 'to_dict'):
                pathway_data['best_hits_high'] = pathway_data['best_hits_high']
            if 'best_hits_low' in pathway_data and hasattr(pathway_data['best_hits_low'], 'to_dict'):
                pathway_data['best_hits_low'] = pathway_data['best_hits_low']

        # Process DataFrame objects in raw_results
        if 'raw_results' in serializable_results:
            for db_name, raw_data in serializable_results['raw_results'].items():
                if raw_data:
                    # Convert DataFrame to list of dicts
                    for key in ['high_threshold_hits', 'low_threshold_hits', 'cultivability_hits', 'salt_tolerance_hits']:
                        if key in raw_data and hasattr(raw_data[key], 'to_dict'):
                            raw_data[key] = raw_data[key].to_dict('records')
        
        # Save JSON results
        json_file = os.path.join(self.results_dir, f"{output_prefix}_diamond_results.json")
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(serializable_results, f, indent=2, ensure_ascii=False)
        
        # Save CSV summary with new columns
        csv_file = os.path.join(self.results_dir, f"{output_prefix}_diamond_summary.csv")
        summary_data = []
        
        for db_name, pathway_data in results['pathway_results'].items():
            if db_name == 'CULTIVATION':
                # Cultivability database CSV format
                summary_data.append({
                    'Database': db_name,
                    'Pathway_Name': pathway_data['pathway_name'],
                    'Reference_Sequences': pathway_data['reference_sequences'],
                    'Cultivability_Hits': pathway_data.get('cultivability_hits', 0),
                    'Cultivability_Percentage': pathway_data.get('cultivability_percentage', 0),
                    'Cultivability_Status': pathway_data.get('cultivability_status', ''),
                    'Average_Identity': pathway_data.get('average_identity', 0)
                })
            elif db_name == 'SALT_TOLERANCE':
                # Salt tolerance database CSV format
                summary_data.append({
                    'Database': db_name,
                    'Pathway_Name': pathway_data['pathway_name'],
                    'Reference_Sequences': pathway_data['reference_sequences'],
                    'Salt_Tolerance_Hits': pathway_data.get('salt_tolerance_hits', 0),
                    'Salt_Tolerance_Level': pathway_data.get('salt_tolerance_level', ''),
                    'Average_Identity': pathway_data.get('average_identity', 0)
                })
            else:
                # Metabolic pathway database CSV format
                summary_data.append({
                    'Database': db_name,
                    'Pathway_Name': pathway_data['pathway_name'],
                    'Reference_Sequences': pathway_data['reference_sequences'],
                    'High_Threshold_Hits': pathway_data.get('high_threshold_hits', 0),
                    'Low_Threshold_Hits': pathway_data.get('low_threshold_hits', 0),
                    'High_Completeness_Percentage': pathway_data.get('high_completeness_percentage', 0),
                    'Low_Completeness_Percentage': pathway_data.get('low_completeness_percentage', 0),
                    'Average_Identity_High': pathway_data.get('average_identity_high', 0),
                    'Average_Identity_Low': pathway_data.get('average_identity_low', 0),
                    'Metabolic_Pathway_Completeness': pathway_data.get('metabolic_pathway_completeness', '')
                })
        
        if summary_data:
            df_summary = pd.DataFrame(summary_data)
            df_summary.to_csv(csv_file, index=False)
        
        self.logger.info(f"Results saved to {json_file} and {csv_file}")
    
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
        
        # Categorize pathways
        methane_dbs = [db for db in results['pathway_results'].keys() if 'CH4' in db]
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
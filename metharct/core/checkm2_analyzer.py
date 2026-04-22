#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CheckM2 genome quality analyzer for MethArCT

Assesses genome quality and completeness using CheckM2 for cultivability evaluation.
"""

import os
import subprocess
import json
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Union
from ..utils.config import Config
from ..utils.logger import get_logger
from ..utils.file_utils import FileUtils

class CheckM2Analyzer:
    """CheckM2 genome quality analyzer"""
    
    def __init__(self, config: Optional[Config] = None):
        self.config = config or Config()
        self.logger = get_logger("checkm2_analyzer")
        
        # Tool configuration - WSL support
        self.use_wsl = self.config.get('tools.checkm2.use_wsl', False)
        
        if self.use_wsl:
            self.checkm2_path = self.config.get('tools.checkm2.wsl_path', 'wsl checkm2')
        else:
            self.checkm2_path = self.config.get('tools.checkm2.path', 'checkm2')
            
        self.threads = self.config.get('tools.checkm2.threads', 4)
        
        # Results directory
        self.results_dir = self.config.get('output.base_dir', 'results')
        FileUtils.ensure_dir(self.results_dir)
        
        # Quality thresholds
        self.high_quality_threshold = {
            'completeness': 90.0,
            'contamination': 5.0
        }
        self.medium_quality_threshold = {
            'completeness': 50.0,
            'contamination': 10.0
        }
        
        # Check tool availability
        self._check_checkm2_availability()
    
    def _check_checkm2_availability(self) -> bool:
        """
        Check if CheckM2 tool is available

        Returns:
            True if CheckM2 is available, False otherwise
        """
        try:
            if self.use_wsl:
                # WSL environment check
                result = subprocess.run(
                    ['wsl', 'checkm2', '--version'],
                    capture_output=True,
                    text=True,
                    timeout=30
                )
            else:
                # Local environment check
                result = subprocess.run(
                    [self.checkm2_path, '--version'],
                    capture_output=True,
                    text=True,
                    timeout=30
                )

            if result.returncode == 0:
                version = result.stdout.strip()
                env_type = "WSL" if self.use_wsl else "Local"
                self.logger.info(f"CheckM2 available ({env_type}): {version}")
                return True
            else:
                # Try fallback help command
                if self.use_wsl:
                    result = subprocess.run(
                        ['wsl', 'checkm2', '--help'],
                        capture_output=True,
                        text=True,
                        timeout=30
                    )
                else:
                    result = subprocess.run(
                        [self.checkm2_path, '--help'],
                        capture_output=True,
                        text=True,
                        timeout=30
                    )

                if result.returncode == 0 and 'checkm2' in result.stdout.lower():
                    env_type = "WSL" if self.use_wsl else "Local"
                    self.logger.info(f"CheckM2 available ({env_type}) (version unknown)")
                    return True
                else:
                    self.logger.error(f"CheckM2 not available: {result.stderr}")
                    return False
                
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            self.logger.error(f"Error checking CheckM2 availability: {str(e)}")
            return False
    
    def analyze_genome_quality(self,
                              input_path: Union[str, Path],
                              output_prefix: Optional[str] = None,
                              input_type: str = 'fasta') -> Dict[str, any]:
        """
        Analyze genome quality using CheckM2

        Args:
            input_path: Path to input file or directory
            output_prefix: Prefix for output files
            input_type: Input type ('fasta' for single file, 'directory' for multiple genomes)

        Returns:
            Analysis results dictionary
        """
        input_path = Path(input_path)

        if not input_path.exists():
            raise FileNotFoundError(f"Input path not found: {input_path}")

        # Set output prefix
        if output_prefix is None:
            output_prefix = input_path.stem if input_path.is_file() else input_path.name

        # Check if input looks like protein sequences (many short sequences)
        # CheckM2 with --genes flag supports protein sequences
        input_is_protein = False
        if input_path.is_file():
            try:
                with open(input_path, 'r') as f:
                    first_lines = []
                    seq_count = 0
                    total_length = 0
                    for line in f:
                        if line.startswith('>'):
                            seq_count += 1
                            if seq_count > 10:
                                break
                        else:
                            total_length += len(line.strip())
                    avg_seq_length = total_length / max(seq_count, 1)

                    # Protein sequences are typically shorter than genome sequences
                    # If average sequence length is < 500 aa, it's likely protein sequences
                    if avg_seq_length < 500:
                        self.logger.info(
                            f"Input appears to be protein sequences (avg length: {avg_seq_length:.1f} aa). "
                            f"Will use CheckM2 --genes flag for protein sequence analysis."
                        )
                        input_is_protein = True
            except Exception as e:
                self.logger.warning(f"Could not determine input type: {e}")

        self.logger.info(f"Starting CheckM2 analysis for {input_path.name}")

        # Run CheckM2 analysis (pass input_is_protein flag)
        checkm2_result = self._run_checkm2_analysis(input_path, output_prefix, input_type, input_is_protein)

        if checkm2_result['status'] == 'success':
            # Process and save results
            analysis_results = self._process_results(checkm2_result, input_path)
            self._save_results(analysis_results, output_prefix)
            return analysis_results
        else:
            raise RuntimeError(f"CheckM2 analysis failed: {checkm2_result['message']}")
    
    def _run_checkm2_analysis(self, input_path: Path, output_prefix: str, input_type: str, input_is_protein: bool = False) -> Dict:
        """
        Run CheckM2 analysis

        Args:
            input_path: Input path
            output_prefix: Output file prefix
            input_type: Input type
            input_is_protein: Whether input is protein sequences (use --genes flag)

        Returns:
            Analysis results dictionary
        """
        # Create output directory
        output_dir = os.path.join(self.results_dir, f"{output_prefix}_checkm2")
        FileUtils.ensure_dir(output_dir)

        # Prepare input for CheckM2 - use a separate input directory
        input_dir = os.path.join(self.results_dir, f"{output_prefix}_checkm2_input")
        FileUtils.ensure_dir(input_dir)

        if input_type == 'fasta' and input_path.is_file():
            # Single genome file - copy to input directory
            if input_is_protein:
                genome_file = os.path.join(input_dir, f"{output_prefix}.faa")
            else:
                genome_file = os.path.join(input_dir, f"{output_prefix}.fna")

            import shutil
            shutil.copy2(input_path, genome_file)
            input_for_checkm2 = os.path.abspath(input_dir)
        elif input_type == 'directory' and input_path.is_dir():
            input_for_checkm2 = str(input_path)
        else:
            raise ValueError(f"Invalid input type {input_type} for path {input_path}")

        # Build CheckM2 command
        cmd = [
            self.checkm2_path, 'predict',
            '--threads', str(self.threads),
            '--input', os.path.abspath(input_for_checkm2),
            '--output-directory', os.path.abspath(output_dir),
            '--force'  # Overwrite existing results
        ]

        # Add --genes flag if input is protein sequences
        if input_is_protein:
            cmd.append('--genes')
            cmd.extend(['-x', 'faa'])  # Protein files use .faa extension
            self.logger.info("Using CheckM2 --genes flag for protein sequence analysis")

        # Add database path if configured - use environment variable
        checkm2_db_path = self.config.get('tools.checkm2.database_path', None)
        if checkm2_db_path:
            os.environ['CHECKM2DB'] = checkm2_db_path
            self.logger.info(f"Using CheckM2 database via env: {checkm2_db_path}")

        try:
            self.logger.info(f"Running CheckM2 command: {' '.join(cmd)}")
            self.logger.info(f"Input directory: {os.path.abspath(input_for_checkm2)}")
            self.logger.info(f"Output directory: {os.path.abspath(output_dir)}")

            # Verify input directory exists and has files
            if not os.path.exists(input_for_checkm2):
                raise FileNotFoundError(f"Input directory does not exist: {input_for_checkm2}")

            input_files = os.listdir(input_for_checkm2)
            if not input_files:
                raise ValueError(f"Input directory is empty: {input_for_checkm2}")

            self.logger.info(f"Found {len(input_files)} files in input directory: {input_files}")

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600  # 1 hour timeout
            )

            # Log output for debugging
            if result.stdout:
                self.logger.debug(f"CheckM2 stdout: {result.stdout}")
            if result.stderr:
                self.logger.debug(f"CheckM2 stderr: {result.stderr}")

            # Check if database is missing
            if result.returncode != 0 and 'database' in result.stderr.lower():
                self.logger.warning("CheckM2 database not found. Please download using: checkm2 database --download")
                return {
                    'status': 'skipped',
                    'message': 'CheckM2 database not found. Please download using: checkm2 database --download',
                    'reason': 'missing_database'
                }

            if result.returncode == 0:
                # Look for results file
                results_file = os.path.join(output_dir, 'quality_report.tsv')
                if not os.path.exists(results_file):
                    # Try alternative file names
                    for alt_name in ['checkm2_results.tsv', 'results.tsv', 'output.tsv']:
                        alt_file = os.path.join(output_dir, alt_name)
                        if os.path.exists(alt_file):
                            results_file = alt_file
                            break
                    else:
                        # Look for any TSV file in the output directory
                        tsv_files = list(Path(output_dir).glob('*.tsv'))
                        if tsv_files:
                            results_file = str(tsv_files[0])
                            self.logger.info(f"Using alternative results file: {results_file}")
                        else:
                            return {
                                'status': 'error',
                                'message': 'CheckM2 completed but no results file found'
                            }
                
                return {
                    'status': 'success',
                    'output_dir': output_dir,
                    'results_file': results_file,
                    'stdout': result.stdout,
                    'stderr': result.stderr
                }
            else:
                return {
                    'status': 'error',
                    'message': f"CheckM2 failed with return code {result.returncode}: {result.stderr}"
                }
                
        except subprocess.TimeoutExpired:
            self.logger.error("CheckM2 analysis timeout")
            return {
                'status': 'error',
                'message': 'CheckM2 analysis timeout (1 hour)'
            }
        except Exception as e:
            self.logger.error(f"Error running CheckM2 analysis: {str(e)}")
            return {
                'status': 'error',
                'message': f'Error running CheckM2: {str(e)}'
            }
    
    def _process_results(self, checkm2_result: Dict, input_path: Path) -> Dict:
        """
        Process CheckM2 analysis results
        
        Args:
            checkm2_result: Raw CheckM2 results
            input_path: Input path
            
        Returns:
            Processed analysis results
        """
        results_file = checkm2_result['results_file']
        
        try:
            # Parse results file
            parsed_results = self._parse_results_file(results_file)
            
            # Build analysis results
            analysis_results = {
                'input_path': str(input_path),
                'analysis_timestamp': pd.Timestamp.now().isoformat(),
                'tool': 'CheckM2',
                'results': parsed_results,
                'output_dir': checkm2_result['output_dir'],
                'summary': self._generate_summary(parsed_results)
            }
            
            return analysis_results
            
        except Exception as e:
            self.logger.error(f"Error processing CheckM2 results: {str(e)}")
            raise
    
    def _parse_results_file(self, results_file: str) -> List[Dict]:
        """
        Parse CheckM2 results file
        
        Args:
            results_file: Path to results file
            
        Returns:
            List of genome quality results
        """
        try:
            # Read TSV file
            df = pd.read_csv(results_file, sep='\t')
            
            self.logger.debug(f"CheckM2 results columns: {df.columns.tolist()}")
            
            # Convert to list of dictionaries
            results = []
            for _, row in df.iterrows():
                genome_result = {
                    'genome_name': row.get('Name', row.get('Genome', 'Unknown')),
                    'completeness': float(row.get('Completeness', row.get('completeness', 0))),
                    'contamination': float(row.get('Contamination', row.get('contamination', 0))),
                    'quality_score': self._calculate_quality_score(
                        float(row.get('Completeness', row.get('completeness', 0))),
                        float(row.get('Contamination', row.get('contamination', 0)))
                    ),
                    'quality_category': self._categorize_quality(
                        float(row.get('Completeness', row.get('completeness', 0))),
                        float(row.get('Contamination', row.get('contamination', 0)))
                    )
                }
                
                # Add additional fields if available
                for col in df.columns:
                    if col.lower() not in ['name', 'genome', 'completeness', 'contamination']:
                        genome_result[col.lower()] = row[col]
                
                results.append(genome_result)
            
            return results
            
        except Exception as e:
            self.logger.error(f"Error parsing CheckM2 results file: {str(e)}")
            raise
    
    def _calculate_quality_score(self, completeness: float, contamination: float) -> float:
        """
        Calculate genome quality score
        
        Args:
            completeness: Genome completeness percentage
            contamination: Genome contamination percentage
            
        Returns:
            Quality score (0-100)
        """
        # Simple quality score: completeness - 5 * contamination
        score = completeness - (5 * contamination)
        return max(0, min(100, score))
    
    def _categorize_quality(self, completeness: float, contamination: float) -> str:
        """
        Categorize genome quality
        
        Args:
            completeness: Genome completeness percentage
            contamination: Genome contamination percentage
            
        Returns:
            Quality category string
        """
        if (completeness >= self.high_quality_threshold['completeness'] and 
            contamination <= self.high_quality_threshold['contamination']):
            return 'high_quality'
        elif (completeness >= self.medium_quality_threshold['completeness'] and 
              contamination <= self.medium_quality_threshold['contamination']):
            return 'medium_quality'
        else:
            return 'low_quality'
    
    def _generate_summary(self, results: List[Dict]) -> Dict:
        """
        Generate summary statistics
        
        Args:
            results: List of genome results
            
        Returns:
            Summary statistics dictionary
        """
        if not results:
            return {
                'total_genomes': 0,
                'high_quality_genomes': 0,
                'medium_quality_genomes': 0,
                'low_quality_genomes': 0,
                'average_completeness': 0,
                'average_contamination': 0,
                'average_quality_score': 0
            }
        
        # Count quality categories
        quality_counts = {'high_quality': 0, 'medium_quality': 0, 'low_quality': 0}
        for result in results:
            quality_counts[result['quality_category']] += 1
        
        # Calculate averages
        avg_completeness = sum(r['completeness'] for r in results) / len(results)
        avg_contamination = sum(r['contamination'] for r in results) / len(results)
        avg_quality_score = sum(r['quality_score'] for r in results) / len(results)
        
        return {
            'total_genomes': len(results),
            'high_quality_genomes': quality_counts['high_quality'],
            'medium_quality_genomes': quality_counts['medium_quality'],
            'low_quality_genomes': quality_counts['low_quality'],
            'average_completeness': round(avg_completeness, 2),
            'average_contamination': round(avg_contamination, 2),
            'average_quality_score': round(avg_quality_score, 2),
            'cultivability_assessment': self._assess_cultivability(results)
        }
    
    def _assess_cultivability(self, results: List[Dict]) -> Dict:
        """
        Assess cultivability based on genome quality
        
        Args:
            results: List of genome results
            
        Returns:
            Cultivability assessment dictionary
        """
        if not results:
            return {
                'overall_cultivability': 'unknown',
                'confidence': 0,
                'recommendation': 'No genome data available'
            }
        
        # Calculate cultivability score based on quality
        high_quality_count = sum(1 for r in results if r['quality_category'] == 'high_quality')
        medium_quality_count = sum(1 for r in results if r['quality_category'] == 'medium_quality')
        total_genomes = len(results)
        
        high_quality_ratio = high_quality_count / total_genomes
        medium_quality_ratio = medium_quality_count / total_genomes
        
        if high_quality_ratio >= 0.7:
            cultivability = 'high'
            confidence = 0.9
            recommendation = 'Excellent candidates for cultivation'
        elif high_quality_ratio >= 0.3 or medium_quality_ratio >= 0.7:
            cultivability = 'medium'
            confidence = 0.7
            recommendation = 'Good candidates for cultivation with optimization'
        else:
            cultivability = 'low'
            confidence = 0.5
            recommendation = 'May require significant optimization for cultivation'
        
        return {
            'overall_cultivability': cultivability,
            'confidence': confidence,
            'recommendation': recommendation,
            'high_quality_ratio': round(high_quality_ratio, 2),
            'medium_quality_ratio': round(medium_quality_ratio, 2)
        }
    
    def _save_results(self, results: Dict, output_prefix: str):
        """
        Save analysis results to files
        
        Args:
            results: Analysis results
            output_prefix: Output file prefix
        """
        # Save JSON results
        json_file = os.path.join(self.results_dir, f"{output_prefix}_checkm2_results.json")
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        
        # Save CSV summary
        csv_file = os.path.join(self.results_dir, f"{output_prefix}_checkm2_summary.csv")
        
        if results['results']:
            df_results = pd.DataFrame(results['results'])
            df_results.to_csv(csv_file, index=False)
        
        self.logger.info(f"Results saved to {json_file} and {csv_file}")
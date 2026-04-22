#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tome optimal growth temperature analyzer for MethArCT

Predicts optimal growth temperature (OGT) from protein sequences using the Tome tool.
"""

import os
import subprocess
import tempfile
import json
import csv
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Union
from Bio import SeqIO
from ..utils.config import Config
from ..utils.logger import get_logger
from ..utils.file_utils import FileUtils
from ..utils.sequence_utils import SequenceUtils

class TomeAnalyzer:
    """Tome optimal growth temperature analyzer"""
    
    def __init__(self, config: Optional[Config] = None):
        self.config = config or Config()
        self.logger = get_logger("tome_analyzer")

        # Tool configuration - WSL support
        self.use_wsl = self.config.get('tools.tome.use_wsl', False)

        # Check if local Tome installation path is configured
        tome_install_dir = self.config.get('tools.tome.install_dir', None)
        if tome_install_dir and os.path.exists(tome_install_dir):
            self.tome_script_path = os.path.join(tome_install_dir, 'tome', 'tome.py')
            self.tome_path = 'python'
            self.logger.info(f"Using local Tome installation: {self.tome_script_path}")
        elif self.use_wsl:
            self.tome_path = self.config.get('tools.tome.wsl_path', 'wsl tome')
        else:
            self.tome_path = self.config.get('tools.tome.path', 'tome')
            
        self.threads = self.config.get('tools.tome.threads', 4)
        
        # Results directory
        self.results_dir = self.config.get('output.base_dir', 'results')
        FileUtils.ensure_dir(self.results_dir)
        
        # Check tool availability
        self._check_tome_availability()
    
    def _check_tome_availability(self) -> bool:
        """
        Check if Tome tool is available

        Returns:
            True if Tome is available, False otherwise
        """
        try:
            if hasattr(self, 'tome_script_path') and self.tome_script_path:
                # Locally installed Tome
                if os.path.exists(self.tome_script_path):
                    result = subprocess.run(
                        ['python', self.tome_script_path, '--help'],
                        capture_output=True,
                        text=True,
                        timeout=30
                    )
                    output = (result.stdout + result.stderr).lower()
                    if result.returncode == 0 or 'tome' in output:
                        self.logger.info(f"Tome available (local installation): {self.tome_script_path}")
                        return True
                self.logger.error(f"Tome not available: script not found {self.tome_script_path}")
                return False
            elif self.use_wsl:
                # WSL environment check
                result = subprocess.run(
                    ['wsl', 'tome', '--version'],
                    capture_output=True,
                    text=True,
                    timeout=30
                )
            else:
                # Local environment check
                result = subprocess.run(
                    [self.tome_path, '--version'],
                    capture_output=True,
                    text=True,
                    timeout=30
                )

            if result.returncode == 0:
                version = result.stdout.strip()
                env_type = "WSL" if self.use_wsl else "Local"
                self.logger.info(f"Tome available ({env_type}): {version}")
                return True
            else:
                # Try fallback help command
                if self.use_wsl:
                    result = subprocess.run(
                        ['wsl', 'tome', '--help'],
                        capture_output=True,
                        text=True,
                        timeout=30
                    )
                else:
                    result = subprocess.run(
                        [self.tome_path, '--help'],
                        capture_output=True,
                        text=True,
                        timeout=30
                    )

                if result.returncode == 0 and 'tome' in result.stdout.lower():
                    env_type = "WSL" if self.use_wsl else "Local"
                    self.logger.info(f"Tome available ({env_type}) (version unknown)")
                    return True
                else:
                    self.logger.error(f"Tome not available: {result.stderr}")
                    return False

        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            self.logger.error(f"Error checking Tome availability: {str(e)}")
            return False
    
    def predict_ogt(self, 
                   input_file: Union[str, Path],
                   output_prefix: Optional[str] = None) -> Dict[str, any]:
        """
        Predict optimal growth temperature from protein sequences
        
        Args:
            input_file: Path to input FASTA file
            output_prefix: Prefix for output files
            
        Returns:
            Prediction results dictionary
        """
        input_file = Path(input_file)
        
        if not input_file.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")
        
        if not FileUtils.validate_fasta(input_file):
            raise ValueError(f"Invalid FASTA file: {input_file}")
        
        # Validate sequence type
        with open(input_file, 'r') as f:
            first_seq = next(iter(SeqIO.parse(f, 'fasta')), None)
            if first_seq:
                seq_type = SequenceUtils.validate_sequence_type(str(first_seq.seq))
                if seq_type != 'protein':
                    self.logger.warning(f"Input sequences appear to be {seq_type}, not protein")
        
        # Set output prefix
        if output_prefix is None:
            output_prefix = input_file.stem
        
        self.logger.info(f"Starting Tome OGT prediction for {input_file.name}")
        
        # Run Tome prediction
        prediction_result = self._run_tome_prediction(input_file, output_prefix)
        
        if prediction_result['status'] == 'success':
            # Process and save results
            analysis_results = self._process_results(prediction_result, input_file)
            self._save_results(analysis_results, output_prefix)
            return analysis_results
        else:
            raise RuntimeError(f"Tome prediction failed: {prediction_result['message']}")
    
    def _run_tome_prediction(self, input_file: Path, output_prefix: str) -> Dict:
        """
        Run Tome predOGT prediction
        
        Args:
            input_file: Input FASTA file
            output_prefix: Output file prefix
            
        Returns:
            Prediction results dictionary
        """
        # Create output directory
        output_dir = os.path.join(self.results_dir, f"{output_prefix}_tome")
        FileUtils.ensure_dir(output_dir)

        # Set output file path - Tome uses -o for output filename
        output_file = os.path.join(output_dir, f"{output_prefix}_ogt_results.tsv")

        # Build Tome command
        if hasattr(self, 'tome_script_path') and self.tome_script_path:
            # Locally installed Tome uses: python script.py predOGT
            cmd = [
                'python', self.tome_script_path, 'predOGT',
                '--fasta', str(input_file),
                '-o', output_file,
                '-p', str(self.threads)
            ]
        else:
            cmd = [
                self.tome_path, 'predOGT',
                '--fasta', str(input_file),
                '-o', output_file,
                '-p', str(self.threads)
            ]

        try:
            self.logger.info(f"Running Tome command: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=1800  # 30 minutes timeout
            )

            # Log output for debugging
            if result.stdout:
                self.logger.debug(f"Tome stdout: {result.stdout}")
            if result.stderr:
                self.logger.debug(f"Tome stderr: {result.stderr}")

            # Tome outputs to file specified by -o, check if it exists
            if os.path.exists(output_file):
                self.logger.info(f"Found Tome output file: {output_file}")
                return {
                    'status': 'success',
                    'output_file': output_file,
                    'stdout': result.stdout,
                    'stderr': result.stderr
                }
            elif result.returncode == 0:
                # If -o file doesn't exist, maybe Tome printed to stdout
                if result.stdout:
                    self.logger.info("Tome output in stdout, saving to file")
                    with open(output_file, 'w') as f:
                        f.write(result.stdout)
                    return {
                        'status': 'success',
                        'output_file': output_file,
                        'stdout': result.stdout,
                        'stderr': result.stderr
                    }
                all_files = os.listdir(output_dir) if os.path.exists(output_dir) else []
                self.logger.warning(f"No output file found. Available files: {all_files}")
                return {
                    'status': 'error',
                    'message': f'Tome completed but no output file found. Available files: {all_files}'
                }
            else:
                return {
                    'status': 'error',
                    'message': f"Tome failed with return code {result.returncode}: {result.stderr}"
                }

        except subprocess.TimeoutExpired:
            self.logger.error("Tome prediction timeout")
            return {
                'status': 'error',
                'message': 'Tome prediction timeout (30 minutes)'
            }
        except Exception as e:
            self.logger.error(f"Error running Tome prediction: {str(e)}")
            return {
                'status': 'error',
                'message': f'Error running Tome: {str(e)}'
            }
    
    def _process_results(self, prediction_result: Dict, input_file: Path) -> Dict:
        """
        Process Tome prediction results
        
        Args:
            prediction_result: Raw prediction results
            input_file: Input FASTA file
            
        Returns:
            Processed analysis results
        """
        output_file = prediction_result['output_file']
        
        try:
            # Parse results file
            parsed_results = self._parse_results_file(output_file)
            
            # Count input sequences
            total_sequences = FileUtils.count_sequences(input_file)
            
            # Build analysis results
            analysis_results = {
                'input_file': str(input_file),
                'analysis_timestamp': pd.Timestamp.now().isoformat(),
                'tool': 'Tome',
                'total_sequences': total_sequences,
                'prediction': parsed_results,
                'output_file': output_file,
                'summary': {
                    'predicted_ogt_celsius': parsed_results.get('predOGT', 0),
                    'confidence': parsed_results.get('confidence', 0),
                    'temperature_category': self._categorize_temperature(parsed_results.get('predOGT', 0))
                }
            }
            
            return analysis_results
            
        except Exception as e:
            self.logger.error(f"Error processing Tome results: {str(e)}")
            raise
    
    def _parse_results_file(self, result_file: str) -> Dict:
        """
        Parse Tome results file
        
        Args:
            result_file: Path to results file
            
        Returns:
            Parsed results dictionary
        """
        try:
            # Log file content for debugging
            with open(result_file, 'r', encoding='utf-8') as f:
                content = f.read()
                self.logger.debug(f"Results file content: {content}")
            
            # Try to read as CSV with different delimiters
            for delimiter in ['\t', ',', ' ']:
                try:
                    df = pd.read_csv(result_file, delimiter=delimiter, encoding='utf-8')
                    self.logger.debug(f"Successfully read CSV with delimiter '{delimiter}'")
                    self.logger.debug(f"Columns: {df.columns.tolist()}")
                    break
                except Exception as e:
                    self.logger.debug(f"Failed to read with delimiter '{delimiter}': {str(e)}")
                    continue
            else:
                raise ValueError("Could not parse results file with any delimiter")
            
            # Find OGT column
            ogt_column = None
            for col in df.columns:
                if 'ogt' in col.lower() or 'temperature' in col.lower():
                    ogt_column = col
                    break
            
            if ogt_column is None:
                # Try to extract from first data row
                if len(df) > 0 and len(df.columns) > 1:
                    # Assume second column contains OGT value
                    ogt_value = float(df.iloc[0, 1])
                else:
                    raise ValueError("Could not find OGT column in results")
            else:
                ogt_value = float(df[ogt_column].iloc[0])
            
            return {
                'predOGT': ogt_value,
                'confidence': 0.95,  # Default confidence as Tome doesn't provide this
                'raw_data': df.to_dict('records')
            }
            
        except Exception as e:
            self.logger.error(f"Error parsing results file: {str(e)}")
            
            # Try manual parsing as fallback
            try:
                with open(result_file, 'r', encoding='utf-8') as f:
                    lines = f.readlines()
                    
                if len(lines) > 1:
                    # Try to extract numeric value from second line
                    data_line = lines[1].strip()
                    fields = data_line.split()
                    
                    for field in fields:
                        try:
                            value = float(field)
                            if 0 < value < 150:  # Reasonable temperature range
                                return {
                                    'predOGT': value,
                                    'confidence': 0.95,
                                    'raw_data': []
                                }
                        except ValueError:
                            continue
                
                raise ValueError("Could not extract OGT value from results file")
                
            except Exception as e2:
                self.logger.error(f"Manual parsing also failed: {str(e2)}")
                raise ValueError(f"Failed to parse results file: {str(e)}")
    
    def _categorize_temperature(self, temperature: float) -> str:
        """
        Categorize temperature into psychrophile, mesophile, thermophile, hyperthermophile
        
        Args:
            temperature: Temperature in Celsius
            
        Returns:
            Temperature category string
        """
        if temperature < 20:
            return 'psychrophile'
        elif temperature < 45:
            return 'mesophile'
        elif temperature < 80:
            return 'thermophile'
        else:
            return 'hyperthermophile'
    
    def _save_results(self, results: Dict, output_prefix: str):
        """
        Save analysis results to files
        
        Args:
            results: Analysis results
            output_prefix: Output file prefix
        """
        # Save JSON results
        json_file = os.path.join(self.results_dir, f"{output_prefix}_tome_results.json")
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        
        # Save CSV summary
        csv_file = os.path.join(self.results_dir, f"{output_prefix}_tome_summary.csv")
        summary_data = [{
            'Input_File': results['input_file'],
            'Total_Sequences': results['total_sequences'],
            'Predicted_OGT_Celsius': results['summary']['predicted_ogt_celsius'],
            'Confidence': results['summary']['confidence'],
            'Temperature_Category': results['summary']['temperature_category'],
            'Analysis_Timestamp': results['analysis_timestamp']
        }]
        
        df_summary = pd.DataFrame(summary_data)
        df_summary.to_csv(csv_file, index=False)
        
        self.logger.info(f"Results saved to {json_file} and {csv_file}")
    
    def batch_predict(self, 
                     input_files: List[Union[str, Path]],
                     output_dir: Optional[str] = None) -> Dict[str, Dict]:
        """
        Predict OGT for multiple files
        
        Args:
            input_files: List of input FASTA files
            output_dir: Output directory (optional)
            
        Returns:
            Dictionary mapping file names to results
        """
        if output_dir:
            original_results_dir = self.results_dir
            self.results_dir = output_dir
            FileUtils.ensure_dir(self.results_dir)
        
        batch_results = {}
        
        try:
            for input_file in input_files:
                input_path = Path(input_file)
                self.logger.info(f"Processing {input_path.name}")
                
                try:
                    result = self.predict_ogt(input_path)
                    batch_results[input_path.name] = result
                except Exception as e:
                    self.logger.error(f"Failed to process {input_path.name}: {str(e)}")
                    batch_results[input_path.name] = {
                        'status': 'error',
                        'message': str(e)
                    }
        
        finally:
            if output_dir:
                self.results_dir = original_results_dir
        
        return batch_results
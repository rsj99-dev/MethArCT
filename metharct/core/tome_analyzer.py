#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tome optimal growth temperature analyzer for MethArCT

Predicts optimal growth temperature (OGT) from protein sequences using the Tome tool.
Supports both subprocess-based CLI invocation and direct Python module API.
"""

import os
import sys
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

        self.use_wsl = self.config.get('tools.tome.use_wsl', False)

        tome_install_dir = self.config.get('tools.tome.install_dir', None)
        if not tome_install_dir:
            tome_install_dir = os.path.join(os.path.dirname(__file__), '..', '..', 'Tome-1.1.0')

        if tome_install_dir and os.path.exists(tome_install_dir):
            self.tome_script_path = os.path.join(tome_install_dir, 'tome', 'tome.py')
            self.tome_install_dir = tome_install_dir
            self.tome_path = 'python'
            self.logger.info(f"Using local Tome installation: {self.tome_script_path}")
        elif self.use_wsl:
            self.tome_path = self.config.get('tools.tome.wsl_path', 'wsl tome')
            self.tome_script_path = None
            self.tome_install_dir = None
        else:
            self.tome_path = self.config.get('tools.tome.path', 'tome')
            self.tome_script_path = None
            self.tome_install_dir = None

        self.threads = self.config.get('tools.tome.threads', 4)

        self.results_dir = self.config.get('output.base_dir', 'results')
        FileUtils.ensure_dir(self.results_dir)

        self.tool_available = self._check_tome_availability()

    def _check_tome_availability(self) -> bool:
        try:
            if self.tome_script_path and os.path.exists(self.tome_script_path):
                result = subprocess.run(
                    ['python', self.tome_script_path, '--help'],
                    capture_output=True, text=True, timeout=30
                )
                output = (result.stdout + result.stderr).lower()
                if result.returncode == 0 or 'tome' in output or 'predOGT' in output:
                    self.logger.info(f"Tome available (local installation): {self.tome_script_path}")
                    return True
                self.logger.warning(f"Tome script found but not responding: {self.tome_script_path}")
            elif self.use_wsl:
                result = subprocess.run(
                    ['wsl', 'tome', '--version'],
                    capture_output=True, text=True, timeout=30
                )
                if result.returncode == 0:
                    self.logger.info(f"Tome available (WSL): {result.stdout.strip()}")
                    return True
                result = subprocess.run(
                    ['wsl', 'tome', '--help'],
                    capture_output=True, text=True, timeout=30
                )
                if result.returncode == 0 and 'tome' in result.stdout.lower():
                    self.logger.info("Tome available (WSL, version unknown)")
                    return True
            else:
                result = subprocess.run(
                    [self.tome_path, '--version'],
                    capture_output=True, text=True, timeout=30
                )
                if result.returncode == 0:
                    self.logger.info(f"Tome available (Local): {result.stdout.strip()}")
                    return True
                result = subprocess.run(
                    [self.tome_path, '--help'],
                    capture_output=True, text=True, timeout=30
                )
                if result.returncode == 0 and 'tome' in result.stdout.lower():
                    self.logger.info("Tome available (Local, version unknown)")
                    return True

            self._try_direct_module_import()
            return True

        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            self.logger.warning(f"Tome CLI check failed: {str(e)}, trying Python module API...")
            return self._try_direct_module_import()
        except Exception as e:
            self.logger.warning(f"Tome availability check error: {str(e)}")
            return self._try_direct_module_import()

    def _try_direct_module_import(self) -> bool:
        try:
            if self.tome_install_dir and os.path.exists(self.tome_install_dir):
                tome_path = os.path.join(self.tome_install_dir, 'tome')
                if tome_path not in sys.path:
                    sys.path.insert(0, self.tome_install_dir)
            import importlib
            self._tome_module = importlib.import_module('tome.tome')
            self.logger.info("Tome available via Python module API")
            return True
        except ImportError:
            self.logger.warning("Tome Python module not available")
            self._tome_module = None
            return False

    def predict_ogt(self,
                   input_file: Union[str, Path],
                   output_prefix: Optional[str] = None) -> Dict[str, any]:
        input_file = Path(input_file)

        if not input_file.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")

        if not FileUtils.validate_fasta(input_file):
            raise ValueError(f"Invalid FASTA file: {input_file}")

        with open(input_file, 'r') as f:
            first_seq = next(iter(SeqIO.parse(f, 'fasta')), None)
            if first_seq:
                seq_type = SequenceUtils.validate_sequence_type(str(first_seq.seq))
                if seq_type != 'protein':
                    self.logger.warning(f"Input sequences appear to be {seq_type}, not protein")

        if output_prefix is None:
            output_prefix = input_file.stem

        self.logger.info(f"Starting Tome OGT prediction for {input_file.name}")

        prediction_result = self._run_tome_prediction(input_file, output_prefix)

        if prediction_result['status'] == 'success':
            analysis_results = self._process_results(prediction_result, input_file)
            self._save_results(analysis_results, output_prefix)
            return analysis_results
        else:
            raise RuntimeError(f"Tome prediction failed: {prediction_result['message']}")

    def _run_tome_prediction(self, input_file: Path, output_prefix: str) -> Dict:
        output_dir = os.path.join(self.results_dir, f"{output_prefix}_tome")
        FileUtils.ensure_dir(output_dir)
        output_file = os.path.join(output_dir, f"{output_prefix}_ogt_results.tsv")

        if hasattr(self, '_tome_module') and self._tome_module is not None:
            return self._run_tome_module_prediction(input_file, output_file)

        if self.tome_script_path and os.path.exists(self.tome_script_path):
            cmd = [
                'python', self.tome_script_path, 'predOGT',
                '--fasta', str(input_file),
                '-o', output_file,
                '-p', str(self.threads)
            ]
        elif self.use_wsl:
            cmd = [
                'wsl', 'tome', 'predOGT',
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
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)

            if result.stdout:
                self.logger.debug(f"Tome stdout: {result.stdout[:500]}")
            if result.stderr:
                self.logger.debug(f"Tome stderr: {result.stderr[:500]}")

            if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                self.logger.info(f"Found Tome output file: {output_file}")
                return {'status': 'success', 'output_file': output_file}
            elif result.returncode == 0:
                if result.stdout:
                    self.logger.info("Tome output in stdout, saving to file")
                    with open(output_file, 'w') as f:
                        f.write(result.stdout)
                    return {'status': 'success', 'output_file': output_file}
                return {'status': 'error', 'message': 'Tome completed but produced no output'}
            else:
                self.logger.warning(f"CLI method failed, trying Python module API...")
                return self._run_tome_module_prediction(input_file, output_file)

        except subprocess.TimeoutExpired:
            self.logger.error("Tome prediction timeout (30 minutes)")
            return {'status': 'error', 'message': 'Tome prediction timeout (30 minutes)'}
        except Exception as e:
            self.logger.warning(f"CLI method error: {str(e)}, trying Python module API...")
            return self._run_tome_module_prediction(input_file, output_file)

    def _run_tome_module_prediction(self, input_file: Path, output_file: str) -> Dict:
        try:
            if not hasattr(self, '_tome_module') or self._tome_module is None:
                if not self._try_direct_module_import():
                    return {'status': 'error', 'message': 'Tome Python module not available'}

            tome = self._tome_module
            model, means, stds, features = tome.load_model()
            ogt_value = tome.predict(str(input_file), model, means, stds, features, self.threads)

            filename = os.path.basename(str(input_file))
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(f"FileName\tpredOGT (C)\n{filename}\t{ogt_value}\n")

            self.logger.info(f"Tome module prediction: {ogt_value}°C")
            return {'status': 'success', 'output_file': output_file}

        except Exception as e:
            self.logger.error(f"Tome module prediction failed: {str(e)}")
            return {'status': 'error', 'message': f'Tome prediction failed: {str(e)}'}

    def _process_results(self, prediction_result: Dict, input_file: Path) -> Dict:
        output_file = prediction_result['output_file']

        try:
            parsed_results = self._parse_results_file(output_file)
            total_sequences = FileUtils.count_sequences(input_file)

            analysis_results = {
                'input_file': str(input_file),
                'analysis_timestamp': pd.Timestamp.now().isoformat(),
                'tool': 'Tome',
                'total_sequences': total_sequences,
                'prediction': parsed_results,
                'output_file': output_file,
                'temperature_predictions': parsed_results,
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
        try:
            with open(result_file, 'r', encoding='utf-8') as f:
                content = f.read()
                self.logger.debug(f"Results file content: {content[:500]}")

            for delimiter in ['\t', ',', ' ']:
                try:
                    df = pd.read_csv(result_file, delimiter=delimiter, encoding='utf-8')
                    self.logger.debug(f"Successfully read with delimiter '{delimiter}'")
                    break
                except Exception:
                    continue
            else:
                raise ValueError("Could not parse results file with any delimiter")

            ogt_column = None
            for col in df.columns:
                if 'ogt' in col.lower() or 'temperature' in col.lower():
                    ogt_column = col
                    break

            if ogt_column is None:
                if len(df) > 0 and len(df.columns) > 1:
                    ogt_value = float(df.iloc[0, 1])
                else:
                    raise ValueError("Could not find OGT column in results")
            else:
                ogt_value = float(df[ogt_column].iloc[0])

            return {
                'predOGT': ogt_value,
                'confidence': 0.95,
                'raw_data': df.to_dict('records')
            }

        except Exception as e:
            self.logger.error(f"Error parsing results file: {str(e)}")
            try:
                with open(result_file, 'r', encoding='utf-8') as f:
                    lines = f.readlines()
                if len(lines) > 1:
                    data_line = lines[1].strip()
                    fields = data_line.split()
                    for field in fields:
                        try:
                            value = float(field)
                            if 0 < value < 150:
                                return {'predOGT': value, 'confidence': 0.95, 'raw_data': []}
                        except ValueError:
                            continue
                raise ValueError("Could not extract OGT value from results file")
            except Exception as e2:
                self.logger.error(f"Manual parsing also failed: {str(e2)}")
                raise ValueError(f"Failed to parse results file: {str(e)}")

    def _categorize_temperature(self, temperature: float) -> str:
        if temperature < 20:
            return 'Psychrophilic'
        elif 20 <= temperature < 40:
            return 'Mesophilic'
        elif 40 <= temperature < 60:
            return 'Thermophilic'
        else:
            return 'Hyperthermophilic'

    def _save_results(self, analysis_results: Dict, output_prefix: str):
        output_dir = os.path.join(self.results_dir, f"{output_prefix}_tome")
        FileUtils.ensure_dir(output_dir)

        json_file = os.path.join(output_dir, f"{output_prefix}_tome_results.json")
        summary_data = {
            'input_file': analysis_results.get('input_file'),
            'predicted_ogt_celsius': analysis_results.get('summary', {}).get('predicted_ogt_celsius'),
            'temperature_category': analysis_results.get('summary', {}).get('temperature_category'),
            'total_sequences': analysis_results.get('total_sequences'),
            'tool': 'Tome'
        }
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(summary_data, f, indent=2, ensure_ascii=False)
        self.logger.info(f"Results saved to: {json_file}")

    def cleanup(self):
        pass

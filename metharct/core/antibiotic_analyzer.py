#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Antibiotic resistance prediction module for MethArCT.

Compares a query protein FASTA file against methanogen reference databases
in data/databases/kangshengsu/ and recommends appropriate antibiotics
based on AAI (Average Amino Acid Identity) thresholds.

Rules:
    - Bacitracin: AAI > 50% with Methanobacterium.faa or Methanopyrus.faa
    - Bacitracin + Tunicamycin: AAI > 50% with Methanothrix.faa or Methanomassiliicoccus.faa
    - Vanadate: AAI > 65% with Methanosarcina.faa
    - Tunicamycin: AAI < 50% with ALL kangshengsu reference files
"""

import os
import sys
import tempfile
import shutil
import subprocess
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO

from ..utils.config import Config
from ..utils.logger import get_logger


class AntibioticAnalyzer:
    """Select antibiotics based on AAI comparison with methanogen reference databases.

    Uses DIAMOND blastp to compare query sequences against reference genomes
    of methanogenic archaea and recommends antibiotics based on AAI thresholds.
    """

    # Reference database directory (relative to project data dir)
    DEFAULT_DB_SUBDIR = 'kangshengsu'

    # Antibiotic selection rules:
    ANTIBIOTIC_RULES = [
        {
            'name': 'Tunicamycin',
            'conditions': [
                ('Methanobacterium.faa', 50.0, '>'),
                ('Methanopyrus.faa', 50.0, '>'),
            ],
            'logic': 'any',
        },
        {
            'name': 'Bacitracin + Tunicamycin',
            'conditions': [
                ('Methanothrix.faa', 50.0, '>'),
                ('Methanomassiliicoccus.faa', 50.0, '>'),
            ],
            'logic': 'any',
        },
        {
            'name': 'Vanadate',
            'conditions': [
                ('Methanosarcina.faa', 65.0, '>'),
            ],
            'logic': 'any',
        },
        {
            'name': 'Bacitracin',
            'conditions': [
                ('Methanobacterium.faa', 50.0, '<'),
                ('Methanopyrus.faa', 50.0, '<'),
                ('Methanothrix.faa', 50.0, '<'),
                ('Methanomassiliicoccus.faa', 50.0, '<'),
                ('Methanosarcina.faa', 50.0, '<'),
            ],
            'logic': 'all',
        },
    ]

    def __init__(
        self,
        config: Optional[Config] = None,
        db_dir: Optional[str] = None,
        cpus: int = 1,
        evalue: float = 1e-5,
        tmp_dir: Optional[str] = None,
    ):
        """Initialize AntibioticAnalyzer.

        Parameters
        ----------
        config : Config or None
            MethArCT configuration object. If provided, uses database paths
            and settings from config.
        db_dir : str or None
            Path to the directory containing reference .faa files.
            If None, auto-detects from config or project structure.
        cpus : int
            Number of CPUs for DIAMOND.
        evalue : float
            E-value threshold for DIAMOND search.
        tmp_dir : str or None
            Directory for temporary files.
        """
        self.logger = get_logger("antibiotic_analyzer")
        self.config = config or Config()

        if db_dir is None:
            # Resolve kangshengsu database directory (same approach as CultivationAnalyzer)
            data_base = self.config.get(
                'databases.base_dir', 'data/databases'
            )
            cfg_subdir = self.config.get('databases.kangshengsu', self.DEFAULT_DB_SUBDIR)
            if not os.path.isabs(data_base):
                data_base = os.path.abspath(data_base)
            self.db_dir = os.path.join(data_base, cfg_subdir)
        else:
            self.db_dir = db_dir

        self.cpus = cpus
        self.evalue = evalue
        self.tmp_dir = tmp_dir

    def _validate_inputs(self, query_faa: str):
        """Validate that query file and database directory exist.

        Parameters
        ----------
        query_faa : str
            Path to the query protein FASTA file.

        Raises
        ------
        FileNotFoundError
            If query file or database directory doesn't exist.
        """
        if not os.path.isfile(query_faa):
            raise FileNotFoundError(
                f"Query protein file not found: {query_faa}"
            )

        if not os.path.isdir(self.db_dir):
            raise FileNotFoundError(
                f"Reference database directory not found: {self.db_dir}"
            )

        for rule in self.ANTIBIOTIC_RULES:
            for ref_file, _, _ in rule['conditions']:
                ref_path = os.path.join(self.db_dir, ref_file)
                if not os.path.isfile(ref_path):
                    raise FileNotFoundError(
                        f"Reference file not found: {ref_path}"
                    )

    def _create_diamond_db(self, ref_faa: str, output_dir: str) -> str:
        """Create a DIAMOND protein database from a reference .faa file.

        Parameters
        ----------
        ref_faa : str
            Path to the reference FASTA file.
        output_dir : str
            Directory to store the DIAMOND database.

        Returns
        -------
        str
            Path to the DIAMOND database (without .dmnd extension).
        """
        db_name = os.path.splitext(os.path.basename(ref_faa))[0]
        db_path = os.path.join(output_dir, db_name)

        cmd = [
            'diamond', 'makedb',
            '--in', ref_faa,
            '-d', db_path,
            '--threads', str(self.cpus),
        ]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except FileNotFoundError:
            raise RuntimeError(
                "DIAMOND not found. Please install DIAMOND "
                "and ensure it is in your PATH."
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"diamond makedb failed for {ref_faa}: {e.stderr}"
            )

        return db_path

    def _run_diamond_blastp(
        self, query_faa: str, db_path: str, output_file: str
    ) -> str:
        """Run DIAMOND blastp search of query against a reference database.

        Parameters
        ----------
        query_faa : str
            Path to the query protein FASTA file.
        db_path : str
            Path to the DIAMOND database.
        output_file : str
            Path for the DIAMOND tabular output.

        Returns
        -------
        str
            Path to the output file.
        """
        cmd = [
            'diamond', 'blastp',
            '-q', query_faa,
            '-d', db_path,
            '-o', output_file,
            '-e', str(self.evalue),
            '--max-target-seqs', '500',
            '--threads', str(self.cpus),
            '--outfmt', '6',
            'qseqid', 'sseqid', 'pident', 'length',
            'mismatch', 'gapopen', 'qstart', 'qend',
            'sstart', 'send', 'evalue', 'bitscore',
        ]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except FileNotFoundError:
            raise RuntimeError(
                "DIAMOND not found. Please install DIAMOND "
                "and ensure it is in your PATH."
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"diamond blastp failed: {e.stderr}"
            )

        return output_file

    def _parse_diamond_output(self, diamond_output: str) -> Dict[str, float]:
        """Parse DIAMOND tabular output and extract best hit per query.

        For each query sequence, keeps the hit with the highest bitscore.

        Parameters
        ----------
        diamond_output : str
            Path to DIAMOND tabular output file (outfmt 6).

        Returns
        -------
        dict
            d[query_id] -> best percent identity (float)
        """
        best_hits = {}

        try:
            with open(diamond_output, 'r') as f:
                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) < 12:
                        continue

                    query_id = fields[0]
                    pident = float(fields[2])
                    bitscore = float(fields[11])

                    prev = best_hits.get(query_id)
                    if prev is None or bitscore > prev[1]:
                        best_hits[query_id] = (pident, bitscore)
        except Exception as e:
            self.logger.warning(
                f"Failed to parse DIAMOND output {diamond_output}: {e}"
            )

        return {qid: pident for qid, (pident, _) in best_hits.items()}

    def _calculate_aai(
        self, query_faa: str, ref_faa: str, work_dir: str
    ) -> Tuple[float, int, int]:
        """Calculate AAI between query protein file and a reference file.

        AAI is the mean percent identity of the best DIAMOND hit
        (by bitscore) for each query protein.

        Parameters
        ----------
        query_faa : str
            Path to the query protein FASTA file.
        ref_faa : str
            Path to the reference protein FASTA file.
        work_dir : str
            Working directory for intermediate files.

        Returns
        -------
        tuple
            (aai_value: float, num_hits: int, num_queries: int)
        """
        ref_name = os.path.basename(ref_faa)

        db_path = self._create_diamond_db(ref_faa, work_dir)

        diamond_output = os.path.join(work_dir, f"diamond_{ref_name}.tsv")
        self._run_diamond_blastp(query_faa, db_path, diamond_output)

        best_hits = self._parse_diamond_output(diamond_output)

        num_queries = sum(1 for _ in SeqIO.parse(query_faa, 'fasta'))

        if best_hits:
            aai = sum(best_hits.values()) / len(best_hits)
            return aai, len(best_hits), num_queries
        else:
            return 0.0, 0, num_queries

    def _check_rule(
        self, rule: Dict, aai_results: Dict[str, float]
    ) -> Tuple[bool, List[str]]:
        """Check if an antibiotic rule is satisfied.

        Parameters
        ----------
        rule : dict
            Antibiotic rule definition.
        aai_results : dict
            d[reference_filename] -> aai_value

        Returns
        -------
        tuple
            (satisfied: bool, details: list of str)
        """
        details = []
        conditions_met = []

        for ref_file, threshold, operator in rule['conditions']:
            aai = aai_results.get(ref_file, 0.0)

            if operator == '>':
                met = aai > threshold
                symbol = '>'
            else:
                met = aai < threshold
                symbol = '<'

            conditions_met.append(met)
            status = 'PASS' if met else 'FAIL'
            details.append(
                f"  {status} {ref_file}: AAI = {aai:.2f}% "
                f"({symbol} {threshold}%)"
            )

        if rule['logic'] == 'any':
            satisfied = any(conditions_met)
        else:
            satisfied = all(conditions_met)

        return satisfied, details

    def predict_antibiotics(
        self,
        input_file: str,
        output_prefix: Optional[str] = None,
    ) -> Dict:
        """Run antibiotic resistance prediction.

        Parameters
        ----------
        input_file : str
            Path to the query protein FASTA file.
        output_prefix : str or None
            Prefix for output files. If None, derived from input filename.

        Returns
        -------
        dict
            Results containing:
            - 'aai_results': dict of reference -> AAI values
            - 'recommended_antibiotics': list of recommended antibiotics
            - 'all_rules': detailed results for all rules
            - 'output_file': path to results TSV
            - 'status': 'success' or 'failed'
        """
        self._validate_inputs(input_file)

        if output_prefix is None:
            output_prefix = os.path.splitext(os.path.basename(input_file))[0]

        output_dir = os.path.dirname(os.path.abspath(output_prefix))
        if not output_dir:
            output_dir = '.'
        os.makedirs(output_dir, exist_ok=True)

        if self.tmp_dir:
            work_dir = tempfile.mkdtemp(
                prefix='antibiotic_', dir=self.tmp_dir
            )
        else:
            work_dir = tempfile.mkdtemp(prefix='antibiotic_')

        try:
            self.logger.info(f"Query file: {input_file}")
            self.logger.info(f"Reference database: {self.db_dir}")
            self.logger.info(
                "Calculating AAI against reference databases..."
            )

            unique_refs = set()
            for rule in self.ANTIBIOTIC_RULES:
                for ref_file, _, _ in rule['conditions']:
                    unique_refs.add(ref_file)

            aai_results = {}
            aai_details = {}
            for ref_file in sorted(unique_refs):
                ref_path = os.path.join(self.db_dir, ref_file)
                self.logger.info(f"  Comparing against {ref_file}...")

                aai, num_hits, num_queries = self._calculate_aai(
                    input_file, ref_path, work_dir
                )
                aai_results[ref_file] = aai
                aai_details[ref_file] = {
                    'aai': aai,
                    'num_hits': num_hits,
                    'num_queries': num_queries,
                }

                self.logger.info(
                    f"    AAI = {aai:.2f}% "
                    f"({num_hits}/{num_queries} query proteins with hits)"
                )

            self.logger.info("\nEvaluating antibiotic selection rules...")
            recommended = []
            all_rule_results = []

            for rule in self.ANTIBIOTIC_RULES:
                satisfied, details = self._check_rule(rule, aai_results)
                rule_result = {
                    'name': rule['name'],
                    'satisfied': satisfied,
                    'details': details,
                    'logic': rule['logic'],
                }
                all_rule_results.append(rule_result)

                status = 'MATCHED' if satisfied else 'NOT MATCHED'
                self.logger.info(f"  [{status}] {rule['name']}")
                for detail in details:
                    self.logger.info(f"    {detail}")

                if satisfied:
                    recommended.append(rule['name'])

            output_file = os.path.join(
                output_dir, 'antibiotic_selection_results.tsv'
            )
            self._write_results(
                output_file, aai_results, all_rule_results, recommended
            )
            self.logger.info(f"\nResults written to: {output_file}")

            return {
                'status': 'success',
                'aai_results': aai_results,
                'aai_details': aai_details,
                'recommended_antibiotics': recommended,
                'all_rules': all_rule_results,
                'output_file': output_file,
                'summary': {
                    'total_references': len(aai_results),
                    'recommended_count': len(recommended),
                    'recommended_antibiotics': recommended,
                },
            }

        except Exception as e:
            self.logger.error(f"Antibiotic prediction failed: {str(e)}")
            return {
                'status': 'failed',
                'error': str(e),
            }

        finally:
            if os.path.exists(work_dir):
                shutil.rmtree(work_dir)

    def _write_results(
        self,
        output_file: str,
        aai_results: Dict[str, float],
        all_rule_results: List[Dict],
        recommended: List[str],
    ):
        """Write analysis results to a TSV file.

        Parameters
        ----------
        output_file : str
            Path to the output file.
        aai_results : dict
            AAI values per reference file.
        all_rule_results : list
            Detailed results for all rules.
        recommended : list
            List of recommended antibiotics.
        """
        with open(output_file, 'w', encoding='utf-8') as fout:
            fout.write("Reference\tAAI (%)\n")
            for ref_file in sorted(aai_results.keys()):
                fout.write(f"{ref_file}\t{aai_results[ref_file]:.2f}\n")
            fout.write("\n")

            fout.write("Rule\tStatus\n")
            for rule_result in all_rule_results:
                status = (
                    'MATCHED' if rule_result['satisfied']
                    else 'NOT MATCHED'
                )
                fout.write(f"{rule_result['name']}\t{status}\n")
            fout.write("\n")

            fout.write("Recommended Antibiotics\n")
            if recommended:
                for antibiotic in recommended:
                    fout.write(f"{antibiotic}\n")
            else:
                fout.write("None\n")

    def print_summary(self, results: Dict):
        """Print a concise summary to stdout.

        Parameters
        ----------
        results : dict
            Results from the predict_antibiotics() method.
        """
        print("")
        print("  [Antibiotic Resistance Prediction Result]")
        if results.get('recommended_antibiotics'):
            for antibiotic in results['recommended_antibiotics']:
                print(f"    -> {antibiotic}")
        else:
            print("    No matching antibiotic recommendations")

        print("")
        print(f"  Output: {results.get('output_file', 'N/A')}")
        print("")

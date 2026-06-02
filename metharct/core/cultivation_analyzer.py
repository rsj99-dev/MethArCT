# -*- coding: utf-8 -*-
"""
Cultivation Analyzer Module

This module provides cultivation analysis functionality for MethArCT.
It integrates the new cultivation analysis package to assess microbial
cultivability based on metabolic pathway analysis.
"""

import os
import subprocess
import tempfile
import pandas as pd
from typing import Dict, List, Any, Tuple, Optional
import re
import json
import concurrent.futures
import threading
from datetime import datetime
import time
import glob

class CultivationAnalyzer:
    """
    Cultivation Analyzer - Microbial cultivability assessment based on metabolic pathway analysis
    
    This class uses Diamond sequence alignment and amino acid metabolic pathway analysis to assess microbial cultivability,
    focusing on amino acid, vitamin, coenzyme, nucleotide synthesis and ATPase metabolic pathway analysis.
    """
    
    # Chinese to English pathway name mapping
    PATHWAY_NAME_MAP = {
        # Amino acid biosynthesis (path.aa)
        '半胱氨酸生物合成-蛋氨酸': 'Cysteine biosynthesis-methionine',
        '半胱氨酸生物合成-丝氨酸': 'Cysteine biosynthesis-serine',
        '苯丙氨酸生物合成-chorisate': 'Phenylalanine biosynthesis-chorismate',
        '苯丙氨酸生物合成-chorismate': 'Phenylalanine biosynthesis-chorismate',
        '蛋氨酸生物合成-天冬氨酸': 'Methionine biosynthesis-aspartate',
        '脯氨酸生物合成-谷氨酸': 'Proline biosynthesis-glutamate',
        '精氨酸生物合成-谷氨酸': 'Arginine biosynthesis-glutamate',
        '精氨酸生物合成-鸟氨酸': 'Arginine biosynthesis-ornithine',
        '赖氨酸生物合成-琥珀酰DAP途径': 'Lysine biosynthesis-succinyl-DAP pathway',
        '赖氨酸生物合成-氧代戊二酸': 'Lysine biosynthesis-oxoglutarate',
        '赖氨酸生物合成-乙酰DAP途径': 'Lysine biosynthesis-acetyl-DAP pathway',
        '赖氨酸生物合成-AAA途径': 'Lysine biosynthesis-AAA pathway',
        '赖氨酸生物合成-DAP脱氢酶途径': 'Lysine biosynthesis-DAP dehydrogenase pathway',
        '赖氨酸生物合成-DAP转氨酶途径': 'Lysine biosynthesis-DAP transaminase pathway',
        '酪氨酸生物合成-胆碱': 'Tyrosine biosynthesis-choline',
        '酪氨酸生物合成-chorismate': 'Tyrosine biosynthesis-chorismate',
        '亮氨酸生物合成-氧代异戊酸酯': 'Leucine biosynthesis-oxoisovalerate',
        '鸟氨酸生物合成-谷氨酸': 'Ornithine biosynthesis-glutamate',
        '色氨酸生物合成-chorismate': 'Tryptophan biosynthesis-chorismate',
        '丝氨酸生物合成-甘油酸': 'Serine biosynthesis-glycerate',
        '苏氨酸生物合成-天冬氨酸': 'Threonine biosynthesis-aspartate',
        '缬氨酸和异亮氨酸生物合成-丙酮酸': 'Valine and isoleucine biosynthesis-pyruvate',
        '异亮氨酸生物合成-丙酮酸': 'Isoleucine biosynthesis-pyruvate',
        '异亮氨酸生物合成-苏氨酸': 'Isoleucine biosynthesis-threonine',
        '组氨酸生物合成-PRPP': 'Histidine biosynthesis-PRPP',
        # Vitamin and coenzyme biosynthesis (path.vc)
        '吡哆醛生物合成-赤藓糖': 'Pyridoxal biosynthesis-erythrose',
        '吡哆醛生物合成-R5P': 'Pyridoxal biosynthesis-R5P',
        '泛酸生物合成-1': 'Pantothenate biosynthesis-1',
        '泛酸生物合成-2': 'Pantothenate biosynthesis-2',
        '泛酸生物合成-3': 'Pantothenate biosynthesis-3',
        '泛酸生物合成-4': 'Pantothenate biosynthesis-4',
        '泛酸生物合成-5': 'Pantothenate biosynthesis-5',
        '泛酸生物合成-6': 'Pantothenate biosynthesis-6',
        '泛酸酯生物合成-精胺': 'Pantothenate biosynthesis-spermine',
        '辅酶A生物合成-1': 'Coenzyme A biosynthesis-1',
        '辅酶A生物合成-古菌': 'Coenzyme A biosynthesis-archaeal',
        '辅酶F430生物合成-硅盐酸': 'Coenzyme F430 biosynthesis-sirohydrochlorin',
        '钴胺素生物合成-需氧': 'Cobalamin biosynthesis-aerobic',
        '钴胺素生物合成-厌氧': 'Cobalamin biosynthesis-anaerobic',
        '钴胺素生物合成': 'Cobalamin biosynthesis',
        '核黄素生物合成-真菌': 'Riboflavin biosynthesis-fungal',
        '核黄素生物合成-GTP': 'Riboflavin biosynthesis-GTP',
        '硫胺素生物合成-甘氨酸': 'Thiamine biosynthesis-glycine',
        '硫胺素生物合成-古菌': 'Thiamine biosynthesis-archaeal',
        '硫胺素生物合成-酪氨酸': 'Thiamine biosynthesis-tyrosine',
        '硫胺素生物合成-植物': 'Thiamine biosynthesis-plant',
        '硫辛酸生物合成-动物和细菌': 'Lipoic acid biosynthesis-animal and bacterial',
        '硫辛酸生物合成-辛酰辅酶A': 'Lipoic acid biosynthesis-octanoyl-CoA',
        '硫辛酸生物合成-真核': 'Lipoic acid biosynthesis-eukaryotic',
        '硫辛酸生物合成-植物与原核': 'Lipoic acid biosynthesis-plant and prokaryotic',
        '钼辅因子生物合成-GTP': 'Molybdenum cofactor biosynthesis-GTP',
        '生物素生物合成-吡美酰ACP或辅酶A': 'Biotin biosynthesis-pimeloyl-ACP or CoA',
        '生物素生物合成-BioI途径': 'Biotin biosynthesis-BioI pathway',
        '生物素生物合成-BioU': 'Biotin biosynthesis-BioU',
        '生物素生物合成-BioW通路': 'Biotin biosynthesis-BioW pathway',
        '四氢生物蝶呤生物合成-GTP': 'Tetrahydrobiopterin biosynthesis-GTP',
        '四氢叶酸生物合成-GTP': 'Tetrahydrofolate biosynthesis-GTP',
        '四氢叶酸生物合成-PTPS介导': 'Tetrahydrofolate biosynthesis-PTPS mediated',
        '四氢叶酸生物合成-ribA介导': 'Tetrahydrofolate biosynthesis-ribA mediated',
        '苏氏四氢生物蝶呤生物合成-GTP': 'Sulfo-tetrahydrobiopterin biosynthesis-GTP',
        'NAD生物合成-1': 'NAD biosynthesis-1',
        'NAD生物合成-2': 'NAD biosynthesis-2',
        'NAD生物合成-色氨酸': 'NAD biosynthesis-tryptophan',
        # Nucleotide biosynthesis (path_hesuan)
        '嘧啶生物合成-1': 'Pyrimidine biosynthesis-1',
        '嘧啶生物合成-2': 'Pyrimidine biosynthesis-2',
        '嘧啶生物合成-3': 'Pyrimidine biosynthesis-3',
        '嘧啶生物合成-4': 'Pyrimidine biosynthesis-4',
        '嘧啶生物合成-5': 'Pyrimidine biosynthesis-5',
        '嘌呤生物合成-1': 'Purine biosynthesis-1',
        '嘌呤生物合成-2': 'Purine biosynthesis-2',
        '嘌呤生物合成-3': 'Purine biosynthesis-3',
        '嘌呤生物合成-4': 'Purine biosynthesis-4',
        '嘌呤生物合成-5': 'Purine biosynthesis-5',
        # ATPase (path_atp)
        'ATP酶-F1': 'ATP synthase-F1',
        'ATP酶-F2': 'ATP synthase-F2',
        'ATP酶-F3': 'ATP synthase-F3',
        'ATP酶-F4': 'ATP synthase-F4',
        'ATP酶-V1': 'ATP synthase-V1',
        'ATP酶-V2': 'ATP synthase-V2',
        'ATP酶-VA1': 'ATP synthase-VA1',
    }
    
    # Mapping of 20 standard amino acids to their biosynthesis pathway names (English)
    # An amino acid is considered "present" if ANY of its pathways has completeness >= threshold
    AMINO_ACID_PATHWAY_MAP = {
        'Alanine (Ala)': {
            'pathways': [],
            'note': 'Synthesized by simple transamination from pyruvate; no dedicated pathway in database'
        },
        'Arginine (Arg)': {
            'pathways': ['Arginine biosynthesis-glutamate', 'Arginine biosynthesis-ornithine'],
            'note': ''
        },
        'Asparagine (Asn)': {
            'pathways': [],
            'note': 'Synthesized from aspartate by asparagine synthetase; no dedicated pathway in database'
        },
        'Aspartic acid (Asp)': {
            'pathways': [],
            'note': 'Central metabolite derived from oxaloacetate; no dedicated pathway in database'
        },
        'Cysteine (Cys)': {
            'pathways': ['Cysteine biosynthesis-serine', 'Cysteine biosynthesis-methionine'],
            'note': ''
        },
        'Glutamic acid (Glu)': {
            'pathways': [],
            'note': 'Central metabolite from alpha-ketoglutarate; no dedicated pathway in database'
        },
        'Glutamine (Gln)': {
            'pathways': [],
            'note': 'Synthesized from glutamate by glutamine synthetase; no dedicated pathway in database'
        },
        'Glycine (Gly)': {
            'pathways': [],
            'note': 'Derived from serine via serine hydroxymethyltransferase; no dedicated pathway in database'
        },
        'Histidine (His)': {
            'pathways': ['Histidine biosynthesis-PRPP'],
            'note': ''
        },
        'Isoleucine (Ile)': {
            'pathways': ['Isoleucine biosynthesis-threonine', 'Isoleucine biosynthesis-pyruvate',
                         'Valine and isoleucine biosynthesis-pyruvate'],
            'note': ''
        },
        'Leucine (Leu)': {
            'pathways': ['Leucine biosynthesis-oxoisovalerate'],
            'note': ''
        },
        'Lysine (Lys)': {
            'pathways': ['Lysine biosynthesis-succinyl-DAP pathway', 'Lysine biosynthesis-oxoglutarate',
                         'Lysine biosynthesis-acetyl-DAP pathway', 'Lysine biosynthesis-AAA pathway',
                         'Lysine biosynthesis-DAP dehydrogenase pathway', 'Lysine biosynthesis-DAP transaminase pathway'],
            'note': ''
        },
        'Methionine (Met)': {
            'pathways': ['Methionine biosynthesis-aspartate'],
            'note': ''
        },
        'Phenylalanine (Phe)': {
            'pathways': ['Phenylalanine biosynthesis-chorismate'],
            'note': ''
        },
        'Proline (Pro)': {
            'pathways': ['Proline biosynthesis-glutamate'],
            'note': ''
        },
        'Serine (Ser)': {
            'pathways': ['Serine biosynthesis-glycerate'],
            'note': ''
        },
        'Threonine (Thr)': {
            'pathways': ['Threonine biosynthesis-aspartate'],
            'note': ''
        },
        'Tryptophan (Trp)': {
            'pathways': ['Tryptophan biosynthesis-chorismate'],
            'note': ''
        },
        'Tyrosine (Tyr)': {
            'pathways': ['Tyrosine biosynthesis-choline', 'Tyrosine biosynthesis-chorismate'],
            'note': ''
        },
        'Valine (Val)': {
            'pathways': ['Valine and isoleucine biosynthesis-pyruvate'],
            'note': ''
        },
    }
    
    def __init__(self, config: Dict[str, Any]):
        """
        Initialize cultivation analyzer
        
        Args:
            config: Configuration dictionary containing tool paths, database paths and other parameters
        """
        # Get parameters from configuration
        self.diamond_path = config.get('tools', {}).get('diamond', {}).get('path', 'diamond.exe')
        self.use_wsl = config.get('tools', {}).get('diamond', {}).get('use_wsl', False)
        self.threads = config.get('tools', {}).get('diamond', {}).get('threads', 4)
        self.evalue = config.get('tools', {}).get('diamond', {}).get('evalue', 1e-5)
        
        # Database paths - support absolute and relative paths
        base_dir = config.get('databases', {}).get('base_dir', 'data/databases')
        
        # If base_dir is relative path, make it relative to current working directory
        if not os.path.isabs(base_dir):
            base_dir = os.path.abspath(base_dir)
            
        self.cultivation_dir = os.path.join(base_dir, 'cultivation')
        
        # Metabolic pathway directories - updated to new database structure
        self.amino_acid_dir = os.path.join(self.cultivation_dir, 'path.aa')
        self.vitamin_dir = os.path.join(self.cultivation_dir, 'path.vc')
        self.nucleotide_dir = os.path.join(self.cultivation_dir, 'path_hesuan')
        self.atp_dir = os.path.join(self.cultivation_dir, 'path_atp')
        
        # Automatically get all metabolic pathway files from directory
        self.amino_acid_pathways = self._discover_pathways(self.amino_acid_dir, "Amino acid")
        self.vitamin_pathways = self._discover_pathways(self.vitamin_dir, "Vitamin and coenzyme")
        self.nucleotide_pathways = self._discover_pathways(self.nucleotide_dir, "Nucleotide synthesis")
        self.atp_pathways = self._discover_pathways(self.atp_dir, "ATPase")
        
        total_pathways = len(self.amino_acid_pathways) + len(self.vitamin_pathways) + len(self.nucleotide_pathways) + len(self.atp_pathways)
        print(f"Initialize cultivation analyzer, will analyze {total_pathways} metabolic pathways")
        print(f"  - Amino acid metabolic pathways: {len(self.amino_acid_pathways)} ")
        print(f"  - Vitamin and coenzyme metabolic pathways: {len(self.vitamin_pathways)} ")
        print(f"  - Nucleotide synthesis pathways: {len(self.nucleotide_pathways)} ")
        print(f"  - ATPase metabolic pathways: {len(self.atp_pathways)} ")
    
    def _discover_pathways(self, directory: str, pathway_type: str) -> Dict[str, str]:
        """
        Discover all metabolic pathways from specified directory
        
        Args:
            directory: Directory path containing FASTA files
            pathway_type: Pathway type description (for log output)
            
        Returns:
            Dictionary {pathway_name: file_path}
        """
        pathways = {}
        
        if not os.path.exists(directory):
            print(f"Warning: {pathway_type} FASTA directory does not exist: {directory}")
            return pathways
            
        # Find all .fasta files
        fasta_files = glob.glob(os.path.join(directory, "*.fasta"))
        
        for file_path in fasta_files:
            file_name = os.path.basename(file_path)
            pathway_name = file_name.replace('.fasta', '')
            # Translate Chinese pathway names to English
            english_name = self.PATHWAY_NAME_MAP.get(pathway_name, pathway_name)
            pathways[english_name] = file_path
            
        return pathways
    
    def _parse_fasta_file(self, fasta_file: str) -> Dict[str, Dict[str, str]]:
        """
        Parse FASTA file, extract genes and sequences
        
        Args:
            fasta_file: FASTA file path
            
        Returns:
            Gene information dictionary {gene_name: {sequence: str, description: str}}
        """
        genes = {}
        
        if not os.path.exists(fasta_file):
            print(f"Warning: FASTA file does not exist: {fasta_file}")
            return genes
            
        with open(fasta_file, 'r') as f:
            current_gene = None
            current_seq = []
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous gene
                    if current_gene:
                        genes[current_gene] = {
                            'sequence': ''.join(current_seq),
                            'description': description
                        }
                    
                    # Parse new gene header
                    header = line[1:]  # Remove '>'
                    
                    # Extract gene name from new format header information
                    # Format like: >eco:b0907 K00831 phosphoserine aminotransferase [EC:2.6.1.52] | (RefSeq) serC; phosphoserine/phosphohydroxythreonine aminotransferase (A)
                    # Extract gene name (e.g., serC)
                    gene_match = re.search(r'\((RefSeq|GenBank)\)\s*([^;]+);', header)
                    if gene_match:
                        current_gene = gene_match.group(2).strip()
                    else:
                        # If not found, try to extract gene part from species:gene format
                        species_gene_match = re.match(r'([^:]+):([^:]+)', header.split()[0])
                        if species_gene_match:
                            current_gene = species_gene_match.group(2)
                        else:
                            # Finally try to extract first word as gene name
                            current_gene = header.split()[0]
                    
                    description = header
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Save last gene
            if current_gene:
                genes[current_gene] = {
                    'sequence': ''.join(current_seq),
                    'description': description
                }
        
        return genes
    
    def _run_diamond_search(self, query_seq: str, genome_file: str, evalue: float = 1e-5) -> Tuple[bool, Dict[str, Any]]:
        """
        Use Diamond for sequence alignment
        
        Args:
            query_seq: Query sequence
            genome_file: Genome protein sequence file path
            evalue: E-value threshold
            
        Returns:
            (Whether match found, match details)
        """
        try:
            # Create temporary query file
            with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.faa') as temp_query:
                temp_query.write(f">query\n{query_seq}\n")
                temp_query_path = temp_query.name
            
            # Create temporary output file
            with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.txt') as temp_output:
                temp_output_path = temp_output.name
            
            # Create temporary database file
            temp_db_path = temp_query_path.replace('.faa', '.dmnd')
            
            # Build Diamond command
            diamond_cmd = self.diamond_path
            if self.use_wsl:
                diamond_cmd = f"wsl {diamond_cmd}"
            
            # Create Diamond database for genome file
            make_db_cmd = [
                diamond_cmd, 'makedb',
                '--in', genome_file,
                '--db', temp_db_path
            ]
            
            db_result = subprocess.run(make_db_cmd, capture_output=True, text=True)
            
            if db_result.returncode != 0:
                # Clean up temporary files
                self._cleanup_temp_files([temp_query_path, temp_output_path, temp_db_path])
                return False, {'found': False, 'error': 'Database creation failed'}
            
            # Run diamond blastp command
            blast_cmd = [
                diamond_cmd, 'blastp',
                '--query', temp_query_path,
                '--db', temp_db_path,
                '--outfmt', '6',  # Table format
                '--evalue', str(evalue),
                '--max-target-seqs', '1',  # Only return best match
                '--out', temp_output_path
            ]
            
            result = subprocess.run(blast_cmd, capture_output=True, text=True)
            
            match_details = {'found': False, 'hits': []}
            
            if result.returncode == 0 and os.path.exists(temp_output_path):
                with open(temp_output_path, 'r') as f:
                    for line in f:
                        if line.strip():
                            parts = line.strip().split('\t')
                            if len(parts) >= 12:
                                hit_info = {
                                    'subject_id': parts[1],
                                    'identity': float(parts[2]),
                                    'alignment_length': int(parts[3]),
                                    'evalue': float(parts[10]),
                                    'bitscore': float(parts[11])
                                }
                                
                                # Set more relaxed conditions because distant species may have lower homology
                                if hit_info['evalue'] <= evalue and hit_info['identity'] >= 20 and hit_info['alignment_length'] >= 30:
                                    match_details['found'] = True
                                    match_details['hits'].append(hit_info)
            
            # Clean up temporary files
            self._cleanup_temp_files([temp_query_path, temp_output_path, temp_db_path])
            
            return match_details['found'], match_details
            
        except Exception as e:
            print(f"Diamond search error: {e}")
            return False, {'found': False, 'error': str(e)}
    
    def _cleanup_temp_files(self, file_paths: List[str]):
        """Clean up temporary files"""
        for file_path in file_paths:
            try:
                if os.path.exists(file_path):
                    os.unlink(file_path)
            except Exception as e:
                print(f"Clean up temporary files failed {file_path}: {e}")
    
    def analyze_genome(self, genome_file: str) -> Dict[str, Dict[str, Any]]:
        """
        Analyze metabolic pathways in genome (amino acid, vitamin, coenzyme, nucleotide synthesis and ATPase metabolism)
        
        Args:
            genome_file: genome file path
            
        Returns:
            Analysis result dictionary
        """
        if not os.path.exists(genome_file):
            print(f"Error: Genome file does not exist: {genome_file}")
            return {}
        
        results = {}
        
        # Merge all pathways
        all_pathways = {**self.amino_acid_pathways, **self.vitamin_pathways, **self.nucleotide_pathways, **self.atp_pathways}
        
        # Analyze each pathway
        for pathway_name, fasta_file in all_pathways.items():
            if not os.path.exists(fasta_file):
                print(f"Warning: FASTA file does not exist: {fasta_file}")
                continue
            
            # Parse FASTA file to get genes
            genes = self._parse_fasta_file(fasta_file)
            
            if not genes:
                print(f"Warning: Cannot parse from {fasta_file} parse genes")
                continue
            
            # Analyze pathway
            pathway_result = {
                'pathway': pathway_name,
                'total_genes': len(genes),
                'found_genes': 0,
                'missing_genes': [],
                'gene_details': {},
                'completeness': 0.0
            }
            
            # Record found genes to avoid duplicate alignment
            found_genes = set()
            
            # Align each gene
            for gene_name, gene_info in genes.items():
                # If gene already found, skip
                if gene_name in found_genes:
                    pathway_result['gene_details'][gene_name] = {
                        'found': True,
                        'match_type': 'duplicate',
                        'description': gene_info['description']
                    }
                    pathway_result['found_genes'] += 1
                    continue
                
                # Use Diamond for sequence alignment
                found, match_details = self._run_diamond_search(
                    gene_info['sequence'], 
                    genome_file
                )
                
                pathway_result['gene_details'][gene_name] = {
                    'found': found,
                    'match_details': match_details,
                    'description': gene_info['description']
                }
                
                if found:
                    found_genes.add(gene_name)
                    pathway_result['found_genes'] += 1
                    print(f"Found {pathway_name} pathway {gene_name} gene")
                else:
                    pathway_result['missing_genes'].append(gene_name)
            
            # Calculate completeness
            if pathway_result['total_genes'] > 0:
                pathway_result['completeness'] = pathway_result['found_genes'] / pathway_result['total_genes']
            
            results[pathway_name] = pathway_result
        
        return results
    
    def generate_report(self, results: Dict[str, Dict[str, Any]]) -> str:
        """
        Generate analysis report
        
        Args:
            results: analysis results
            
        Returns:
            Report text
        """
        report = []
        report.append("Metabolic Pathway Analysis Report")
        report.append("=" * 50)
        report.append("")
        
        # Overall statistics
        total_pathways = len(results)
        complete_pathways = sum(1 for r in results.values() if r['completeness'] == 1.0)
        
        report.append(f"Total pathways: {total_pathways}")
        report.append(f"Complete pathways: {complete_pathways}")
        report.append(f"Incomplete pathways: {total_pathways - complete_pathways}")
        report.append("")
        
        # Classification statistics
        aa_pathways = {k: v for k, v in results.items() if k in self.amino_acid_pathways}
        vc_pathways = {k: v for k, v in results.items() if k in self.vitamin_pathways}
        nt_pathways = {k: v for k, v in results.items() if k in self.nucleotide_pathways}
        atp_pathways = {k: v for k, v in results.items() if k in self.atp_pathways}
        
        aa_complete = sum(1 for r in aa_pathways.values() if r['completeness'] == 1.0)
        vc_complete = sum(1 for r in vc_pathways.values() if r['completeness'] == 1.0)
        nt_complete = sum(1 for r in nt_pathways.values() if r['completeness'] == 1.0)
        atp_complete = sum(1 for r in atp_pathways.values() if r['completeness'] == 1.0)
        
        report.append("Pathway classification statistics:")
        report.append(f"  Amino acid metabolic pathways: {len(aa_pathways)}  (complete: {aa_complete})")
        report.append(f"  Vitamin and coenzyme metabolic pathways: {len(vc_pathways)}  (complete: {vc_complete})")
        report.append(f"  Nucleotide synthesis pathways: {len(nt_pathways)}  (complete: {nt_complete})")
        report.append(f"  ATPase metabolic pathways: {len(atp_pathways)}  (complete: {atp_complete})")
        report.append("")
        
        # Detailed information for each pathway
        report.append("Detailed information for each pathway:")
        report.append("-" * 50)
        
        for pathway, result in results.items():
            if pathway in self.amino_acid_pathways:
                pathway_type = "Amino acid"
            elif pathway in self.vitamin_pathways:
                pathway_type = "Vitamin/coenzyme"
            elif pathway in self.nucleotide_pathways:
                pathway_type = "Nucleotide synthesis"
            else:
                pathway_type = "ATPase"
                
            report.append(f"[{pathway_type}] {pathway} pathway:")
            report.append(f"  Completeness: {result['completeness']:.2%} ({result['found_genes']}/{result['total_genes']})")
            
            if result['missing_genes']:
                report.append(f"  Missing genes: {', '.join(result['missing_genes'])}")
            else:
                report.append("  All genes present")
            
            report.append("")
        
        # Complete pathway summary
        complete_aa_pathways = [pathway for pathway, result in aa_pathways.items() if result['completeness'] == 1.0]
        complete_vc_pathways = [pathway for pathway, result in vc_pathways.items() if result['completeness'] == 1.0]
        complete_nt_pathways = [pathway for pathway, result in nt_pathways.items() if result['completeness'] == 1.0]
        complete_atp_pathways = [pathway for pathway, result in atp_pathways.items() if result['completeness'] == 1.0]
        
        if complete_aa_pathways:
            report.append("Complete amino acid metabolic pathways:")
            report.append(f"  {', '.join(complete_aa_pathways)}")
            report.append("")
        
        if complete_vc_pathways:
            report.append("Complete vitamin and coenzyme metabolic pathways:")
            report.append(f"  {', '.join(complete_vc_pathways)}")
            report.append("")
        
        if complete_nt_pathways:
            report.append("Complete nucleotide synthesis pathways:")
            report.append(f"  {', '.join(complete_nt_pathways)}")
            report.append("")
        
        if complete_atp_pathways:
            report.append("Complete ATPase metabolic pathways:")
            report.append(f"  {', '.join(complete_atp_pathways)}")
            report.append("")
        
        return '\n'.join(report)
    
    def generate_summary_recommendations(self, results: Dict[str, Dict[str, Any]]) -> str:
        """
        Generate summary recommendations
        
        Args:
            results: analysis results
            
        Returns:
            Recommendation text
        """
        # Classification statistics
        aa_pathways = {k: v for k, v in results.items() if k in self.amino_acid_pathways}
        vc_pathways = {k: v for k, v in results.items() if k in self.vitamin_pathways}
        nt_pathways = {k: v for k, v in results.items() if k in self.nucleotide_pathways}
        atp_pathways = {k: v for k, v in results.items() if k in self.atp_pathways}
        
        # Calculate completeness for each pathway type
        aa_complete = sum(1 for r in aa_pathways.values() if r['completeness'] == 1.0)
        vc_complete = sum(1 for r in vc_pathways.values() if r['completeness'] == 1.0)
        nt_complete = sum(1 for r in nt_pathways.values() if r['completeness'] == 1.0)
        atp_complete = sum(1 for r in atp_pathways.values() if r['completeness'] == 1.0)
        
        aa_avg_completeness = sum(r['completeness'] for r in aa_pathways.values()) / len(aa_pathways) if aa_pathways else 0
        vc_avg_completeness = sum(r['completeness'] for r in vc_pathways.values()) / len(vc_pathways) if vc_pathways else 0
        nt_avg_completeness = sum(r['completeness'] for r in nt_pathways.values()) / len(nt_pathways) if nt_pathways else 0
        atp_avg_completeness = sum(r['completeness'] for r in atp_pathways.values()) / len(atp_pathways) if atp_pathways else 0
        
        # Generate recommendations
        recommendations = []
        recommendations.append("Cultivability Analysis Recommendations")
        recommendations.append("=" * 50)
        recommendations.append("")
        
        # Amino acid metabolism capability assessment
        if aa_avg_completeness >= 0.8:
            recommendations.append("This microorganism has strong amino acid synthesis capability, may be able to synthesize most essential amino acids independently, low cultivation difficulty.")
        elif aa_avg_completeness >= 0.5:
            recommendations.append("This microorganism has some amino acid synthesis capability, but may need supplementation with some essential amino acids.")
        else:
            recommendations.append("This microorganism has weak amino acid synthesis capability, it is recommended to supplement multiple amino acids in the culture medium.")
        
        # Vitamin and coenzyme metabolism capability assessment
        if vc_avg_completeness >= 0.8:
            recommendations.append("This microorganism has strong vitamin and coenzyme synthesis capability, may not need additional vitamin supplementation.")
        elif vc_avg_completeness >= 0.5:
            recommendations.append("This microorganism has some vitamin and coenzyme synthesis capability, but may need supplementation with some vitamins.")
        else:
            recommendations.append("This microorganism has weak vitamin and coenzyme synthesis capability, it is recommended to supplement multiple vitamins in the culture medium.")
        
        # Nucleotide synthesis capability assessment
        if nt_avg_completeness >= 0.8:
            recommendations.append("This microorganism has strong nucleotide synthesis capability, can independently synthesize purines and pyrimidines.")
        elif nt_avg_completeness >= 0.5:
            recommendations.append("This microorganism has some nucleotide synthesis capability, but may need supplementation with some nucleotide precursors.")
        else:
            recommendations.append("This microorganism has weak nucleotide synthesis capability, it is recommended to supplement nucleotides or their precursors in the culture medium.")
        
        # ATPase metabolism capability assessment
        if atp_avg_completeness >= 0.8:
            recommendations.append("This microorganism has complete ATPase system, strong energy metabolism capability.")
        elif atp_avg_completeness >= 0.5:
            recommendations.append("This microorganism has partially complete ATPase system, moderate energy metabolism capability.")
        else:
            recommendations.append("This microorganism has incomplete ATPase system, energy metabolism may be limited, it is recommended to optimize culture conditions.")
        
        recommendations.append("")
        
        # Overall cultivability assessment
        total_avg_completeness = (aa_avg_completeness + vc_avg_completeness + nt_avg_completeness + atp_avg_completeness) / 4
        
        if total_avg_completeness >= 0.8:
            recommendations.append("Overall assessment: This microorganism is relatively easy to culture, may have broad nutritional adaptability.")
        elif total_avg_completeness >= 0.5:
            recommendations.append("Overall assessment: This microorganism has moderate culture difficulty, requires appropriate nutritional supplementation.")
        else:
            recommendations.append("Overall assessment: This microorganism is difficult to culture, requires complex nutritional formulation and optimized culture conditions.")
        
        return '\n'.join(recommendations)
    
    def get_summary(self, results: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """
        Get analysis result summary
        
        Args:
            results: analysis results
            
        Returns:
            Summary dictionary
        """
        # Classification statistics
        aa_pathways = {k: v for k, v in results.items() if k in self.amino_acid_pathways}
        vc_pathways = {k: v for k, v in results.items() if k in self.vitamin_pathways}
        nt_pathways = {k: v for k, v in results.items() if k in self.nucleotide_pathways}
        atp_pathways = {k: v for k, v in results.items() if k in self.atp_pathways}
        
        # Calculate completeness for each pathway type
        aa_complete = sum(1 for r in aa_pathways.values() if r['completeness'] == 1.0)
        vc_complete = sum(1 for r in vc_pathways.values() if r['completeness'] == 1.0)
        nt_complete = sum(1 for r in nt_pathways.values() if r['completeness'] == 1.0)
        atp_complete = sum(1 for r in atp_pathways.values() if r['completeness'] == 1.0)
        
        aa_avg_completeness = sum(r['completeness'] for r in aa_pathways.values()) / len(aa_pathways) if aa_pathways else 0
        vc_avg_completeness = sum(r['completeness'] for r in vc_pathways.values()) / len(vc_pathways) if vc_pathways else 0
        nt_avg_completeness = sum(r['completeness'] for r in nt_pathways.values()) / len(nt_pathways) if nt_pathways else 0
        atp_avg_completeness = sum(r['completeness'] for r in atp_pathways.values()) / len(atp_pathways) if atp_pathways else 0
        
        total_pathways = len(results)
        complete_pathways = sum(1 for r in results.values() if r['completeness'] == 1.0)
        avg_completeness = sum(r['completeness'] for r in results.values()) / total_pathways if results else 0
        
        return {
            'total_pathways': total_pathways,
            'complete_pathways': complete_pathways,
            'avg_completeness': avg_completeness,
            'amino_acid': {
                'total': len(aa_pathways),
                'complete': aa_complete,
                'avg_completeness': aa_avg_completeness
            },
            'vitamin_coenzyme': {
                'total': len(vc_pathways),
                'complete': vc_complete,
                'avg_completeness': vc_avg_completeness
            },
            'nucleotide': {
                'total': len(nt_pathways),
                'complete': nt_complete,
                'avg_completeness': nt_avg_completeness
            },
            'atp': {
                'total': len(atp_pathways),
                'complete': atp_complete,
                'avg_completeness': atp_avg_completeness
            }
        }
    
    def export_to_csv(self, results: Dict[str, Dict[str, Any]], output_file: str):
        """
        Export analysis results to CSV file
        
        Args:
            results: analysis results
            output_file: output file path
        """
        rows = []
        
        for pathway_name, pathway_result in results.items():
            # Determine pathway type
            if pathway_name in self.amino_acid_pathways:
                pathway_type = "Amino acid"
            elif pathway_name in self.vitamin_pathways:
                pathway_type = "Vitamin/coenzyme"
            elif pathway_name in self.nucleotide_pathways:
                pathway_type = "Nucleotide synthesis"
            else:
                pathway_type = "ATPase"
            
            # Add pathway overall information
            rows.append({
                'Pathway_Type': pathway_type,
                'Pathway_Name': pathway_name,
                'Gene_Name': 'Overall',
                'Found': pathway_result['found_genes'],
                'Total': pathway_result['total_genes'],
                'Completeness': pathway_result['completeness'],
                'Description': f"Completeness: {pathway_result['completeness']:.2%}"
            })
            
            # Add detailed information for each gene
            for gene_name, gene_detail in pathway_result['gene_details'].items():
                rows.append({
                    'Pathway_Type': pathway_type,
                    'Pathway_Name': pathway_name,
                    'Gene_Name': gene_name,
                    'Found': 1 if gene_detail['found'] else 0,
                    'Total': 1,
                    'Completeness': 1.0 if gene_detail['found'] else 0.0,
                    'Description': gene_detail['description']
                })
        
        # Create DataFrame and save as CSV
        df = pd.DataFrame(rows)
        df.to_csv(output_file, index=False)
        print(f"Results exported to: {output_file}")
    
    def generate_cultivability_assessment(self, results: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """
        Generate cultivability assessment based on metabolic pathway analysis
        
        Args:
            results: analysis results
            
        Returns:
            Cultivability assessment dictionary
        """
        summary = self.get_summary(results)
        
        total_avg = (
            summary['amino_acid']['avg_completeness'] +
            summary['vitamin_coenzyme']['avg_completeness'] +
            summary['nucleotide']['avg_completeness'] +
            summary['atp']['avg_completeness']
        ) / 4
        
        if total_avg >= 0.8:
            cultivability = 'high'
            description = 'Organism has strong metabolic synthesis capability'
        elif total_avg >= 0.5:
            cultivability = 'medium'
            description = 'Organism has moderate metabolic synthesis capability'
        else:
            cultivability = 'low'
            description = 'Organism has limited metabolic synthesis capability'
        
        return {
            'overall_cultivability': cultivability,
            'description': description,
            'average_pathway_completeness': round(total_avg, 4),
            'amino_acid_completeness': round(summary['amino_acid']['avg_completeness'], 4),
            'vitamin_coenzyme_completeness': round(summary['vitamin_coenzyme']['avg_completeness'], 4),
            'nucleotide_completeness': round(summary['nucleotide']['avg_completeness'], 4),
            'atp_completeness': round(summary['atp']['avg_completeness'], 4),
            'confidence': 0.85 if summary['total_pathways'] > 10 else 0.6
        }
    
    def export_results(self, results: Dict[str, Dict[str, Any]], 
                       cultivability_assessment: Dict[str, Any],
                       output_file: str, format: str = 'csv'):
        """
        Export cultivation analysis results and cultivability assessment
        
        Args:
            results: analysis results
            cultivability_assessment: cultivability assessment
            output_file: output file path
            format: output format ('csv' or 'json')
        """
        if format == 'csv':
            self.export_to_csv(results, output_file)
        elif format == 'json':
            export_data = {
                'cultivability_assessment': cultivability_assessment,
                'pathway_results': {}
            }
            for pathway_name, pathway_result in results.items():
                export_data['pathway_results'][pathway_name] = {
                    'total_genes': pathway_result['total_genes'],
                    'found_genes': pathway_result['found_genes'],
                    'completeness': pathway_result['completeness'],
                    'missing_genes': pathway_result['missing_genes']
                }
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(export_data, f, indent=2, ensure_ascii=False)
    
    def analyze_amino_acid_biosynthesis(self, results: Dict[str, Dict[str, Any]],
                                         completeness_threshold: float = 0.5) -> Dict[str, Any]:
        """
        Analyze which of the 20 standard amino acid biosynthesis pathways are present or missing
        
        Args:
            results: Cultivation analysis results from analyze_genome()
            completeness_threshold: Minimum completeness to consider a pathway as "present" (default 0.5)
            
        Returns:
            Dictionary with amino acid biosynthesis status for all 20 standard amino acids
        """
        aa_status = {}
        present_list = []
        missing_list = []
        no_pathway_list = []
        
        for aa_name, aa_info in self.AMINO_ACID_PATHWAY_MAP.items():
            pathways = aa_info['pathways']
            note = aa_info['note']
            
            if not pathways:
                # No dedicated pathway in database
                aa_status[aa_name] = {
                    'status': 'no_dedicated_pathway',
                    'pathways_checked': [],
                    'best_completeness': None,
                    'note': note
                }
                no_pathway_list.append(aa_name)
                continue
            
            # Check each pathway for this amino acid
            best_completeness = 0.0
            best_pathway = None
            pathway_details = []
            
            for pathway_name in pathways:
                if pathway_name in results:
                    completeness = results[pathway_name]['completeness']
                    found_genes = results[pathway_name]['found_genes']
                    total_genes = results[pathway_name]['total_genes']
                    pathway_details.append({
                        'name': pathway_name,
                        'completeness': completeness,
                        'found_genes': found_genes,
                        'total_genes': total_genes
                    })
                    if completeness > best_completeness:
                        best_completeness = completeness
                        best_pathway = pathway_name
                else:
                    pathway_details.append({
                        'name': pathway_name,
                        'completeness': 0.0,
                        'found_genes': 0,
                        'total_genes': 0,
                        'note': 'Pathway not detected'
                    })
            
            if best_completeness >= completeness_threshold:
                status = 'present'
                present_list.append(aa_name)
            elif best_completeness > 0:
                status = 'partial'
                missing_list.append(aa_name)
            else:
                status = 'missing'
                missing_list.append(aa_name)
            
            aa_status[aa_name] = {
                'status': status,
                'pathways_checked': pathway_details,
                'best_completeness': round(best_completeness, 4),
                'best_pathway': best_pathway,
                'note': note
            }
        
        return {
            'amino_acid_status': aa_status,
            'summary': {
                'total_amino_acids': 20,
                'present': len(present_list),
                'missing_or_partial': len(missing_list),
                'no_dedicated_pathway': len(no_pathway_list),
                'present_list': present_list,
                'missing_list': missing_list,
                'no_pathway_list': no_pathway_list
            },
            'completeness_threshold': completeness_threshold
        }
    
    def print_amino_acid_report(self, aa_analysis: Dict[str, Any]):
        """
        Print a formatted amino acid biosynthesis report
        
        Args:
            aa_analysis: Results from analyze_amino_acid_biosynthesis()
        """
        aa_status = aa_analysis['amino_acid_status']
        summary = aa_analysis['summary']
        threshold = aa_analysis['completeness_threshold']
        
        print("\n" + "=" * 70)
        print("Amino Acid Biosynthesis Pathway Analysis (20 Standard Amino Acids)")
        print("=" * 70)
        print(f"Completeness threshold: {threshold:.0%}")
        print(f"Present: {summary['present']}/20  |  Missing/Partial: {summary['missing_or_partial']}/20  |  No dedicated pathway: {summary['no_dedicated_pathway']}/20")
        print("-" * 70)
        
        # Print detailed table
        print(f"{'Amino Acid':<28} {'Status':<12} {'Best Completeness':<20} {'Best Pathway'}")
        print("-" * 70)
        
        for aa_name, info in sorted(aa_status.items()):
            status = info['status']
            if status == 'present':
                status_str = 'PRESENT'
            elif status == 'partial':
                status_str = 'PARTIAL'
            elif status == 'missing':
                status_str = 'MISSING'
            else:
                status_str = 'N/A'
            
            completeness = info['best_completeness']
            comp_str = f"{completeness:.1%}" if completeness is not None else 'N/A'
            
            best_pathway = info.get('best_pathway', '') or ''
            if status == 'no_dedicated_pathway':
                best_pathway = info.get('note', '')
            
            print(f"{aa_name:<28} {status_str:<12} {comp_str:<20} {best_pathway}")
        
        print("-" * 70)
        
        # Print missing amino acids summary
        if summary['missing_list']:
            print("\nMissing or incomplete amino acid biosynthesis pathways:")
            for aa in summary['missing_list']:
                info = aa_status[aa]
                pathways_detail = ', '.join(
                    f"{p['name']} ({p['completeness']:.1%})"
                    for p in info.get('pathways_checked', [])
                )
                print(f"  - {aa}: {info['status'].upper()} [{pathways_detail}]")
        
        if summary['no_pathway_list']:
            print("\nAmino acids without dedicated pathway in database:")
            for aa in summary['no_pathway_list']:
                note = aa_status[aa].get('note', '')
                print(f"  - {aa}: {note}")
        
        print("=" * 70 + "\n")
    
    def export_amino_acid_csv(self, aa_analysis: Dict[str, Any], output_file: str):
        """
        Export amino acid biosynthesis analysis to CSV
        
        Args:
            aa_analysis: Results from analyze_amino_acid_biosynthesis()
            output_file: Output CSV file path
        """
        rows = []
        aa_status = aa_analysis['amino_acid_status']
        
        for aa_name, info in sorted(aa_status.items()):
            # Short amino acid name (e.g., "Alanine" from "Alanine (Ala)")
            short_name = aa_name.split(' (')[0]
            abbr = aa_name.split('(')[1].rstrip(')') if '(' in aa_name else ''
            
            if info['status'] == 'no_dedicated_pathway':
                rows.append({
                    'Amino_Acid': short_name,
                    'Abbreviation': abbr,
                    'Status': 'no_dedicated_pathway',
                    'Best_Completeness': '',
                    'Best_Pathway': '',
                    'Pathway_Details': '',
                    'Note': info.get('note', '')
                })
            else:
                pathway_details = '; '.join(
                    f"{p['name']} ({p['completeness']:.1%}, {p.get('found_genes', 0)}/{p.get('total_genes', 0)} genes)"
                    for p in info.get('pathways_checked', [])
                )
                rows.append({
                    'Amino_Acid': short_name,
                    'Abbreviation': abbr,
                    'Status': info['status'],
                    'Best_Completeness': info['best_completeness'],
                    'Best_Pathway': info.get('best_pathway', ''),
                    'Pathway_Details': pathway_details,
                    'Note': info.get('note', '')
                })
        
        df = pd.DataFrame(rows)
        df.to_csv(output_file, index=False)
        print(f"Amino acid biosynthesis report exported to: {output_file}")
    
    def analyze_genomes_batch(self, faa_files: List[str], output_dir: str = None, parallel: bool = True, threads: int = 4) -> Dict[str, Dict[str, Any]]:
        """
        Batch analyze multiple genomes
        
        Args:
            faa_files: list of genome file paths
            output_dir: output directory path (optional)
            parallel: whether to use parallel processing
            threads: number of threads (only used in parallel mode)
            
        Returns:
            Analysis result dictionary for all genomes
        """
        if not faa_files:
            print("Error: No genome files provided")
            return {}
        
        # Create output directory
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        start_time = time.time()
        all_results = {}
        
        print(f"Starting batch analysis of {len(faa_files)} genomes")
        
        if parallel:
            # Parallel processing
            all_results = self._analyze_genomes_parallel(faa_files, output_dir, threads)
        else:
            # Serial processing
            for i, genome_file in enumerate(faa_files):
                genome_name = os.path.basename(genome_file).replace('.faa', '')
                print(f"\nAnalyzing {i+1}/{len(faa_files)}: {genome_name}")
                
                # Analyze genome
                results = self.analyze_genome(genome_file)
                
                if results:
                    all_results[genome_name] = results
                    
                    # Save single genome results
                    if output_dir:
                        self._save_genome_results(genome_name, results, output_dir)
        
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"\nBatch analysis completed, time elapsed: {elapsed_time:.2f} seconds")
        
        # Generate batch analysis summary report
        if output_dir and all_results:
            self._generate_batch_summary(all_results, output_dir)
        
        return all_results
    
    def _analyze_genomes_parallel(self, faa_files: List[str], output_dir: str, threads: int) -> Dict[str, Dict[str, Any]]:
        """
        Use multithreading to parallel analyze multiple genomes
        
        Args:
            faa_files: list of genome file paths
            output_dir: output directory path
            threads: number of threads
            
        Returns:
            Analysis result dictionary for all genomes
        """
        all_results = {}
        
        # Create thread lock to protect shared resources
        lock = threading.Lock()
        
        # Define single genome analysis task
        def analyze_single_genome(genome_file):
            genome_name = os.path.basename(genome_file).replace('.faa', '')
            
            # Analyze genome
            results = self.analyze_genome(genome_file)
            
            if results:
                # Use thread lock to protect shared resources
                with lock:
                    all_results[genome_name] = results
                    
                    # Save single genome results
                    if output_dir:
                        self._save_genome_results(genome_name, results, output_dir)
            
            return genome_name, results is not None
        
        # Use thread pool to execute parallel analysis
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            # Submit all tasks
            future_to_file = {executor.submit(analyze_single_genome, genome_file): genome_file for genome_file in faa_files}
            
            # Process completed tasks
            completed = 0
            for future in concurrent.futures.as_completed(future_to_file):
                genome_file = future_to_file[future]
                try:
                    genome_name, success = future.result()
                    completed += 1
                    status = "Success" if success else "Failed"
                    print(f"Progress: {completed}/{len(faa_files)} - {genome_name} analysis {status}")
                except Exception as exc:
                    print(f"Exception occurred while analyzing {genome_file}: {exc}")
        
        return all_results
    
    def _save_genome_results(self, genome_name: str, results: Dict[str, Any], output_dir: str):
        """
        Save single genome analysis results
        
        Args:
            genome_name: genome name
            results: analysis results
            output_dir: output directory path
        """
        # Export CSV
        csv_file = os.path.join(output_dir, f"{genome_name}_pathways.csv")
        self.export_to_csv(results, csv_file)
        
        # Generate recommendations report
        recommendations = self.generate_summary_recommendations(results)
        recommendations_file = os.path.join(output_dir, f"{genome_name}_recommendations.txt")
        with open(recommendations_file, 'w', encoding='utf-8') as f:
            f.write(recommendations)
        
        # Save complete results as JSON
        json_file = os.path.join(output_dir, f"{genome_name}_full_results.json")
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
    
    def _generate_batch_summary(self, all_results: Dict[str, Dict[str, Any]], output_dir: str):
        """
        Generate batch analysis summary report
        
        Args:
            all_results: analysis results for all genomes
            output_dir: output directory path
        """
        # Create summary table
        summary_rows = []
        
        for genome_name, results in all_results.items():
            summary = self.get_summary(results)
            row = {
                'Genome': genome_name,
                'Total_Pathways': summary['total_pathways'],
                'Complete_Pathways': summary['complete_pathways'],
                'Avg_Completeness': f"{summary['avg_completeness']:.2%}",
                'AA_Total': summary['amino_acid']['total'],
                'AA_Complete': summary['amino_acid']['complete'],
                'AA_Avg_Completeness': f"{summary['amino_acid']['avg_completeness']:.2%}",
                'VC_Total': summary['vitamin_coenzyme']['total'],
                'VC_Complete': summary['vitamin_coenzyme']['complete'],
                'VC_Avg_Completeness': f"{summary['vitamin_coenzyme']['avg_completeness']:.2%}",
                'NT_Total': summary['nucleotide']['total'],
                'NT_Complete': summary['nucleotide']['complete'],
                'NT_Avg_Completeness': f"{summary['nucleotide']['avg_completeness']:.2%}",
                'ATP_Total': summary['atp']['total'],
                'ATP_Complete': summary['atp']['complete'],
                'ATP_Avg_Completeness': f"{summary['atp']['avg_completeness']:.2%}"
            }
            summary_rows.append(row)
        
        # Save summary table as CSV
        summary_df = pd.DataFrame(summary_rows)
        summary_csv = os.path.join(output_dir, "batch_summary.csv")
        summary_df.to_csv(summary_csv, index=False)
        
        # Generate text summary report
        report = []
        report.append("Metabolic pathway batch analysis summary report")
        report.append("=" * 50)
        report.append(f"Number of genomes analyzed: {len(all_results)}")
        report.append("")
        
        # Sort by amino acid completeness
        summary_rows_sorted = sorted(summary_rows, 
                                    key=lambda x: (int(x['AA_Complete']), float(x['AA_Avg_Completeness'].rstrip('%'))/100), 
                                    reverse=True)
        
        report.append("Genome amino acid synthesis capability ranking:")
        report.append("-" * 50)
        for i, row in enumerate(summary_rows_sorted):
            report.append(f"{i+1}. {row['Genome']}: {row['AA_Complete']}/{row['AA_Total']} complete amino acid pathways (avg completeness: {row['AA_Avg_Completeness']})")
        
        report.append("")
        report.append("Vitamin and coenzyme metabolism capability ranking:")
        report.append("-" * 50)
        summary_rows_vc_sorted = sorted(summary_rows, 
                                      key=lambda x: (int(x['VC_Complete']), float(x['VC_Avg_Completeness'].rstrip('%'))/100), 
                                      reverse=True)
        for i, row in enumerate(summary_rows_vc_sorted):
            report.append(f"{i+1}. {row['Genome']}: {row['VC_Complete']}/{row['VC_Total']} complete vitamin/coenzyme pathways (avg completeness: {row['VC_Avg_Completeness']})")
        
        report.append("")
        report.append("Nucleotide synthesis capability ranking:")
        report.append("-" * 50)
        summary_rows_nt_sorted = sorted(summary_rows, 
                                      key=lambda x: (int(x['NT_Complete']), float(x['NT_Avg_Completeness'].rstrip('%'))/100), 
                                      reverse=True)
        for i, row in enumerate(summary_rows_nt_sorted):
            report.append(f"{i+1}. {row['Genome']}: {row['NT_Complete']}/{row['NT_Total']} complete nucleotide synthesis pathways (avg completeness: {row['NT_Avg_Completeness']})")
        
        report.append("")
        report.append("ATPase metabolic pathway capability ranking:")
        report.append("-" * 50)
        summary_rows_atp_sorted = sorted(summary_rows, 
                                       key=lambda x: (int(x['ATP_Complete']), float(x['ATP_Avg_Completeness'].rstrip('%'))/100), 
                                       reverse=True)
        for i, row in enumerate(summary_rows_atp_sorted):
            report.append(f"{i+1}. {row['Genome']}: {row['ATP_Complete']}/{row['ATP_Total']} complete ATPase metabolic pathways (avg completeness: {row['ATP_Avg_Completeness']})")
        
        report.append("")
        report.append("Detailed results:")
        report.append("-" * 50)
        for row in summary_rows_sorted:
            report.append(f"{row['Genome']}:")
            report.append(f"  Total pathways: {row['Total_Pathways']} (complete: {row['Complete_Pathways']})")
            report.append(f"  Amino acid metabolism: {row['AA_Total']} pathways (complete: {row['AA_Complete']}, avg completeness: {row['AA_Avg_Completeness']})")
            report.append(f"  Vitamin/coenzyme: {row['VC_Total']} pathways (complete: {row['VC_Complete']}, avg completeness: {row['VC_Avg_Completeness']})")
            report.append(f"  Nucleotide synthesis: {row['NT_Total']} pathways (complete: {row['NT_Complete']}, avg completeness: {row['NT_Avg_Completeness']})")
            report.append(f"  ATPase metabolism: {row['ATP_Total']} pathways (complete: {row['ATP_Complete']}, avg completeness: {row['ATP_Avg_Completeness']})")
            report.append("")
        
        # Save text report
        summary_report = os.path.join(output_dir, "batch_summary_report.txt")
        with open(summary_report, 'w', encoding='utf-8') as f:
            f.write('\n'.join(report))
        
        print(f"Summary table saved to: {summary_csv}")
        print(f"Summary report saved to: {summary_report}")
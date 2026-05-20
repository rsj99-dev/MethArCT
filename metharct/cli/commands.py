#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MethArCT CLI Commands

Implements the individual command functions for the MethArCT CLI.
"""

import sys
import time
from pathlib import Path
from typing import Optional

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from metharct.core import (
    DiamondAnalyzer,
    TomeAnalyzer,
    CheckM2Analyzer,
    PathwayPredictor,
    CultivationAnalyzer
)
from metharct.utils.config import Config
from metharct.utils.logger import get_logger
from metharct.utils.file_utils import FileUtils

def comprehensive_command(input_path: str,
                         output_prefix: str,
                         config: Config,
                         skip_tome: bool = False,
                         skip_checkm2: bool = False):
    """
    Run comprehensive analysis command
    
    Args:
        input_path: Path to input FASTA file
        output_prefix: Output file prefix
        config: Configuration object
        skip_tome: Skip Tome analysis
        skip_checkm2: Skip CheckM2 analysis
    """
    logger = get_logger("comprehensive_command")
    
    print("\n" + "=" * 60)
    print("MethArCT Comprehensive Analysis")
    print("=" * 60)
    print(f"Input file: {input_path}")
    print(f"Output prefix: {output_prefix}")
    print(f"Skip Tome: {skip_tome}")
    print(f"Skip CheckM2: {skip_checkm2}")
    print("=" * 60 + "\n")
    
    try:
        # Validate input file
        if not FileUtils.validate_fasta(input_path):
            raise ValueError(f"Invalid FASTA file: {input_path}")
        
        # Get sequence count for progress tracking
        seq_count = FileUtils.count_sequences(input_path)
        print(f"Processing {seq_count} sequences...\n")
        
        # Initialize pathway predictor
        predictor = PathwayPredictor(config)
        
        # Start analysis
        start_time = time.time()
        
        results = predictor.predict_comprehensive(
            input_path=input_path,
            output_prefix=output_prefix,
            include_tome=not skip_tome,
            include_checkm2=not skip_checkm2
        )
        
        end_time = time.time()
        analysis_time = end_time - start_time
        
        # Print summary
        print("\n" + "=" * 60)
        print("Analysis Summary")
        print("=" * 60)
        
        analyses_performed = results.get('analyses_performed', [])
        print(f"Analyses completed: {', '.join(analyses_performed)}")
        print(f"Total analysis time: {analysis_time:.2f} seconds")
        
        # Print integrated results summary if available
        if 'integrated_analysis' in results:
            integrated = results['integrated_analysis']
            
            if 'overall_assessment' in integrated:
                overall = integrated['overall_assessment']
                print(f"\nOrganism type: {overall.get('organism_type', 'Unknown')}")
                print(f"Metabolic complexity: {overall.get('metabolic_complexity', 'Unknown')}")
                print(f"Cultivation potential: {overall.get('cultivation_potential', 'Unknown')}")
                print(f"Overall confidence: {overall.get('confidence', 0):.2f}")
            
            if 'recommendations' in integrated and integrated['recommendations']:
                print("\nKey Recommendations:")
                for i, rec in enumerate(integrated['recommendations'][:5], 1):
                    print(f"  {i}. {rec}")
                
                if len(integrated['recommendations']) > 5:
                    print(f"  ... and {len(integrated['recommendations']) - 5} more (see detailed report)")
        
        print("\n" + "=" * 60)
        
        logger.info(f"Comprehensive analysis completed in {analysis_time:.2f} seconds")
        
    except Exception as e:
        logger.error(f"Comprehensive analysis failed: {str(e)}")
        raise

def diamond_command(input_path: str,
                   output_prefix: str,
                   config: Config):
    """
    Run Diamond analysis command
    
    Args:
        input_path: Path to input FASTA file
        output_prefix: Output file prefix
        config: Configuration object
    """
    logger = get_logger("diamond_command")
    
    print("\n" + "=" * 60)
    print("MethArCT Diamond Analysis")
    print("=" * 60)
    print(f"Input file: {input_path}")
    print(f"Output prefix: {output_prefix}")
    print("=" * 60 + "\n")
    
    try:
        # Validate input file
        if not FileUtils.validate_fasta(input_path):
            raise ValueError(f"Invalid FASTA file: {input_path}")
        
        # Get sequence count
        seq_count = FileUtils.count_sequences(input_path)
        print(f"Processing {seq_count} sequences...\n")
        
        # Initialize analyzer
        analyzer = DiamondAnalyzer(config)
        
        # Start analysis
        start_time = time.time()
        print("Starting Diamond analysis...")
        
        results = analyzer.analyze_sequence(
            input_file=input_path,
            output_prefix=output_prefix
        )
        
        end_time = time.time()
        analysis_time = end_time - start_time
        
        # Print summary
        print("\n" + "=" * 60)
        print("Diamond Analysis Summary")
        print("=" * 60)
        
        if 'summary' in results:
            summary = results['summary']
            print(f"Total pathways detected: {summary.get('total_pathways', 0)}")
            print(f"Methane pathways: {summary.get('methane_pathways', 0)}")
            print(f"Sulfur pathways: {summary.get('sulfur_pathways', 0)}")
            print(f"Nitrogen pathways: {summary.get('nitrogen_pathways', 0)}")
            print(f"Salt tolerance hits: {summary.get('salt_tolerance_hits', 0)}")
            print(f"Cultivability hits: {summary.get('cultivability_hits', 0)}")
        
        print(f"\nAnalysis time: {analysis_time:.2f} seconds")
        print("=" * 60)
        
        logger.info(f"Diamond analysis completed in {analysis_time:.2f} seconds")
        
    except Exception as e:
        logger.error(f"Diamond analysis failed: {str(e)}")
        raise

def cultivation_command(input_path: str,
                       output_prefix: str,
                       config: Config):
    """
    Run cultivation analysis command
    
    Args:
        input_path: Path to input FASTA file
        output_prefix: Output file prefix
        config: Configuration object
    """
    logger = get_logger("cultivation_command")
    
    print("\n" + "=" * 60)
    print("MethArCT Cultivation Analysis")
    print("=" * 60)
    print(f"Input file: {input_path}")
    print(f"Output prefix: {output_prefix}")
    print("=" * 60 + "\n")
    
    try:
        # Validate input file
        if not FileUtils.validate_fasta(input_path):
            raise ValueError(f"Invalid FASTA file: {input_path}")
        
        # Get sequence count
        seq_count = FileUtils.count_sequences(input_path)
        print(f"Processing {seq_count} sequences...\n")
        
        # Initialize analyzer
        analyzer = CultivationAnalyzer(config.config)
        
        # Start analysis
        start_time = time.time()
        print("Starting cultivation analysis...")
        print("Note: This analysis evaluates cultivability based on metabolic pathways...\n")
        
        # Analyze genome
        results = analyzer.analyze_genome(
            genome_file=input_path
        )
        
        # Generate report
        report = analyzer.generate_report(results=results)
        
        # Save report to file
        report_file = f"{output_prefix}_cultivation_report.txt"
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(report)
        
        end_time = time.time()
        analysis_time = end_time - start_time
        
        # Print summary
        print("\n" + "=" * 60)
        print("Cultivation Analysis Summary")
        print("=" * 60)
        
        if 'summary' in results:
            summary = results['summary']
            print(f"Total pathways analyzed: {summary.get('total_pathways', 0)}")
            print(f"Amino acid pathways: {summary.get('amino_acid_pathways', 0)}")
            print(f"Vitamin pathways: {summary.get('vitamin_pathways', 0)}")
            print(f"Nucleotide pathways: {summary.get('nucleotide_pathways', 0)}")
            print(f"ATP pathways: {summary.get('atp_pathways', 0)}")
            print(f"Cultivability score: {summary.get('cultivability_score', 0):.2f}")
            print(f"Cultivability level: {summary.get('cultivability_level', 'Unknown')}")
        
        if 'recommendations' in results and results['recommendations']:
            print("\nCultivation Recommendations:")
            for i, rec in enumerate(results['recommendations'][:5], 1):
                print(f"  {i}. {rec}")
            
            if len(results['recommendations']) > 5:
                print(f"  ... and {len(results['recommendations']) - 5} more (see detailed report)")
        
        print(f"\nAnalysis time: {analysis_time:.2f} seconds")
        print(f"Report saved to: {output_prefix}_cultivation_report.txt")
        print("=" * 60)
        
        logger.info(f"Cultivation analysis completed in {analysis_time:.2f} seconds")
        
    except Exception as e:
        logger.error(f"Cultivation analysis failed: {str(e)}")
        raise

def tome_command(input_path: str,
                output_prefix: str,
                config: Config,
                batch_size: int = 100):
    """
    Run Tome analysis command
    
    Args:
        input_path: Path to input FASTA file
        output_prefix: Output file prefix
        config: Configuration object
        batch_size: Batch size for processing
    """
    logger = get_logger("tome_command")
    
    print("\n" + "=" * 60)
    print("MethArCT Tome Analysis")
    print("=" * 60)
    print(f"Input file: {input_path}")
    print(f"Output prefix: {output_prefix}")
    print(f"Batch size: {batch_size}")
    print("=" * 60 + "\n")
    
    try:
        # Validate input file
        if not FileUtils.validate_fasta(input_path):
            raise ValueError(f"Invalid FASTA file: {input_path}")
        
        # Get sequence count
        seq_count = FileUtils.count_sequences(input_path)
        print(f"Processing {seq_count} protein sequences...\n")
        
        # Initialize analyzer
        analyzer = TomeAnalyzer(config)
        
        # Update batch size in config
        config.set('tools.tome.batch_size', batch_size)
        
        # Start analysis
        start_time = time.time()
        print("Starting Tome analysis...")
        print("Note: This may take a while for large files...\n")
        
        results = analyzer.predict_ogt(
            protein_file=input_path,
            output_dir=output_prefix
        )
        
        end_time = time.time()
        analysis_time = end_time - start_time
        
        # Print summary
        print("\n" + "=" * 60)
        print("Tome Analysis Summary")
        print("=" * 60)
        
        if 'summary' in results:
            summary = results['summary']
            print(f"Sequences processed: {summary.get('total_sequences', 0)}")
            print(f"Successful predictions: {summary.get('successful_predictions', 0)}")
            print(f"Average OGT: {summary.get('average_ogt', 0):.1f}°C")
            print(f"Temperature category: {summary.get('temperature_category', 'Unknown')}")
            print(f"Confidence: {summary.get('confidence', 0):.2f}")
            
            # Show temperature distribution if available
            if 'temperature_distribution' in summary:
                dist = summary['temperature_distribution']
                print("\nTemperature distribution:")
                for category, count in dist.items():
                    print(f"  {category}: {count} sequences")
        
        print(f"\nAnalysis time: {analysis_time:.2f} seconds")
        print("=" * 60)
        
        logger.info(f"Tome analysis completed in {analysis_time:.2f} seconds")
        
    except Exception as e:
        logger.error(f"Tome analysis failed: {str(e)}")
        raise

def checkm2_command(input_path: str,
                   output_prefix: str,
                   config: Config,
                   input_type: str = 'fasta'):
    """
    Run CheckM2 analysis command
    
    Args:
        input_path: Path to input file or directory
        output_prefix: Output file prefix
        config: Configuration object
        input_type: Type of input ('fasta' or 'directory')
    """
    logger = get_logger("checkm2_command")
    
    print("\n" + "=" * 60)
    print("MethArCT CheckM2 Analysis")
    print("=" * 60)
    print(f"Input path: {input_path}")
    print(f"Input type: {input_type}")
    print(f"Output prefix: {output_prefix}")
    print("=" * 60 + "\n")
    
    try:
        # Validate input
        input_path_obj = Path(input_path)
        if not input_path_obj.exists():
            raise FileNotFoundError(f"Input path not found: {input_path}")
        
        if input_type == 'fasta' and not FileUtils.validate_fasta(input_path):
            raise ValueError(f"Invalid FASTA file: {input_path}")
        
        # Initialize analyzer
        analyzer = CheckM2Analyzer(config)
        
        # Start analysis
        start_time = time.time()
        print("Starting CheckM2 analysis...")
        print("Note: This analysis may take significant time...\n")
        
        results = analyzer.analyze_genome_quality(
            input_path=input_path,
            output_prefix=output_prefix,
            input_type=input_type
        )
        
        end_time = time.time()
        analysis_time = end_time - start_time
        
        # Print summary
        print("\n" + "=" * 60)
        print("CheckM2 Analysis Summary")
        print("=" * 60)
        
        if 'summary' in results:
            summary = results['summary']
            print(f"Total genomes analyzed: {summary.get('total_genomes', 0)}")
            print(f"High quality genomes: {summary.get('high_quality_genomes', 0)}")
            print(f"Medium quality genomes: {summary.get('medium_quality_genomes', 0)}")
            print(f"Low quality genomes: {summary.get('low_quality_genomes', 0)}")
            print(f"Average completeness: {summary.get('average_completeness', 0):.1f}%")
            print(f"Average contamination: {summary.get('average_contamination', 0):.1f}%")
            print(f"Average quality score: {summary.get('average_quality_score', 0):.1f}")
            
            # Cultivability assessment
            if 'cultivability_assessment' in summary:
                cult_assess = summary['cultivability_assessment']
                print(f"\nCultivability assessment: {cult_assess.get('overall_cultivability', 'Unknown')}")
                print(f"Confidence: {cult_assess.get('confidence', 0):.2f}")
                print(f"Recommendation: {cult_assess.get('recommendation', 'No recommendation')}")
        
        print(f"\nAnalysis time: {analysis_time:.2f} seconds")
        print("=" * 60)
        
        logger.info(f"CheckM2 analysis completed in {analysis_time:.2f} seconds")
        
    except Exception as e:
        logger.error(f"CheckM2 analysis failed: {str(e)}")
        raise

def print_progress(current: int, total: int, prefix: str = "Progress"):
    """
    Print progress bar
    
    Args:
        current: Current progress
        total: Total items
        prefix: Progress prefix text
    """
    if total == 0:
        return
    
    percent = (current / total) * 100
    bar_length = 50
    filled_length = int(bar_length * current // total)
    
    bar = '█' * filled_length + '-' * (bar_length - filled_length)
    
    print(f'\r{prefix}: |{bar}| {percent:.1f}% ({current}/{total})', end='', flush=True)
    
    if current == total:
        print()  # New line when complete

def format_time(seconds: float) -> str:
    """
    Format time duration
    
    Args:
        seconds: Time in seconds
        
    Returns:
        Formatted time string
    """
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        minutes = seconds / 60
        return f"{minutes:.1f}m"
    else:
        hours = seconds / 3600
        return f"{hours:.1f}h"

def print_file_info(file_path: str):
    """
    Print information about input file
    
    Args:
        file_path: Path to file
    """
    try:
        path = Path(file_path)
        file_size = FileUtils.get_file_size(file_path)
        
        print(f"File: {path.name}")
        print(f"Size: {FileUtils.format_file_size(file_size)}")
        
        if FileUtils.validate_fasta(file_path):
            seq_count = FileUtils.count_sequences(file_path)
            print(f"Sequences: {seq_count}")
        
    except Exception as e:
        print(f"Could not read file info: {str(e)}")
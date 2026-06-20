#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MethArCT CLI Commands

Implements the individual command functions for the MethArCT CLI.
"""

import csv
import sys
import time
from pathlib import Path
from typing import Optional, List, Dict, Any

from metharct.core import (
    DiamondAnalyzer,
    TomeAnalyzer,
    CheckM2Analyzer,
    PathwayPredictor,
    CultivationAnalyzer,
    SuShaAnalyzer,
    PHAnalyzer,
    AntibioticAnalyzer,
)
from metharct.utils.config import Config
from metharct.utils.logger import get_logger
from metharct.utils.file_utils import FileUtils

def comprehensive_command(input_path: str,
                         output_prefix: str,
                         config: Config,
                         skip_tome: bool = False,
                         skip_checkm2: bool = False,
                         skip_susha: bool = False,
                         skip_ph: bool = False,
                         skip_antibiotic: bool = False):
    """
    Run comprehensive analysis command
    
    Args:
        input_path: Path to input FASTA file
        output_prefix: Output file prefix
        config: Configuration object
        skip_tome: Skip Tome analysis
        skip_checkm2: Skip CheckM2 analysis
        skip_susha: Skip SuSha salinity prediction
        skip_ph: Skip pH preference prediction
        skip_antibiotic: Skip antibiotic resistance prediction
    """
    logger = get_logger("comprehensive_command")
    
    print("\n" + "=" * 60)
    print("MethArCT Comprehensive Analysis")
    print("=" * 60)
    print(f"Input file: {input_path}")
    print(f"Output prefix: {output_prefix}")
    print(f"Skip Tome: {skip_tome}")
    print(f"Skip CheckM2: {skip_checkm2}")
    print(f"Skip SuSha: {skip_susha}")
    print(f"Skip pH: {skip_ph}")
    print(f"Skip Antibiotic: {skip_antibiotic}")
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
            include_checkm2=not skip_checkm2,
            include_susha=not skip_susha,
            include_ph=not skip_ph,
            include_antibiotic=not skip_antibiotic
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

                # Print key characteristics if available
                key_chars = overall.get('key_characteristics', [])
                if key_chars:
                    print("\nKey Characteristics:")
                    for char in key_chars:
                        print(f"  - {char}")

            # Print pH prediction summary if available
            env_adaptation = integrated.get('environmental_adaptation', {})
            ph_data = env_adaptation.get('ph', {})
            if ph_data:
                print("\npH Preference Prediction:")
                ph_opt = ph_data.get('ph_optimum')
                ph_max = ph_data.get('ph_max')
                ph_min = ph_data.get('ph_min')
                is_novel = ph_data.get('is_novel', False)
                if ph_opt is not None:
                    print(f"  pH optimum: {ph_opt:.2f}")
                if ph_max is not None:
                    print(f"  pH maximum: {ph_max:.2f}")
                if ph_min is not None:
                    print(f"  pH minimum: {ph_min:.2f}")
                if is_novel:
                    print("  Note: Input genome differs from training set; predictions may be less reliable.")
            elif 'ph' in results.get('results', {}):
                ph_result = results['results']['ph']
                if ph_result.get('status') == 'failed':
                    print(f"\npH Preference Prediction: FAILED ({ph_result.get('error', 'Unknown error')})")

            # Print antibiotic resistance summary
            antibiotic_data = integrated.get('antibiotic_resistance', {})
            if antibiotic_data:
                print("\nAntibiotic Resistance Prediction:")
                if antibiotic_data.get('status') == 'failed':
                    print(f"  FAILED ({antibiotic_data.get('error', 'Unknown error')})")
                else:
                    recommended_abs = antibiotic_data.get('recommended_antibiotics', [])
                    if recommended_abs:
                        print(f"  Recommended: {', '.join(recommended_abs)}")
                    else:
                        print("  No matching antibiotic recommendations")
                    aai = antibiotic_data.get('aai_results', {})
                    if aai:
                        print("  AAI values:")
                        for ref, val in sorted(aai.items()):
                            print(f"    {ref}: {val:.2f}%")
            elif 'antibiotic' in results.get('results', {}):
                ab_result = results['results']['antibiotic']
                if ab_result.get('status') == 'failed':
                    print(f"\nAntibiotic Resistance Prediction: FAILED ({ab_result.get('error', 'Unknown error')})")

            # Print energy metabolism completeness summary
            energy_metabolism = integrated.get('energy_metabolism', {})
            if energy_metabolism:
                print("\nEnergy Metabolism Pathway Completeness:")
                category_labels = {
                    'methane': 'Methane metabolism',
                    'sulfur': 'Sulfur metabolism',
                    'nitrogen': 'Nitrogen metabolism'
                }
                for cat_key, label in category_labels.items():
                    cat_data = energy_metabolism.get(cat_key, {})
                    complete_list = cat_data.get('complete', [])
                    total = cat_data.get('total_count', 0)
                    complete_count = cat_data.get('complete_count', 0)
                    if complete_list:
                        names = [p['name'] for p in complete_list]
                        print(f"  {label} ({complete_count}/{total} complete): {', '.join(names)}")
                    else:
                        print(f"  {label} (0/{total} complete): no complete pathway detected")

            if 'recommendations' in integrated and integrated['recommendations']:
                print("\nKey Recommendations:")
                for i, rec in enumerate(integrated['recommendations'][:10], 1):
                    print(f"  {i}. {rec}")

                if len(integrated['recommendations']) > 10:
                    print(f"  ... and {len(integrated['recommendations']) - 10} more (see detailed report)")
            
            # Print amino acid biosynthesis summary
            aa_data = integrated.get('amino_acid_biosynthesis', {})
            if aa_data:
                aa_summary = aa_data.get('summary', {})
                print(f"\nAmino Acid Biosynthesis (20 standard amino acids):")
                print(f"  Present: {aa_summary.get('present', 0)}/20")
                print(f"  Missing/Partial: {aa_summary.get('missing_or_partial', 0)}/20")
                print(f"  No dedicated pathway in DB: {aa_summary.get('no_dedicated_pathway', 0)}/20")
                
                missing_list = aa_summary.get('missing_list', [])
                if missing_list:
                    missing_names = [aa.split(' (')[0] for aa in missing_list]
                    print(f"  Missing: {', '.join(missing_names)}")
        
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
            pathway_results = results.get('pathway_results', {})
            
            # Count pathways by type
            methane_count = sum(1 for db in pathway_results.keys() if db in ['CO2-CH4', 'JIAAN-CH4', 'JIACHUN-CH4', 'JIALIUCHUN-CH4', 'YISUAN-CH4', 'C16-CH4', 'CO-CH4', 'JIASUAN-CH4', 'JIAYANGJI-CH4', 'ZHIFANGSUAN-CH4', '2JIAAN-CH4', '3JIAAN-CH4', 'Glycine betaine methanogenesis', 'Methylthiopropionate methanogenesis', 'Tetramethylammonium methanogenesis', 'Methanol dismutation methanogenesis'])
            sulfur_count = sum(1 for db in pathway_results.keys() if db in ['ASR', 'SO', 'SOX', 'S4I', 'SR', 'DSR'])
            nitrogen_count = sum(1 for db in pathway_results.keys() if db in ['ANR', 'DEN', 'DNR', 'NIT'])
            
            cult_hits = pathway_results.get('CULTIVATION', {}).get('cultivability_hits', 0)
            
            print(f"Total pathways detected: {summary.get('pathways_detected', 0)}")
            print(f"Methane pathways: {methane_count}")
            print(f"Sulfur pathways: {sulfur_count}")
            print(f"Nitrogen pathways: {nitrogen_count}")
            print(f"Cultivability hits: {cult_hits}")
        
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


def susha_command(input_path: str,
                  output_prefix: str,
                  config: Config):
    """
    Run SuSha salinity adaptation prediction
    
    Args:
        input_path: Path to input FASTA file
        output_prefix: Output file prefix
        config: Configuration object
    """
    logger = get_logger("susha_command")
    
    print("\n" + "=" * 60)
    print("MethArCT SuSha Salinity Prediction")
    print("=" * 60)
    print(f"Input file: {input_path}")
    print(f"Output prefix: {output_prefix}")
    print("=" * 60 + "\n")
    
    try:
        # Validate input file
        if not FileUtils.validate_fasta(input_path):
            raise ValueError(f"Invalid FASTA file: {input_path}")
        
        # Initialize analyzer
        analyzer = SuShaAnalyzer(config)
        
        if not analyzer.tool_available:
            raise RuntimeError("SuSha module is not available. Please ensure 'shap' and 'openpyxl' are installed.")
        
        # Start analysis
        start_time = time.time()
        print("Starting SuSha salinity prediction...")
        
        results = analyzer.predict_salinity(
            input_file=input_path,
            output_prefix=output_prefix
        )
        
        end_time = time.time()
        analysis_time = end_time - start_time
        
        # Print summary
        print("\n" + "=" * 60)
        print("SuSha Prediction Summary")
        print("=" * 60)
        
        if results.get('status') == 'success':
            prediction = results.get('prediction', {})
            print(f"Predicted salinity: {prediction.get('salinity_label', 'Unknown')}")
            print(f"Confidence: {prediction.get('confidence', 0):.2%}")
            
            top3 = results.get('top3_predictions', [])
            if top3:
                print("\nTop 3 predictions:")
                for t in top3:
                    print(f"  {t['rank']}. {t['label']}: {t['probability']:.2%}")
            
            output_files = results.get('output_files', {})
            if output_files.get('tsv'):
                print(f"\nResults saved to: {output_files['tsv']}")
        else:
            print(f"Prediction failed: {results.get('error', 'Unknown error')}")
        
        print(f"\nAnalysis time: {analysis_time:.2f} seconds")
        print("=" * 60)
        
        logger.info(f"SuSha prediction completed in {analysis_time:.2f} seconds")
        
    except Exception as e:
        logger.error(f"SuSha prediction failed: {str(e)}")
        raise


def ph_command(input_path: str,
               output_prefix: str,
               config: Config):
    """
    Run pH preference prediction

    Args:
        input_path: Path to input FASTA file
        output_prefix: Output file prefix
        config: Configuration object
    """
    logger = get_logger("ph_command")

    print("\n" + "=" * 60)
    print("MethArCT pH Preference Prediction")
    print("=" * 60)
    print(f"Input file: {input_path}")
    print(f"Output prefix: {output_prefix}")
    print("=" * 60 + "\n")

    try:
        # Validate input file
        if not FileUtils.validate_fasta(input_path):
            raise ValueError(f"Invalid FASTA file: {input_path}")

        # Initialize analyzer
        analyzer = PHAnalyzer(config)

        if not analyzer.tool_available:
            raise RuntimeError(
                "pH predictor module is not available. "
                "Please ensure 'hmmlearn' is installed."
            )

        # Start analysis
        start_time = time.time()
        print("Starting pH preference prediction...")

        results = analyzer.predict_ph(
            input_file=input_path,
            output_prefix=output_prefix
        )

        end_time = time.time()
        analysis_time = end_time - start_time

        # Print summary
        print("\n" + "=" * 60)
        print("pH Prediction Summary")
        print("=" * 60)

        if results.get('status') == 'success':
            prediction = results.get('prediction', {})
            ph_opt = prediction.get('ph_optimum', {}).get('value')
            ph_max = prediction.get('ph_max', {}).get('value')
            ph_min = prediction.get('ph_min', {}).get('value')
            is_novel = results.get('is_novel', False)

            print(f"pH optimum: {ph_opt:.2f}" if ph_opt is not None else "pH optimum: N/A")
            print(f"pH maximum: {ph_max:.2f}" if ph_max is not None else "pH maximum: N/A")
            print(f"pH minimum: {ph_min:.2f}" if ph_min is not None else "pH minimum: N/A")
            if is_novel:
                print("\nNote: Input genome differs from training set; predictions may be less reliable.")

            output_files = results.get('output_files', {})
            if output_files.get('tsv'):
                print(f"\nResults saved to: {output_files['tsv']}")
        else:
            print(f"Prediction failed: {results.get('error', 'Unknown error')}")

        print(f"\nAnalysis time: {analysis_time:.2f} seconds")
        print("=" * 60)

        logger.info(f"pH prediction completed in {analysis_time:.2f} seconds")

    except Exception as e:
        logger.error(f"pH prediction failed: {str(e)}")
        raise


def antibiotic_command(input_path: str,
                      output_prefix: str,
                      config: Config):
    """
    Run antibiotic resistance prediction

    Args:
        input_path: Path to input FASTA file
        output_prefix: Output file prefix
        config: Configuration object
    """
    logger = get_logger("antibiotic_command")

    print("\n" + "=" * 60)
    print("MethArCT Antibiotic Resistance Prediction")
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
        threads = config.get('tools.diamond.threads', 4)
        evalue = config.get('tools.diamond.evalue', 1e-5)
        analyzer = AntibioticAnalyzer(
            config=config,
            cpus=threads,
            evalue=evalue,
        )

        # Start analysis
        start_time = time.time()
        print("Starting antibiotic resistance prediction...")
        print("Note: This runs DIAMOND blastp against 5 methanogen reference genomes\n")

        results = analyzer.predict_antibiotics(
            input_file=input_path,
            output_prefix=output_prefix,
        )

        end_time = time.time()
        analysis_time = end_time - start_time

        # Print summary
        print("\n" + "=" * 60)
        print("Antibiotic Resistance Prediction Summary")
        print("=" * 60)

        if results.get('status') == 'success':
            # Print AAI results
            aai_results = results.get('aai_results', {})
            if aai_results:
                print("\nAAI against reference genomes:")
                for ref_file in sorted(aai_results.keys()):
                    aai = aai_results[ref_file]
                    details = results.get('aai_details', {}).get(ref_file, {})
                    nhits = details.get('num_hits', 0)
                    nquery = details.get('num_queries', 0)
                    print(f"  {ref_file}: {aai:.2f}% ({nhits}/{nquery} hits)")

            # Print recommendations
            recommended = results.get('recommended_antibiotics', [])
            print("\nRecommended Antibiotics:")
            if recommended:
                for ab in recommended:
                    print(f"  -> {ab}")
            else:
                print("  No matching antibiotic recommendations")

            # Print rule details
            all_rules = results.get('all_rules', [])
            if all_rules:
                print("\nRule Evaluation:")
                for rule_result in all_rules:
                    status = 'MATCHED' if rule_result['satisfied'] else 'NOT MATCHED'
                    print(f"  [{status}] {rule_result['name']}")
                    for detail in rule_result['details']:
                        print(f"    {detail}")

            output_file = results.get('output_file')
            if output_file:
                print(f"\nResults saved to: {output_file}")
        else:
            print(f"Prediction failed: {results.get('error', 'Unknown error')}")

        print(f"\nAnalysis time: {analysis_time:.2f} seconds")
        print("=" * 60)

        logger.info(f"Antibiotic prediction completed in {analysis_time:.2f} seconds")

    except Exception as e:
        logger.error(f"Antibiotic prediction failed: {str(e)}")
        raise


# ============================================================
# Batch comprehensive analysis command
# ============================================================

# Salinity label -> Chinese range mapping
_SALINITY_RANGE_MAP = {
    'Salt-sensitive': '0%-1%',
    'Halotolerant': '1%-3%',
    'Slight halophilic': '3%-5%',
    'Moderate halophilic': '5%-15%',
    'Extreme halophilic': '15%-30%',
}

# Methane pathway db key -> Chinese substrate name
_SUBSTRATE_CN_MAP = {
    # Methane metabolism
    'CO2-CH4': 'CO2/H2',
    'JIAAN-CH4': 'Methylamine',
    'JIACHUN-CH4': 'Methanol',
    'JIALIUCHUN-CH4': 'Methanethiol',
    'YISUAN-CH4': 'Acetate',
    'C16-CH4': 'Long-chain fatty acids',
    'CO-CH4': 'CO',
    'JIASUAN-CH4': 'Formate',
    'JIAYANGJI-CH4': 'Methoxy compounds',
    'ZHIFANGSUAN-CH4': 'Fatty acids',
    '2JIAAN-CH4': 'Dimethylamine',
    '3JIAAN-CH4': 'Trimethylamine',
    'Glycine betaine methanogenesis': 'Glycine betaine',
    'Methylthiopropionate methanogenesis': 'Methylthiopropionate',
    'Tetramethylammonium methanogenesis': 'Tetramethylammonium',
    'Methanol dismutation methanogenesis': 'Methanol dismutation',
    # Sulfur metabolism
    'ASR': 'Assimilatory sulfate reduction',
    'SO': 'Sulfide oxidation',
    'SOX': 'SOX sulfur oxidation',
    'S4I': 'S4I sulfur oxidation',
    'SR': 'Sulfur reduction',
    'DSR': 'Dissimilatory sulfate reduction',
    # Nitrogen metabolism
    'ANR': 'Assimilatory nitrate reduction',
    'DEN': 'Denitrification',
    'DNR': 'Dissimilatory nitrate reduction',
    'NIT': 'Nitrification',
}

# Batch output CSV columns
_BATCH_CSV_COLUMNS = [
    'FAA_filename', 'Temperature_range', 'Optimal_temperature',
    'pH_range', 'Optimal_pH', 'Salinity_range',
    'Substrate_metabolism', 'Additional_amino_acids_required', 'Recommended_antibiotics',
]


def _estimate_temp_range(ogt: float) -> str:
    """Estimate growth temperature range from OGT."""
    if ogt < 20:
        spread = 8
    elif ogt < 40:
        spread = 12
    elif ogt < 60:
        spread = 15
    else:
        spread = 10
    t_min = round(max(ogt - spread, 0.0), 1)
    t_max = round(ogt + spread, 1)
    return f"{t_min}-{t_max}"


# Chinese amino acid name -> English translation
_AA_CN_TO_EN = {
    '半胱氨酸': 'Cysteine',
    '苯丙氨酸': 'Phenylalanine',
    '蛋氨酸': 'Methionine',
    '脯氨酸': 'Proline',
    '精氨酸': 'Arginine',
    '赖氨酸': 'Lysine',
    '酪氨酸': 'Tyrosine',
    '亮氨酸': 'Leucine',
    '异亮氨酸': 'Isoleucine',
    '色氨酸': 'Tryptophan',
    '苏氨酸': 'Threonine',
    '缬氨酸': 'Valine',
    '丝氨酸': 'Serine',
    '组氨酸': 'Histidine',
}


def _extract_missing_amino_acids(cult_results: Dict) -> str:
    """
    Extract amino acids whose biosynthesis pathways are all incomplete.

    Returns a comma-separated string of English amino acid names.
    """
    # Build amino acid -> [completeness values] mapping
    aa_completeness: Dict[str, List[float]] = {}
    for pathway_name, pathway_data in cult_results.items():
        if '生物合成' not in pathway_name:
            continue
        aa_name = pathway_name.split('生物合成')[0]
        if aa_name not in aa_completeness:
            aa_completeness[aa_name] = []
        aa_completeness[aa_name].append(pathway_data.get('completeness', 0.0))

    missing = []
    for aa_name, completeness_list in aa_completeness.items():
        # All pathways for this amino acid are incomplete (< 100%)
        if all(c < 1.0 for c in completeness_list):
            en_name = _AA_CN_TO_EN.get(aa_name, aa_name)
            missing.append(en_name)

    return ', '.join(missing) if missing else ''


def _run_single_batch_analysis(
    faa_path: Path,
    output_prefix: str,
    config: Config,
    analyzers: Dict[str, Any],
    logger,
    mags: bool = False,
) -> Dict[str, str]:
    """
    Run all analyses on a single .faa file and return a row dict
    matching the batch CSV columns.

    Args:
        mags: If True, show pathways with completeness >= 70%;
              otherwise only show pathways with completeness == 100%.
    """
    row = {col: '' for col in _BATCH_CSV_COLUMNS}
    row['FAA_filename'] = faa_path.name

    input_str = str(faa_path)

    # ---- Tome (temperature) ----
    try:
        tome = analyzers['tome']
        tome_results = tome.predict_ogt(input_file=input_str, output_prefix=output_prefix)
        ogt = tome_results.get('summary', {}).get('predicted_ogt_celsius')
        if ogt is not None:
            row['Optimal_temperature'] = round(ogt, 1)
            row['Temperature_range'] = _estimate_temp_range(ogt)
    except Exception as e:
        logger.warning(f"[{faa_path.name}] Tome failed: {e}")

    # ---- pH ----
    try:
        ph = analyzers['ph']
        ph_results = ph.predict_ph(input_file=input_str, output_prefix=output_prefix)
        if ph_results.get('status') == 'success':
            ph_min = ph_results.get('summary', {}).get('ph_min')
            ph_max = ph_results.get('summary', {}).get('ph_max')
            ph_opt = ph_results.get('summary', {}).get('ph_optimum')
            if ph_opt is not None:
                row['Optimal_pH'] = round(ph_opt, 1)
            if ph_min is not None and ph_max is not None:
                row['pH_range'] = f"{round(ph_min, 1)}-{round(ph_max, 1)}"
    except Exception as e:
        logger.warning(f"[{faa_path.name}] pH failed: {e}")

    # ---- SuSha (salinity) ----
    try:
        susha = analyzers['susha']
        susha_results = susha.predict_salinity(input_file=input_str, output_prefix=output_prefix)
        if susha_results.get('status') == 'success':
            label = susha_results.get('prediction', {}).get('salinity_label', '')
            row['Salinity_range'] = _SALINITY_RANGE_MAP.get(label, label)
    except Exception as e:
        logger.warning(f"[{faa_path.name}] SuSha failed: {e}")

    # ---- Diamond (energy metabolism pathways) ----
    try:
        diamond = analyzers['diamond']
        diamond_results = diamond.analyze_sequence(input_file=input_str, output_prefix=output_prefix)
        pathway_results = diamond_results.get('pathway_results', {})
        # Completeness threshold: 100% by default, 70% for MAGs
        completeness_threshold = 70.0 if mags else 100.0
        substrates = []
        for db_key, cn_name in _SUBSTRATE_CN_MAP.items():
            if db_key in pathway_results:
                pw_data = pathway_results[db_key]
                completeness = pw_data.get('low_completeness_percentage', 0.0)
                if completeness >= completeness_threshold:
                    substrates.append(cn_name)
        row['Substrate_metabolism'] = ', '.join(substrates) if substrates else ''
    except Exception as e:
        logger.warning(f"[{faa_path.name}] Diamond failed: {e}")

    # ---- Cultivation (amino acids needed) ----
    try:
        cultivation = analyzers['cultivation']
        cult_results = cultivation.analyze_genome(genome_file=input_str)
        row['Additional_amino_acids_required'] = _extract_missing_amino_acids(cult_results)
    except Exception as e:
        logger.warning(f"[{faa_path.name}] Cultivation failed: {e}")

    # ---- Antibiotic ----
    try:
        antibiotic = analyzers['antibiotic']
        ab_results = antibiotic.predict_antibiotics(
            input_file=input_str, output_prefix=output_prefix
        )
        if ab_results.get('status') == 'success':
            recommended = ab_results.get('recommended_antibiotics', [])
            row['Recommended_antibiotics'] = ', '.join(recommended) if recommended else ''
    except Exception as e:
        logger.warning(f"[{faa_path.name}] Antibiotic failed: {e}")

    return row


def batch_command(input_dir: str,
                  output_csv: str,
                  config: Config,
                  mags: bool = False):
    """
    Batch comprehensive analysis of all .faa files in a directory.

    Args:
        input_dir: Path to directory containing .faa files
        output_csv: Path to output CSV file
        config: Configuration object
        mags: If True, show pathways with completeness >= 70% (for MAGs);
              otherwise only show pathways with completeness == 100%.
    """
    logger = get_logger("batch_command")

    # Scan for .faa files
    input_path = Path(input_dir)
    if not input_path.is_dir():
        raise ValueError(f"Input is not a directory: {input_dir}")

    faa_files = sorted(
        p for p in input_path.iterdir()
        if p.suffix.lower() in ('.faa', '.fa', '.fasta') and p.is_file()
    )

    if not faa_files:
        raise FileNotFoundError(f"No FASTA files found in: {input_dir}")

    print("\n" + "=" * 60)
    print("MethArCT Batch Comprehensive Analysis")
    print("=" * 60)
    print(f"Input directory: {input_dir}")
    print(f"FASTA files found: {len(faa_files)}")
    print(f"Output CSV: {output_csv}")
    print(f"Mode: {'MAGs (pathway completeness >= 70%)' if mags else 'Isolate genome (pathway completeness == 100%)'}")
    print("=" * 60 + "\n")

    # Ensure output directory exists
    output_path = Path(output_csv)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Set analyzer output base directory
    config.set('output.base_dir', str(output_path.parent))

    # Initialize all analyzers once
    print("Initializing analyzers...")
    analyzers: Dict[str, Any] = {}

    try:
        analyzers['tome'] = TomeAnalyzer(config)
        print("  [OK] Tome (temperature)")
    except Exception as e:
        logger.warning(f"Tome unavailable: {e}")
        print(f"  [SKIP] Tome: {e}")

    try:
        ph_analyzer = PHAnalyzer(config)
        if ph_analyzer.tool_available:
            analyzers['ph'] = ph_analyzer
            print("  [OK] pH predictor")
        else:
            print("  [SKIP] pH predictor not available")
    except Exception as e:
        logger.warning(f"pH unavailable: {e}")
        print(f"  [SKIP] pH: {e}")

    try:
        susha_analyzer = SuShaAnalyzer(config)
        if susha_analyzer.tool_available:
            analyzers['susha'] = susha_analyzer
            print("  [OK] SuSha (salinity)")
        else:
            print("  [SKIP] SuSha not available")
    except Exception as e:
        logger.warning(f"SuSha unavailable: {e}")
        print(f"  [SKIP] SuSha: {e}")

    try:
        analyzers['diamond'] = DiamondAnalyzer(config)
        print("  [OK] Diamond (metabolic pathways)")
    except Exception as e:
        logger.warning(f"Diamond unavailable: {e}")
        print(f"  [SKIP] Diamond: {e}")

    try:
        analyzers['cultivation'] = CultivationAnalyzer(config.config)
        print("  [OK] Cultivation (amino acid biosynthesis)")
    except Exception as e:
        logger.warning(f"Cultivation unavailable: {e}")
        print(f"  [SKIP] Cultivation: {e}")

    try:
        threads = config.get('tools.diamond.threads', 4)
        evalue = config.get('tools.diamond.evalue', 1e-5)
        analyzers['antibiotic'] = AntibioticAnalyzer(config=config, cpus=threads, evalue=evalue)
        print("  [OK] Antibiotic resistance")
    except Exception as e:
        logger.warning(f"Antibiotic unavailable: {e}")
        print(f"  [SKIP] Antibiotic: {e}")

    if not analyzers:
        raise RuntimeError("No analyzers available. Please check your installation.")

    print(f"\nProcessing {len(faa_files)} FASTA files...\n")

    # Process each file
    all_rows: List[Dict[str, str]] = []
    total_start = time.time()

    for idx, faa_file in enumerate(faa_files, 1):
        file_start = time.time()
        print(f"[{idx}/{len(faa_files)}] {faa_file.name}")

        prefix = faa_file.stem
        row = _run_single_batch_analysis(
            faa_path=faa_file,
            output_prefix=prefix,
            config=config,
            analyzers=analyzers,
            logger=logger,
            mags=mags,
        )
        all_rows.append(row)

        elapsed = time.time() - file_start
        print(f"    Done in {elapsed:.1f}s")

    # Write CSV
    with open(output_csv, 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=_BATCH_CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(all_rows)

    total_elapsed = time.time() - total_start
    print(f"\n{'=' * 60}")
    print(f"Batch analysis completed in {total_elapsed:.1f}s")
    print(f"Results saved to: {output_csv}")
    print(f"{'=' * 60}")

    logger.info(f"Batch analysis completed: {len(faa_files)} files in {total_elapsed:.1f}s")
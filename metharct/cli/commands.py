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
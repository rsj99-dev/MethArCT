#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MethArCT Example Usage Script

This script demonstrates how to use MethArCT for various analyses:
- Core functionality: Diamond metabolic pathway analysis (required)
- Optional functionality: Tome temperature prediction, CheckM2 genome quality assessment
"""

import os
import sys
from pathlib import Path

project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from metharct.core.diamond_analyzer import DiamondAnalyzer
from metharct.core.tome_analyzer import TomeAnalyzer
from metharct.core.checkm2_analyzer import CheckM2Analyzer
from metharct.core.pathway_predictor import PathwayPredictor
from metharct.utils.config import Config
from metharct.utils.logger import setup_logger

def main():
    """
    Main function: Demonstrates basic usage of MethArCT

    Note:
    - Diamond analysis is the core functionality and must be installed
    - Tome and CheckM2 are optional; if not installed, corresponding analyses will be skipped
    """
    logger = setup_logger("example", level="INFO")
    logger.info("Starting MethArCT example analysis")

    config = Config()

    input_protein_file = "example_protein.fasta"
    input_genome_file = "example_genome.fasta"
    output_dir = "example_results"

    os.makedirs(output_dir, exist_ok=True)

    if not os.path.exists(input_protein_file):
        logger.warning(f"Protein file {input_protein_file} not found, skipping related analysis")
        print("Please rename your protein sequence file to 'example_protein.fasta' or modify the file path in this script")

    if not os.path.exists(input_genome_file):
        logger.warning(f"Genome file {input_genome_file} not found, skipping related analysis")
        print("Please rename your genome sequence file to 'example_genome.fasta' or modify the file path in this script")

    try:
        if os.path.exists(input_protein_file):
            logger.info("Starting Diamond metabolic pathway analysis (core functionality)")
            diamond_analyzer = DiamondAnalyzer(config)

            if diamond_analyzer.check_availability():
                diamond_results = diamond_analyzer.analyze_file(
                    input_file=input_protein_file,
                    output_dir=os.path.join(output_dir, "diamond"),
                    threads=4,
                    evalue=1e-5
                )
                logger.info(f"Diamond analysis completed, results saved to: {diamond_results['output_file']}")
                print(f"Detected metabolic pathways: {diamond_results['pathways']}")
            else:
                logger.error("Diamond tool not available, please check installation (this is a required core functionality)")

        if os.path.exists(input_protein_file):
            logger.info("Starting Tome temperature prediction analysis (optional functionality)")
            tome_analyzer = TomeAnalyzer(config)

            if tome_analyzer.check_availability():
                tome_results = tome_analyzer.predict_file(
                    input_file=input_protein_file,
                    output_dir=os.path.join(output_dir, "tome")
                )
                logger.info(f"Tome analysis completed, results saved to: {tome_results['output_file']}")
                print(f"Predicted optimal growth temperature: {tome_results['predicted_ogt']}°C")
                print(f"Temperature classification: {tome_results['temperature_class']}")
            else:
                logger.warning("Tome tool not available, skipping temperature prediction (optional functionality)")

        if os.path.exists(input_genome_file):
            logger.info("Starting CheckM2 genome quality assessment (optional functionality)")
            checkm2_analyzer = CheckM2Analyzer(config)

            if checkm2_analyzer.check_availability():
                checkm2_results = checkm2_analyzer.analyze_file(
                    input_file=input_genome_file,
                    output_dir=os.path.join(output_dir, "checkm2"),
                    threads=4
                )
                logger.info(f"CheckM2 analysis completed, results saved to: {checkm2_results['output_file']}")
                print(f"Genome completeness: {checkm2_results['completeness']:.2f}%")
                print(f"Contamination: {checkm2_results['contamination']:.2f}%")
                print(f"Quality grade: {checkm2_results['quality_grade']}")
            else:
                logger.warning("CheckM2 tool not available, skipping genome quality assessment (optional functionality)")

        if os.path.exists(input_protein_file) and os.path.exists(input_genome_file):
            logger.info("Starting comprehensive analysis (requires all tools to be available)")
            pathway_predictor = PathwayPredictor(config)

            try:
                comprehensive_results = pathway_predictor.comprehensive_analysis(
                    protein_file=input_protein_file,
                    genome_file=input_genome_file,
                    output_dir=os.path.join(output_dir, "comprehensive"),
                    threads=4
                )

                logger.info(f"Comprehensive analysis completed, results saved to: {comprehensive_results['output_dir']}")
                print("\n=== Comprehensive Analysis Results Summary ===")
                print(f"Metabolic pathways: {comprehensive_results.get('pathways', [])}")
                print(f"Salt tolerance: {comprehensive_results.get('salt_tolerance', 'Unknown')}")
                print(f"Optimal growth temperature: {comprehensive_results.get('optimal_temperature', 'Unknown')}°C")
                print(f"Culturability: {comprehensive_results.get('culturability', 'Unknown')}")
            except Exception as e:
                logger.warning(f"Comprehensive analysis failed, possibly due to missing optional tools: {str(e)}")
                print("\nNote: Comprehensive analysis requires all tools (Diamond, Tome, CheckM2) to be available")

        logger.info("All analyses completed!")
        print(f"\nAll results have been saved in the '{output_dir}' directory")

    except Exception as e:
        logger.error(f"Error during analysis: {str(e)}")
        print(f"Error: {str(e)}")
        return 1

    return 0

def show_help():
    """
    Display help information
    """
    help_text = """
    MethArCT Example Usage Script

    Usage:
        python example_usage.py

    Prerequisites:
    1. Required tools:
       - Diamond (core functionality, must be installed)
    2. Optional tools (for extended functionality):
       - Tome (temperature prediction)
       - CheckM2 (genome quality assessment)
    3. Prepare input files:
       - example_protein.fasta: Protein sequence file
       - example_genome.fasta: Genome sequence file (required for CheckM2 analysis)
    4. Run the script for analysis

    Output Results:
    - example_results/diamond/: Diamond analysis results (core functionality)
    - example_results/tome/: Tome analysis results (optional functionality)
    - example_results/checkm2/: CheckM2 analysis results (optional functionality)
    - example_results/comprehensive/: Comprehensive analysis results (requires all tools)

    Notes:
    - Only Diamond is required, other tools are optional
    - If optional tools are missing, corresponding analyses will be skipped
    - Please modify the input file names in the script according to your actual file paths
    - Ensure sufficient disk space for results
    - Analysis time depends on input file size and system performance
    """
    print(help_text)

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] in ["-h", "--help", "help"]:
        show_help()
    else:
        exit_code = main()
        sys.exit(exit_code)

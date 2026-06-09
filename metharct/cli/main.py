#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MethArCT Command Line Interface Main Module

Provides the main entry point for the MethArCT command-line tool.
"""

import argparse
import sys
import os
from pathlib import Path
from typing import Optional

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

try:
    from metharct.utils.config import Config
    from metharct.utils.logger import setup_logger
    from .commands import (
        diamond_command,
        tome_command,
        checkm2_command,
        comprehensive_command,
        susha_command,
        ph_command,
        antibiotic_command
    )
except ImportError as e:
    print(f"Error importing modules: {e}")
    print("Please make sure MethArCT is properly installed.")
    sys.exit(1)

def create_parser() -> argparse.ArgumentParser:
    """
    Create the main argument parser
    
    Returns:
        Configured ArgumentParser instance
    """
    parser = argparse.ArgumentParser(
        prog='metharct',
        description='MethArCT: Metabolic pathway prediction and cultivability assessment tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Comprehensive analysis (recommended)
  metharct comprehensive input.fasta -o results
  
  # Individual analyses
  metharct diamond input.fasta -o diamond_results
  metharct tome input.fasta -o tome_results
  metharct susha input.fasta -o susha_results
  metharct ph input.fasta -o ph_results
  metharct antibiotic input.fasta -o antibiotic_results
  metharct checkm2 input.fasta -o checkm2_results
  
  # With custom configuration
  metharct comprehensive input.fasta -c config.yaml -o results
  
  # Verbose output
  metharct comprehensive input.fasta -o results -v

For more information, visit: https://github.com/rsj99-dev/MethArCT
"""
    )
    
    # Global options
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    
    parser.add_argument(
        '-c', '--config',
        type=str,
        help='Path to configuration file (YAML format)'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='MethArCT 0.5.5'
    )
    
    # Create subparsers for different commands
    subparsers = parser.add_subparsers(
        dest='command',
        help='Available commands',
        metavar='COMMAND'
    )
    
    # Comprehensive analysis command
    comprehensive_parser = subparsers.add_parser(
        'comprehensive',
        help='Run comprehensive analysis (Diamond + Tome + SuSha + pH + CheckM2)',
        description='Perform comprehensive metabolic pathway prediction, temperature analysis, salinity prediction, pH prediction, and genome quality assessment'
    )
    _add_comprehensive_args(comprehensive_parser)
    
    # Diamond analysis command
    diamond_parser = subparsers.add_parser(
        'diamond',
        help='Run Diamond analysis for metabolic pathways and cultivability',
        description='Analyze metabolic pathways and cultivability using Diamond BLAST'
    )
    _add_diamond_args(diamond_parser)
    
    # Tome analysis command
    tome_parser = subparsers.add_parser(
        'tome',
        help='Run Tome analysis for optimal growth temperature prediction',
        description='Predict optimal growth temperature using Tome tool'
    )
    _add_tome_args(tome_parser)
    
    # CheckM2 analysis command
    checkm2_parser = subparsers.add_parser(
        'checkm2',
        help='Run CheckM2 analysis for genome quality assessment',
        description='Assess genome quality and completeness using CheckM2'
    )
    _add_checkm2_args(checkm2_parser)
    
    # Cultivation analysis command
    cultivation_parser = subparsers.add_parser(
        'cultivation',
        help='Run cultivation analysis for cultivability assessment',
        description='Assess cultivability based on metabolic pathways'
    )
    _add_cultivation_args(cultivation_parser)
    
    # SuSha salinity prediction command
    susha_parser = subparsers.add_parser(
        'susha',
        help='Run SuSha salinity adaptation prediction',
        description='Predict microbial salinity adaptation using SuSha ensemble learning model'
    )
    _add_susha_args(susha_parser)

    # pH preference prediction command
    ph_parser = subparsers.add_parser(
        'ph',
        help='Run pH preference prediction',
        description='Predict microbial growth pH preference (optimum, max, min) using GenomeSpot models'
    )
    _add_ph_args(ph_parser)

    # Antibiotic resistance prediction command
    antibiotic_parser = subparsers.add_parser(
        'antibiotic',
        help='Run antibiotic resistance prediction based on AAI',
        description='Predict antibiotic resistance by comparing query proteins against methanogen reference genomes using AAI thresholds'
    )
    _add_antibiotic_args(antibiotic_parser)

    return parser

def _add_comprehensive_args(parser: argparse.ArgumentParser):
    """Add arguments for comprehensive analysis command"""
    parser.add_argument(
        'input',
        type=str,
        help='Input FASTA file path'
    )
    
    parser.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help='Output directory or file prefix'
    )
    
    parser.add_argument(
        '--skip-tome',
        action='store_true',
        help='Skip Tome temperature analysis'
    )
    
    parser.add_argument(
        '--skip-checkm2',
        action='store_true',
        help='Skip CheckM2 genome quality analysis'
    )
    
    parser.add_argument(
        '--skip-susha',
        action='store_true',
        help='Skip SuSha salinity prediction'
    )

    parser.add_argument(
        '--skip-ph',
        action='store_true',
        help='Skip pH preference prediction'
    )

    parser.add_argument(
        '--skip-antibiotic',
        action='store_true',
        help='Skip antibiotic resistance prediction'
    )
    
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of threads to use (default: 4)'
    )
    
    parser.add_argument(
        '--evalue',
        type=float,
        default=1e-5,
        help='E-value threshold for Diamond (default: 1e-5)'
    )

def _add_diamond_args(parser: argparse.ArgumentParser):
    """Add arguments for Diamond analysis command"""
    parser.add_argument(
        'input',
        type=str,
        help='Input FASTA file path'
    )
    
    parser.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help='Output directory or file prefix'
    )
    
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of threads to use (default: 4)'
    )
    
    parser.add_argument(
        '--evalue',
        type=float,
        default=1e-5,
        help='E-value threshold (default: 1e-5)'
    )
    
    parser.add_argument(
        '--max-target-seqs',
        type=int,
        default=10,
        help='Maximum number of target sequences (default: 10)'
    )

def _add_tome_args(parser: argparse.ArgumentParser):
    """Add arguments for Tome analysis command"""
    parser.add_argument(
        'input',
        type=str,
        help='Input FASTA file path (protein sequences)'
    )
    
    parser.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help='Output directory or file prefix'
    )
    
    parser.add_argument(
        '--batch-size',
        type=int,
        default=100,
        help='Batch size for processing sequences (default: 100)'
    )

def _add_cultivation_args(parser: argparse.ArgumentParser):
    """
    Add arguments for cultivation command
    
    Args:
        parser: Parser to add arguments to
    """
    parser.add_argument(
        "input",
        help="Path to input FASTA file"
    )
    
    parser.add_argument(
        "-o", "--output",
        default="cultivation_output",
        help="Output prefix for result files (default: cultivation_output)"
    )
    
    parser.add_argument(
        "--diamond-path",
        help="Path to Diamond executable (overrides config)"
    )
    
    parser.add_argument(
        "--threads",
        type=int,
        help="Number of threads for Diamond (overrides config)"
    )
    
    parser.add_argument(
        "--evalue",
        type=float,
        help="E-value threshold for Diamond (overrides config)"
    )
    
    parser.add_argument(
        "--use-wsl",
        action="store_true",
        help="Use WSL for Diamond (overrides config)"
    )
    
    parser.add_argument(
        "--pathway-weights",
        help="JSON file with pathway weights (overrides config)"
    )
    
    parser.add_argument(
        "--cultivability-threshold",
        type=float,
        help="Cultivability score threshold (overrides config)"
    )

def _add_susha_args(parser: argparse.ArgumentParser):
    """Add arguments for SuSha salinity prediction command"""
    parser.add_argument(
        'input',
        type=str,
        help='Input FASTA file path (protein sequences)'
    )

    parser.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help='Output file prefix'
    )

def _add_ph_args(parser: argparse.ArgumentParser):
    """Add arguments for pH preference prediction command"""
    parser.add_argument(
        'input',
        type=str,
        help='Input FASTA file path (protein sequences)'
    )

    parser.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help='Output file prefix'
    )

def _add_antibiotic_args(parser: argparse.ArgumentParser):
    """Add arguments for antibiotic resistance prediction command"""
    parser.add_argument(
        'input',
        type=str,
        help='Input FASTA file path (protein sequences)'
    )

    parser.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help='Output file prefix'
    )

    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of threads for DIAMOND (default: 4)'
    )

    parser.add_argument(
        '--evalue',
        type=float,
        default=1e-5,
        help='E-value threshold for DIAMOND (default: 1e-5)'
    )

def _add_checkm2_args(parser: argparse.ArgumentParser):
    """
    Add arguments for checkm2 command
    
    Args:
        parser: Parser to add arguments to
    """
    parser.add_argument(
        "input",
        help="Path to input FASTA file or directory"
    )
    
    parser.add_argument(
        "-o", "--output",
        default="checkm2_output",
        help="Output prefix for result files (default: checkm2_output)"
    )
    
    parser.add_argument(
        "-t", "--type",
        choices=["fasta", "directory"],
        default="fasta",
        help="Input type: fasta file or directory (default: fasta)"
    )
    
    parser.add_argument(
        "--checkm2-path",
        help="Path to CheckM2 executable (overrides config)"
    )
    
    parser.add_argument(
        "--threads",
        type=int,
        help="Number of threads for CheckM2 (overrides config)"
    )
    
    parser.add_argument(
        "--use-wsl",
        action="store_true",
        help="Use WSL for CheckM2 (overrides config)"
    )
    
    parser.add_argument(
        "--input-type",
        choices=['fasta', 'directory'],
        default='fasta',
        help='Input type: single fasta file or directory of genomes (default: fasta)'
    )

def setup_logging(verbose: bool = False):
    """
    Setup logging configuration
    
    Args:
        verbose: Enable verbose logging
    """
    log_level = 'DEBUG' if verbose else 'INFO'
    setup_logger(
        name='metharct',
        level=log_level,
        log_file='metharct.log',
        console=True
    )

def load_config(config_path: Optional[str] = None) -> Config:
    """
    Load configuration
    
    Args:
        config_path: Path to configuration file
        
    Returns:
        Configuration object
    """
    config = Config(config_path)
    
    if config_path:
        if os.path.exists(config_path):
            print(f"Loaded configuration from {config_path}")
        else:
            print(f"Warning: Configuration file {config_path} not found, using defaults")
    
    return config

def validate_input_file(input_path: str) -> Path:
    """
    Validate input file exists and is readable
    
    Args:
        input_path: Path to input file
        
    Returns:
        Validated Path object
        
    Raises:
        FileNotFoundError: If file doesn't exist
        PermissionError: If file is not readable
    """
    path = Path(input_path)
    
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    if not path.is_file():
        raise ValueError(f"Input path is not a file: {input_path}")
    
    if not os.access(path, os.R_OK):
        raise PermissionError(f"Input file is not readable: {input_path}")
    
    return path

def create_output_directory(output_path: str) -> Path:
    """
    Create output directory if it doesn't exist
    
    Args:
        output_path: Output directory path
        
    Returns:
        Path object for output directory
    """
    path = Path(output_path)
    
    # If output_path looks like a file (has extension), use its parent directory
    if path.suffix:
        output_dir = path.parent
        output_prefix = path.stem
    else:
        output_dir = path
        output_prefix = None
    
    # Create directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    
    return output_dir, output_prefix

def main():
    """
    Main entry point for the MethArCT CLI
    """
    try:
        parser = create_parser()
        
        # Parse arguments
        args = parser.parse_args()
        
        # If no command specified, show help
        if not args.command:
            parser.print_help()
            sys.exit(1)
        
        try:
            # Setup logging
            setup_logging(args.verbose)
            
            # Load configuration
            config = load_config(args.config)
            
            # Validate input file (except for some commands that might not need it)
            input_path = None
            if hasattr(args, 'input'):
                input_path = validate_input_file(args.input)
            
            # Create output directory
            output_dir = None
            output_prefix = None
            if hasattr(args, 'output'):
                output_dir, output_prefix = create_output_directory(args.output)
                
                # Use the original output argument as prefix if it looks like a filename
                if output_prefix is None:
                    output_prefix = Path(args.input).stem if hasattr(args, 'input') else 'metharct_results'
            
            # Update config with command-line arguments
            if hasattr(args, 'threads'):
                config.set('tools.diamond.threads', args.threads)
                config.set('tools.checkm2.threads', args.threads)
            
            if hasattr(args, 'evalue'):
                config.set('tools.diamond.evalue', args.evalue)
            
            if hasattr(args, 'max_target_seqs'):
                config.set('tools.diamond.max_target_seqs', args.max_target_seqs)
            
            # Set output directory in config
            if hasattr(args, 'output') and output_dir:
                config.set('output.base_dir', str(output_dir))
            
            # Execute the appropriate command
            if args.command == 'comprehensive':
                comprehensive_command(
                    input_path=str(input_path),
                    output_prefix=output_prefix,
                    config=config,
                    skip_tome=args.skip_tome,
                    skip_checkm2=args.skip_checkm2,
                    skip_susha=getattr(args, 'skip_susha', False),
                    skip_ph=getattr(args, 'skip_ph', False),
                    skip_antibiotic=getattr(args, 'skip_antibiotic', False)
                )
            
            elif args.command == 'diamond':
                diamond_command(
                    input_path=str(input_path),
                    output_prefix=output_prefix,
                    config=config
                )
            
            elif args.command == 'tome':
                tome_command(
                    input_path=str(input_path),
                    output_prefix=output_prefix,
                    config=config,
                    batch_size=getattr(args, 'batch_size', 100)
                )
            
            elif args.command == 'checkm2':
                checkm2_command(
                    input_path=str(input_path),
                    output_prefix=output_prefix,
                    config=config,
                    input_type=getattr(args, 'input_type', 'fasta')
                )
            
            elif args.command == 'cultivation':
                from .commands import cultivation_command
                
                # Override config with command line arguments
                if hasattr(args, 'diamond_path') and args.diamond_path:
                    config.set('cultivation.diamond_path', args.diamond_path)
                
                if hasattr(args, 'threads') and args.threads:
                    config.set('cultivation.threads', args.threads)
                
                if hasattr(args, 'evalue') and args.evalue:
                    config.set('cultivation.evalue', args.evalue)
                
                if hasattr(args, 'use_wsl') and args.use_wsl:
                    config.set('cultivation.use_wsl', args.use_wsl)
                
                if hasattr(args, 'pathway_weights') and args.pathway_weights:
                    config.set('cultivation.pathway_weights', args.pathway_weights)
                
                if hasattr(args, 'cultivability_threshold') and args.cultivability_threshold:
                    config.set('cultivation.cultivability_threshold', args.cultivability_threshold)
                
                cultivation_command(
                    input_path=str(input_path),
                    output_prefix=output_prefix,
                    config=config
                )
            
            elif args.command == 'susha':
                susha_command(
                    input_path=str(input_path),
                    output_prefix=output_prefix,
                    config=config
                )

            elif args.command == 'ph':
                ph_command(
                    input_path=str(input_path),
                    output_prefix=output_prefix,
                    config=config
                )

            elif args.command == 'antibiotic':
                # Override config with command line arguments
                if hasattr(args, 'threads') and args.threads:
                    config.set('tools.diamond.threads', args.threads)
                if hasattr(args, 'evalue') and args.evalue:
                    config.set('tools.diamond.evalue', args.evalue)

                antibiotic_command(
                    input_path=str(input_path),
                    output_prefix=output_prefix,
                    config=config
                )
            
            print("\nAnalysis completed successfully!")
            if output_dir:
                print(f"Results saved to: {output_dir}")
            
        except KeyboardInterrupt:
            print("\nAnalysis interrupted by user")
            sys.exit(1)
        
        except Exception as e:
            print(f"\nError: {str(e)}")
            if hasattr(args, 'verbose') and args.verbose:
                import traceback
                traceback.print_exc()
            sys.exit(1)
    except Exception as e:
        print(f"\nCritical error in command line interface: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()

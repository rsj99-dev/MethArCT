# -*- coding: utf-8 -*-
"""
MethArCT Main Analyzer

Primary analysis class integrating Diamond, CheckM2, Antibiotic and other analysis tools
"""

import os
import json
from pathlib import Path
from typing import Dict, List, Optional, Union

from .core.diamond_analyzer import DiamondAnalyzer
from .core.checkm2_analyzer import CheckM2Analyzer
from .core.antibiotic_analyzer import AntibioticAnalyzer
from .core.huainanzi_analyzer import HuainanziAnalyzer
from .core.pathway_predictor import PathwayPredictor
from .utils.config import Config
from .utils.logger import get_logger
from .utils.file_utils import FileUtils, NumpyEncoder


class MethArCTAnalyzer:
    """
    MethArCT Main Analyzer Class
    
    Main interface class integrating all analysis functionalities
    """
    
    def __init__(self, config: Optional[Config] = None):
        """
        Initialize MethArCT analyzer
        
        Args:
            config: Configuration object, uses default if None
        """
        self.config = config or Config()
        self.logger = get_logger("metharct_analyzer")
        
        # Initialize each analyzer
        self.diamond_analyzer = DiamondAnalyzer(self.config)
        self.checkm2_analyzer = CheckM2Analyzer(self.config)
        self.antibiotic_analyzer = AntibioticAnalyzer(self.config)
        self.huainanzi_analyzer = HuainanziAnalyzer(self.config)
        self.pathway_predictor = PathwayPredictor(self.config)
        
        # Results output directory
        self.results_dir = self.config.get('output.base_dir', 'results')
        FileUtils.ensure_dir(self.results_dir)
        
        self.logger.info("MethArCT analyzer initialized")
    
    def analyze_comprehensive(self, 
                            input_file: str,
                            output_dir: Optional[str] = None,
                            run_diamond: bool = True,
                            run_checkm2: bool = True,
                            run_antibiotic: bool = True,
                            **kwargs) -> Dict:
        """
        Execute comprehensive analysis

        The analysis is orchestrated by ``PathwayPredictor`` only once; the
        previous implementation ran DIAMOND and the antibiotic prediction twice
        and called a non-existent CheckM2 method.
        
        Args:
            input_file: Input file path (FASTA format)
            output_dir: Output directory, uses default if None
            run_diamond: Whether to run Diamond analysis
            run_checkm2: Whether to run CheckM2 analysis
            run_antibiotic: Whether to run antibiotic resistance prediction
            **kwargs: Other parameters forwarded to PathwayPredictor
                (run_cultivation, run_susha, run_huainanzi, run_ph)
            
        Returns:
            Dictionary containing all analysis results
        """
        if output_dir is None:
            output_dir = self.results_dir
            
        FileUtils.ensure_dir(output_dir)
        output_prefix = os.path.join(output_dir, Path(input_file).stem)
        
        results = {
            'input_file': input_file,
            'output_dir': output_dir,
            'analysis_results': {}
        }
        
        self.logger.info(f"Starting comprehensive analysis: {input_file}")
        
        try:
            comprehensive_results = self.pathway_predictor.predict_comprehensive(
                input_path=input_file,
                output_prefix=output_prefix,
                include_diamond=run_diamond,
                include_checkm2=run_checkm2,
                include_cultivation=kwargs.get('run_cultivation', True),
                include_susha=kwargs.get('run_susha', True),
                include_huainanzi=kwargs.get('run_huainanzi', True),
                include_ph=kwargs.get('run_ph', True),
                include_antibiotic=run_antibiotic,
            )

            results['analysis_results'] = comprehensive_results.get('results', {})
            results['integrated_analysis'] = comprehensive_results.get(
                'integrated_analysis', {}
            )
            results['analyses_performed'] = comprehensive_results.get(
                'analyses_performed', []
            )
                
            # Save comprehensive results
            results_file = os.path.join(output_dir, 'comprehensive_results.json')
            with open(results_file, 'w', encoding='utf-8') as f:
                json.dump(results, f, indent=2, ensure_ascii=False, cls=NumpyEncoder)
                
            self.logger.info(f"Comprehensive analysis completed, results saved to: {results_file}")
            
        except Exception as e:
            self.logger.error(f"Comprehensive analysis failed: {str(e)}")
            results['error'] = str(e)
            raise
            
        return results
    
    def analyze_diamond(self, input_file: str, output_dir: Optional[str] = None) -> Dict:
        """
        Execute Diamond analysis only
        
        Args:
            input_file: Input file path
            output_dir: Output directory
            
        Returns:
            Diamond analysis results
        """
        if output_dir is None:
            output_dir = os.path.join(self.results_dir, 'diamond')
            
        output_prefix = os.path.join(output_dir, Path(input_file).stem)
        FileUtils.ensure_dir(output_dir)
        return self.diamond_analyzer.analyze_sequence(input_file, output_prefix)
    
    def analyze_checkm2(self, input_file: str, output_dir: Optional[str] = None) -> Dict:
        """
        Execute CheckM2 analysis only
        
        Args:
            input_file: Input file path
            output_dir: Output directory
            
        Returns:
            CheckM2 analysis results
        """
        if output_dir is None:
            output_dir = os.path.join(self.results_dir, 'checkm2')
            
        FileUtils.ensure_dir(output_dir)
        output_prefix = os.path.join(output_dir, Path(input_file).stem)
        return self.checkm2_analyzer.analyze_genome_quality(
            input_file, output_prefix=output_prefix
        )
    
    def analyze_antibiotic(self, input_file: str, output_dir: Optional[str] = None) -> Dict:
        """
        Execute antibiotic resistance prediction only
    
        Args:
            input_file: Input file path
            output_dir: Output directory
    
        Returns:
            Antibiotic prediction results
        """
        if output_dir is None:
            output_dir = os.path.join(self.results_dir, 'antibiotic')
    
        FileUtils.ensure_dir(output_dir)
        return self.antibiotic_analyzer.predict_antibiotics(
            input_file=input_file,
            output_prefix=os.path.join(output_dir, Path(input_file).stem),
        )
    
    def get_analysis_summary(self, results: Dict) -> Dict:
        """
        Generate analysis results summary
        
        Args:
            results: Analysis results dictionary
            
        Returns:
            Results summary dictionary
        """
        summary = {
            'input_file': results.get('input_file'),
            'total_analyses': len(results.get('analysis_results', {})),
            'completed_analyses': [],
            'failed_analyses': []
        }
        
        for analysis_type, result in results.get('analysis_results', {}).items():
            if 'error' in result:
                summary['failed_analyses'].append(analysis_type)
            else:
                summary['completed_analyses'].append(analysis_type)
                
        return summary
    
    def validate_input(self, input_file: str) -> bool:
        """
        Validate input file
        
        Args:
            input_file: Input file path
            
        Returns:
            Whether validation passed
        """
        if not os.path.exists(input_file):
            self.logger.error(f"Input file does not exist: {input_file}")
            return False
            
        # Check file format (simple check)
        try:
            with open(input_file, 'r') as f:
                first_line = f.readline().strip()
                if not first_line.startswith('>'):
                    self.logger.error(f"Input file is not valid FASTA format: {input_file}")
                    return False
        except Exception as e:
            self.logger.error(f"Failed to read input file: {str(e)}")
            return False
            
        return True
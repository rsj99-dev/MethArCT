#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pathway predictor for MethArCT

Integrates Diamond, Cultivation, Tome and CheckM2 analysis results to provide comprehensive
metabolic pathway predictions and cultivability assessments.
"""

import os
import json
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple
from ..utils.config import Config
from ..utils.logger import get_logger
from ..utils.file_utils import FileUtils
from .diamond_analyzer import DiamondAnalyzer
from .cultivation_analyzer import CultivationAnalyzer
from .tome_analyzer import TomeAnalyzer
from .checkm2_analyzer import CheckM2Analyzer

class PathwayPredictor:
    """Comprehensive pathway predictor integrating multiple analysis tools"""
    
    def __init__(self, config: Optional[Config] = None):
        self.config = config or Config()
        self.logger = get_logger("pathway_predictor")
        
        # Initialize analyzers
        self.diamond_analyzer = DiamondAnalyzer(self.config)
        self.cultivation_analyzer = CultivationAnalyzer(self.config.config)
        self.tome_analyzer = TomeAnalyzer(self.config)
        self.checkm2_analyzer = CheckM2Analyzer(self.config)
        
        # Results directory
        self.results_dir = self.config.get('output.base_dir', 'results')
        FileUtils.ensure_dir(self.results_dir)
        
        # Pathway confidence thresholds
        self.confidence_thresholds = {
            'high': 0.8,
            'medium': 0.5,
            'low': 0.2
        }
    
    def predict_comprehensive(self, 
                            input_path: Union[str, Path],
                            output_prefix: Optional[str] = None,
                            include_tome: bool = True,
                            include_checkm2: bool = True,
                            include_cultivation: bool = True) -> Dict[str, any]:
        """
        Perform comprehensive pathway prediction analysis
        
        Args:
            input_path: Path to input FASTA file
            output_prefix: Prefix for output files
            include_tome: Whether to include Tome analysis
            include_checkm2: Whether to include CheckM2 analysis
            include_cultivation: Whether to include Cultivation analysis
            
        Returns:
            Comprehensive analysis results
        """
        input_path = Path(input_path)
        
        if not input_path.exists():
            raise FileNotFoundError(f"Input file not found: {input_path}")
        
        # Set output prefix
        if output_prefix is None:
            output_prefix = input_path.stem
        
        self.logger.info(f"Starting comprehensive analysis for {input_path.name}")
        
        # Initialize results dictionary
        comprehensive_results = {
            'input_path': str(input_path),
            'analysis_timestamp': pd.Timestamp.now().isoformat(),
            'output_prefix': output_prefix,
            'analyses_performed': [],
            'results': {}
        }
        
        # 1. Diamond analysis (metabolic pathways, salt tolerance, cultivability)
        try:
            self.logger.info("Running Diamond analysis...")
            diamond_results = self.diamond_analyzer.analyze_sequence(
                input_file=input_path, output_prefix=output_prefix
            )
            comprehensive_results['analyses_performed'].append('diamond')
            comprehensive_results['results']['diamond'] = diamond_results
            self.logger.info("Diamond analysis completed")
        except Exception as e:
            self.logger.error(f"Diamond analysis failed: {str(e)}")
            comprehensive_results['results']['diamond'] = {
                'status': 'failed',
                'error': str(e)
            }
        
        # 1.5. Cultivation analysis (metabolic pathway-based cultivability assessment)
        if include_cultivation:
            try:
                self.logger.info("Running Cultivation analysis...")
                cultivation_results = self.cultivation_analyzer.analyze_genome(str(input_path))
                cultivability_assessment = self.cultivation_analyzer.generate_cultivability_assessment(cultivation_results)
                
                # Export cultivation results
                cultivation_output = os.path.join(self.results_dir, f"{output_prefix}_cultivation")
                self.cultivation_analyzer.export_results(
                    cultivation_results, 
                    cultivability_assessment, 
                    cultivation_output + ".csv", 
                    'csv'
                )
                
                # Amino acid biosynthesis analysis (20 standard amino acids)
                aa_analysis = self.cultivation_analyzer.analyze_amino_acid_biosynthesis(cultivation_results)
                self.cultivation_analyzer.print_amino_acid_report(aa_analysis)
                
                # Export amino acid CSV
                aa_csv_output = os.path.join(self.results_dir, f"{output_prefix}_amino_acid_biosynthesis.csv")
                self.cultivation_analyzer.export_amino_acid_csv(aa_analysis, aa_csv_output)
                
                comprehensive_results['analyses_performed'].append('cultivation')
                comprehensive_results['results']['cultivation'] = {
                    'pathway_results': cultivation_results,
                    'cultivability_assessment': cultivability_assessment,
                    'amino_acid_biosynthesis': aa_analysis
                }
                self.logger.info("Cultivation analysis completed")
            except Exception as e:
                self.logger.error(f"Cultivation analysis failed: {str(e)}")
                comprehensive_results['results']['cultivation'] = {
                    'status': 'failed',
                    'error': str(e)
                }
        
        # 2. Tome analysis (optimal growth temperature)
        if include_tome:
            try:
                self.logger.info("Running Tome analysis...")
                tome_results = self.tome_analyzer.predict_ogt(
                    input_file=input_path, output_prefix=output_prefix
                )
                comprehensive_results['analyses_performed'].append('tome')
                comprehensive_results['results']['tome'] = tome_results
                self.logger.info("Tome analysis completed")
            except Exception as e:
                self.logger.error(f"Tome analysis failed: {str(e)}")
                comprehensive_results['results']['tome'] = {
                    'status': 'failed',
                    'error': str(e)
                }
        
        # 3. CheckM2 analysis (genome quality and cultivability)
        if include_checkm2:
            try:
                self.logger.info("Running CheckM2 analysis...")
                checkm2_results = self.checkm2_analyzer.analyze_genome_quality(
                    input_path, output_prefix
                )
                comprehensive_results['analyses_performed'].append('checkm2')
                comprehensive_results['results']['checkm2'] = checkm2_results
                self.logger.info("CheckM2 analysis completed")
            except Exception as e:
                self.logger.error(f"CheckM2 analysis failed: {str(e)}")
                comprehensive_results['results']['checkm2'] = {
                    'status': 'failed',
                    'error': str(e)
                }
        
        # 4. Integrate results
        try:
            integrated_results = self._integrate_results(comprehensive_results['results'])
            comprehensive_results['integrated_analysis'] = integrated_results
            self.logger.info("Results integration completed")
        except Exception as e:
            self.logger.error(f"Results integration failed: {str(e)}")
            comprehensive_results['integrated_analysis'] = {
                'status': 'failed',
                'error': str(e)
            }
        
        # 5. Save comprehensive results
        self._save_comprehensive_results(comprehensive_results, output_prefix)
        
        return comprehensive_results
    
    def _integrate_results(self, analysis_results: Dict) -> Dict:
        """
        Integrate results from multiple analyses
        
        Args:
            analysis_results: Dictionary containing results from different analyses
            
        Returns:
            Integrated analysis results
        """
        integrated = {
            'metabolic_profile': {},
            'environmental_adaptation': {},
            'cultivability_assessment': {},
            'overall_assessment': {},
            'recommendations': []
        }
        
        # Process Diamond results
        if 'diamond' in analysis_results and 'status' not in analysis_results['diamond']:
            diamond_data = analysis_results['diamond']
            
            # Metabolic pathways
            if 'pathway_analysis' in diamond_data:
                integrated['metabolic_profile'] = self._process_metabolic_pathways(
                    diamond_data['pathway_analysis']
                )
            
            # Salt tolerance
            if 'salt_tolerance' in diamond_data:
                integrated['environmental_adaptation']['salt_tolerance'] = \
                    diamond_data['salt_tolerance']
        
        # Process Tome results
        if 'tome' in analysis_results and 'status' not in analysis_results['tome']:
            tome_data = analysis_results['tome']
            
            if 'temperature_predictions' in tome_data:
                integrated['environmental_adaptation']['temperature'] = \
                    self._process_temperature_data(tome_data['temperature_predictions'])
        
        # Process CheckM2 results
        if 'checkm2' in analysis_results and 'status' not in analysis_results['checkm2']:
            checkm2_data = analysis_results['checkm2']
            
            if 'summary' in checkm2_data:
                integrated['cultivability_assessment']['genome_quality'] = \
                    checkm2_data['summary']
        
        # Process Cultivation results
        if 'cultivation' in analysis_results and 'status' not in analysis_results['cultivation']:
            cultivation_data = analysis_results['cultivation']
            
            if 'cultivability_assessment' in cultivation_data:
                integrated['cultivability_assessment']['metabolic_pathways'] = \
                    cultivation_data['cultivability_assessment']
            
            # Amino acid biosynthesis analysis
            if 'amino_acid_biosynthesis' in cultivation_data:
                aa_data = cultivation_data['amino_acid_biosynthesis']
                integrated['amino_acid_biosynthesis'] = {
                    'summary': aa_data.get('summary', {}),
                    'amino_acid_status': aa_data.get('amino_acid_status', {})
                }
        
        # Generate overall assessment
        integrated['overall_assessment'] = self._generate_overall_assessment(integrated)
        
        # Generate recommendations
        integrated['recommendations'] = self._generate_recommendations(integrated)
        
        return integrated
    
    def _process_metabolic_pathways(self, pathway_data: Dict) -> Dict:
        """
        Process metabolic pathway data
        
        Args:
            pathway_data: Raw pathway analysis data
            
        Returns:
            Processed metabolic profile
        """
        metabolic_profile = {
            'methane_metabolism': {},
            'sulfur_metabolism': {},
            'nitrogen_metabolism': {},
            'overall_metabolic_diversity': 0,
            'key_pathways': []
        }
        
        # Process each pathway category
        for category, pathways in pathway_data.items():
            if category in ['methane_pathways', 'sulfur_pathways', 'nitrogen_pathways']:
                category_key = category.replace('_pathways', '_metabolism')
                
                pathway_summary = {
                    'detected_pathways': [],
                    'pathway_count': 0,
                    'confidence_scores': [],
                    'dominant_pathway': None
                }
                
                if isinstance(pathways, dict):
                    for pathway_name, pathway_info in pathways.items():
                        if isinstance(pathway_info, dict) and pathway_info.get('hits', 0) > 0:
                            confidence = self._calculate_pathway_confidence(pathway_info)
                            
                            pathway_summary['detected_pathways'].append({
                                'name': pathway_name,
                                'hits': pathway_info.get('hits', 0),
                                'confidence': confidence,
                                'coverage': pathway_info.get('coverage', 0)
                            })
                            
                            pathway_summary['confidence_scores'].append(confidence)
                    
                    # Sort by confidence and get dominant pathway
                    if pathway_summary['detected_pathways']:
                        pathway_summary['detected_pathways'].sort(
                            key=lambda x: x['confidence'], reverse=True
                        )
                        pathway_summary['dominant_pathway'] = \
                            pathway_summary['detected_pathways'][0]['name']
                        pathway_summary['pathway_count'] = len(pathway_summary['detected_pathways'])
                
                metabolic_profile[category_key] = pathway_summary
        
        # Calculate overall metabolic diversity
        total_pathways = sum(
            profile['pathway_count'] 
            for profile in metabolic_profile.values() 
            if isinstance(profile, dict) and 'pathway_count' in profile
        )
        metabolic_profile['overall_metabolic_diversity'] = total_pathways
        
        # Identify key pathways
        for category_key in ['methane_metabolism', 'sulfur_metabolism', 'nitrogen_metabolism']:
            category_data = metabolic_profile[category_key]
            if category_data.get('detected_pathways'):
                # Add high-confidence pathways to key pathways
                for pathway in category_data['detected_pathways']:
                    if pathway['confidence'] >= self.confidence_thresholds['high']:
                        metabolic_profile['key_pathways'].append({
                            'category': category_key,
                            'pathway': pathway['name'],
                            'confidence': pathway['confidence']
                        })
        
        return metabolic_profile
    
    def _calculate_pathway_confidence(self, pathway_info: Dict) -> float:
        """
        Calculate confidence score for a pathway
        
        Args:
            pathway_info: Pathway information dictionary
            
        Returns:
            Confidence score (0-1)
        """
        hits = pathway_info.get('hits', 0)
        coverage = pathway_info.get('coverage', 0)
        total_sequences = pathway_info.get('total_sequences', 1)
        
        # Simple confidence calculation based on hits and coverage
        hit_ratio = min(hits / max(total_sequences * 0.1, 1), 1.0)  # Normalize by expected hits
        coverage_score = coverage / 100.0 if coverage <= 100 else 1.0
        
        # Weighted average
        confidence = (hit_ratio * 0.6) + (coverage_score * 0.4)
        
        return min(confidence, 1.0)
    
    def _process_temperature_data(self, temperature_data: Dict) -> Dict:
        """
        Process temperature prediction data
        
        Args:
            temperature_data: Raw temperature prediction data
            
        Returns:
            Processed temperature profile
        """
        if not temperature_data:
            return {'status': 'no_data'}
        
        # Extract temperature information
        avg_ogt = temperature_data.get('average_ogt', 0)
        temperature_category = temperature_data.get('temperature_category', 'unknown')
        confidence = temperature_data.get('confidence', 0)
        
        return {
            'optimal_growth_temperature': avg_ogt,
            'temperature_category': temperature_category,
            'confidence': confidence,
            'adaptation_assessment': self._assess_temperature_adaptation(avg_ogt)
        }
    
    def _assess_temperature_adaptation(self, ogt: float) -> Dict:
        """
        Assess temperature adaptation based on OGT
        
        Args:
            ogt: Optimal growth temperature
            
        Returns:
            Temperature adaptation assessment
        """
        if ogt < 15:
            return {
                'category': 'psychrophile',
                'description': 'Cold-adapted organism',
                'cultivation_temperature': '4-15°C',
                'special_requirements': 'Cold cultivation conditions required'
            }
        elif ogt < 45:
            return {
                'category': 'mesophile',
                'description': 'Moderate temperature organism',
                'cultivation_temperature': '20-40°C',
                'special_requirements': 'Standard cultivation conditions'
            }
        elif ogt < 80:
            return {
                'category': 'thermophile',
                'description': 'Heat-loving organism',
                'cultivation_temperature': '45-75°C',
                'special_requirements': 'Elevated temperature cultivation required'
            }
        else:
            return {
                'category': 'hyperthermophile',
                'description': 'Extreme heat-loving organism',
                'cultivation_temperature': '>80°C',
                'special_requirements': 'High-temperature cultivation systems required'
            }
    
    def _generate_overall_assessment(self, integrated_data: Dict) -> Dict:
        """
        Generate overall assessment based on integrated data
        
        Args:
            integrated_data: Integrated analysis data
            
        Returns:
            Overall assessment dictionary
        """
        assessment = {
            'organism_type': 'unknown',
            'metabolic_complexity': 'low',
            'environmental_adaptation': 'standard',
            'cultivation_potential': 'unknown',
            'confidence': 0.0,
            'key_characteristics': []
        }
        
        # Assess organism type based on metabolic profile
        metabolic_profile = integrated_data.get('metabolic_profile', {})
        
        # Check for methanogenic characteristics
        methane_metabolism = metabolic_profile.get('methane_metabolism', {})
        if methane_metabolism.get('pathway_count', 0) > 0:
            assessment['organism_type'] = 'methanogen'
            assessment['key_characteristics'].append('Methane production capability')
        
        # Assess metabolic complexity
        total_pathways = metabolic_profile.get('overall_metabolic_diversity', 0)
        if total_pathways >= 5:
            assessment['metabolic_complexity'] = 'high'
        elif total_pathways >= 2:
            assessment['metabolic_complexity'] = 'medium'
        
        # Assess environmental adaptation
        env_adaptation = integrated_data.get('environmental_adaptation', {})
        
        # Temperature adaptation
        temp_data = env_adaptation.get('temperature', {})
        if temp_data.get('temperature_category') in ['thermophile', 'hyperthermophile']:
            assessment['environmental_adaptation'] = 'thermophilic'
            assessment['key_characteristics'].append('High temperature adaptation')
        elif temp_data.get('temperature_category') == 'psychrophile':
            assessment['environmental_adaptation'] = 'psychrophilic'
            assessment['key_characteristics'].append('Cold adaptation')
        
        # Salt tolerance
        salt_tolerance = env_adaptation.get('salt_tolerance', {})
        if salt_tolerance.get('tolerance_level') == 'high':
            assessment['environmental_adaptation'] = 'halophilic'
            assessment['key_characteristics'].append('High salt tolerance')
        
        # Cultivation potential
        cultivability = integrated_data.get('cultivability_assessment', {})
        genome_quality = cultivability.get('genome_quality', {})
        
        if genome_quality.get('overall_cultivability') == 'high':
            assessment['cultivation_potential'] = 'high'
        elif genome_quality.get('overall_cultivability') == 'medium':
            assessment['cultivation_potential'] = 'medium'
        else:
            assessment['cultivation_potential'] = 'challenging'
        
        # Calculate overall confidence
        confidence_scores = []
        
        if metabolic_profile.get('key_pathways'):
            avg_pathway_confidence = sum(
                p['confidence'] for p in metabolic_profile['key_pathways']
            ) / len(metabolic_profile['key_pathways'])
            confidence_scores.append(avg_pathway_confidence)
        
        if temp_data.get('confidence'):
            confidence_scores.append(temp_data['confidence'])
        
        if genome_quality.get('confidence'):
            confidence_scores.append(genome_quality['confidence'])
        
        if confidence_scores:
            assessment['confidence'] = sum(confidence_scores) / len(confidence_scores)
        
        return assessment
    
    def _generate_recommendations(self, integrated_data: Dict) -> List[str]:
        """
        Generate cultivation and research recommendations
        
        Args:
            integrated_data: Integrated analysis data
            
        Returns:
            List of recommendations
        """
        recommendations = []
        
        # Temperature recommendations
        env_adaptation = integrated_data.get('environmental_adaptation', {})
        temp_data = env_adaptation.get('temperature', {})
        
        if 'adaptation_assessment' in temp_data:
            temp_assessment = temp_data['adaptation_assessment']
            recommendations.append(
                f"Cultivation temperature: {temp_assessment.get('cultivation_temperature', 'Standard')}"
            )
            
            if temp_assessment.get('special_requirements'):
                recommendations.append(temp_assessment['special_requirements'])
        
        # Salt tolerance recommendations
        salt_tolerance = env_adaptation.get('salt_tolerance', {})
        if salt_tolerance.get('tolerance_level') == 'high':
            recommendations.append("Use high-salt media for cultivation")
        elif salt_tolerance.get('tolerance_level') == 'low':
            recommendations.append("Use low-salt or freshwater media")
        
        # Metabolic recommendations
        metabolic_profile = integrated_data.get('metabolic_profile', {})
        
        # Methane metabolism
        methane_metabolism = metabolic_profile.get('methane_metabolism', {})
        if methane_metabolism.get('pathway_count', 0) > 0:
            recommendations.append("Consider anaerobic cultivation conditions for methanogenesis")
            recommendations.append("Monitor methane production during cultivation")
        
        # Sulfur metabolism
        sulfur_metabolism = metabolic_profile.get('sulfur_metabolism', {})
        if sulfur_metabolism.get('pathway_count', 0) > 0:
            recommendations.append("Consider sulfur-containing compounds in media")
        
        # Cultivation potential recommendations
        cultivability = integrated_data.get('cultivability_assessment', {})
        genome_quality = cultivability.get('genome_quality', {})
        
        if genome_quality.get('overall_cultivability') == 'high':
            recommendations.append("High cultivation success probability - proceed with standard protocols")
        elif genome_quality.get('overall_cultivability') == 'medium':
            recommendations.append("Medium cultivation potential - consider media optimization")
        else:
            recommendations.append("Challenging cultivation - may require specialized conditions")
        
        # Amino acid biosynthesis recommendations
        aa_biosynthesis = integrated_data.get('amino_acid_biosynthesis', {})
        if aa_biosynthesis:
            aa_summary = aa_biosynthesis.get('summary', {})
            missing_aa = aa_summary.get('missing_list', [])
            if missing_aa:
                missing_names = [aa.split(' (')[0] for aa in missing_aa]
                recommendations.append(
                    f"Missing or incomplete amino acid biosynthesis: {', '.join(missing_names)} - consider supplementing culture medium"
                )
        
        # Overall assessment recommendations
        overall = integrated_data.get('overall_assessment', {})
        if overall.get('confidence', 0) < 0.5:
            recommendations.append("Low confidence predictions - consider additional sequencing or analysis")
        
        return recommendations
    
    def _save_comprehensive_results(self, results: Dict, output_prefix: str):
        """
        Save comprehensive analysis results
        
        Args:
            results: Comprehensive analysis results
            output_prefix: Output file prefix
        """
        # Save complete JSON results
        json_file = os.path.join(self.results_dir, f"{output_prefix}_comprehensive_analysis.json")
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        
        # Save summary report
        summary_file = os.path.join(self.results_dir, f"{output_prefix}_analysis_summary.txt")
        self._generate_summary_report(results, summary_file)
        
        # Save integrated results as CSV if possible
        if 'integrated_analysis' in results and 'status' not in results['integrated_analysis']:
            csv_file = os.path.join(self.results_dir, f"{output_prefix}_integrated_summary.csv")
            self._save_integrated_csv(results['integrated_analysis'], csv_file)
        
        self.logger.info(f"Comprehensive results saved to {json_file}, {summary_file}")
    
    def _generate_summary_report(self, results: Dict, output_file: str):
        """
        Generate human-readable summary report
        
        Args:
            results: Analysis results
            output_file: Output file path
        """
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("MethArCT Comprehensive Analysis Report\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"Input file: {results.get('input_path', 'Unknown')}\n")
            f.write(f"Analysis date: {results.get('analysis_timestamp', 'Unknown')}\n")
            f.write(f"Analyses performed: {', '.join(results.get('analyses_performed', []))}\n\n")
            
            # Overall assessment
            if 'integrated_analysis' in results and 'overall_assessment' in results['integrated_analysis']:
                overall = results['integrated_analysis']['overall_assessment']
                f.write("Overall Assessment:\n")
                f.write("-" * 20 + "\n")
                f.write(f"Organism type: {overall.get('organism_type', 'Unknown')}\n")
                f.write(f"Metabolic complexity: {overall.get('metabolic_complexity', 'Unknown')}\n")
                f.write(f"Environmental adaptation: {overall.get('environmental_adaptation', 'Unknown')}\n")
                f.write(f"Cultivation potential: {overall.get('cultivation_potential', 'Unknown')}\n")
                f.write(f"Confidence: {overall.get('confidence', 0):.2f}\n\n")
                
                if overall.get('key_characteristics'):
                    f.write("Key characteristics:\n")
                    for char in overall['key_characteristics']:
                        f.write(f"  - {char}\n")
                    f.write("\n")
            
            # Recommendations
            if 'integrated_analysis' in results and 'recommendations' in results['integrated_analysis']:
                recommendations = results['integrated_analysis']['recommendations']
                if recommendations:
                    f.write("Recommendations:\n")
                    f.write("-" * 15 + "\n")
                    for i, rec in enumerate(recommendations, 1):
                        f.write(f"{i}. {rec}\n")
                    f.write("\n")
            
            # Individual analysis summaries
            f.write("Individual Analysis Results:\n")
            f.write("-" * 30 + "\n")
            
            for analysis_type in results.get('analyses_performed', []):
                if analysis_type in results.get('results', {}):
                    analysis_data = results['results'][analysis_type]
                    f.write(f"\n{analysis_type.upper()} Analysis:\n")
                    
                    if 'status' in analysis_data:
                        f.write(f"  Status: {analysis_data['status']}\n")
                        if analysis_data['status'] == 'failed':
                            f.write(f"  Error: {analysis_data.get('error', 'Unknown error')}\n")
                    else:
                        f.write("  Status: Completed successfully\n")
                        
                        # Add specific summaries for each analysis type
                        if analysis_type == 'diamond' and 'summary' in analysis_data:
                            summary = analysis_data['summary']
                            f.write(f"  Total pathways detected: {summary.get('total_pathways', 0)}\n")
                            f.write(f"  Methane pathways: {summary.get('methane_pathways', 0)}\n")
                            f.write(f"  Salt tolerance hits: {summary.get('salt_tolerance_hits', 0)}\n")
                        
                        elif analysis_type == 'tome' and 'summary' in analysis_data:
                            summary = analysis_data['summary']
                            f.write(f"  Average OGT: {summary.get('average_ogt', 0):.1f}°C\n")
                            f.write(f"  Temperature category: {summary.get('temperature_category', 'Unknown')}\n")
                        
                        elif analysis_type == 'checkm2' and 'summary' in analysis_data:
                            summary = analysis_data['summary']
                            f.write(f"  Average completeness: {summary.get('average_completeness', 0):.1f}%\n")
                            f.write(f"  Average contamination: {summary.get('average_contamination', 0):.1f}%\n")
                            f.write(f"  Cultivability: {summary.get('cultivability_assessment', {}).get('overall_cultivability', 'Unknown')}\n")
                        
                        # Add amino acid biosynthesis summary for cultivation analysis
                        if analysis_type == 'cultivation' and 'amino_acid_biosynthesis' in analysis_data:
                            aa_data = analysis_data['amino_acid_biosynthesis']
                            aa_summary = aa_data.get('summary', {})
                            f.write(f"\n  Amino Acid Biosynthesis (20 standard amino acids):\n")
                            f.write(f"  Present: {aa_summary.get('present', 0)}/20\n")
                            f.write(f"  Missing/Partial: {aa_summary.get('missing_or_partial', 0)}/20\n")
                            
                            # List missing amino acids
                            missing_list = aa_summary.get('missing_list', [])
                            if missing_list:
                                missing_names = [aa.split(' (')[0] for aa in missing_list]
                                f.write(f"  Missing pathways: {', '.join(missing_names)}\n")
                            
                            present_list = aa_summary.get('present_list', [])
                            if present_list:
                                present_names = [aa.split(' (')[0] for aa in present_list]
                                f.write(f"  Present pathways: {', '.join(present_names)}\n")
    
    def _save_integrated_csv(self, integrated_data: Dict, output_file: str):
        """
        Save integrated analysis data as CSV
        
        Args:
            integrated_data: Integrated analysis data
            output_file: Output CSV file path
        """
        try:
            # Create a flattened summary for CSV export
            csv_data = []
            
            # Overall assessment
            overall = integrated_data.get('overall_assessment', {})
            row = {
                'analysis_type': 'overall_assessment',
                'organism_type': overall.get('organism_type', ''),
                'metabolic_complexity': overall.get('metabolic_complexity', ''),
                'environmental_adaptation': overall.get('environmental_adaptation', ''),
                'cultivation_potential': overall.get('cultivation_potential', ''),
                'confidence': overall.get('confidence', 0)
            }
            csv_data.append(row)
            
            # Metabolic pathways
            metabolic_profile = integrated_data.get('metabolic_profile', {})
            for pathway_type in ['methane_metabolism', 'sulfur_metabolism', 'nitrogen_metabolism']:
                if pathway_type in metabolic_profile:
                    pathway_data = metabolic_profile[pathway_type]
                    row = {
                        'analysis_type': pathway_type,
                        'pathway_count': pathway_data.get('pathway_count', 0),
                        'dominant_pathway': pathway_data.get('dominant_pathway', ''),
                        'confidence': sum(pathway_data.get('confidence_scores', [])) / len(pathway_data.get('confidence_scores', [1])) if pathway_data.get('confidence_scores') else 0
                    }
                    csv_data.append(row)
            
            # Save to CSV
            if csv_data:
                df = pd.DataFrame(csv_data)
                df.to_csv(output_file, index=False)
                
        except Exception as e:
            self.logger.warning(f"Could not save integrated CSV: {str(e)}")
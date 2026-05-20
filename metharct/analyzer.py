# -*- coding: utf-8 -*-
"""
MethArCT主分析器

整合Diamond、Tome、CheckM2等分析工具的主要分析类
"""

import os
import json
from pathlib import Path
from typing import Dict, List, Optional, Union

from .core.diamond_analyzer import DiamondAnalyzer
from .core.tome_analyzer import TomeAnalyzer
from .core.checkm2_analyzer import CheckM2Analyzer
from .core.pathway_predictor import PathwayPredictor
from .utils.config import Config
from .utils.logger import get_logger
from .utils.file_utils import FileUtils


class MethArCTAnalyzer:
    """
    MethArCT主分析器类
    
    整合所有分析功能的主要接口类
    """
    
    def __init__(self, config: Optional[Config] = None):
        """
        初始化MethArCT分析器
        
        Args:
            config: 配置对象，如果为None则使用默认配置
        """
        self.config = config or Config()
        self.logger = get_logger("metharct_analyzer")
        
        # 初始化各个分析器
        self.diamond_analyzer = DiamondAnalyzer(self.config)
        self.tome_analyzer = TomeAnalyzer(self.config)
        self.checkm2_analyzer = CheckM2Analyzer(self.config)
        self.pathway_predictor = PathwayPredictor(self.config)
        
        # 结果存储目录
        self.results_dir = self.config.get('output.base_dir', 'results')
        FileUtils.ensure_dir(self.results_dir)
        
        self.logger.info("MethArCT分析器初始化完成")
    
    def analyze_comprehensive(self, 
                            input_file: str,
                            output_dir: Optional[str] = None,
                            run_diamond: bool = True,
                            run_tome: bool = True,
                            run_checkm2: bool = True,
                            **kwargs) -> Dict:
        """
        执行综合分析
        
        Args:
            input_file: 输入文件路径（FASTA格式）
            output_dir: 输出目录，如果为None则使用默认目录
            run_diamond: 是否运行Diamond分析
            run_tome: 是否运行Tome分析
            run_checkm2: 是否运行CheckM2分析
            **kwargs: 其他参数
            
        Returns:
            包含所有分析结果的字典
        """
        if output_dir is None:
            output_dir = self.results_dir
            
        FileUtils.ensure_dir(output_dir)
        
        results = {
            'input_file': input_file,
            'output_dir': output_dir,
            'analysis_results': {}
        }
        
        self.logger.info(f"开始综合分析: {input_file}")
        
        try:
            # Diamond代谢通路分析
            if run_diamond:
                self.logger.info("运行Diamond代谢通路分析...")
                diamond_results = self.diamond_analyzer.analyze_sequence(
                    input_file, 
                    os.path.join(output_dir, 'diamond')
                )
                results['analysis_results']['diamond'] = diamond_results
                
            # Tome温度预测分析
            if run_tome:
                self.logger.info("运行Tome温度预测分析...")
                tome_results = self.tome_analyzer.analyze(
                    input_file,
                    os.path.join(output_dir, 'tome')
                )
                results['analysis_results']['tome'] = tome_results
                
            # CheckM2基因组质量评估
            if run_checkm2:
                self.logger.info("运行CheckM2基因组质量评估...")
                checkm2_results = self.checkm2_analyzer.analyze(
                    input_file,
                    os.path.join(output_dir, 'checkm2')
                )
                results['analysis_results']['checkm2'] = checkm2_results
                
            # 代谢通路预测（基于Diamond结果）
            if run_diamond:
                self.logger.info("运行代谢通路预测...")
                pathway_results = self.pathway_predictor.predict_comprehensive(
                    input_path=input_file,
                    output_prefix=os.path.basename(input_file).replace('.faa', ''),
                    include_tome=False,
                    include_checkm2=False
                )
                results['analysis_results']['pathways'] = pathway_results
                
            # 保存综合结果
            results_file = os.path.join(output_dir, 'comprehensive_results.json')
            with open(results_file, 'w', encoding='utf-8') as f:
                json.dump(results, f, indent=2, ensure_ascii=False)
                
            self.logger.info(f"综合分析完成，结果保存至: {results_file}")
            
        except Exception as e:
            self.logger.error(f"综合分析失败: {str(e)}")
            results['error'] = str(e)
            raise
            
        return results
    
    def analyze_diamond(self, input_file: str, output_dir: Optional[str] = None) -> Dict:
        """
        仅执行Diamond分析
        
        Args:
            input_file: 输入文件路径
            output_dir: 输出目录
            
        Returns:
            Diamond分析结果
        """
        if output_dir is None:
            output_dir = os.path.join(self.results_dir, 'diamond')
            
        return self.diamond_analyzer.analyze_sequence(input_file, output_dir)
    
    def analyze_tome(self, input_file: str, output_dir: Optional[str] = None) -> Dict:
        """
        仅执行Tome分析
        
        Args:
            input_file: 输入文件路径
            output_dir: 输出目录
            
        Returns:
            Tome分析结果
        """
        if output_dir is None:
            output_dir = os.path.join(self.results_dir, 'tome')
            
        return self.tome_analyzer.analyze(input_file, output_dir)
    
    def analyze_checkm2(self, input_file: str, output_dir: Optional[str] = None) -> Dict:
        """
        仅执行CheckM2分析
        
        Args:
            input_file: 输入文件路径
            output_dir: 输出目录
            
        Returns:
            CheckM2分析结果
        """
        if output_dir is None:
            output_dir = os.path.join(self.results_dir, 'checkm2')
            
        return self.checkm2_analyzer.analyze(input_file, output_dir)
    
    def get_analysis_summary(self, results: Dict) -> Dict:
        """
        生成分析结果摘要
        
        Args:
            results: 分析结果字典
            
        Returns:
            结果摘要字典
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
        验证输入文件
        
        Args:
            input_file: 输入文件路径
            
        Returns:
            验证是否通过
        """
        if not os.path.exists(input_file):
            self.logger.error(f"输入文件不存在: {input_file}")
            return False
            
        # 检查文件格式（简单检查）
        try:
            with open(input_file, 'r') as f:
                first_line = f.readline().strip()
                if not first_line.startswith('>'):
                    self.logger.error(f"输入文件不是有效的FASTA格式: {input_file}")
                    return False
        except Exception as e:
            self.logger.error(f"读取输入文件失败: {str(e)}")
            return False
            
        return True
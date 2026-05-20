#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MethArCT 主分析脚本
Methanogenic Archaea Metabolic Pathway Analysis Tool

使用方法:
    python metharct_analysis.py <输入文件> [选项]
    
示例:
    python metharct_analysis.py protein.fasta
    python metharct_analysis.py protein.fasta -g genome.fasta -o results
    python metharct_analysis.py protein.fasta --diamond-only
"""

import os
import sys
import argparse
import logging
from pathlib import Path
from datetime import datetime

# 添加项目路径到Python路径
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

try:
    from metharct.core.diamond_analyzer import DiamondAnalyzer
    from metharct.core.tome_analyzer import TomeAnalyzer
    from metharct.core.checkm2_analyzer import CheckM2Analyzer
    from metharct.core.pathway_predictor import PathwayPredictor
    from metharct.utils.config import Config
    from metharct.utils.logger import setup_logger
except ImportError as e:
    print(f"❌ 导入MethArCT模块失败: {e}")
    print("请确保MethArCT已正确安装: pip install -e .")
    sys.exit(1)

def setup_logging(log_level="DEBUG"):
    """设置日志系统"""
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(f'metharct_analysis_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
        ]
    )
    # 设置所有相关模块的日志级别为DEBUG
    logging.getLogger("metharct").setLevel(logging.DEBUG)
    logging.getLogger("diamond_analyzer").setLevel(logging.DEBUG)
    return logging.getLogger("metharct_analysis")

def check_input_file(file_path):
    """检查输入文件是否存在且格式正确"""
    if not os.path.exists(file_path):
        return False, f"文件不存在: {file_path}"
    
    # 检查文件大小
    file_size = os.path.getsize(file_path)
    if file_size == 0:
        return False, f"文件为空: {file_path}"
    
    # 检查文件格式（简单检查）
    try:
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
            if not first_line.startswith('>'):
                return False, f"文件格式不正确，应该是FASTA格式: {file_path}"
    except Exception as e:
        return False, f"无法读取文件: {e}"
    
    return True, "文件检查通过"

def run_diamond_analysis(config, input_file, output_dir, logger):
    """运行Diamond代谢通路分析"""
    try:
        logger.info("开始Diamond代谢通路分析")
        diamond_analyzer = DiamondAnalyzer(config)
        
        results = diamond_analyzer.analyze_sequence(
            input_file=input_file,
            output_prefix=os.path.join(output_dir, "diamond")
        )
        
        logger.info(f"Diamond分析完成")
        return results
        
    except Exception as e:
        logger.error(f"Diamond分析失败: {str(e)}")
        return None

def run_tome_analysis(config, input_file, output_dir, logger):
    """运行Tome温度预测分析"""
    try:
        logger.info("开始Tome温度预测分析")
        tome_analyzer = TomeAnalyzer(config)
        
        # TomeAnalyzer在初始化时已经检查了可用性
        # 如果工具不可用，初始化时会抛出异常
        results = tome_analyzer.predict_ogt(
            protein_file=input_file,
            output_dir=output_dir
        )
        
        logger.info(f"Tome分析完成")
        return results
        
    except Exception as e:
        logger.warning(f"Tome分析失败: {str(e)}")
        return None

def run_checkm2_analysis(config, input_file, output_dir, logger):
    """运行CheckM2基因组质量评估"""
    try:
        logger.info("开始CheckM2基因组质量评估")
        checkm2_analyzer = CheckM2Analyzer(config)
        
        # 使用output_dir作为output_prefix
        output_prefix = output_dir
        
        results = checkm2_analyzer.analyze_genome_quality(
            input_path=input_file,
            output_prefix=output_prefix,
            input_type='fasta'
        )
        
        logger.info(f"CheckM2分析完成")
        return results
        
    except Exception as e:
        logger.warning(f"CheckM2分析失败: {str(e)}")
        return None

def run_comprehensive_analysis(config, protein_file, genome_file, output_dir, logger):
    """运行综合分析"""
    try:
        logger.info("开始综合分析")
        pathway_predictor = PathwayPredictor(config)
        
        results = pathway_predictor.comprehensive_analysis(
            protein_file=protein_file,
            genome_file=genome_file,
            output_dir=os.path.join(output_dir, "comprehensive"),
            threads=4
        )
        
        logger.info(f"综合分析完成: {results.get('output_dir', 'Unknown')}")
        return results
        
    except Exception as e:
        logger.warning(f"综合分析失败: {str(e)}")
        return None

def print_results_summary(diamond_results, tome_results, checkm2_results, comprehensive_results, logger):
    """打印分析结果摘要"""
    print("\n" + "="*60)
    print("MethArCT 分析结果摘要")
    print("="*60)
    
    if diamond_results:
        print(f"✓ Diamond代谢通路分析: 完成")
        if 'pathways' in diamond_results:
            print(f"  检测到的代谢通路: {diamond_results['pathways']}")
    else:
        print("✗ Diamond代谢通路分析: 失败")
    
    if tome_results:
        print(f"✓ Tome温度预测分析: 完成")
        if 'predicted_ogt' in tome_results:
            print(f"  预测的最适生长温度: {tome_results['predicted_ogt']}°C")
        if 'temperature_class' in tome_results:
            print(f"  温度分类: {tome_results['temperature_class']}")
    else:
        print("✗ Tome温度预测分析: 跳过或失败")
    
    if checkm2_results:
        print(f"✓ CheckM2基因组质量评估: 完成")
        if 'completeness' in checkm2_results:
            print(f"  基因组完整性: {checkm2_results['completeness']:.2f}%")
        if 'contamination' in checkm2_results:
            print(f"  污染度: {checkm2_results['contamination']:.2f}%")
        if 'quality_grade' in checkm2_results:
            print(f"  质量等级: {checkm2_results['quality_grade']}")
    else:
        print("✗ CheckM2基因组质量评估: 跳过或失败")
    
    if comprehensive_results:
        print(f"✓ 综合分析: 完成")
        if 'pathways' in comprehensive_results:
            print(f"  代谢通路: {comprehensive_results['pathways']}")
        if 'salt_tolerance' in comprehensive_results:
            print(f"  耐盐性: {comprehensive_results['salt_tolerance']}")
        if 'optimal_temperature' in comprehensive_results:
            print(f"  最适生长温度: {comprehensive_results['optimal_temperature']}°C")
        if 'culturability' in comprehensive_results:
            print(f"  可培养性: {comprehensive_results['culturability']}")
    else:
        print("✗ 综合分析: 跳过或失败")
    
    print("="*60)

def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='MethArCT - Methanogenic Archaea Metabolic Pathway Analysis Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  %(prog)s protein.fasta                    # 基本分析
  %(prog)s protein.fasta -g genome.fasta    # 完整分析
  %(prog)s protein.fasta --diamond-only     # 仅核心分析
  %(prog)s protein.fasta --log-level DEBUG # 调试模式
        """
    )
    
    parser.add_argument('input_file', help='输入蛋白质序列文件 (FASTA格式)')
    parser.add_argument('-g', '--genome', help='基因组序列文件 (可选，用于CheckM2分析)')
    parser.add_argument('-o', '--output', default='metharct_results', help='输出目录 (默认: metharct_results)')
    parser.add_argument('--diamond-only', action='store_true', help='仅运行Diamond分析')
    parser.add_argument('--log-level', default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], help='日志级别')
    parser.add_argument('--threads', type=int, default=4, help='线程数 (默认: 4)')
    
    args = parser.parse_args()
    
    # 设置日志
    logger = setup_logging(args.log_level)
    logger.info("MethArCT分析开始")
    
    # 检查输入文件
    is_valid, message = check_input_file(args.input_file)
    if not is_valid:
        logger.error(message)
        print(f"错误: {message}")
        return 1
    
    # 创建输出目录
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"输出目录: {output_dir}")
    
    try:
        # 加载配置
        config = Config()
        logger.info("配置加载成功")
        
        # 运行Diamond分析（必需）
        diamond_results = run_diamond_analysis(config, args.input_file, output_dir, logger)
        
        if not diamond_results:
            logger.error("Diamond分析失败，无法继续")
            return 1
        
        tome_results = None
        checkm2_results = None
        comprehensive_results = None
        
        # 如果不是仅Diamond分析，运行其他分析
        if not args.diamond_only:
            # Tome温度预测分析
            tome_results = run_tome_analysis(config, args.input_file, output_dir, logger)
            
            # CheckM2基因组质量评估（如果有基因组文件）
            if args.genome and os.path.exists(args.genome):
                checkm2_results = run_checkm2_analysis(config, args.genome, output_dir, logger)
            else:
                logger.info("未提供基因组文件，跳过CheckM2分析")
            
            # 综合分析（如果有基因组文件）
            if args.genome and os.path.exists(args.genome):
                comprehensive_results = run_comprehensive_analysis(
                    config, args.input_file, args.genome, output_dir, logger
                )
        
        # 打印结果摘要
        print_results_summary(diamond_results, tome_results, checkm2_results, comprehensive_results, logger)
        
        logger.info("所有分析完成")
        print(f"\n所有结果已保存在 '{output_dir}' 目录中")
        
        return 0
        
    except Exception as e:
        logger.error(f"分析过程中出现错误: {str(e)}")
        print(f"错误: {str(e)}")
        return 1

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Diamond sequence alignment analyzer for MethArCT

Performs protein sequence alignment using Diamond BLAST for metabolic pathway prediction,
salt tolerance assessment, and cultivability evaluation.
"""

import os
import subprocess
import tempfile
import shutil
import json
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from ..utils.config import Config
from ..utils.logger import get_logger
from ..utils.file_utils import FileUtils
from ..utils.sequence_utils import SequenceUtils

class DiamondAnalyzer:
    """Diamond sequence alignment analyzer"""
    
    def __init__(self, config: Optional[Config] = None):
        self.config = config or Config()
        self.logger = get_logger("diamond_analyzer")
        
        # Tool configuration - 支持WSL
        self.use_wsl = self.config.get('tools.diamond.use_wsl', False)
        
        if self.use_wsl:
            self.diamond_path = self.config.get('tools.diamond.wsl_path', 'wsl diamond')
        else:
            self.diamond_path = self.config.get('tools.diamond.path', 'diamond')
            
        self.threads = self.config.get('tools.diamond.threads', 4)
        self.evalue = self.config.get('tools.diamond.evalue', 1e-5)
        self.max_target_seqs = self.config.get('tools.diamond.max_target_seqs', 1)
        self.identity_threshold = self.config.get('tools.diamond.identity_threshold', 30.0)
        
        # 双阈值分析参数
        self.HIGH_E_VALUE_THRESHOLD = 1e-100
        self.HIGH_BITSCORE_THRESHOLD = 400
        self.LOW_E_VALUE_THRESHOLD = 1e-5
        self.LOW_BITSCORE_THRESHOLD = 40
        self.HIGH_QUALITY_THRESHOLD = 60
        
        # 栽培评估阈值
        self.CULTIVABILITY_E_VALUE_THRESHOLD = 1e-3
        self.CULTIVABILITY_BITSCORE_THRESHOLD = 50
        
        # 耐盐性评估阈值
        self.SALT_TOLERANCE_E_VALUE_THRESHOLD = 1e-5  # 使用diamond默认阈值
        self.SALT_TOLERANCE_BITSCORE_THRESHOLD = 100    # 使用diamond默认值（无比特分数过滤）
        
        # Database paths
        self.db_base_dir = self.config.get('databases.base_dir', 'data/databases')
        self.db_dir = os.path.join(self.db_base_dir, 'diamond')
        
        # Results directory
        self.results_dir = self.config.get('output.base_dir', 'results')
        FileUtils.ensure_dir(self.results_dir)
        FileUtils.ensure_dir(self.db_dir)
        
        # Database information
        self.pathway_names = self.config.get('pathway_names', {})
        self.reference_counts = self.config.get('reference_sequence_counts', {})
        
        # Check tool availability
        self._check_diamond_availability()
        
        # Prepare databases
        self._prepare_databases()
    
    def _check_diamond_availability(self) -> bool:
        """
        Check if Diamond tool is available - 支持WSL环境
        
        Returns:
            True if Diamond is available, False otherwise
        """
        try:
            if self.use_wsl:
                # WSL环境下的检查
                result = subprocess.run(
                    ['wsl', 'diamond', 'version'],
                    capture_output=True,
                    text=True,
                    timeout=30
                )
            else:
                # 本地环境下的检查
                result = subprocess.run(
                    [self.diamond_path, 'version'],
                    capture_output=True,
                    text=True,
                    timeout=30
                )
            
            if result.returncode == 0:
                version = result.stdout.strip()
                env_type = "WSL" if self.use_wsl else "本地"
                self.logger.info(f"Diamond工具可用 ({env_type}): {version}")
                return True
            else:
                self.logger.error(f"Diamond工具不可用: {result.stderr}")
                return False
                
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            self.logger.error(f"Error checking Diamond availability: {str(e)}")
            return False
    
    def _prepare_databases(self):
        """
        Prepare Diamond databases from reference sequences
        """
        self.logger.info("Preparing Diamond databases...")
        
        # Get all database paths from config
        db_paths = self.config.get_all_database_paths()
        
        for db_name, fasta_path in db_paths.items():
            db_file = os.path.join(self.db_dir, f"{db_name}.dmnd")
            
            # Check if database exists and is newer than source file
            if (os.path.exists(db_file) and 
                os.path.exists(fasta_path) and
                os.path.getmtime(db_file) > os.path.getmtime(fasta_path)):
                self.logger.debug(f"Database {db_name} is up to date")
                continue
            
            # Create or update database
            if os.path.exists(fasta_path):
                self._create_diamond_database(fasta_path, db_name)
            else:
                self.logger.warning(f"Reference file not found: {fasta_path}")
    
    def _create_diamond_database(self, fasta_path: str, db_name: str) -> bool:
        """
        Create Diamond database from FASTA file
        
        Args:
            fasta_path: Path to input FASTA file
            db_name: Database name
            
        Returns:
            True if successful, False otherwise
        """
        try:
            db_path = os.path.join(self.db_dir, db_name)
            
            cmd = [
                self.diamond_path, 'makedb',
                '--in', fasta_path,
                '--db', db_path,
                '--threads', str(self.threads)
            ]
            
            self.logger.info(f"Creating Diamond database: {db_name}")
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300  # 5 minutes timeout
            )
            
            if result.returncode == 0:
                self.logger.info(f"Successfully created database: {db_name}")
                return True
            else:
                self.logger.error(f"Failed to create database {db_name}: {result.stderr}")
                return False
                
        except subprocess.TimeoutExpired:
            self.logger.error(f"Timeout creating database: {db_name}")
            return False
        except Exception as e:
            self.logger.error(f"Error creating database {db_name}: {str(e)}")
            return False
    
    def analyze_sequence(self, 
                        input_file: Union[str, Path],
                        output_prefix: Optional[str] = None,
                        databases: Optional[List[str]] = None) -> Dict[str, any]:
        """
        Analyze protein sequences against reference databases
        
        Args:
            input_file: Path to input FASTA file
            output_prefix: Prefix for output files
            databases: List of databases to search against (None for all)
            
        Returns:
            Analysis results dictionary
        """
        input_file = Path(input_file)
        
        if not input_file.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")
        
        if not FileUtils.validate_fasta(input_file):
            raise ValueError(f"Invalid FASTA file: {input_file}")
        
        # Set output prefix
        if output_prefix is None:
            output_prefix = input_file.stem
        
        # Get databases to search
        if databases is None:
            databases = list(self.config.get_all_database_paths().keys())
        
        self.logger.info(f"Starting Diamond analysis for {input_file.name}")
        self.logger.info(f"Searching against {len(databases)} databases")
        
        # Run Diamond searches
        search_results = {}
        for db_name in databases:
            result = self._run_diamond_search(input_file, db_name, output_prefix)
            if result:
                search_results[db_name] = result
        
        # Process and summarize results
        analysis_results = self._process_results(search_results, input_file)
        
        # Save results
        self._save_results(analysis_results, output_prefix)
        
        return analysis_results
    
    def _run_diamond_search(self, 
                           input_file: Path, 
                           db_name: str, 
                           output_prefix: str) -> Optional[Dict]:
        """
        Run Diamond search against a specific database with dual threshold analysis
        
        Args:
            input_file: Input FASTA file
            db_name: Database name
            output_prefix: Output file prefix
            
        Returns:
            Search results dictionary or None if failed
        """
        db_path = os.path.join(self.db_dir, f"{db_name}.dmnd")
        
        if not os.path.exists(db_path):
            self.logger.warning(f"Database not found: {db_path}")
            return None
        
        # Create output file in the same directory as output_prefix
        output_dir = os.path.dirname(output_prefix) if os.path.dirname(output_prefix) else self.results_dir
        output_file = os.path.join(
            output_dir, 
            f"{os.path.basename(output_prefix)}_{db_name}_diamond.tsv"
        )
        
        # 根据数据库类型选择E值阈值
        if db_name == 'CULTIVATION':
            search_evalue = self.CULTIVABILITY_E_VALUE_THRESHOLD
        elif db_name == 'NAIYAN':
            search_evalue = self.SALT_TOLERANCE_E_VALUE_THRESHOLD
        else:
            search_evalue = self.LOW_E_VALUE_THRESHOLD
        
        cmd = [
            self.diamond_path, 'blastp',
            '--query', str(input_file),
            '--db', db_path,
            '--out', output_file,
            '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 
                       'mismatch', 'gapopen', 'qstart', 'qend', 
                       'sstart', 'send', 'evalue', 'bitscore',
            '--evalue', str(search_evalue),
            '--max-target-seqs', str(self.max_target_seqs),
            '--threads', str(self.threads),
            '--sensitive'
        ]
        
        try:
            self.logger.debug(f"Running Diamond search: {db_name}")
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600  # 10 minutes timeout
            )
            
            if result.returncode == 0:
                # Parse results
                if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                    df = pd.read_csv(
                        output_file,
                        sep='\t',
                        header=None,
                        names=['qseqid', 'sseqid', 'pident', 'length',
                               'mismatch', 'gapopen', 'qstart', 'qend',
                               'sstart', 'send', 'evalue', 'bitscore']
                    )
                    
                    # 根据数据库类型返回不同的结果结构
                    if db_name == 'CULTIVATION':
                        # 栽培评估数据库只进行栽培评估过滤
                        cultivability_hits = df[
                            (df['evalue'] <= self.CULTIVABILITY_E_VALUE_THRESHOLD) & 
                            (df['bitscore'] >= self.CULTIVABILITY_BITSCORE_THRESHOLD)
                        ]
                        
                        return {
                            'database': db_name,
                            'total_hits': len(df),
                            'cultivability_hits': cultivability_hits,
                            'output_file': output_file
                        }
                    
                    elif db_name == 'NAIYAN':
                        # 耐盐性数据库只进行耐盐性评估过滤
                        salt_tolerance_hits = df[
                            (df['evalue'] <= self.SALT_TOLERANCE_E_VALUE_THRESHOLD) & 
                            (df['bitscore'] >= self.SALT_TOLERANCE_BITSCORE_THRESHOLD)
                        ]
                        
                        return {
                            'database': db_name,
                            'total_hits': len(df),
                            'salt_tolerance_hits': salt_tolerance_hits,
                            'output_file': output_file
                        }
                    
                    else:
                        # 代谢途径数据库进行双阈值过滤
                        # 检查是否为甲烷代谢途径，如果是则使用单阈值分析
                        if db_name in ['CO2-CH4', 'JIAAN-CH4', 'JIACHUN-CH4', 'JIALIUCHUN-CH4', 
                                     'YISUAN-CH4', 'C16-CH4', 'CO-CH4', 'JIASUAN-CH4', 
                                     'JIAYANGJI-CH4', 'ZHIFANGSUAN-CH4', '2JIAAN-CH4', '3JIAAN-CH4',
                                     '甘氨酸甜菜碱产甲烷', '硫代丙酸甲酯产甲烷', '四甲基铵产甲烷',
                                     '甲醇歧化产甲烷']:
                            # 甲烷代谢途径：单阈值分析（只保留低阈值）
                            low_threshold_hits = df[
                                (df['evalue'] <= self.LOW_E_VALUE_THRESHOLD) & 
                                (df['bitscore'] >= self.LOW_BITSCORE_THRESHOLD)
                            ]
                            
                            return {
                                'database': db_name,
                                'total_hits': len(df),
                                'low_threshold_hits': low_threshold_hits,
                                'output_file': output_file
                            }
                        else:
                            # 其他代谢途径（硫、氮）保持双阈值分析
                            high_threshold_hits = df[
                                (df['evalue'] <= self.HIGH_E_VALUE_THRESHOLD) & 
                                (df['bitscore'] >= self.HIGH_BITSCORE_THRESHOLD)
                            ]
                            
                            low_threshold_hits = df[
                                (df['evalue'] <= self.LOW_E_VALUE_THRESHOLD) & 
                                (df['bitscore'] >= self.LOW_BITSCORE_THRESHOLD)
                            ]
                            
                            return {
                                'database': db_name,
                                'total_hits': len(df),
                                'high_threshold_hits': high_threshold_hits,
                                'low_threshold_hits': low_threshold_hits,
                                'output_file': output_file
                            }
                        
                else:
                    self.logger.info(f"No hits found for database: {db_name}")
                    # 根据数据库类型返回相应的空DataFrame结构
                    if db_name == 'CULTIVATION':
                        return {
                            'database': db_name,
                            'total_hits': 0,
                            'cultivability_hits': pd.DataFrame(),
                            'output_file': output_file
                        }
                    elif db_name == 'NAIYAN':
                        return {
                            'database': db_name,
                            'total_hits': 0,
                            'salt_tolerance_hits': pd.DataFrame(),
                            'output_file': output_file
                        }
                    else:
                        return {
                            'database': db_name,
                            'total_hits': 0,
                            'high_threshold_hits': pd.DataFrame(),
                            'low_threshold_hits': pd.DataFrame(),
                            'output_file': output_file
                        }
            else:
                self.logger.error(f"Diamond search failed for {db_name}: {result.stderr}")
                return None
                
        except subprocess.TimeoutExpired:
            self.logger.error(f"Diamond search timeout for database: {db_name}")
            return None
        except Exception as e:
            self.logger.error(f"Error in Diamond search for {db_name}: {str(e)}")
            return None
    
    def _process_results(self, search_results: Dict, input_file: Path) -> Dict:
        """
        Process and summarize Diamond search results with dual threshold analysis
        
        Args:
            search_results: Raw search results
            input_file: Input FASTA file
            
        Returns:
            Processed analysis results
        """
        # Count input sequences
        total_sequences = FileUtils.count_sequences(input_file)
        
        # Initialize results
        pathway_results = {}
        summary_stats = {
            'total_input_sequences': total_sequences,
            'databases_searched': len(search_results),
            'pathways_detected': 0,
            'total_hits': 0
        }
        
        # Process each database result
        for db_name, result in search_results.items():
            if result and result['total_hits'] > 0:
                reference_count = self.reference_counts.get(db_name, 1)
                
                if db_name == 'CULTIVATION':
                    # 栽培评估数据库处理
                    # 计算有多少个不同的查询序列匹配到了栽培相关基因（去重）
                    cultivability_matched = len(result['cultivability_hits']['qseqid'].drop_duplicates())
                    # 修改计算逻辑：匹配项/自身的基因数乘以100%
                    cultivability_percentage = (cultivability_matched / total_sequences) * 100
                    
                    if cultivability_percentage <= 11:
                        cultivability_status = "Cultivable"
                    elif cultivability_percentage <= 13:
                        cultivability_status = "Co-cultivation dependent"
                    else:
                        cultivability_status = "Obligate symbiont"
                    
                    pathway_results[db_name] = {
                        'pathway_name': self.pathway_names.get(db_name, db_name),
                        'reference_sequences': reference_count,
                        'cultivability_hits': cultivability_matched,
                        'cultivability_percentage': round(cultivability_percentage, 2),
                        'cultivability_status': cultivability_status,
                        'average_identity': round(
                            result['cultivability_hits']['pident'].mean(), 2
                        ) if len(result['cultivability_hits']) > 0 else 0,
                        'best_hits': result['cultivability_hits'].nlargest(5, 'bitscore').to_dict('records')
                    }
                    
                elif db_name == 'NAIYAN':
                    # 耐盐性数据库处理
                    # 计算有多少个不同的查询序列匹配到了耐盐性相关基因（去重）
                    salt_tolerance_matched = len(result['salt_tolerance_hits']['qseqid'].drop_duplicates())
                    
                    if 0 <= salt_tolerance_matched <= 2:
                        salt_tolerance_level = "No Resistance"
                    elif 3 <= salt_tolerance_matched <= 5:
                        salt_tolerance_level = "Moderate"
                    elif 6 <= salt_tolerance_matched <= 8:
                        salt_tolerance_level = "Strong"
                    elif salt_tolerance_matched >= 9:
                        salt_tolerance_level = "Very Strong"
                    else:
                        salt_tolerance_level = "Unknown"
                    
                    pathway_results[db_name] = {
                        'pathway_name': self.pathway_names.get(db_name, db_name),
                        'reference_sequences': reference_count,
                        'salt_tolerance_hits': salt_tolerance_matched,
                        'salt_tolerance_level': salt_tolerance_level,
                        'average_identity': round(
                            result['salt_tolerance_hits']['pident'].mean(), 2
                        ) if len(result['salt_tolerance_hits']) > 0 else 0,
                        'best_hits': result['salt_tolerance_hits'].nlargest(5, 'bitscore').to_dict('records')
                    }
                    
                else:
                    # 代谢途径数据库处理（双阈值分析）
                    # 检查是否为甲烷代谢途径，如果是则使用单阈值分析
                    if db_name in ['CO2-CH4', 'JIAAN-CH4', 'JIACHUN-CH4', 'JIALIUCHUN-CH4', 
                                 'YISUAN-CH4', 'C16-CH4', 'CO-CH4', 'JIASUAN-CH4', 
                                 'JIAYANGJI-CH4', 'ZHIFANGSUAN-CH4', '2JIAAN-CH4', '3JIAAN-CH4',
                                 '甘氨酸甜菜碱产甲烷', '硫代丙酸甲酯产甲烷', '四甲基铵产甲烷',
                                 '甲醇歧化产甲烷']:
                        # 甲烷代谢途径：单阈值分析（只保留低阈值）
                        # 注意：同源基因（如FwdF/FwdG/FwdH、Hmd/Mtd）虽然序列相似，但代表不同基因功能，应该分别统计
                        low_threshold_matched = len(result['low_threshold_hits'].drop_duplicates(subset=['sseqid']))
                        
                        # 对于甲酸产甲烷途径，添加特殊化处理
                        if db_name == 'JIASUAN-CH4':
                            # 调试信息：确认代码被执行
                            self.logger.debug("=== JIASUAN-CH4 通路评判逻辑开始执行 ===")
                            
                            # 获取所有比对到的唯一基因名称（第二列）
                            unique_hits = result['low_threshold_hits'].iloc[:, 1].unique()
                            
                            # 检查关键基因的比对情况
                            fwd_fgh_present = any(gene in unique_hits for gene in ['FwdF', 'FwdG', 'FwdH'])
                            hmd_mtd_present = any(gene in unique_hits for gene in ['Hmd', 'Mtd'])
                            
                            # 检查FdhA和FdhB基因的高严格比对条件
                            fdhA_hits = result['low_threshold_hits'][
                                (result['low_threshold_hits'].iloc[:, 1] == 'FdhA') &
                                (result['low_threshold_hits']['bitscore'] > 380) &
                                (result['low_threshold_hits']['evalue'] <= 1e-100)
                            ]
                            fdhB_hits = result['low_threshold_hits'][
                                (result['low_threshold_hits'].iloc[:, 1] == 'FdhB') &
                                (result['low_threshold_hits']['bitscore'] > 380) &
                                (result['low_threshold_hits']['evalue'] <= 1e-100)
                            ]
                            fdhA_present = len(fdhA_hits) > 0
                            fdhB_present = len(fdhB_hits) > 0
                            
                            # 调试信息：输出关键基因检测状态
                            self.logger.debug(f"JIASUAN-CH4 关键基因检测状态:")
                            self.logger.debug(f"  FwdF/FwdG/FwdH: {fwd_fgh_present}")
                            self.logger.debug(f"  Hmd/Mtd: {hmd_mtd_present}")
                            self.logger.debug(f"  FdhA: {fdhA_present} (命中数: {len(fdhA_hits)})")
                            self.logger.debug(f"  FdhB: {fdhB_present} (命中数: {len(fdhB_hits)})")
                            self.logger.debug(f"  基础基因覆盖: {len(unique_hits)}/{reference_count}")
                            
                            # 优化评判逻辑：采用加权评分法
                            # 基础基因覆盖（权重40%）：检测到的基因数量
                            base_score = min(len(unique_hits) / reference_count * 40, 40)
                            
                            # 关键功能基因（权重60%）：
                            # - FwdF/FwdG/FwdH（权重15%）
                            # - Hmd/Mtd（权重15%）
                            # - FdhA（权重15%）
                            # - FdhB（权重15%）
                            fwd_fgh_score = 15 if fwd_fgh_present else 0
                            hmd_mtd_score = 15 if hmd_mtd_present else 0
                            fdhA_score = 15 if fdhA_present else 0
                            fdhB_score = 15 if fdhB_present else 0
                            
                            # 计算总得分
                            total_score = base_score + fwd_fgh_score + hmd_mtd_score + fdhA_score + fdhB_score
                            
                            # 根据得分确定匹配基因数
                            # 如果所有关键功能基因都满足，直接认定为100%完整度
                            if fwd_fgh_present and hmd_mtd_present and fdhA_present and fdhB_present:
                                self.logger.debug(f"JIASUAN-CH4: 所有关键功能基因满足，使用参考序列数: {reference_count}")
                                low_threshold_matched = reference_count
                            elif total_score >= 80:  # 80分以上认为途径完整
                                self.logger.debug(f"JIASUAN-CH4: 总得分{total_score}≥80，使用参考序列数: {reference_count}")
                                low_threshold_matched = reference_count
                            elif total_score >= 60:  # 60-79分使用实际检测数
                                self.logger.debug(f"JIASUAN-CH4: 总得分{total_score}在60-79之间，使用实际检测数: {len(unique_hits)}")
                                low_threshold_matched = len(unique_hits)
                            else:  # 60分以下使用实际检测数
                                self.logger.debug(f"JIASUAN-CH4: 总得分{total_score}<60，使用实际检测数: {len(unique_hits)}")
                                low_threshold_matched = len(unique_hits)
                        
                        # 对于CO2还原产甲烷途径，添加特殊化处理
                        elif db_name == 'CO2-CH4':
                            # 调试信息：确认代码被执行
                            self.logger.debug("=== CO2-CH4 通路评判逻辑开始执行 ===")
                            
                            # 获取所有比对到的唯一基因名称（第二列）
                            unique_hits = result['low_threshold_hits'].iloc[:, 1].unique()
                            
                            # 检查关键基因的比对情况
                            # CO2-CH4通路的关键基因是FwdF/FwdG/FwdH和Hmd/Mtd，不需要FdhA/FdhB
                            fwd_fgh_present = any(gene in unique_hits for gene in ['FwdF', 'FwdG', 'FwdH'])
                            hmd_mtd_present = any(gene in unique_hits for gene in ['Hmd', 'Mtd'])
                            
                            # 调试信息：输出关键基因检测状态
                            self.logger.debug(f"CO2-CH4 关键基因检测状态:")
                            self.logger.debug(f"  FwdF/FwdG/FwdH: {fwd_fgh_present}")
                            self.logger.debug(f"  Hmd/Mtd: {hmd_mtd_present}")
                            self.logger.debug(f"  基础基因覆盖: {len(unique_hits)}/{reference_count}")
                            
                            # 优化评判逻辑：采用加权评分法
                            # 基础基因覆盖（权重60%）：检测到的基因数量
                            base_score = min(len(unique_hits) / reference_count * 60, 60)
                            
                            # 关键功能基因（权重40%）：
                            # - FwdF/FwdG/FwdH（权重20%）
                            # - Hmd/Mtd（权重20%）
                            fwd_fgh_score = 20 if fwd_fgh_present else 0
                            hmd_mtd_score = 20 if hmd_mtd_present else 0
                            
                            # 计算总得分
                            total_score = base_score + fwd_fgh_score + hmd_mtd_score
                            
                            # 根据得分确定匹配基因数
                            # 如果所有关键功能基因都满足，直接认定为100%完整度
                            if fwd_fgh_present and hmd_mtd_present:
                                self.logger.debug(f"CO2-CH4: 所有关键功能基因满足，使用参考序列数: {reference_count}")
                                low_threshold_matched = reference_count
                            elif total_score >= 80:  # 80分以上认为途径完整
                                self.logger.debug(f"CO2-CH4: 总得分{total_score}≥80，使用参考序列数: {reference_count}")
                                low_threshold_matched = reference_count
                            elif total_score >= 60:  # 60-79分使用实际检测数
                                self.logger.debug(f"CO2-CH4: 总得分{total_score}在60-79之间，使用实际检测数: {len(unique_hits)}")
                                low_threshold_matched = len(unique_hits)
                            else:  # 60分以下使用实际检测数
                                self.logger.debug(f"CO2-CH4: 总得分{total_score}<60，使用实际检测数: {len(unique_hits)}")
                                low_threshold_matched = len(unique_hits)
                        
                        # 对于甲硫醇产甲烷途径，添加特殊化处理
                        elif db_name == 'JIALIUCHUN-CH4':
                            # 调试信息：确认代码被执行
                            self.logger.debug("=== JIALIUCHUN-CH4 通路评判逻辑开始执行 ===")
                            
                            # 获取所有比对到的唯一基因名称（第二列）
                            unique_hits = result['low_threshold_hits'].iloc[:, 1].unique()
                            
                            # 检查关键基因的比对情况
                            # 甲硫醇特异性基因：MtsA1, MtsA2（添加严格的比对条件：bitscore > 200 且 evalue <= 1e-100）
                            mtsA1_hits = result['low_threshold_hits'][
                                (result['low_threshold_hits'].iloc[:, 1] == 'MtsA1') &
                                (result['low_threshold_hits']['bitscore'] > 200) &
                                (result['low_threshold_hits']['evalue'] <= 1e-100)
                            ]
                            mtsA2_hits = result['low_threshold_hits'][
                                (result['low_threshold_hits'].iloc[:, 1] == 'MtsA2') &
                                (result['low_threshold_hits']['bitscore'] > 200) &
                                (result['low_threshold_hits']['evalue'] <= 1e-100)
                            ]
                            mtsA1_present = len(mtsA1_hits) > 0
                            mtsA2_present = len(mtsA2_hits) > 0
                            
                            # 甲基辅酶M还原酶相关基因：KYC55281.1, KYC55283.1, KYC55284.1, KYC55314.1
                            mcr_genes = ['KYC55281.1', 'KYC55283.1', 'KYC55284.1', 'KYC55314.1']
                            mcr_genes_present = [gene in unique_hits for gene in mcr_genes]
                            mcr_genes_count = sum(mcr_genes_present)
                            
                            # 检查甲基辅酶M还原酶相关基因中任意3个是否比对到
                            mcr_3_present = mcr_genes_count >= 3
                            
                            # 调试信息：输出关键基因检测状态
                            self.logger.debug(f"JIALIUCHUN-CH4 关键基因检测状态:")
                            self.logger.debug(f"  MtsA1: {mtsA1_present} (命中数: {len(mtsA1_hits)})")
                            self.logger.debug(f"  MtsA2: {mtsA2_present} (命中数: {len(mtsA2_hits)})")
                            self.logger.debug(f"  甲基辅酶M还原酶基因比对情况: {mcr_genes_count}/4")
                            self.logger.debug(f"  甲基辅酶M还原酶基因任意3个比对到: {mcr_3_present}")
                            self.logger.debug(f"  基础基因覆盖: {len(unique_hits)}/{reference_count}")
                            
                            # 输出详细的MtsA1和MtsA2比对分数信息用于调试
                            if len(mtsA1_hits) > 0:
                                self.logger.debug(f"  MtsA1最高bitscore: {mtsA1_hits['bitscore'].max()}")
                            if len(mtsA2_hits) > 0:
                                self.logger.debug(f"  MtsA2最高bitscore: {mtsA2_hits['bitscore'].max()}")
                            
                            # 根据新策略确定匹配基因数
                            # 如果MtsA1和MtsA2都比对到，且甲基辅酶M还原酶相关基因中任意3个比对到，则认定为100%完整度
                            if mtsA1_present and mtsA2_present and mcr_3_present:
                                self.logger.debug(f"JIALIUCHUN-CH4: 关键基因条件满足，使用参考序列数: {reference_count}")
                                low_threshold_matched = reference_count
                            else:
                                # 不满足条件时，如果MtsA1未满足严格条件，从实际检测数中排除MtsA1基因
                                # 因为MtsA1是甲硫醇通路的关键基因，必须满足严格条件
                                adjusted_hits = len(unique_hits)
                                if not mtsA1_present:
                                    # 如果MtsA1未满足条件，从检测数中减去1（因为MtsA1是必需基因）
                                    adjusted_hits = max(0, len(unique_hits) - 1)
                                    self.logger.debug(f"JIALIUCHUN-CH4: MtsA1未满足严格条件，调整检测数: {len(unique_hits)} -> {adjusted_hits}")
                                else:
                                    self.logger.debug(f"JIALIUCHUN-CH4: 关键基因条件不满足，使用实际检测数: {len(unique_hits)}")
                                low_threshold_matched = adjusted_hits
                        
                        # 对于甘氨酸甜菜碱产甲烷途径，添加特殊化处理
                        elif db_name == '甘氨酸甜菜碱产甲烷':
                            # 调试信息：确认代码被执行
                            self.logger.debug("=== 甘氨酸甜菜碱产甲烷通路评判逻辑开始执行 ===")
                            
                            # 获取所有比对到的唯一基因名称（第二列）
                            unique_hits = result['low_threshold_hits'].iloc[:, 1].unique()
                            
                            # 检查关键基因的比对情况
                            # 甘氨酸甜菜碱产甲烷途径的关键基因：MtgB, dimethylamine_corrinoid_protein_3, MV10360
                            mtgB_present = 'MtgB' in unique_hits
                            dimethylamine_corrinoid_present = 'dimethylamine_corrinoid_protein_3' in unique_hits
                            mv10360_present = 'MV10360' in unique_hits
                            
                            # 检查所有关键基因是否都比对到
                            all_key_genes_present = mtgB_present and dimethylamine_corrinoid_present and mv10360_present
                            
                            # 调试信息：输出关键基因检测状态
                            self.logger.debug(f"甘氨酸甜菜碱产甲烷 关键基因检测状态:")
                            self.logger.debug(f"  MtgB: {mtgB_present}")
                            self.logger.debug(f"  dimethylamine_corrinoid_protein_3: {dimethylamine_corrinoid_present}")
                            self.logger.debug(f"  MV10360: {mv10360_present}")
                            self.logger.debug(f"  所有关键基因比对到: {all_key_genes_present}")
                            self.logger.debug(f"  基础基因覆盖: {len(unique_hits)}/{reference_count}")
                            
                            # 根据新策略确定匹配基因数
                            # 如果所有3个关键基因都比对到，则认定为100%完整度
                            if all_key_genes_present:
                                self.logger.debug(f"甘氨酸甜菜碱产甲烷: 所有关键基因条件满足，使用参考序列数: {reference_count}")
                                low_threshold_matched = reference_count
                            else:
                                # 不满足条件时使用实际检测数
                                self.logger.debug(f"甘氨酸甜菜碱产甲烷: 关键基因条件不满足，使用实际检测数: {len(unique_hits)}")
                                low_threshold_matched = len(unique_hits)
                        
                        # 对于硫代丙酸甲酯产甲烷途径，添加特殊化处理
                        elif db_name == '硫代丙酸甲酯产甲烷':
                            # 调试信息：确认代码被执行
                            self.logger.debug("=== 硫代丙酸甲酯产甲烷通路评判逻辑开始执行 ===")
                            
                            # 获取所有比对到的唯一基因名称（第二列）
                            unique_hits = result['low_threshold_hits'].iloc[:, 1].unique()
                            
                            # 检查关键基因的比对情况，要求比对分数大于100且E值小于1e-5
                            # 硫代丙酸甲酯产甲烷途径的关键基因：mtpA1, mtsA1, mtpA2, mtsA2
                            mtpA1_hits = result['low_threshold_hits'][
                                (result['low_threshold_hits'].iloc[:, 1] == 'mtpA1') &
                                (result['low_threshold_hits']['bitscore'] > 100) &
                                (result['low_threshold_hits']['evalue'] <= 1e-5)
                            ]
                            mtsA1_hits = result['low_threshold_hits'][
                                (result['low_threshold_hits'].iloc[:, 1] == 'mtsA1') &
                                (result['low_threshold_hits']['bitscore'] > 100) &
                                (result['low_threshold_hits']['evalue'] <= 1e-5)
                            ]
                            mtpA2_hits = result['low_threshold_hits'][
                                (result['low_threshold_hits'].iloc[:, 1] == 'mtpA2') &
                                (result['low_threshold_hits']['bitscore'] > 100) &
                                (result['low_threshold_hits']['evalue'] <= 1e-5)
                            ]
                            mtsA2_hits = result['low_threshold_hits'][
                                (result['low_threshold_hits'].iloc[:, 1] == 'mtsA2') &
                                (result['low_threshold_hits']['bitscore'] > 100) &
                                (result['low_threshold_hits']['evalue'] <= 1e-5)
                            ]
                            
                            mtpA1_present = len(mtpA1_hits) > 0
                            mtsA1_present = len(mtsA1_hits) > 0
                            mtpA2_present = len(mtpA2_hits) > 0
                            mtsA2_present = len(mtsA2_hits) > 0
                            
                            # 检查所有关键基因是否都比对到
                            all_key_genes_present = mtpA1_present and mtsA1_present and mtpA2_present and mtsA2_present
                            
                            # 调试信息：输出关键基因检测状态
                            self.logger.debug(f"硫代丙酸甲酯产甲烷 关键基因检测状态:")
                            self.logger.debug(f"  mtpA1: {mtpA1_present} (命中数: {len(mtpA1_hits)})")
                            self.logger.debug(f"  mtsA1: {mtsA1_present} (命中数: {len(mtsA1_hits)})")
                            self.logger.debug(f"  mtpA2: {mtpA2_present} (命中数: {len(mtpA2_hits)})")
                            self.logger.debug(f"  mtsA2: {mtsA2_present} (命中数: {len(mtsA2_hits)})")
                            self.logger.debug(f"  所有关键基因比对到: {all_key_genes_present}")
                            self.logger.debug(f"  基础基因覆盖: {len(unique_hits)}/{reference_count}")
                            
                            # 根据新策略确定匹配基因数
                            # 如果所有4个关键基因都比对到（比对分数>100且E值<1e-5），则认定为100%完整度
                            if all_key_genes_present:
                                self.logger.debug(f"硫代丙酸甲酯产甲烷: 所有关键基因条件满足，使用参考序列数: {reference_count}")
                                low_threshold_matched = reference_count
                            else:
                                # 不满足条件时使用实际检测数
                                self.logger.debug(f"硫代丙酸甲酯产甲烷: 关键基因条件不满足，使用实际检测数: {len(unique_hits)}")
                                low_threshold_matched = len(unique_hits)
                        
                        # 对于四甲基铵产甲烷途径，添加特殊化处理
                        elif db_name == '四甲基铵产甲烷':
                            # 调试信息：确认代码被执行
                            self.logger.debug("=== 四甲基铵产甲烷通路评判逻辑开始执行 ===")
                            
                            # 获取所有比对到的唯一基因名称（第二列）
                            unique_hits = result['low_threshold_hits'].iloc[:, 1].unique()
                            
                            # 检查关键基因的比对情况，要求比对分数大于200且E值小于1e-100
                            # 四甲基铵产甲烷途径的关键基因：MtqA/MT2, MtqB, MtqC
                            mtqA_hits = result['low_threshold_hits'][
                                (result['low_threshold_hits'].iloc[:, 1] == 'MtqA/MT2') &
                                (result['low_threshold_hits']['bitscore'] > 200) &
                                (result['low_threshold_hits']['evalue'] <= 1e-100)
                            ]
                            mtqB_hits = result['low_threshold_hits'][
                                (result['low_threshold_hits'].iloc[:, 1] == 'MtqB') &
                                (result['low_threshold_hits']['bitscore'] > 200) &
                                (result['low_threshold_hits']['evalue'] <= 1e-100)
                            ]
                            mtqC_hits = result['low_threshold_hits'][
                                (result['low_threshold_hits'].iloc[:, 1] == 'MtqC') &
                                (result['low_threshold_hits']['bitscore'] > 200) &
                                (result['low_threshold_hits']['evalue'] <= 1e-100)
                            ]
                            
                            mtqA_present = len(mtqA_hits) > 0
                            mtqB_present = len(mtqB_hits) > 0
                            mtqC_present = len(mtqC_hits) > 0
                            
                            # 检查所有关键基因是否都比对到
                            all_key_genes_present = mtqA_present and mtqB_present and mtqC_present
                            
                            # 调试信息：输出关键基因检测状态
                            self.logger.debug(f"四甲基铵产甲烷 关键基因检测状态:")
                            self.logger.debug(f"  MtqA/MT2: {mtqA_present} (命中数: {len(mtqA_hits)})")
                            self.logger.debug(f"  MtqB: {mtqB_present} (命中数: {len(mtqB_hits)})")
                            self.logger.debug(f"  MtqC: {mtqC_present} (命中数: {len(mtqC_hits)})")
                            self.logger.debug(f"  所有关键基因比对到: {all_key_genes_present}")
                            self.logger.debug(f"  基础基因覆盖: {len(unique_hits)}/{reference_count}")
                            
                            # 根据新策略确定匹配基因数
                            # 如果所有3个关键基因都比对到（比对分数>200且E值<1e-100），则认定为100%完整度
                            if all_key_genes_present:
                                self.logger.debug(f"四甲基铵产甲烷: 所有关键基因条件满足，使用参考序列数: {reference_count}")
                                low_threshold_matched = reference_count
                            else:
                                # 不满足条件时使用实际检测数
                                self.logger.debug(f"四甲基铵产甲烷: 关键基因条件不满足，使用实际检测数: {len(unique_hits)}")
                                low_threshold_matched = len(unique_hits)
                        
                        # 对于甲醇歧化产甲烷途径，添加特殊化处理
                        elif db_name == '甲醇歧化产甲烷':
                            # 调试信息：确认代码被执行
                            self.logger.debug("=== 甲醇歧化产甲烷通路评判逻辑开始执行 ===")
                            
                            # 获取所有比对到的唯一基因名称（第二列）
                            unique_hits = result['low_threshold_hits'].iloc[:, 1].unique()
                            
                            # 检查关键基因的比对情况
                            # 甲醇歧化产甲烷途径的关键基因：MvhA和elpA（两个基因中有一个比对到即可）
                            mvhA_present = 'MvhA' in unique_hits
                            elpA_present = 'elpA' in unique_hits
                            
                            # 检查其他基因的比对情况（除了MvhA和elpA之外的所有基因）
                            # 甲醇歧化产甲烷途径的其他关键基因：FwdA-FwdH, Ftr, Mch, Hmd, Mtd, elpB, elpC等
                            other_genes = ['FwdA', 'FwdB', 'FwdC', 'FwdD', 'FwdE', 'FwdF', 'FwdG', 'FwdH', 
                                         'Ftr', 'Mch', 'Hmd', 'Mtd', 'elpB', 'elpC']
                            other_genes_present = [gene in unique_hits for gene in other_genes]
                            other_genes_count = sum(other_genes_present)
                            
                            # 检查其他基因是否全部比对到（除了MvhA和elpA）
                            all_other_genes_present = other_genes_count == len(other_genes)
                            
                            # 调试信息：输出关键基因检测状态
                            self.logger.debug(f"甲醇歧化产甲烷 关键基因检测状态:")
                            self.logger.debug(f"  MvhA: {mvhA_present}")
                            self.logger.debug(f"  elpA: {elpA_present}")
                            self.logger.debug(f"  MvhA或elpA比对到: {mvhA_present or elpA_present}")
                            self.logger.debug(f"  其他基因比对情况: {other_genes_count}/{len(other_genes)}")
                            self.logger.debug(f"  所有其他基因比对到: {all_other_genes_present}")
                            self.logger.debug(f"  基础基因覆盖: {len(unique_hits)}/{reference_count}")
                            
                            # 根据新策略确定匹配基因数
                            # 如果其他基因全部比对到，且MvhA和elpA中有一个比对到，则认定为100%完整度
                            if all_other_genes_present and (mvhA_present or elpA_present):
                                self.logger.debug(f"甲醇歧化产甲烷: 其他基因全部比对到且MvhA/elpA中有一个比对到，使用参考序列数: {reference_count}")
                                low_threshold_matched = reference_count
                            else:
                                # 不满足条件时使用实际检测数
                                self.logger.debug(f"甲醇歧化产甲烷: 条件不满足，使用实际检测数: {len(unique_hits)}")
                                low_threshold_matched = len(unique_hits)
                            
                            # 确保完整度不超过100%
                            if low_threshold_matched > reference_count:
                                self.logger.debug(f"甲醇歧化产甲烷: 检测基因数({low_threshold_matched})超过参考序列数({reference_count})，限制为100%")
                                low_threshold_matched = reference_count
                        
                        low_completeness = (low_threshold_matched / reference_count) * 100
                        
                        # 确保完整度百分比不超过100%
                        if low_completeness > 100:
                            self.logger.debug(f"完整度百分比超过100%，限制为100%")
                            low_completeness = 100
                        
                        pathway_results[db_name] = {
                            'pathway_name': self.pathway_names.get(db_name, db_name),
                            'low_threshold_hits': low_threshold_matched,
                            'reference_sequences': reference_count,
                            'low_completeness_percentage': round(low_completeness, 2),
                            'average_identity_low': round(
                                result['low_threshold_hits']['pident'].mean(), 2
                            ) if len(result['low_threshold_hits']) > 0 else 0,
                            'best_hits_low': result['low_threshold_hits'].nlargest(5, 'bitscore').to_dict('records')
                        }
                    else:
                        # 其他代谢途径（硫、氮）保持双阈值分析
                        high_threshold_matched = len(result['high_threshold_hits'].drop_duplicates(subset=['sseqid']))
                        low_threshold_matched = len(result['low_threshold_hits'].drop_duplicates(subset=['sseqid']))
                        
                        # 确保检测基因数不超过参考序列数
                        if high_threshold_matched > reference_count:
                            self.logger.debug(f"{db_name}: 高阈值检测基因数({high_threshold_matched})超过参考序列数({reference_count})，限制为100%")
                            high_threshold_matched = reference_count
                        if low_threshold_matched > reference_count:
                            self.logger.debug(f"{db_name}: 低阈值检测基因数({low_threshold_matched})超过参考序列数({reference_count})，限制为100%")
                            low_threshold_matched = reference_count
                        
                        high_completeness = (high_threshold_matched / reference_count) * 100
                        low_completeness = (low_threshold_matched / reference_count) * 100
                        
                        # 确保完整度百分比不超过100%
                        if high_completeness > 100:
                            self.logger.debug(f"{db_name}: 高阈值完整度超过100%，限制为100%")
                            high_completeness = 100
                        if low_completeness > 100:
                            self.logger.debug(f"{db_name}: 低阈值完整度超过100%，限制为100%")
                            low_completeness = 100
                        
                        # 代谢途径完整度格式："高阈值%~低阈值%"
                        metabolic_pathway_completeness = f"{high_completeness:.2f}%~{low_completeness:.2f}%"
                        
                        pathway_results[db_name] = {
                            'pathway_name': self.pathway_names.get(db_name, db_name),
                            'high_threshold_hits': high_threshold_matched,
                            'low_threshold_hits': low_threshold_matched,
                            'reference_sequences': reference_count,
                            'high_completeness_percentage': round(high_completeness, 2),
                            'low_completeness_percentage': round(low_completeness, 2),
                            'metabolic_pathway_completeness': metabolic_pathway_completeness,
                            'average_identity_high': round(
                                result['high_threshold_hits']['pident'].mean(), 2
                            ) if len(result['high_threshold_hits']) > 0 else 0,
                            'average_identity_low': round(
                                result['low_threshold_hits']['pident'].mean(), 2
                            ) if len(result['low_threshold_hits']) > 0 else 0,
                            'best_hits_high': result['high_threshold_hits'].nlargest(5, 'bitscore').to_dict('records'),
                            'best_hits_low': result['low_threshold_hits'].nlargest(5, 'bitscore').to_dict('records')
                        }
                
                summary_stats['pathways_detected'] += 1
                if 'low_threshold_hits' in result:
                    summary_stats['total_hits'] += len(result['low_threshold_hits'])
                elif 'cultivability_hits' in result:
                    summary_stats['total_hits'] += len(result['cultivability_hits'])
                elif 'salt_tolerance_hits' in result:
                    summary_stats['total_hits'] += len(result['salt_tolerance_hits'])
        
        return {
            'input_file': str(input_file),
            'analysis_timestamp': pd.Timestamp.now().isoformat(),
            'summary': summary_stats,
            'pathway_results': pathway_results,
            'raw_results': search_results
        }
    
    def _save_results(self, results: Dict, output_prefix: str):
        """
        Save analysis results to files with simplified output
        
        Args:
            results: Analysis results
            output_prefix: Output file prefix
        """
        # 创建可序列化的结果副本，将DataFrame转换为字典列表
        serializable_results = results.copy()
        
        # 处理pathway_results中的DataFrame对象
        for db_name, pathway_data in serializable_results['pathway_results'].items():
            # 将best_hits中的DataFrame转换为字典列表
            if 'best_hits_high' in pathway_data and hasattr(pathway_data['best_hits_high'], 'to_dict'):
                pathway_data['best_hits_high'] = pathway_data['best_hits_high'].to_dict('records')
            if 'best_hits_low' in pathway_data and hasattr(pathway_data['best_hits_low'], 'to_dict'):
                pathway_data['best_hits_low'] = pathway_data['best_hits_low'].to_dict('records')
        
        # 处理raw_results中的DataFrame对象
        if 'raw_results' in serializable_results:
            for db_name, raw_data in serializable_results['raw_results'].items():
                if raw_data:
                    # 将DataFrame转换为字典列表
                    for key in ['high_threshold_hits', 'low_threshold_hits', 'cultivability_hits', 'salt_tolerance_hits']:
                        if key in raw_data and hasattr(raw_data[key], 'to_dict'):
                            raw_data[key] = raw_data[key].to_dict('records')
        
        # 只保存一个汇总CSV文件，删除重复的JSON和CSV文件输出
        output_dir = os.path.dirname(output_prefix) if os.path.dirname(output_prefix) else self.results_dir
        summary_csv_file = os.path.join(output_dir, f"{os.path.basename(output_prefix)}_metabolic_summary.csv")
        
        # 准备汇总数据
        summary_data = []
        
        # 获取输入文件信息
        input_file = results.get('input_file', 'Unknown')
        analysis_timestamp = results.get('analysis_timestamp', '')
        
        # 处理每个代谢途径的结果
        for db_name, pathway_data in results['pathway_results'].items():
            # 确定代谢途径类型
            pathway_type = 'Other'
            if db_name in ['CO2-CH4', 'JIAAN-CH4', 'JIACHUN-CH4', 'JIALIUCHUN-CH4', 
                         'YISUAN-CH4', 'C16-CH4', 'CO-CH4', 'JIASUAN-CH4', 
                         'JIAYANGJI-CH4', 'ZHIFANGSUAN-CH4', '2JIAAN-CH4', '3JIAAN-CH4',
                         '甘氨酸甜菜碱产甲烷', '硫代丙酸甲酯产甲烷', '四甲基铵产甲烷', '甲醇歧化产甲烷']:
                pathway_type = 'Methane'
            elif db_name in ['ASR', 'SO', 'SOX', 'S4I', 'SR', 'DSR']:
                pathway_type = 'Sulfur'
            elif db_name in ['ANR', 'DEN', 'DNR', 'NIT']:
                pathway_type = 'Nitrogen'
            elif db_name == 'CULTIVATION':
                pathway_type = 'Cultivation'
            elif db_name == 'NAIYAN':
                pathway_type = 'Salt_Tolerance'
            
            # 根据代谢途径类型提取相应的数据
            if pathway_type == 'Methane':
                # 甲烷代谢途径
                summary_data.append({
                    'Input_File': input_file,
                    'Analysis_Timestamp': analysis_timestamp,
                    'Pathway_Type': pathway_type,
                    'Database': db_name,
                    'Pathway_Name': pathway_data['pathway_name'],
                    'Reference_Sequences': pathway_data['reference_sequences'],
                    'Low_Threshold_Hits': pathway_data.get('low_threshold_hits', 0),
                    'Low_Completeness_Percentage': pathway_data.get('low_completeness_percentage', 0),
                    'Average_Identity_Low': pathway_data.get('average_identity_low', 0),
                    'Detection_Status': 'Detected' if pathway_data.get('low_threshold_hits', 0) > 0 else 'Not Detected'
                })
            elif pathway_type in ['Sulfur', 'Nitrogen']:
                # 硫和氮代谢途径
                summary_data.append({
                    'Input_File': input_file,
                    'Analysis_Timestamp': analysis_timestamp,
                    'Pathway_Type': pathway_type,
                    'Database': db_name,
                    'Pathway_Name': pathway_data['pathway_name'],
                    'Reference_Sequences': pathway_data['reference_sequences'],
                    'High_Threshold_Hits': pathway_data.get('high_threshold_hits', 0),
                    'Low_Threshold_Hits': pathway_data.get('low_threshold_hits', 0),
                    'High_Completeness_Percentage': pathway_data.get('high_completeness_percentage', 0),
                    'Low_Completeness_Percentage': pathway_data.get('low_completeness_percentage', 0),
                    'Average_Identity_High': pathway_data.get('average_identity_high', 0),
                    'Average_Identity_Low': pathway_data.get('average_identity_low', 0),
                    'Metabolic_Pathway_Completeness': pathway_data.get('metabolic_pathway_completeness', ''),
                    'Detection_Status': 'Detected' if pathway_data.get('low_threshold_hits', 0) > 0 else 'Not Detected'
                })
            elif pathway_type == 'Cultivation':
                # 栽培评估
                summary_data.append({
                    'Input_File': input_file,
                    'Analysis_Timestamp': analysis_timestamp,
                    'Pathway_Type': pathway_type,
                    'Database': db_name,
                    'Pathway_Name': pathway_data['pathway_name'],
                    'Reference_Sequences': pathway_data['reference_sequences'],
                    'Cultivability_Hits': pathway_data.get('cultivability_hits', 0),
                    'Cultivability_Percentage': pathway_data.get('cultivability_percentage', 0),
                    'Cultivability_Status': pathway_data.get('cultivability_status', ''),
                    'Average_Identity': pathway_data.get('average_identity', 0),
                    'Detection_Status': 'Detected' if pathway_data.get('cultivability_hits', 0) > 0 else 'Not Detected'
                })
            elif pathway_type == 'Salt_Tolerance':
                # 耐盐性评估
                summary_data.append({
                    'Input_File': input_file,
                    'Analysis_Timestamp': analysis_timestamp,
                    'Pathway_Type': pathway_type,
                    'Database': db_name,
                    'Pathway_Name': pathway_data['pathway_name'],
                    'Reference_Sequences': pathway_data['reference_sequences'],
                    'Salt_Tolerance_Hits': pathway_data.get('salt_tolerance_hits', 0),
                    'Salt_Tolerance_Level': pathway_data.get('salt_tolerance_level', ''),
                    'Average_Identity': pathway_data.get('average_identity', 0),
                    'Detection_Status': 'Detected' if pathway_data.get('salt_tolerance_hits', 0) > 0 else 'Not Detected'
                })
        
        if summary_data:
            df_summary = pd.DataFrame(summary_data)
            df_summary.to_csv(summary_csv_file, index=False)
        
        self.logger.info(f"Results saved to {summary_csv_file}")
    
    def get_pathway_summary(self, results: Dict) -> Dict[str, Dict]:
        """
        Get pathway analysis summary
        
        Args:
            results: Analysis results
            
        Returns:
            Pathway summary dictionary
        """
        pathway_summary = {
            'methane_pathways': {},
            'sulfur_pathways': {},
            'nitrogen_pathways': {},
            'other_features': {}
        }
        
        # Categorize pathways - use consistent logic with _save_results method
        methane_dbs = ['CO2-CH4', 'JIAAN-CH4', 'JIACHUN-CH4', 'JIALIUCHUN-CH4', 
                     'YISUAN-CH4', 'C16-CH4', 'CO-CH4', 'JIASUAN-CH4', 
                     'JIAYANGJI-CH4', 'ZHIFANGSUAN-CH4', '2JIAAN-CH4', '3JIAAN-CH4',
                     '甘氨酸甜菜碱产甲烷', '硫代丙酸甲酯产甲烷', '四甲基铵产甲烷', '甲醇歧化产甲烷']
        sulfur_dbs = ['ASR', 'SO', 'SOX', 'S4I', 'SR', 'DSR']
        nitrogen_dbs = ['ANR', 'DEN', 'DNR', 'NIT']
        
        for db_name, pathway_data in results['pathway_results'].items():
            if db_name in methane_dbs:
                pathway_summary['methane_pathways'][db_name] = pathway_data
            elif db_name in sulfur_dbs:
                pathway_summary['sulfur_pathways'][db_name] = pathway_data
            elif db_name in nitrogen_dbs:
                pathway_summary['nitrogen_pathways'][db_name] = pathway_data
            else:
                pathway_summary['other_features'][db_name] = pathway_data
        
        return pathway_summary
    
    def cleanup_temp_files(self, output_prefix: str):
        """
        Clean up temporary files
        
        Args:
            output_prefix: Output file prefix
        """
        temp_files = []
        for file in os.listdir(self.results_dir):
            if file.startswith(output_prefix) and file.endswith('_diamond.tsv'):
                temp_files.append(os.path.join(self.results_dir, file))
        
        FileUtils.cleanup_temp_files(temp_files)
        self.logger.debug(f"Cleaned up {len(temp_files)} temporary files")
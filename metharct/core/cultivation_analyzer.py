# -*- coding: utf-8 -*-
"""
Cultivation Analyzer Module

This module provides cultivation analysis functionality for MethArCT.
It integrates the new cultivation analysis package to assess microbial
cultivability based on metabolic pathway analysis.
"""

import os
import subprocess
import tempfile
import pandas as pd
from typing import Dict, List, Any, Tuple, Optional
import re
import json
import concurrent.futures
import threading
from datetime import datetime
import time
import glob

class CultivationAnalyzer:
    """
    Cultivation Analyzer - 基于代谢通路分析的微生物培养性评估
    
    该类使用Diamond序列比对和氨基酸代谢通路分析来评估微生物的可培养性，
    专注于氨基酸、维生素、辅酶、核酸合成和ATP酶代谢途径分析。
    """
    
    def __init__(self, config: Dict[str, Any]):
        """
        初始化培养性分析器
        
        Args:
            config: 配置字典，包含工具路径、数据库路径等参数
        """
        # 从配置中获取参数
        self.diamond_path = config.get('tools', {}).get('diamond', {}).get('path', 'diamond.exe')
        self.use_wsl = config.get('tools', {}).get('diamond', {}).get('use_wsl', False)
        self.threads = config.get('tools', {}).get('diamond', {}).get('threads', 4)
        self.evalue = config.get('tools', {}).get('diamond', {}).get('evalue', 1e-5)
        
        # 数据库路径 - 支持绝对路径和相对路径
        base_dir = config.get('databases', {}).get('base_dir', 'data/databases')
        
        # 如果base_dir是相对路径，则相对于当前工作目录
        if not os.path.isabs(base_dir):
            base_dir = os.path.abspath(base_dir)
            
        self.cultivation_dir = os.path.join(base_dir, 'cultivation')
        
        # 各类代谢通路目录 - 更新为新的数据库结构
        self.amino_acid_dir = os.path.join(self.cultivation_dir, 'path.aa')
        self.vitamin_dir = os.path.join(self.cultivation_dir, 'path.vc')
        self.nucleotide_dir = os.path.join(self.cultivation_dir, 'path_hesuan')
        self.atp_dir = os.path.join(self.cultivation_dir, 'path_atp')
        
        # 从目录中自动获取所有代谢途径文件
        self.amino_acid_pathways = self._discover_pathways(self.amino_acid_dir, "氨基酸")
        self.vitamin_pathways = self._discover_pathways(self.vitamin_dir, "维生素和辅酶")
        self.nucleotide_pathways = self._discover_pathways(self.nucleotide_dir, "核酸合成")
        self.atp_pathways = self._discover_pathways(self.atp_dir, "ATP酶")
        
        total_pathways = len(self.amino_acid_pathways) + len(self.vitamin_pathways) + len(self.nucleotide_pathways) + len(self.atp_pathways)
        print(f"初始化培养性分析器，将分析 {total_pathways} 个代谢通路")
        print(f"  - 氨基酸代谢通路: {len(self.amino_acid_pathways)} 个")
        print(f"  - 维生素和辅酶代谢通路: {len(self.vitamin_pathways)} 个")
        print(f"  - 核酸合成途径: {len(self.nucleotide_pathways)} 个")
        print(f"  - ATP酶代谢途径: {len(self.atp_pathways)} 个")
    
    def _discover_pathways(self, directory: str, pathway_type: str) -> Dict[str, str]:
        """
        从指定目录中发现所有代谢途径
        
        Args:
            directory: 包含FASTA文件的目录路径
            pathway_type: 通路类型描述（用于日志输出）
            
        Returns:
            字典 {pathway_name: file_path}
        """
        pathways = {}
        
        if not os.path.exists(directory):
            print(f"警告: {pathway_type}FASTA目录不存在: {directory}")
            return pathways
            
        # 查找所有.fasta文件
        fasta_files = glob.glob(os.path.join(directory, "*.fasta"))
        
        for file_path in fasta_files:
            file_name = os.path.basename(file_path)
            pathway_name = file_name.replace('.fasta', '')
            pathways[pathway_name] = file_path
            
        return pathways
    
    def _parse_fasta_file(self, fasta_file: str) -> Dict[str, Dict[str, str]]:
        """
        解析FASTA文件，提取基因和序列
        
        Args:
            fasta_file: FASTA文件路径
            
        Returns:
            基因信息字典 {gene_name: {sequence: str, description: str}}
        """
        genes = {}
        
        if not os.path.exists(fasta_file):
            print(f"警告: FASTA文件不存在: {fasta_file}")
            return genes
            
        with open(fasta_file, 'r') as f:
            current_gene = None
            current_seq = []
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # 保存前一个基因
                    if current_gene:
                        genes[current_gene] = {
                            'sequence': ''.join(current_seq),
                            'description': description
                        }
                    
                    # 解析新基因头部
                    header = line[1:]  # 移除'>'
                    
                    # 从新格式的头部信息中提取基因名
                    # 格式如: >eco:b0907 K00831 phosphoserine aminotransferase [EC:2.6.1.52] | (RefSeq) serC; phosphoserine/phosphohydroxythreonine aminotransferase (A)
                    # 提取基因名 (如 serC)
                    gene_match = re.search(r'\((RefSeq|GenBank)\)\s*([^;]+);', header)
                    if gene_match:
                        current_gene = gene_match.group(2).strip()
                    else:
                        # 如果没有找到，尝试提取物种:基因格式中的基因部分
                        species_gene_match = re.match(r'([^:]+):([^:]+)', header.split()[0])
                        if species_gene_match:
                            current_gene = species_gene_match.group(2)
                        else:
                            # 最后尝试提取第一个单词作为基因名
                            current_gene = header.split()[0]
                    
                    description = header
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # 保存最后一个基因
            if current_gene:
                genes[current_gene] = {
                    'sequence': ''.join(current_seq),
                    'description': description
                }
        
        return genes
    
    def _run_diamond_search(self, query_seq: str, genome_file: str, evalue: float = 1e-5) -> Tuple[bool, Dict[str, Any]]:
        """
        使用Diamond进行序列比对
        
        Args:
            query_seq: 查询序列
            genome_file: 基因组蛋白质序列文件路径
            evalue: E值阈值
            
        Returns:
            (是否找到匹配, 匹配详情)
        """
        try:
            # 创建临时查询文件
            with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.faa') as temp_query:
                temp_query.write(f">query\n{query_seq}\n")
                temp_query_path = temp_query.name
            
            # 创建临时输出文件
            with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.txt') as temp_output:
                temp_output_path = temp_output.name
            
            # 创建临时数据库文件
            temp_db_path = temp_query_path.replace('.faa', '.dmnd')
            
            # 构建Diamond命令
            diamond_cmd = self.diamond_path
            if self.use_wsl:
                diamond_cmd = f"wsl {diamond_cmd}"
            
            # 为基因组文件创建Diamond数据库
            make_db_cmd = [
                diamond_cmd, 'makedb',
                '--in', genome_file,
                '--db', temp_db_path
            ]
            
            db_result = subprocess.run(make_db_cmd, capture_output=True, text=True)
            
            if db_result.returncode != 0:
                # 清理临时文件
                self._cleanup_temp_files([temp_query_path, temp_output_path, temp_db_path])
                return False, {'found': False, 'error': '数据库创建失败'}
            
            # 运行diamond blastp命令
            blast_cmd = [
                diamond_cmd, 'blastp',
                '--query', temp_query_path,
                '--db', temp_db_path,
                '--outfmt', '6',  # 表格格式
                '--evalue', str(evalue),
                '--max-target-seqs', '1',  # 只返回最佳匹配
                '--out', temp_output_path
            ]
            
            result = subprocess.run(blast_cmd, capture_output=True, text=True)
            
            match_details = {'found': False, 'hits': []}
            
            if result.returncode == 0 and os.path.exists(temp_output_path):
                with open(temp_output_path, 'r') as f:
                    for line in f:
                        if line.strip():
                            parts = line.strip().split('\t')
                            if len(parts) >= 12:
                                hit_info = {
                                    'subject_id': parts[1],
                                    'identity': float(parts[2]),
                                    'alignment_length': int(parts[3]),
                                    'evalue': float(parts[10]),
                                    'bitscore': float(parts[11])
                                }
                                
                                # 设置更宽松的条件，因为远缘物种可能同源性较低
                                if hit_info['evalue'] <= evalue and hit_info['identity'] >= 20 and hit_info['alignment_length'] >= 30:
                                    match_details['found'] = True
                                    match_details['hits'].append(hit_info)
            
            # 清理临时文件
            self._cleanup_temp_files([temp_query_path, temp_output_path, temp_db_path])
            
            return match_details['found'], match_details
            
        except Exception as e:
            print(f"Diamond搜索出错: {e}")
            return False, {'found': False, 'error': str(e)}
    
    def _cleanup_temp_files(self, file_paths: List[str]):
        """清理临时文件"""
        for file_path in file_paths:
            try:
                if os.path.exists(file_path):
                    os.unlink(file_path)
            except Exception as e:
                print(f"清理临时文件失败 {file_path}: {e}")
    
    def analyze_genome(self, genome_file: str) -> Dict[str, Dict[str, Any]]:
        """
        分析基因组中的代谢通路（氨基酸、维生素、辅酶、核酸合成和ATP酶代谢）
        
        Args:
            genome_file: 基因组文件路径
            
        Returns:
            分析结果字典
        """
        if not os.path.exists(genome_file):
            print(f"错误: 基因组文件不存在: {genome_file}")
            return {}
        
        results = {}
        
        # 合并所有通路
        all_pathways = {**self.amino_acid_pathways, **self.vitamin_pathways, **self.nucleotide_pathways, **self.atp_pathways}
        
        # 分析每个通路
        for pathway_name, fasta_file in all_pathways.items():
            if not os.path.exists(fasta_file):
                print(f"警告: FASTA文件不存在: {fasta_file}")
                continue
            
            # 解析FASTA文件获取基因
            genes = self._parse_fasta_file(fasta_file)
            
            if not genes:
                print(f"警告: 无法从 {fasta_file} 解析基因")
                continue
            
            # 分析通路
            pathway_result = {
                'pathway': pathway_name,
                'total_genes': len(genes),
                'found_genes': 0,
                'missing_genes': [],
                'gene_details': {},
                'completeness': 0.0
            }
            
            # 记录已找到的基因，避免重复比对
            found_genes = set()
            
            # 对每个基因进行比对
            for gene_name, gene_info in genes.items():
                # 如果已经找到该基因，跳过
                if gene_name in found_genes:
                    pathway_result['gene_details'][gene_name] = {
                        'found': True,
                        'match_type': 'duplicate',
                        'description': gene_info['description']
                    }
                    pathway_result['found_genes'] += 1
                    continue
                
                # 使用Diamond进行序列比对
                found, match_details = self._run_diamond_search(
                    gene_info['sequence'], 
                    genome_file
                )
                
                pathway_result['gene_details'][gene_name] = {
                    'found': found,
                    'match_details': match_details,
                    'description': gene_info['description']
                }
                
                if found:
                    found_genes.add(gene_name)
                    pathway_result['found_genes'] += 1
                    print(f"找到 {pathway_name} 通路的 {gene_name} 基因")
                else:
                    pathway_result['missing_genes'].append(gene_name)
            
            # 计算完整度
            if pathway_result['total_genes'] > 0:
                pathway_result['completeness'] = pathway_result['found_genes'] / pathway_result['total_genes']
            
            results[pathway_name] = pathway_result
        
        return results
    
    def generate_report(self, results: Dict[str, Dict[str, Any]]) -> str:
        """
        生成分析报告
        
        Args:
            results: 分析结果
            
        Returns:
            报告文本
        """
        report = []
        report.append("代谢通路分析报告")
        report.append("=" * 50)
        report.append("")
        
        # 总体统计
        total_pathways = len(results)
        complete_pathways = sum(1 for r in results.values() if r['completeness'] == 1.0)
        
        report.append(f"总通路数: {total_pathways}")
        report.append(f"完整通路数: {complete_pathways}")
        report.append(f"不完整通路数: {total_pathways - complete_pathways}")
        report.append("")
        
        # 分类统计
        aa_pathways = {k: v for k, v in results.items() if k in self.amino_acid_pathways}
        vc_pathways = {k: v for k, v in results.items() if k in self.vitamin_pathways}
        nt_pathways = {k: v for k, v in results.items() if k in self.nucleotide_pathways}
        atp_pathways = {k: v for k, v in results.items() if k in self.atp_pathways}
        
        aa_complete = sum(1 for r in aa_pathways.values() if r['completeness'] == 1.0)
        vc_complete = sum(1 for r in vc_pathways.values() if r['completeness'] == 1.0)
        nt_complete = sum(1 for r in nt_pathways.values() if r['completeness'] == 1.0)
        atp_complete = sum(1 for r in atp_pathways.values() if r['completeness'] == 1.0)
        
        report.append("通路分类统计:")
        report.append(f"  氨基酸代谢通路: {len(aa_pathways)} 个 (完整: {aa_complete})")
        report.append(f"  维生素和辅酶代谢通路: {len(vc_pathways)} 个 (完整: {vc_complete})")
        report.append(f"  核酸合成途径: {len(nt_pathways)} 个 (完整: {nt_complete})")
        report.append(f"  ATP酶代谢途径: {len(atp_pathways)} 个 (完整: {atp_complete})")
        report.append("")
        
        # 各通路详细情况
        report.append("各通路详细情况:")
        report.append("-" * 50)
        
        for pathway, result in results.items():
            if pathway in self.amino_acid_pathways:
                pathway_type = "氨基酸"
            elif pathway in self.vitamin_pathways:
                pathway_type = "维生素/辅酶"
            elif pathway in self.nucleotide_pathways:
                pathway_type = "核酸合成"
            else:
                pathway_type = "ATP酶"
                
            report.append(f"[{pathway_type}] {pathway} 通路:")
            report.append(f"  完整度: {result['completeness']:.2%} ({result['found_genes']}/{result['total_genes']})")
            
            if result['missing_genes']:
                report.append(f"  缺失基因: {', '.join(result['missing_genes'])}")
            else:
                report.append("  所有基因均存在")
            
            report.append("")
        
        # 完整通路总结
        complete_aa_pathways = [pathway for pathway, result in aa_pathways.items() if result['completeness'] == 1.0]
        complete_vc_pathways = [pathway for pathway, result in vc_pathways.items() if result['completeness'] == 1.0]
        complete_nt_pathways = [pathway for pathway, result in nt_pathways.items() if result['completeness'] == 1.0]
        complete_atp_pathways = [pathway for pathway, result in atp_pathways.items() if result['completeness'] == 1.0]
        
        if complete_aa_pathways:
            report.append("完整的氨基酸代谢通路:")
            report.append(f"  {', '.join(complete_aa_pathways)}")
            report.append("")
        
        if complete_vc_pathways:
            report.append("完整的维生素和辅酶代谢通路:")
            report.append(f"  {', '.join(complete_vc_pathways)}")
            report.append("")
        
        if complete_nt_pathways:
            report.append("完整的核酸合成途径:")
            report.append(f"  {', '.join(complete_nt_pathways)}")
            report.append("")
        
        if complete_atp_pathways:
            report.append("完整的ATP酶代谢途径:")
            report.append(f"  {', '.join(complete_atp_pathways)}")
            report.append("")
        
        return '\n'.join(report)
    
    def generate_summary_recommendations(self, results: Dict[str, Dict[str, Any]]) -> str:
        """
        生成汇总建议
        
        Args:
            results: 分析结果
            
        Returns:
            建议文本
        """
        # 分类统计
        aa_pathways = {k: v for k, v in results.items() if k in self.amino_acid_pathways}
        vc_pathways = {k: v for k, v in results.items() if k in self.vitamin_pathways}
        nt_pathways = {k: v for k, v in results.items() if k in self.nucleotide_pathways}
        atp_pathways = {k: v for k, v in results.items() if k in self.atp_pathways}
        
        # 计算各类通路的完整度
        aa_complete = sum(1 for r in aa_pathways.values() if r['completeness'] == 1.0)
        vc_complete = sum(1 for r in vc_pathways.values() if r['completeness'] == 1.0)
        nt_complete = sum(1 for r in nt_pathways.values() if r['completeness'] == 1.0)
        atp_complete = sum(1 for r in atp_pathways.values() if r['completeness'] == 1.0)
        
        aa_avg_completeness = sum(r['completeness'] for r in aa_pathways.values()) / len(aa_pathways) if aa_pathways else 0
        vc_avg_completeness = sum(r['completeness'] for r in vc_pathways.values()) / len(vc_pathways) if vc_pathways else 0
        nt_avg_completeness = sum(r['completeness'] for r in nt_pathways.values()) / len(nt_pathways) if nt_pathways else 0
        atp_avg_completeness = sum(r['completeness'] for r in atp_pathways.values()) / len(atp_pathways) if atp_pathways else 0
        
        # 生成建议
        recommendations = []
        recommendations.append("培养性分析建议")
        recommendations.append("=" * 50)
        recommendations.append("")
        
        # 氨基酸代谢能力评估
        if aa_avg_completeness >= 0.8:
            recommendations.append("该微生物具有强大的氨基酸合成能力，可能能够自主合成大部分必需氨基酸，培养难度较低。")
        elif aa_avg_completeness >= 0.5:
            recommendations.append("该微生物具有一定的氨基酸合成能力，但可能需要补充部分必需氨基酸。")
        else:
            recommendations.append("该微生物氨基酸合成能力较弱，建议在培养基中补充多种氨基酸。")
        
        # 维生素和辅酶代谢能力评估
        if vc_avg_completeness >= 0.8:
            recommendations.append("该微生物具有强大的维生素和辅酶合成能力，可能不需要额外补充维生素。")
        elif vc_avg_completeness >= 0.5:
            recommendations.append("该微生物具有一定的维生素和辅酶合成能力，但可能需要补充部分维生素。")
        else:
            recommendations.append("该微生物维生素和辅酶合成能力较弱，建议在培养基中补充多种维生素。")
        
        # 核酸合成能力评估
        if nt_avg_completeness >= 0.8:
            recommendations.append("该微生物具有强大的核酸合成能力，能够自主合成嘌呤和嘧啶。")
        elif nt_avg_completeness >= 0.5:
            recommendations.append("该微生物具有一定的核酸合成能力，但可能需要补充部分核苷酸前体。")
        else:
            recommendations.append("该微生物核酸合成能力较弱，建议在培养基中补充核苷酸或其前体。")
        
        # ATP酶代谢能力评估
        if atp_avg_completeness >= 0.8:
            recommendations.append("该微生物具有完整的ATP酶系统，能量代谢能力强。")
        elif atp_avg_completeness >= 0.5:
            recommendations.append("该微生物ATP酶系统部分完整，能量代谢能力一般。")
        else:
            recommendations.append("该微生物ATP酶系统不完整，能量代谢可能受限，建议优化培养条件。")
        
        recommendations.append("")
        
        # 总体培养性评估
        total_avg_completeness = (aa_avg_completeness + vc_avg_completeness + nt_avg_completeness + atp_avg_completeness) / 4
        
        if total_avg_completeness >= 0.8:
            recommendations.append("总体评估：该微生物培养难度较低，可能具有广泛的营养适应性。")
        elif total_avg_completeness >= 0.5:
            recommendations.append("总体评估：该微生物培养难度中等，需要适当的营养补充。")
        else:
            recommendations.append("总体评估：该微生物培养难度较高，需要复杂的营养配方和优化的培养条件。")
        
        return '\n'.join(recommendations)
    
    def get_summary(self, results: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """
        获取分析结果摘要
        
        Args:
            results: 分析结果
            
        Returns:
            摘要字典
        """
        # 分类统计
        aa_pathways = {k: v for k, v in results.items() if k in self.amino_acid_pathways}
        vc_pathways = {k: v for k, v in results.items() if k in self.vitamin_pathways}
        nt_pathways = {k: v for k, v in results.items() if k in self.nucleotide_pathways}
        atp_pathways = {k: v for k, v in results.items() if k in self.atp_pathways}
        
        # 计算各类通路的完整度
        aa_complete = sum(1 for r in aa_pathways.values() if r['completeness'] == 1.0)
        vc_complete = sum(1 for r in vc_pathways.values() if r['completeness'] == 1.0)
        nt_complete = sum(1 for r in nt_pathways.values() if r['completeness'] == 1.0)
        atp_complete = sum(1 for r in atp_pathways.values() if r['completeness'] == 1.0)
        
        aa_avg_completeness = sum(r['completeness'] for r in aa_pathways.values()) / len(aa_pathways) if aa_pathways else 0
        vc_avg_completeness = sum(r['completeness'] for r in vc_pathways.values()) / len(vc_pathways) if vc_pathways else 0
        nt_avg_completeness = sum(r['completeness'] for r in nt_pathways.values()) / len(nt_pathways) if nt_pathways else 0
        atp_avg_completeness = sum(r['completeness'] for r in atp_pathways.values()) / len(atp_pathways) if atp_pathways else 0
        
        total_pathways = len(results)
        complete_pathways = sum(1 for r in results.values() if r['completeness'] == 1.0)
        avg_completeness = sum(r['completeness'] for r in results.values()) / total_pathways if results else 0
        
        return {
            'total_pathways': total_pathways,
            'complete_pathways': complete_pathways,
            'avg_completeness': avg_completeness,
            'amino_acid': {
                'total': len(aa_pathways),
                'complete': aa_complete,
                'avg_completeness': aa_avg_completeness
            },
            'vitamin_coenzyme': {
                'total': len(vc_pathways),
                'complete': vc_complete,
                'avg_completeness': vc_avg_completeness
            },
            'nucleotide': {
                'total': len(nt_pathways),
                'complete': nt_complete,
                'avg_completeness': nt_avg_completeness
            },
            'atp': {
                'total': len(atp_pathways),
                'complete': atp_complete,
                'avg_completeness': atp_avg_completeness
            }
        }
    
    def export_to_csv(self, results: Dict[str, Dict[str, Any]], output_file: str):
        """
        将分析结果导出为CSV文件
        
        Args:
            results: 分析结果
            output_file: 输出文件路径
        """
        rows = []
        
        for pathway_name, pathway_result in results.items():
            # 确定通路类型
            if pathway_name in self.amino_acid_pathways:
                pathway_type = "氨基酸"
            elif pathway_name in self.vitamin_pathways:
                pathway_type = "维生素/辅酶"
            elif pathway_name in self.nucleotide_pathways:
                pathway_type = "核酸合成"
            else:
                pathway_type = "ATP酶"
            
            # 添加通路总体信息
            rows.append({
                'Pathway_Type': pathway_type,
                'Pathway_Name': pathway_name,
                'Gene_Name': '总体',
                'Found': pathway_result['found_genes'],
                'Total': pathway_result['total_genes'],
                'Completeness': pathway_result['completeness'],
                'Description': f"完整度: {pathway_result['completeness']:.2%}"
            })
            
            # 添加每个基因的详细信息
            for gene_name, gene_detail in pathway_result['gene_details'].items():
                rows.append({
                    'Pathway_Type': pathway_type,
                    'Pathway_Name': pathway_name,
                    'Gene_Name': gene_name,
                    'Found': 1 if gene_detail['found'] else 0,
                    'Total': 1,
                    'Completeness': 1.0 if gene_detail['found'] else 0.0,
                    'Description': gene_detail['description']
                })
        
        # 创建DataFrame并保存为CSV
        df = pd.DataFrame(rows)
        df.to_csv(output_file, index=False)
        print(f"结果已导出到: {output_file}")
    
    def analyze_genomes_batch(self, faa_files: List[str], output_dir: str = None, parallel: bool = True, threads: int = 4) -> Dict[str, Dict[str, Any]]:
        """
        批量分析多个基因组
        
        Args:
            faa_files: 基因组文件路径列表
            output_dir: 输出目录路径（可选）
            parallel: 是否使用并行处理
            threads: 线程数（仅在并行模式下使用）
            
        Returns:
            所有基因组的分析结果字典
        """
        if not faa_files:
            print("错误: 没有提供基因组文件")
            return {}
        
        # 创建输出目录
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        start_time = time.time()
        all_results = {}
        
        print(f"开始批量分析 {len(faa_files)} 个基因组")
        
        if parallel:
            # 并行处理
            all_results = self._analyze_genomes_parallel(faa_files, output_dir, threads)
        else:
            # 串行处理
            for i, genome_file in enumerate(faa_files):
                genome_name = os.path.basename(genome_file).replace('.faa', '')
                print(f"\n正在分析 {i+1}/{len(faa_files)}: {genome_name}")
                
                # 分析基因组
                results = self.analyze_genome(genome_file)
                
                if results:
                    all_results[genome_name] = results
                    
                    # 保存单个基因组的结果
                    if output_dir:
                        self._save_genome_results(genome_name, results, output_dir)
        
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"\n批量分析完成，耗时: {elapsed_time:.2f} 秒")
        
        # 生成批量分析汇总报告
        if output_dir and all_results:
            self._generate_batch_summary(all_results, output_dir)
        
        return all_results
    
    def _analyze_genomes_parallel(self, faa_files: List[str], output_dir: str, threads: int) -> Dict[str, Dict[str, Any]]:
        """
        使用多线程并行分析多个基因组
        
        Args:
            faa_files: 基因组文件路径列表
            output_dir: 输出目录路径
            threads: 线程数
            
        Returns:
            所有基因组的分析结果字典
        """
        all_results = {}
        
        # 创建线程锁，用于保护共享资源
        lock = threading.Lock()
        
        # 定义单个基因组分析任务
        def analyze_single_genome(genome_file):
            genome_name = os.path.basename(genome_file).replace('.faa', '')
            
            # 分析基因组
            results = self.analyze_genome(genome_file)
            
            if results:
                # 使用线程锁保护共享资源
                with lock:
                    all_results[genome_name] = results
                    
                    # 保存单个基因组的结果
                    if output_dir:
                        self._save_genome_results(genome_name, results, output_dir)
            
            return genome_name, results is not None
        
        # 使用线程池执行并行分析
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            # 提交所有任务
            future_to_file = {executor.submit(analyze_single_genome, genome_file): genome_file for genome_file in faa_files}
            
            # 处理完成的任务
            completed = 0
            for future in concurrent.futures.as_completed(future_to_file):
                genome_file = future_to_file[future]
                try:
                    genome_name, success = future.result()
                    completed += 1
                    status = "成功" if success else "失败"
                    print(f"进度: {completed}/{len(faa_files)} - {genome_name} 分析{status}")
                except Exception as exc:
                    print(f"分析 {genome_file} 时发生异常: {exc}")
        
        return all_results
    
    def _save_genome_results(self, genome_name: str, results: Dict[str, Any], output_dir: str):
        """
        保存单个基因组的分析结果
        
        Args:
            genome_name: 基因组名称
            results: 分析结果
            output_dir: 输出目录路径
        """
        # 导出CSV
        csv_file = os.path.join(output_dir, f"{genome_name}_pathways.csv")
        self.export_to_csv(results, csv_file)
        
        # 生成建议报告
        recommendations = self.generate_summary_recommendations(results)
        recommendations_file = os.path.join(output_dir, f"{genome_name}_recommendations.txt")
        with open(recommendations_file, 'w', encoding='utf-8') as f:
            f.write(recommendations)
        
        # 保存完整结果为JSON
        json_file = os.path.join(output_dir, f"{genome_name}_full_results.json")
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
    
    def _generate_batch_summary(self, all_results: Dict[str, Dict[str, Any]], output_dir: str):
        """
        生成批量分析汇总报告
        
        Args:
            all_results: 所有基因组的分析结果
            output_dir: 输出目录路径
        """
        # 创建汇总表格
        summary_rows = []
        
        for genome_name, results in all_results.items():
            summary = self.get_summary(results)
            row = {
                'Genome': genome_name,
                'Total_Pathways': summary['total_pathways'],
                'Complete_Pathways': summary['complete_pathways'],
                'Avg_Completeness': f"{summary['avg_completeness']:.2%}",
                'AA_Total': summary['amino_acid']['total'],
                'AA_Complete': summary['amino_acid']['complete'],
                'AA_Avg_Completeness': f"{summary['amino_acid']['avg_completeness']:.2%}",
                'VC_Total': summary['vitamin_coenzyme']['total'],
                'VC_Complete': summary['vitamin_coenzyme']['complete'],
                'VC_Avg_Completeness': f"{summary['vitamin_coenzyme']['avg_completeness']:.2%}",
                'NT_Total': summary['nucleotide']['total'],
                'NT_Complete': summary['nucleotide']['complete'],
                'NT_Avg_Completeness': f"{summary['nucleotide']['avg_completeness']:.2%}",
                'ATP_Total': summary['atp']['total'],
                'ATP_Complete': summary['atp']['complete'],
                'ATP_Avg_Completeness': f"{summary['atp']['avg_completeness']:.2%}"
            }
            summary_rows.append(row)
        
        # 保存汇总表格为CSV
        summary_df = pd.DataFrame(summary_rows)
        summary_csv = os.path.join(output_dir, "batch_summary.csv")
        summary_df.to_csv(summary_csv, index=False)
        
        # 生成文本汇总报告
        report = []
        report.append("代谢通路批量分析汇总报告")
        report.append("=" * 50)
        report.append(f"分析基因组数量: {len(all_results)}")
        report.append("")
        
        # 按氨基酸完整度排序
        summary_rows_sorted = sorted(summary_rows, 
                                    key=lambda x: (int(x['AA_Complete']), float(x['AA_Avg_Completeness'].rstrip('%'))/100), 
                                    reverse=True)
        
        report.append("基因组氨基酸合成能力排序:")
        report.append("-" * 50)
        for i, row in enumerate(summary_rows_sorted):
            report.append(f"{i+1}. {row['Genome']}: {row['AA_Complete']}/{row['AA_Total']} 氨基酸完整通路 (平均完整度: {row['AA_Avg_Completeness']})")
        
        report.append("")
        report.append("维生素和辅酶代谢能力排序:")
        report.append("-" * 50)
        summary_rows_vc_sorted = sorted(summary_rows, 
                                      key=lambda x: (int(x['VC_Complete']), float(x['VC_Avg_Completeness'].rstrip('%'))/100), 
                                      reverse=True)
        for i, row in enumerate(summary_rows_vc_sorted):
            report.append(f"{i+1}. {row['Genome']}: {row['VC_Complete']}/{row['VC_Total']} 完整维生素/辅酶通路 (平均完整度: {row['VC_Avg_Completeness']})")
        
        report.append("")
        report.append("核酸合成能力排序:")
        report.append("-" * 50)
        summary_rows_nt_sorted = sorted(summary_rows, 
                                      key=lambda x: (int(x['NT_Complete']), float(x['NT_Avg_Completeness'].rstrip('%'))/100), 
                                      reverse=True)
        for i, row in enumerate(summary_rows_nt_sorted):
            report.append(f"{i+1}. {row['Genome']}: {row['NT_Complete']}/{row['NT_Total']} 完整核酸合成通路 (平均完整度: {row['NT_Avg_Completeness']})")
        
        report.append("")
        report.append("ATP酶代谢通路能力排序:")
        report.append("-" * 50)
        summary_rows_atp_sorted = sorted(summary_rows, 
                                       key=lambda x: (int(x['ATP_Complete']), float(x['ATP_Avg_Completeness'].rstrip('%'))/100), 
                                       reverse=True)
        for i, row in enumerate(summary_rows_atp_sorted):
            report.append(f"{i+1}. {row['Genome']}: {row['ATP_Complete']}/{row['ATP_Total']} 完整ATP酶代谢通路 (平均完整度: {row['ATP_Avg_Completeness']})")
        
        report.append("")
        report.append("详细结果:")
        report.append("-" * 50)
        for row in summary_rows_sorted:
            report.append(f"{row['Genome']}:")
            report.append(f"  总通路数: {row['Total_Pathways']} (完整: {row['Complete_Pathways']})")
            report.append(f"  氨基酸代谢: {row['AA_Total']} 个通路 (完整: {row['AA_Complete']}, 平均完整度: {row['AA_Avg_Completeness']})")
            report.append(f"  维生素/辅酶: {row['VC_Total']} 个通路 (完整: {row['VC_Complete']}, 平均完整度: {row['VC_Avg_Completeness']})")
            report.append(f"  核酸合成: {row['NT_Total']} 个通路 (完整: {row['NT_Complete']}, 平均完整度: {row['NT_Avg_Completeness']})")
            report.append(f"  ATP酶代谢: {row['ATP_Total']} 个通路 (完整: {row['ATP_Complete']}, 平均完整度: {row['ATP_Avg_Completeness']})")
            report.append("")
        
        # 保存文本报告
        summary_report = os.path.join(output_dir, "batch_summary_report.txt")
        with open(summary_report, 'w', encoding='utf-8') as f:
            f.write('\n'.join(report))
        
        print(f"汇总表格已保存到: {summary_csv}")
        print(f"汇总报告已保存到: {summary_report}")
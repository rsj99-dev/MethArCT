#!/usr/bin/env python

"""
Tome 批量测试工具

这个脚本提供了增强的批量测试功能，支持：n1. 多个 FASTA 文件的批量预测
2. 目录递归扫描
3. 详细的统计报告
4. 多种输出格式
5. 进度显示

作者: Tome 项目增强版
日期: 2024
"""

import os
import sys
import glob
import time
import argparse
import pandas as pd
from pathlib import Path

# 添加 tome 模块路径
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'tome'))
from tome import load_model, predict, get_dimer_frequency


class TomeBatchTester:
    """Tome 批量测试类"""
    
    def __init__(self, num_processes=1, verbose=True):
        self.num_processes = num_processes
        self.verbose = verbose
        self.model = None
        self.means = None
        self.stds = None
        self.features = None
        self.results = []
        
    def load_tome_model(self):
        """加载 Tome 模型"""
        if self.verbose:
            print("正在加载 Tome 模型...")
        
        try:
            self.model, self.means, self.stds, self.features = load_model()
            if self.verbose:
                print("模型加载成功！")
        except Exception as e:
            print(f"模型加载失败: {e}")
            sys.exit(1)
    
    def find_fasta_files(self, paths, recursive=False):
        """查找 FASTA 文件"""
        fasta_files = []
        fasta_extensions = ['.fasta', '.fa', '.fas', '.faa']
        
        for path in paths:
            path_obj = Path(path)
            
            if path_obj.is_file():
                # 单个文件
                if any(path.lower().endswith(ext) for ext in fasta_extensions):
                    fasta_files.append(str(path_obj.absolute()))
                else:
                    print(f"警告: {path} 不是 FASTA 文件")
            
            elif path_obj.is_dir():
                # 目录
                if recursive:
                    # 递归搜索
                    for ext in fasta_extensions:
                        pattern = f"**/*{ext}"
                        fasta_files.extend([str(p.absolute()) for p in path_obj.glob(pattern)])
                else:
                    # 只搜索当前目录
                    for ext in fasta_extensions:
                        pattern = f"*{ext}"
                        fasta_files.extend([str(p.absolute()) for p in path_obj.glob(pattern)])
            else:
                print(f"警告: {path} 不存在")
        
        # 去重并排序
        fasta_files = sorted(list(set(fasta_files)))
        
        if self.verbose:
            print(f"找到 {len(fasta_files)} 个 FASTA 文件")
        
        return fasta_files
    
    def predict_single_file(self, fasta_file):
        """预测单个文件的 OGT"""
        try:
            start_time = time.time()
            pred_ogt = predict(fasta_file, self.model, self.means, self.stds, self.features, self.num_processes)
            end_time = time.time()
            
            result = {
                'file_path': fasta_file,
                'file_name': os.path.basename(fasta_file),
                'pred_ogt': pred_ogt,
                'processing_time': round(end_time - start_time, 2),
                'status': 'success'
            }
            
            return result
            
        except Exception as e:
            result = {
                'file_path': fasta_file,
                'file_name': os.path.basename(fasta_file),
                'pred_ogt': None,
                'processing_time': None,
                'status': f'error: {str(e)}'
            }
            return result
    
    def run_batch_prediction(self, fasta_files):
        """运行批量预测"""
        if not self.model:
            self.load_tome_model()
        
        self.results = []
        total_files = len(fasta_files)
        
        if self.verbose:
            print(f"\n开始批量预测 {total_files} 个文件...")
            print("=" * 60)
        
        for i, fasta_file in enumerate(fasta_files, 1):
            if self.verbose:
                print(f"[{i}/{total_files}] 处理: {os.path.basename(fasta_file)}")
            
            result = self.predict_single_file(fasta_file)
            self.results.append(result)
            
            if self.verbose:
                if result['status'] == 'success':
                    print(f"  预测 OGT: {result['pred_ogt']}°C (耗时: {result['processing_time']}s)")
                else:
                    print(f"  错误: {result['status']}")
        
        if self.verbose:
            print("=" * 60)
            print("批量预测完成！")
    
    def generate_statistics(self):
        """生成统计信息"""
        successful_results = [r for r in self.results if r['status'] == 'success']
        failed_results = [r for r in self.results if r['status'] != 'success']
        
        if successful_results:
            ogts = [r['pred_ogt'] for r in successful_results]
            times = [r['processing_time'] for r in successful_results]
            
            stats = {
                'total_files': len(self.results),
                'successful': len(successful_results),
                'failed': len(failed_results),
                'success_rate': round(len(successful_results) / len(self.results) * 100, 2),
                'ogt_min': round(min(ogts), 2),
                'ogt_max': round(max(ogts), 2),
                'ogt_mean': round(sum(ogts) / len(ogts), 2),
                'total_time': round(sum(times), 2),
                'avg_time_per_file': round(sum(times) / len(times), 2)
            }
        else:
            stats = {
                'total_files': len(self.results),
                'successful': 0,
                'failed': len(failed_results),
                'success_rate': 0.0,
                'ogt_min': None,
                'ogt_max': None,
                'ogt_mean': None,
                'total_time': 0.0,
                'avg_time_per_file': 0.0
            }
        
        return stats
    
    def save_results(self, output_file, include_stats=True, format='tsv'):
        """保存结果到文件"""
        if not self.results:
            print("没有结果可保存")
            return
        
        # 创建 DataFrame
        df = pd.DataFrame(self.results)
        
        # 重新排列列的顺序
        columns_order = ['file_name', 'pred_ogt', 'processing_time', 'status', 'file_path']
        df = df[columns_order]
        
        # 保存主要结果
        if format.lower() == 'csv':
            df.to_csv(output_file, index=False, encoding='utf-8')
        else:  # 默认 TSV
            df.to_csv(output_file, sep='\t', index=False, encoding='utf-8')
        
        if self.verbose:
            print(f"结果已保存到: {output_file}")
        
        # 保存统计信息
        if include_stats:
            stats = self.generate_statistics()
            stats_file = output_file.replace('.tsv', '_stats.txt').replace('.csv', '_stats.txt')
            
            with open(stats_file, 'w', encoding='utf-8') as f:
                f.write("Tome 批量测试统计报告\n")
                f.write("=" * 40 + "\n")
                f.write(f"总文件数: {stats['total_files']}\n")
                f.write(f"成功处理: {stats['successful']}\n")
                f.write(f"处理失败: {stats['failed']}\n")
                f.write(f"成功率: {stats['success_rate']}%\n")
                f.write("\n")
                
                if stats['successful'] > 0:
                    f.write("OGT 预测统计:\n")
                    f.write(f"  最低 OGT: {stats['ogt_min']}°C\n")
                    f.write(f"  最高 OGT: {stats['ogt_max']}°C\n")
                    f.write(f"  平均 OGT: {stats['ogt_mean']}°C\n")
                    f.write("\n")
                    f.write("性能统计:\n")
                    f.write(f"  总处理时间: {stats['total_time']}秒\n")
                    f.write(f"  平均每文件: {stats['avg_time_per_file']}秒\n")
                
                if stats['failed'] > 0:
                    f.write("\n失败的文件:\n")
                    for result in self.results:
                        if result['status'] != 'success':
                            f.write(f"  {result['file_name']}: {result['status']}\n")
            
            if self.verbose:
                print(f"统计报告已保存到: {stats_file}")
    
    def print_summary(self):
        """打印摘要信息"""
        if not self.results:
            print("没有结果可显示")
            return
        
        stats = self.generate_statistics()
        
        print("\n" + "=" * 50)
        print("批量测试摘要")
        print("=" * 50)
        print(f"总文件数: {stats['total_files']}")
        print(f"成功处理: {stats['successful']}")
        print(f"处理失败: {stats['failed']}")
        print(f"成功率: {stats['success_rate']}%")
        
        if stats['successful'] > 0:
            print(f"\nOGT 预测范围: {stats['ogt_min']}°C - {stats['ogt_max']}°C")
            print(f"平均 OGT: {stats['ogt_mean']}°C")
            print(f"总处理时间: {stats['total_time']}秒")
            print(f"平均每文件: {stats['avg_time_per_file']}秒")
        
        print("=" * 50)


def main():
    parser = argparse.ArgumentParser(
        description='Tome 批量测试工具 - 增强版批量 OGT 预测',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  # 批量测试多个文件
  python batch_test.py file1.fasta file2.fasta file3.fasta -o results.tsv
  
  # 测试目录中的所有 FASTA 文件
  python batch_test.py test/proteomes/ -o results.tsv
  
  # 递归搜索目录
  python batch_test.py data/ -r -o results.tsv
  
  # 使用多进程加速
  python batch_test.py test/proteomes/ -p 4 -o results.tsv
        """
    )
    
    parser.add_argument('inputs', nargs='+', 
                       help='输入的 FASTA 文件或目录路径')
    parser.add_argument('-o', '--output', required=True,
                       help='输出文件路径 (.tsv 或 .csv)')
    parser.add_argument('-r', '--recursive', action='store_true',
                       help='递归搜索子目录中的 FASTA 文件')
    parser.add_argument('-p', '--processes', type=int, default=1,
                       help='并行进程数 (默认: 1)')
    parser.add_argument('--format', choices=['tsv', 'csv'], default='tsv',
                       help='输出格式 (默认: tsv)')
    parser.add_argument('--no-stats', action='store_true',
                       help='不生成统计报告文件')
    parser.add_argument('-q', '--quiet', action='store_true',
                       help='静默模式，减少输出信息')
    
    args = parser.parse_args()
    
    # 创建批量测试器
    tester = TomeBatchTester(
        num_processes=args.processes,
        verbose=not args.quiet
    )
    
    # 查找 FASTA 文件
    fasta_files = tester.find_fasta_files(args.inputs, args.recursive)
    
    if not fasta_files:
        print("错误: 没有找到 FASTA 文件")
        sys.exit(1)
    
    # 运行批量预测
    tester.run_batch_prediction(fasta_files)
    
    # 保存结果
    tester.save_results(
        args.output, 
        include_stats=not args.no_stats,
        format=args.format
    )
    
    # 打印摘要
    if not args.quiet:
        tester.print_summary()


if __name__ == '__main__':
    main()
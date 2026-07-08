"""命令行接口

用法:
    huainanzi -i protein.faa -o result.csv
"""

import argparse
import csv
import sys

from .predict import predict_from_fasta


def main():
    parser = argparse.ArgumentParser(
        prog='huainanzi',
        description='基于蛋白质组序列预测微生物生长温度范围 (T_min, T_opt, T_max)',
        epilog='示例: huainanzi -i protein.faa -o result.csv',
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        metavar='FILE',
        help='输入的蛋白质 FASTA 文件路径 (.faa 或 .fasta)',
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        metavar='FILE',
        help='输出的 CSV 结果文件路径',
    )

    args = parser.parse_args()

    try:
        print(f"正在分析: {args.input}")
        results = predict_from_fasta(args.input)

        # 写入 CSV 输出
        with open(args.output, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(['Metric', 'Temperature (°C)'])
            writer.writerow(['T_min (最低生长温度)', f"{results['T_min']:.2f}"])
            writer.writerow(['T_opt (最适生长温度)', f"{results['T_opt']:.2f}"])
            writer.writerow(['T_max (最高生长温度)', f"{results['T_max']:.2f}"])

        print(f"预测完成，结果已保存到: {args.output}")
        print(f"  最低生长温度 (T_min): {results['T_min']:.2f} °C")
        print(f"  最适生长温度 (T_opt): {results['T_opt']:.2f} °C")
        print(f"  最高生长温度 (T_max): {results['T_max']:.2f} °C")

    except FileNotFoundError as e:
        print(f"错误: 文件未找到 - {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"错误: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"错误: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()

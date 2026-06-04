# MethArCT - Methanogenic Archaeal Culturomics Toolkit 龙芯版本

## 概述

MethArCT（产甲烷古菌培养组学工具包）是一个面向产甲烷古菌宏基因组和基因组分析的综合工具箱。它集成了多种生物信息学分析功能，用于预测微生物代谢通路、耐盐性、最适生长温度和可培养性。

在线访问 MethArCT v0.1：[http://methardb.cn/tools/diamond](http://methardb.cn/tools/diamond)，用于产甲烷古菌的蛋白质功能预测。

## 主要功能

### 核心功能（必需）
- **代谢通路预测**：基于 Diamond 工具的甲烷、硫、氮代谢通路分析（21 条通路）
- **耐盐性评估**：基于基因数据库的耐盐性预测
- **可培养性评估**：基于代谢通路的培养难度评估

### 扩展功能（可选）
- **温度预测**：最适生长温度（OGT）预测 — 需要 Tome 工具
- **基因组质量评估**：基因组完整度估计和污染检测 — 需要 CheckM2 工具

## 系统要求

### 基本要求
- **Python**：>= 3.9
- **必需工具**：Diamond（从源码编译，见安装步骤）
- **目标架构**：LoongArch64（龙芯）

### 可选工具
- **Tome**：用于 OGT 预测（可选）
- **CheckM2**：用于基因组质量评估（可选）

## 快速安装

### 方式一：一键部署（推荐）

```bash
# 获取项目
git clone https://github.com/rsj99-dev/MethArCT.git
cd MethArCT

# 一键部署（自动完成所有步骤）
bash deploy.sh
```

该脚本将自动完成以下步骤：
1. 检查系统编译工具（gcc、g++、cmake、make）
2. 创建 Python 虚拟环境
3. 安装所有 Python 依赖
4. 从源码编译 DIAMOND（使用 `-DX86=OFF` 禁用 x86 SIMD，编译 generic 版本）
5. 安装 Tome（OGT 预测工具）
6. 安装 MethArCT 本体

### 方式二：手动安装

#### Step 1: 获取项目
```bash
git clone https://github.com/rsj99-dev/MethArCT.git
cd MethArCT
```

#### Step 2: 创建 Python 环境
```bash
# 使用 conda（推荐）
conda env create -f environment.yml
conda activate metharct

# 或使用 pip + 虚拟环境
python3 -m venv venv
source venv/bin/activate
```

#### Step 3: 安装依赖
```bash
# 安装 Python 依赖
pip install -r requirements.txt

# 安装项目
pip install -e .
```

#### Step 4: 编译 DIAMOND
```bash
# 使用编译脚本（自动下载源码并编译）
bash build_diamond.sh

# 或手动编译
# 从 https://github.com/bbuchfink/diamond 下载源码
cd diamond_build/diamond-2.1.9/build
cmake .. -DCMAKE_BUILD_TYPE=Release -DX86=OFF
make -j$(nproc)
```

> **注意**：`-DX86=OFF` 是关键参数。DIAMOND 默认启用 x86 SIMD 优化（SSE4.1/AVX2），
> 这些指令集在龙芯上不存在。添加此参数后 CMake 会跳过 SIMD 目标，只编译通用标量版本。

#### Step 5: 验证安装
```bash
metharct --help
```

## 使用方法

### 1. 命令行使用（推荐）

**核心分析（仅 Diamond）**：
```bash
# 代谢通路分析
metharct diamond "protein.faa" -o results/

# 综合分析（跳过可选工具）
metharct comprehensive "protein.faa" -o results/ --skip-tome --skip-checkm2
```

**完整分析（含可选工具）**：
```bash
# 全部功能分析
metharct comprehensive "protein.faa" -o results/

# 单独运行可选分析
metharct tome "protein.faa" -o results/        # 需要 Tome
metharct checkm2 "genome.fasta" -o results/    # 需要 CheckM2
```

### 2. Python API 使用
```python
from metharct import MethArCTAnalyzer

analyzer = MethArCTAnalyzer()

# 核心功能（Diamond 分析）
results = analyzer.analyze_diamond('protein_sequences.faa')

# 完整分析（需安装可选工具）
results = analyzer.comprehensive_analysis('protein_sequences.faa')

print(results)
```

### 3. 输入文件要求

- **支持格式**：FASTA 格式（.fa、.fasta、.faa）
- **文件大小**：建议不超过 50MB
- **字符编码**：UTF-8
- **含特殊字符的文件名**（如括号、空格）需用引号包裹

### 4. 输出说明

**核心功能输出（Diamond 分析）**：
- `[filename]_diamond_summary.csv` — 代谢通路摘要
- `[filename]_diamond_results.json` — 详细分析结果
- `[filename]_*_diamond.tsv` — 各通路的详细比对结果

**综合分析报告**：
- `[filename]_integrated_summary.csv` — 综合评估结果
- `[filename]_comprehensive_analysis.json` — 完整分析数据

**可选功能输出**：
- Tome 结果在 `[filename]_tome/` 目录
- CheckM2 结果在 `[filename]_checkm2/` 目录，含 `quality_report.tsv`

## 故障排查

### 安装问题

**命令行入口不可用**：
```bash
# 替代运行方式
python -m metharct.cli.main --help

# 重新安装
pip uninstall metharct -y
pip install -e .
```

**Diamond 工具未找到**：
```bash
# 检查 DIAMOND 是否编译成功
./loongarch_install/bin/diamond version

# 或添加到 PATH
export PATH="$(pwd)/loongarch_install/bin:$PATH"
```

**DIAMOND 编译失败（x86 SIMD 错误）**：
```bash
# 确保 cmake 命令包含 -DX86=OFF
cd diamond_build/diamond-2.1.9/build
rm -rf *
cmake .. -DCMAKE_BUILD_TYPE=Release -DX86=OFF
make -j$(nproc)
```

### 运行时问题

**文件名含特殊字符**：
```bash
# 错误
metharct diamond protein.faa -o results/

# 正确
metharct diamond "protein.faa" -o results/
```

## 项目结构

```
MethArCT/
├── metharct/                  # 主程序包
│   ├── cli/                   # 命令行接口
│   ├── core/                  # 核心分析模块
│   └── utils/                 # 工具函数
├── data/                      # 参考数据库
│   ├── databases/             # Diamond 数据库
│   ├── methane/               # 甲烷代谢基因
│   ├── nitrogen/              # 氮代谢基因
│   ├── sulfur/                # 硫代谢基因
│   └── salt/                  # 耐盐基因
├── Tome-1.1.0/               # Tome OGT 预测工具
├── deploy.sh                  # 龙芯一键部署脚本
├── build_diamond.sh           # DIAMOND 龙芯编译脚本
├── requirements.txt           # Python 依赖
├── environment.yml            # Conda 环境配置
├── metharct_config.yaml       # 运行配置文件
└── setup.py                   # 打包配置
```

## 许可证

本项目基于 MIT 许可证开源 — 详见 LICENSE 文件。

## 引用

如果您在研究中使用了 MethArCT，请引用：

[即将发布]

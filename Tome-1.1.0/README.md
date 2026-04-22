# Tome 1.1.0: 蛋白质最适生长温度预测工具

基于机器学习的微生物最适生长温度(OGT)预测和酶功能分析工具。

## 安装

```bash
pip install pandas numpy scikit-learn joblib requests biopython
```

## 使用方法

### 单个基因组 OGT 预测
```bash
python tome/tome.py predOGT --fasta 你的文件.fasta
```

### 批量 OGT 预测
```bash
# 多个文件
python batch_test.py file1.fasta file2.fasta file3.fasta -o results.tsv

# 整个目录
python batch_test.py 目录路径/ -o results.tsv

# 多进程加速
python batch_test.py 目录路径/ -p 4 -o results.tsv
```

### 酶功能分析
```bash
python tome/tome.py getEnzymes --fasta 酶文件.fasta
```

## 文档

- 详细说明：[快速开始指南.md](快速开始指南.md)
- English: [Quick_Start_Guide.md](Quick_Start_Guide.md)

---

**模型性能**: RMSE 2.16°C, R² 0.955

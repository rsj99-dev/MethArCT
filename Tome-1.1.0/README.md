Tome 1.1.0: Protein Optimal Growth Temperature Prediction Tool
A machine learning-based tool for predicting microbial optimal growth temperature (OGT) and analyzing enzyme functions.

Installation
bash
pip install pandas numpy scikit-learn joblib requests biopython
Usage
Single Genome OGT Prediction
bash
python tome/tome.py predOGT --fasta your_file.fasta
Batch OGT Prediction
bash
# Multiple files
python batch_test.py file1.fasta file2.fasta file3.fasta -o results.tsv

# Entire directory
python batch_test.py directory_path/ -o results.tsv

# Multi-process acceleration
python batch_test.py directory_path/ -p 4 -o results.tsv
Enzyme Function Analysis
bash
python tome/tome.py getEnzymes --fasta enzyme_file.fasta
Documentation
Detailed guide: Quick_Start_Guide.md (Note: The original Chinese filename is 快速开始指南.md, but they have an English version as Quick_Start_Guide.md. We'll keep the link as is, but the text says "Detailed guide" pointing to the Chinese file? Actually the original lists "详细说明：[快速开始指南.md]" and then "English: [Quick_Start_Guide.md]". So we should translate accordingly: "Detailed guide: [Quick_Start_Guide.md]" might be confusing because the file name is Chinese. Better to keep the Chinese filename as is for the detailed guide link, and then the English link. Alternatively, we can translate the link text but keep the file name. I'll do: "Detailed guide (Chinese): Quick_Start_Guide.md?" No, the original Chinese filename is "快速开始指南.md". In the English translation, we can write "Detailed guide: 快速开始指南.md" and then "English: Quick_Start_Guide.md". That seems accurate.

But the instruction says "全部翻译成英文" - translate everything into English. So we should produce an English version of the entire text, including the links. The links themselves are filenames; we can keep the Chinese filename as is because it's a file name. Or we could change the Chinese filename to an English equivalent? Probably keep as is for accuracy.

Let me produce final output with proper formatting.

markdown
# Tome 1.1.0: Protein Optimal Growth Temperature Prediction Tool

A machine learning-based tool for predicting microbial optimal growth temperature (OGT) and analyzing enzyme functions.

## Installation

```bash
pip install pandas numpy scikit-learn joblib requests biopython
Usage
Single Genome OGT Prediction
bash
python tome/tome.py predOGT --fasta your_file.fasta
Batch OGT Prediction
bash
# Multiple files
python batch_test.py file1.fasta file2.fasta file3.fasta -o results.tsv

# Entire directory
python batch_test.py directory_path/ -o results.tsv

# Multi-process acceleration
python batch_test.py directory_path/ -p 4 -o results.tsv
Enzyme Function Analysis
bash
python tome/tome.py getEnzymes --fasta enzyme_file.fasta

Model Performance: RMSE 2.16°C, R² 0.955

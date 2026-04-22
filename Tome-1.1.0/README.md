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

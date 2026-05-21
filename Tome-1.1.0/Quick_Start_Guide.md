# Tome 1.1.0 Quick Start Guide

## Dependencies

```
pandas
numpy
scikit-learn
joblib
requests
biopython
```

## Install Environment

```bash
pip install pandas numpy scikit-learn joblib requests biopython
```

## Single Genome/Enzyme Testing

### OGT Prediction (Genome)
```bash
python tome/tome.py predOGT --fasta your_file.fasta
```

### Enzyme Function Analysis
```bash
python tome/tome.py getEnzymes --fasta your_file.fasta
```

## Batch Testing

### Batch OGT Prediction (Multiple Genomes)
```bash
# Method 1: Specify multiple files
python batch_test.py file1.fasta file2.fasta file3.fasta -o results.tsv

# Method 2: Entire directory
python batch_test.py directory_path/ -o results.tsv

# Method 3: Use multiprocessing for acceleration
python batch_test.py directory_path/ -p 4 -o results.tsv
```

### Batch Enzyme Function Analysis
```bash
python tome/tome.py getEnzymes --fasta enzyme_file.fasta --outfile enzyme_results.tsv
```

## Output Files

- Results saved in specified `.tsv` file
- Batch testing generates additional `*_stats.txt` statistical report

---

**That's it!** 🚀
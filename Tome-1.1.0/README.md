# Tome 1.1.0: Protein Optimal Growth Temperature Prediction Tool

A machine learning-based tool for predicting microbial optimal growth temperature (OGT) and enzyme function analysis.

**This version: Tome for Windows and MethArCT**

## Download

- **Baidu Pan**: https://pan.baidu.com/s/5SK_-cZ3jWEoWUV-S4638fw#list/path=%2F

## Installation

```bash
pip install pandas numpy scikit-learn joblib requests biopython
```

## Usage

### Single Genome OGT Prediction
```bash
python tome/tome.py predOGT --fasta your_file.fasta
```

### Batch OGT Prediction
```bash
# Multiple files
python batch_test.py file1.fasta file2.fasta file3.fasta -o results.tsv

# Entire directory
python batch_test.py directory_path/ -o results.tsv

# Multi-process acceleration
python batch_test.py directory_path/ -p 4 -o results.tsv
```

### Enzyme Function Analysis
```bash
python tome/tome.py getEnzymes --fasta enzyme_file.fasta
```

## Documentation

- Quick Start Guide: [Quick_Start_Guide.md](Quick_Start_Guide.md)

---

## Changes from Original Repository

This version (Tome for Windows and MethArCT) includes the following modifications compared to the [original repository](https://github.com/EngqvistLab/Tome):

### 1. Cross-Platform Compatibility
- Replaced Unix-specific path splitting (`split('/')`) with `os.path.basename()` for cross-platform path handling
- Added `requests` library for cross-platform file downloading instead of relying on Unix-only tools
- Added `platform` module detection for OS-specific operations
- Added `shutil.which()` for cross-platform executable detection

### 2. Dependency Updates
- Replaced deprecated `pickle` with `joblib` for model serialization (compatible with scikit-learn >= 1.7.0)
- Added `StandardScaler` from sklearn for proper feature standardization
- Updated SVR model parameters for scikit-learn 1.7.0 compatibility (`gamma='scale'`)
- Added `warnings` filter to suppress FutureWarning and UserWarning

### 3. Enhanced Error Handling
- Added graceful fallback when BLAST+ is not installed (enzyme analysis continues without sequence comparison)
- Added proper error handling for model loading with automatic retraining fallback
- Added cross-platform temporary file cleanup with error handling
- Added proper file encoding (`utf-8`) for output files

### 4. New Features
- Added `batch_test.py` - Enhanced batch testing tool with:
  - Multiple file support
  - Directory recursive scanning
  - Detailed statistical reports
  - Progress display
  - Multiple output formats (TSV/CSV)
  - Multi-process support

### 5. Code Improvements
- Fixed file closing issues in output writing
- Improved subprocess handling with `capture_output=True`
- Added proper exception handling throughout the codebase

---

**Model Performance**: RMSE 2.16°C, R² 0.955

## License

This project is licensed under the GNU General Public License v3 or later (GPLv3+).
See the [LICENSE.txt](LICENSE.txt) file for details.

## Author

- **Gang Li** - gangl@chalmers.se

## Repository

Original source: https://github.com/EngqvistLab/Tome

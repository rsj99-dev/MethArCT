# MethArCT - Methanogenic Archaeal Culturomics Toolkit

## MethArCT

## Introduction

MethArCT (Methanogenic Archaeal Culturomics Toolkit) is a comprehensive toolkit specifically designed for the analysis of methanogenic archaeal culturomics. This tool integrates multiple bioinformatics analysis functionalities, enabling the prediction of microbial metabolic pathways, salt tolerance, optimal growth temperature, and culturability.

## Main Features

### Core Features (Required)
- **Metabolic Pathway Prediction**: Analysis of methane, sulfur, and nitrogen metabolic pathways (21 pathways) - Based on the Diamond tool
- **Salt Tolerance Assessment**: Prediction of salt tolerance based on gene databases
- **Culturability Assessment**: Evaluation of cultivation difficulty based on metabolic pathways

### Extended Features (Optional)
- **Temperature Prediction**: Optimal growth temperature (OGT) prediction - Requires installation of the Tome tool
- **Genome Quality Assessment**: Assessment of genome completeness and contamination based on CheckM2 - Requires installation of the CheckM2 tool

## System Requirements

### Basic Requirements
- **Python**: 3.8+
- **Required Tool**: Diamond
- **Operating System**: Windows, Linux, macOS

### Optional Tools
- **Tome**: For optimal growth temperature prediction (optional)
- **CheckM2**: For genome quality assessment (optional)

## Quick Installation Guide

### Step 1: Obtain the Project
```bash
# Clone the project locally
git clone https://github.com/your-username/MethArCT.git
cd MethArCT
```

### Step 2: Create a Python Environment
```bash
# Create an environment using conda (recommended)
conda env create -f environment.yml
conda activate metharct

# Or create a virtual environment using pip
python -m venv metharct_env
# Windows activation:
metharct_env\Scripts\activate
# Linux/macOS activation:
source metharct_env/bin/activate
```

### Step 3: Install Core Dependencies
```bash
# Install Python dependency packages
pip install -r requirements.txt

# Install the project as an executable command
pip install -e .

# Install the required Diamond tool
conda install -c bioconda diamond
```

### Step 4: Verify Installation
```bash
# Test core functionality
python quick_test.py

# Check the command-line tool
metharct --help
```

## Optional Tool Installation (Extended Features)

If you need to use temperature prediction or genome quality assessment features, you can install the following optional tools:

### Tome Tool Installation (Temperature Prediction)
```bash
# Install in the conda environment
conda install -c bioconda tome

# Verify installation
tome --help
```

### CheckM2 Tool Installation (Genome Quality Assessment)
```bash
# Install in the conda environment
conda install -c bioconda checkm2

# Verify installation
checkm2 --help
```

### Special Instructions for Windows Users

**Important Note**: For Windows users, if you encounter compatibility issues when installing Tome and CheckM2, it is recommended to use the WSL (Windows Subsystem for Linux) environment:

1. **Install WSL2**:
   - Run in PowerShell: `wsl --install`
   - Restart the computer

2. **Install Tools in WSL**:
   ```bash
   # Enter the WSL environment
   wsl
   
   # Install conda (if not installed)
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   
   # Create an environment and install tools
   conda create -n bioinfo python=3.8
   conda activate bioinfo
   conda install -c bioconda diamond tome checkm2
   ```

3. **Configure the Project to Use WSL**:
   - Refer to the `example_wsl_usage.py` file in the project
   - Or specify the WSL path in the configuration file

## Detailed Usage Guide

### 1. Basic Usage Methods

#### 1.1 Command-Line Usage (Recommended)

**Core Analysis (Using Diamond Only)**:
```bash
# Metabolic pathway analysis - This is a core feature and requires no additional tools
metharct diamond "protein(1).faa" -o results/

# Comprehensive analysis (skipping optional tools)
metharct comprehensive "protein(1).faa" -o results/ --skip-tome --skip-checkm2
```

**Complete Analysis (Requires Installation of Optional Tools)**:
```bash
# Comprehensive analysis including all features
metharct comprehensive "protein(1).faa" -o results/

# Run optional analyses separately
metharct tome "protein(1).faa" -o results/        # Requires Tome
metharct checkm2 "genome.fasta" -o results/      # Requires CheckM2
```

#### 1.2 Python Script Usage
```python
from metharct import MethArCTAnalyzer

# Create an analyzer instance
analyzer = MethArCTAnalyzer()

# Use only core functionality (Diamond analysis)
results = analyzer.analyze_diamond('protein_sequences.faa')

# Complete analysis (if optional tools are installed)
results = analyzer.comprehensive_analysis('protein_sequences.faa')

print(results)
```

### 2. Input File Requirements

#### 2.1 File Format
- **Supported Formats**: FASTA format (.fa, .fasta, .faa)
- **File Size**: Recommended not to exceed 50MB
- **Character Encoding**: UTF-8

#### 2.2 Sequence Types
- **Protein Sequences**: For Diamond metabolic pathway analysis and Tome temperature prediction
- **Genome Sequences**: For CheckM2 genome quality assessment

#### 2.3 Filename Considerations
- If the filename contains special characters (such as parentheses, spaces), enclose it in quotes
- Example: `"protein(1).faa"` instead of `protein(1).faa`

### 3. Testing the Installation

#### 3.1 Quick Test
```bash
# Test core functionality
python quick_test.py

# Check the command-line tool
metharct --help
```

#### 3.2 Testing with Example Files
```bash
# Test using the example files provided with the project
metharct diamond "protein(1).faa" -o test_results/

# If optional tools are installed, run a complete test
metharct comprehensive "protein(1).faa" -o test_results/
```

### 4. Output Results Description

#### 4.1 Core Function Output (Diamond Analysis)
After the analysis is completed, the following will be generated in the specified output directory:

- **Metabolic Pathway Analysis Results**:
  - `[filename]_diamond_summary.csv` - Metabolic pathway summary table
  - `[filename]_diamond_results.json` - Detailed analysis results
  - `[filename]_*_diamond.tsv` - Detailed alignment results for each pathway

- **Comprehensive Analysis Report**:
  - `[filename]_integrated_summary.csv` - Comprehensive evaluation results
  - `[filename]_comprehensive_analysis.json` - Complete analysis data
  - `[filename]_analysis_summary.txt` - Text format summary

#### 4.2 Optional Function Output

**Tome Analysis Results** (if Tome is installed):
- `[filename]_tome/` - Tome analysis output directory
- Temperature prediction-related files

**CheckM2 Analysis Results** (if CheckM2 is installed):
- `[filename]_checkm2/` - CheckM2 analysis output directory
- `quality_report.tsv` - Genome quality report

### 5. Common Issue Resolution

#### 5.1 Installation-Related Issues

**Command-Line Entry Point Issues**:
```bash
# If the metharct command is unavailable but the package is installed, use the following method
python -m metharct.cli.main --help

# Reinstall the project (fix command-line entry points)
pip uninstall metharct -y
pip install -e .

# Check if the environment is correctly activated
conda activate metharct
```

**Diamond Tool Missing**:
```bash
# Ensure you are in the correct conda environment
conda activate metharct

# Install Diamond (core required tool)
conda install -c bioconda diamond

# Verify installation
diamond --help
```

#### 5.2 Runtime Issues

**Filename Contains Special Characters**:
```bash
# Incorrect way
metharct diamond protein(1).faa -o results/

# Correct way
metharct diamond "protein(1).faa" -o results/
```

**Permission Issues (Windows)**:
- Ensure running in a directory with write permissions
- Avoid creating output files in system directories

#### 5.3 Optional Tool Issues

**Tome or CheckM2 Installation Failure**:
- These are optional tools and do not affect core functionality
- You can skip them using the `--skip-tome` and `--skip-checkm2` parameters
- Windows users are advised to install these tools in a WSL environment

**WSL Configuration (For Optional Tools Only)**:
```bash
# Install WSL2
wsl --install

# Install conda and tools in WSL
wsl
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda install -c bioconda tome checkm2
```

## Quick Start Guide

### Step 1: Verify Installation
```bash
# Check if core tools are successfully installed
metharct --help

# If the above command fails, use the following method (exactly the same functionality)
python -m metharct.cli.main --help

# Check the Diamond tool
diamond --help
```

### Step 2: Prepare Data
- Prepare protein sequence files in FASTA format
- Ensure filenames do not contain special characters, or enclose them in quotes

### Step 3: Run Analysis
```bash
# Core functionality analysis
metharct diamond "your_protein.faa" -o results/

# If the metharct command is unavailable, use the following method (exactly the same functionality) (Recommended)
python -m metharct.cli.main diamond "your_protein.faa" -o results/

# Complete analysis (if optional tools are installed)
metharct comprehensive "your_protein.faa" -o results/
# Or
python -m metharct.cli.main comprehensive "your_protein.faa" -o results/
```

### Step 4: View Results
- Open the `*_integrated_summary.csv` file in the output directory
- Check `*_analysis_summary.txt` for an analysis summary

## Learning Resources

- 📝 **Python Usage Examples**: [example_usage.py](example_usage.py)
- 🖥️ **WSL Usage Examples**: [example_wsl_usage.py](example_wsl_usage.py)
- 🐛 **Issue Reporting**: GitHub Issues

## Beginner Recommendations

1. **Start with Core Features**: Begin with Diamond analysis to familiarize yourself with tool operation
2. **Understand Output Results**: Focus on metabolic complexity and cultivation potential assessment
3. **Expand Gradually**: Install and use optional tools only when needed
4. **Encountering Issues**: Check the common issue resolution section, or submit a GitHub Issue

## Database Description

MethArCT includes the following metabolic pathway databases:

### Methane Metabolic Pathways (8 types)
- Various methanogenesis pathways including formate methane, acetate methane, methanol methane, etc.

### Sulfur Metabolic Pathways (6 types)
- Sulfate reduction, sulfur oxidation, sulfide oxidation, etc.

### Nitrogen Metabolic Pathways (4 types)
- Nitrification, denitrification, ammonia oxidation, etc.

### Other Pathways (3 types)
- Genes related to salt tolerance and culturability assessment

## Technical Support

### Getting Help
- 📖 **Detailed Documentation**: [TESTING.md](TESTING.md)
- 📝 **Usage Examples**: [example_usage.py](example_usage.py)
- 🖥️ **WSL Configuration**: [example_wsl_usage.py](example_wsl_usage.py)
- 🐛 **Issue Reporting**: GitHub Issues

### Common Questions
1. **Installation Issues**: Check the "Common Issue Resolution" section above
2. **Usage Issues**: Refer to the "Quick Start Guide"
3. **Result Interpretation**: Check the "Output Results Description"

## Project Information

### License
MIT License - See the [LICENSE](LICENSE) file for details

### Version Information
- **Current Version**: 1.0.0
- **Python Requirement**: 3.8+
- **Core Dependencies**: Diamond

### Citation
If you use MethArCT in your research, please cite:
```
[Paper information to be published]
```

### Contact
- **GitHub Issues**: Report issues and feature requests
- **Project Homepage**: [[GitHub Repository](https://github.com/rsj99-dev/MethArCT)]

---

**Thank you for using MethArCT!** 🧬

If you find this tool useful, please give us a ⭐️!

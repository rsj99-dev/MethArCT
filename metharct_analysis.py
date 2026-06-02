#!/usr/bin/env python3
import sys
import os
from pathlib import Path

project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from metharct.core.diamond_analyzer import DiamondAnalyzer
from metharct.core.tome_analyzer import TomeAnalyzer
from metharct.utils.config import Config

def main():
    if len(sys.argv) < 2:
        print("Usage: metharct_analysis.py <protein.fasta>")
        return 1
    
    input_file = sys.argv[1]
    if not os.path.exists(input_file):
        print(f"File not found: {input_file}")
        return 1
    
    config = Config()
    da = DiamondAnalyzer(config)
    results = da.analyze_sequence(input_file, output_prefix="diamond_out")
    print("Diamond analysis done.")
    return 0

if __name__ == "__main__":
    sys.exit(main())

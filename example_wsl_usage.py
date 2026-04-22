#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MethArCT WSL Usage Example

This example shows how to configure MethArCT to use tools installed in WSL environment:
- Required tool: Diamond (core functionality)
- Optional tools: CheckM2, Tome (extended functionality)
"""

import os
from metharct.utils.config_wsl import WSLConfig
from metharct.core.metharct_analyzer import MethArCTAnalyzer

def main():
    """
    WSL environment usage example
    """
    print("=== MethArCT WSL Environment Usage Example ===")

    print("\n1. Creating WSL configuration...")
    config = WSLConfig()

    print("\n2. Checking WSL environment...")
    if not config.check_wsl_available():
        print("Error: WSL environment is not available")
        return

    print("WSL environment is available")

    print("\n3. Checking tool availability...")

    print("\nCore tool (required):")
    if config.check_tool_in_wsl('diamond'):
        print("✓ Diamond is available in WSL (core functionality)")
    else:
        print("✗ Diamond is not available in WSL (this is a required core functionality!)")
        print("  Installation: wsl conda install -c bioconda diamond")

    print("\nOptional tools (extended functionality):")
    if config.check_tool_in_wsl('checkm2'):
        print("✓ CheckM2 is available in WSL (genome quality assessment)")
    else:
        print("✗ CheckM2 is not available in WSL (optional functionality)")
        print("  Installation: wsl conda install -c bioconda checkm2")

    if config.check_tool_in_wsl('tome'):
        print("✓ Tome is available in WSL (temperature prediction)")
    else:
        print("✗ Tome is not available in WSL (optional functionality)")
        print("  Installation: wsl conda install -c bioconda tome")

    print("\n4. Configuring WSL usage...")

    config.enable_wsl_mode()

    # config.set_wsl_paths({
    #     'diamond': '/home/username/miniconda3/bin/diamond',
    #     'checkm2': '/home/username/miniconda3/bin/checkm2',
    #     'tome': '/home/username/miniconda3/bin/tome'
    # })

    print("WSL mode enabled")

    print("\n5. Creating analyzer instance...")
    try:
        analyzer = MethArCTAnalyzer(config=config)
        print("✓ MethArCT analyzer created successfully (WSL mode)")

        test_fasta = "test_data/test_sequences.fasta"
        if os.path.exists(test_fasta):
            print(f"\n6. Running test analysis: {test_fasta}")
            results = analyzer.analyze_sequences(test_fasta)
            print(f"Analysis completed, results saved to: {results['output_dir']}")
        else:
            print(f"\n6. Test data file not found: {test_fasta}")
            print("Please prepare a FASTA format sequence file for testing")

    except Exception as e:
        print(f"✗ Failed to create analyzer: {e}")
        print("Please check WSL environment and tool installation")

    print("\n=== WSL Usage Example Complete ===")

def show_wsl_installation_guide():
    """
    Display WSL installation guide
    """
    print("=== WSL Environment Installation Guide ===")
    print()
    print("1. Install WSL:")
    print("   wsl --install")
    print()
    print("2. Install Miniconda in WSL:")
    print("   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh")
    print("   bash Miniconda3-latest-Linux-x86_64.sh")
    print()
    print("3. Install bioinformatics tools in WSL:")
    print("   # Required tool (core functionality)")
    print("   conda install -c bioconda diamond")
    print()
    print("   # Optional tools (extended functionality)")
    print("   conda install -c bioconda checkm2 tome")
    print()
    print("4. Verify installation:")
    print("   # Verify required tool")
    print("   wsl diamond version")
    print()
    print("   # Verify optional tools")
    print("   wsl checkm2 --version")
    print("   wsl tome --version")
    print()
    print("5. Run MethArCT WSL example:")
    print("   python example_wsl_usage.py")
    print()
    print("Note: Only Diamond is required, CheckM2 and Tome are optional tools")
    print()

if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "--install-guide":
        show_wsl_installation_guide()
    else:
        main()

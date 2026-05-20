#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MethArCT Installation Script
Automated installation and setup for MethArCT package
"""

import os
import sys
import subprocess
import platform
from pathlib import Path

def check_python_version():
    """Check if Python version is compatible."""
    if sys.version_info < (3, 8):
        print("❌ Python 3.8 or higher is required")
        print(f"Current version: {sys.version}")
        return False
    print("✅ Python version check passed")
    return True

def check_pip():
    """Check if pip is available."""
    try:
        subprocess.run([sys.executable, "-m", "pip", "--version"], 
                      capture_output=True, check=True)
        print("✅ pip is available")
        return True
    except subprocess.CalledProcessError:
        print("❌ pip is not available")
        return False

def install_package():
    """Install MethArCT package."""
    print("\n📦 Installing MethArCT package...")
    
    try:
        # Install the package in development mode
        result = subprocess.run([
            sys.executable, "-m", "pip", "install", "-e", "."
        ], capture_output=True, text=True)
        
        if result.returncode == 0:
            print("✅ MethArCT package installed successfully")
            return True
        else:
            print(f"❌ Installation failed: {result.stderr}")
            return False
            
    except Exception as e:
        print(f"❌ Installation error: {e}")
        return False

def install_optional_dependencies():
    """Install optional dependencies for full functionality."""
    print("\n🔧 Installing optional dependencies...")
    
    try:
        result = subprocess.run([
            sys.executable, "-m", "pip", "install", ".[full]"
        ], capture_output=True, text=True)
        
        if result.returncode == 0:
            print("✅ Optional dependencies installed")
        else:
            print("⚠️  Optional dependencies installation had issues")
            print(f"   Error: {result.stderr}")
            
    except Exception as e:
        print(f"⚠️  Optional dependencies installation error: {e}")

def check_installation():
    """Verify that MethArCT is properly installed."""
    print("\n🔍 Verifying installation...")
    
    try:
        # Test import
        import metharct
        print("✅ MethArCT package imported successfully")
        
        # Test version
        version = metharct.get_version()
        print(f"✅ MethArCT version: {version}")
        
        # Test configuration
        from metharct.utils.config import Config
        config = Config()
        print("✅ Configuration loaded successfully")
        
        return True
        
    except ImportError as e:
        print(f"❌ Import error: {e}")
        return False
    except Exception as e:
        print(f"❌ Verification error: {e}")
        return False

def create_example_script():
    """Create an example usage script."""
    example_script = '''#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""MethArCT Example Usage Script"""

from metharct.core.pathway_predictor import PathwayPredictor

def main():
    print("MethArCT Example Analysis")
    print("=" * 50)
    
    # Initialize predictor
    predictor = PathwayPredictor()
    
    print("Provide your own protein file for analysis")
    print("Example usage:")
    print("  predictor = PathwayPredictor()")
    print("  results = predictor.comprehensive_analysis('your_protein.faa', 'results/')")
        
    print("Example script completed")

if __name__ == "__main__":
    main()
'''
    
    with open("example_analysis.py", "w", encoding="utf-8") as f:
        f.write(example_script)
    
    print("✅ Created example_analysis.py script")

def main():
    """Main installation function."""
    print("🧬 MethArCT Installation Script")
    print("=" * 50)
    
    # Check prerequisites
    if not check_python_version():
        sys.exit(1)
    
    if not check_pip():
        sys.exit(1)
    
    # Install package
    if not install_package():
        sys.exit(1)
    
    # Install optional dependencies
    install_optional_dependencies()
    
    # Verify installation
    if not check_installation():
        print("\n❌ Installation verification failed")
        sys.exit(1)
    
    # Create example script
    create_example_script()
    
    print("\n" + "=" * 50)
    print("🎉 MethArCT Installation Completed Successfully!")
    print("\n📚 Next Steps:")
    print("   1. Run: python example_analysis.py")
    print("   2. Or use: metharct-analyze your_protein.faa -o results")
    print("   3. See README.md for detailed usage instructions")
    print("\n🔧 For help:")
    print("   • metharct --help")
    print("   • metharct-analyze --help")
    print("   • https://github.com/MethArCT/metharct")

if __name__ == "__main__":
    main()
# -*- coding: utf-8 -*-
"""
WSL environment configuration file
Used to use Linux tools (Diamond, Tome, CheckM2) via WSL on Windows
"""

import os
import subprocess
from .config import Config

class WSLConfig(Config):
    """Configuration class for WSL environment"""
    
    def __init__(self, config_file: str = None):
        super().__init__(config_file)
        self._setup_wsl_paths()
    
    def _setup_wsl_paths(self):
        """Set WSL tool paths"""
        # Check if WSL is available
        if not self._check_wsl_available():
            raise RuntimeError("WSL not available, please install WSL first")
        
        # Update all tool paths to WSL versions
        self.config['tools']['diamond']['path'] = 'wsl diamond'
        self.config['tools']['diamond']['use_wsl'] = True
        
        self.config['tools']['tome']['path'] = 'wsl tome'
        self.config['tools']['tome']['use_wsl'] = True
        
        self.config['tools']['checkm2']['path'] = 'wsl checkm2'
        self.config['tools']['checkm2']['use_wsl'] = True
        
        self.logger.info("All tools configured to use WSL environment")
    
    def _check_wsl_available(self) -> bool:
        """Check if WSL is available"""
        try:
            result = subprocess.run(
                ['wsl', '--version'],
                capture_output=True,
                text=True,
                timeout=10
            )
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return False
    
    def get_wsl_path(self, windows_path: str) -> str:
        """Convert Windows path to WSL path"""
        # Convert Windows path to WSL format
        # Example: D:\data -> /mnt/d/data
        if ':' in windows_path:
            drive, path = windows_path.split(':', 1)
            wsl_path = f"/mnt/{drive.lower()}{path.replace(chr(92), '/')}"
            return wsl_path
        return windows_path
    
    def check_wsl_tools(self) -> dict:
        """Check tool availability in WSL"""
        tools_status = {}
        
        # Check Diamond
        try:
            result = subprocess.run(
                ['wsl', 'diamond', 'version'],
                capture_output=True,
                text=True,
                timeout=30
            )
            tools_status['diamond'] = {
                'available': result.returncode == 0,
                'version': result.stdout.strip() if result.returncode == 0 else None,
                'error': result.stderr if result.returncode != 0 else None
            }
        except Exception as e:
            tools_status['diamond'] = {
                'available': False,
                'version': None,
                'error': str(e)
            }
        
        # Check Tome
        try:
            result = subprocess.run(
                ['wsl', 'tome', '--version'],
                capture_output=True,
                text=True,
                timeout=30
            )
            tools_status['tome'] = {
                'available': result.returncode == 0,
                'version': result.stdout.strip() if result.returncode == 0 else None,
                'error': result.stderr if result.returncode != 0 else None
            }
        except Exception as e:
            tools_status['tome'] = {
                'available': False,
                'version': None,
                'error': str(e)
            }
        
        # Check CheckM2
        try:
            result = subprocess.run(
                ['wsl', 'checkm2', '--version'],
                capture_output=True,
                text=True,
                timeout=30
            )
            tools_status['checkm2'] = {
                'available': result.returncode == 0,
                'version': result.stdout.strip() if result.returncode == 0 else None,
                'error': result.stderr if result.returncode != 0 else None
            }
        except Exception as e:
            tools_status['checkm2'] = {
                'available': False,
                'version': None,
                'error': str(e)
            }
        
        return tools_status
    
    def setup_wsl_environment(self):
        """Set WSL environment variables and paths"""
        # Can add special settings for WSL environment here
        # For example, set database paths to WSL format
        pass
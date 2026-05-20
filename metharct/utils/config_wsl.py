# -*- coding: utf-8 -*-
"""
WSL环境配置文件
用于在Windows上通过WSL使用Linux工具（Diamond、Tome、CheckM2）
"""

import os
import subprocess
from .config import Config

class WSLConfig(Config):
    """WSL环境专用配置类"""
    
    def __init__(self, config_file: str = None):
        super().__init__(config_file)
        self._setup_wsl_paths()
    
    def _setup_wsl_paths(self):
        """设置WSL工具路径"""
        # 检查WSL是否可用
        if not self._check_wsl_available():
            raise RuntimeError("WSL不可用，请先安装WSL")
        
        # 更新所有工具路径为WSL版本
        self.config['tools']['diamond']['path'] = 'wsl diamond'
        self.config['tools']['diamond']['use_wsl'] = True
        
        self.config['tools']['tome']['path'] = 'wsl tome'
        self.config['tools']['tome']['use_wsl'] = True
        
        self.config['tools']['checkm2']['path'] = 'wsl checkm2'
        self.config['tools']['checkm2']['use_wsl'] = True
        
        self.logger.info("已配置所有工具使用WSL环境")
    
    def _check_wsl_available(self) -> bool:
        """检查WSL是否可用"""
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
        """将Windows路径转换为WSL路径"""
        # 将Windows路径转换为WSL格式
        # 例如: D:\data -> /mnt/d/data
        if ':' in windows_path:
            drive, path = windows_path.split(':', 1)
            wsl_path = f"/mnt/{drive.lower()}{path.replace(chr(92), '/')}"
            return wsl_path
        return windows_path
    
    def check_wsl_tools(self) -> dict:
        """检查WSL中工具的可用性"""
        tools_status = {}
        
        # 检查Diamond
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
        
        # 检查Tome
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
        
        # 检查CheckM2
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
        """设置WSL环境变量和路径"""
        # 这里可以添加WSL环境的特殊设置
        # 比如设置数据库路径为WSL格式
        pass
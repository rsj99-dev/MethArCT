#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Configuration management for MethArCT

Handles configuration loading, database paths, and tool settings.
"""

import os
import yaml
import json
from pathlib import Path
from typing import Dict, Any, Optional

class Config:
    """Configuration manager for MethArCT"""
    
    def __init__(self, config_file: Optional[str] = None):
        self.config_file = config_file
        self.config = self._load_default_config()
        
        if config_file and os.path.exists(config_file):
            self._load_config_file(config_file)
    
    def _load_default_config(self) -> Dict[str, Any]:
        """Load default configuration from bundled YAML file."""
        default_yaml = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'data',
            'default_config.yaml',
        )
        if os.path.exists(default_yaml):
            with open(default_yaml, 'r', encoding='utf-8') as f:
                return yaml.safe_load(f)
        # Fallback: minimal config if YAML is missing
        return {
            'tools': {'diamond': {'path': 'diamond', 'threads': 4, 'evalue': 1e-5}},
            'databases': {'base_dir': 'data/databases'},
            'output': {'base_dir': 'results'},
            'logging': {'level': 'INFO'},
        }
    
    def _load_config_file(self, config_file: str):
        """Load configuration from file"""
        try:
            with open(config_file, 'r', encoding='utf-8') as f:
                if config_file.endswith('.yaml') or config_file.endswith('.yml'):
                    user_config = yaml.safe_load(f)
                elif config_file.endswith('.json'):
                    user_config = json.load(f)
                else:
                    raise ValueError(f"Unsupported config file format: {config_file}")
            
            # Merge user config with default config
            self._merge_config(self.config, user_config)
        except Exception as e:
            print(f"Warning: Failed to load config file {config_file}: {e}")
    
    def _merge_config(self, default: Dict, user: Dict):
        """Recursively merge user config into default config"""
        for key, value in user.items():
            if key in default and isinstance(default[key], dict) and isinstance(value, dict):
                self._merge_config(default[key], value)
            else:
                default[key] = value
    
    def get(self, key: str, default: Any = None) -> Any:
        """Get configuration value using dot notation"""
        keys = key.split('.')
        value = self.config
        
        for k in keys:
            if isinstance(value, dict) and k in value:
                value = value[k]
            else:
                return default
        
        return value
    
    def set(self, key: str, value: Any):
        """Set configuration value using dot notation"""
        keys = key.split('.')
        config = self.config
        
        for k in keys[:-1]:
            if k not in config:
                config[k] = {}
            config = config[k]
        
        config[keys[-1]] = value
    
    def get_database_path(self, db_type: str, db_name: str) -> str:
        """Get full path to database file"""
        base_dir = self.get('databases.base_dir', 'data/databases')
        
        if db_type == 'methane':
            rel_path = self.get(f'databases.methane_pathways.{db_name}')
        elif db_type == 'sulfur':
            rel_path = self.get(f'databases.sulfur_pathways.{db_name}')
        elif db_type == 'nitrogen':
            rel_path = self.get(f'databases.nitrogen_pathways.{db_name}')
        elif db_type == 'cultivation':
            rel_path = self.get('databases.cultivation')
        else:
            raise ValueError(f"Unknown database type: {db_type}")
        
        if not rel_path:
            raise ValueError(f"Database {db_name} not found in {db_type} databases")
        
        return os.path.join(base_dir, rel_path)
    
    def get_all_database_paths(self) -> Dict[str, str]:
        """Get all database paths"""
        paths = {}
        base_dir = self.get('databases.base_dir', 'data/databases')
        
        # Methane pathways
        for name, path in self.get('databases.methane_pathways', {}).items():
            paths[name] = os.path.join(base_dir, path)
        
        # Sulfur pathways
        for name, path in self.get('databases.sulfur_pathways', {}).items():
            paths[name] = os.path.join(base_dir, path)
        
        # Nitrogen pathways
        for name, path in self.get('databases.nitrogen_pathways', {}).items():
            paths[name] = os.path.join(base_dir, path)
        
        # Cultivation
        cult_path = self.get('databases.cultivation')
        if cult_path:
            paths['CULTIVATION'] = os.path.join(base_dir, cult_path)
        
        return paths
    
    def get_cultivation_pathway_dir(self, pathway_type: str) -> str:
        """Get directory path for cultivation pathway type"""
        base_dir = self.get('databases.base_dir', 'data/databases')
        rel_path = self.get(f'databases.cultivation_pathways.{pathway_type}')
        
        if not rel_path:
            raise ValueError(f"Cultivation pathway type {pathway_type} not found")
        
        return os.path.join(base_dir, rel_path)
    
    def save_config(self, output_file: str):
        """Save current configuration to file"""
        with open(output_file, 'w', encoding='utf-8') as f:
            if output_file.endswith('.yaml') or output_file.endswith('.yml'):
                yaml.dump(self.config, f, default_flow_style=False, allow_unicode=True)
            elif output_file.endswith('.json'):
                json.dump(self.config, f, indent=2, ensure_ascii=False)
            else:
                raise ValueError(f"Unsupported output format: {output_file}")
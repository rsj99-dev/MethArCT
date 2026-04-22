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
        """Load default configuration"""
        return {
            'tools': {
                'diamond': {
                    'path': 'diamond',
                    'wsl_path': 'wsl diamond',
                    'use_wsl': False,
                    'threads': 4,
                    'evalue': 1e-5,
                    'max_target_seqs': 1,
                    'identity_threshold': 30.0
                },
                'tome': {
                    'path': 'tome',
                    'wsl_path': 'wsl tome',
                    'use_wsl': False,
                    'threads': 4,
                    'install_dir': os.path.join(os.path.dirname(__file__), '..', '..', 'Tome-1.1.0')
                },
                'checkm2': {
                    'path': 'checkm2',
                    'wsl_path': 'wsl checkm2',
                    'use_wsl': False,
                    'threads': 4,
                    'database_path': os.path.join(os.path.dirname(__file__), '..', '..', 'checkm2_db', 'uniref100.KO.1.dmnd')
                }
            },
            'databases': {
                'base_dir': 'data/databases',
                'methane_pathways': {
                    'CO2-CH4': 'methane/processed_CO2-CH4.fasta',
                    'METHYLAMINE-CH4': 'methane/processed_METHYLAMINE-CH4.fasta',
                    'METHANOL-CH4': 'methane/processed_METHANOL-CH4.fasta',
                    'METHANETHIOL-CH4': 'methane/processed_METHANETHIOL-CH4.fasta',
                    'ACETIC_ACID-CH4': 'methane/processed_ACETIC_ACID-CH4.fasta',
                    'C16-CH4': 'methane/processed_C16-CH4.fasta',
                    'CO-CH4': 'methane/processed_CO-CH4.fasta',
                    'FORMIC_ACID-CH4': 'methane/processed_FORMIC_ACID-CH4.fasta',
                    'METHOXY-CH4': 'methane/processed_METHOXY-CH4.fasta',
                    'FATTY_ACID-CH4': 'methane/processed_FATTY_ACID-CH4.fasta',
                    'DIMETHYLAMINE-CH4': 'methane/processed_DIMETHYLAMINE-CH4.fasta',
                    'TRIMETHYLAMINE-CH4': 'methane/processed_TRIMETHYLAMINE-CH4.fasta'
                },
                'sulfur_pathways': {
                    'ASR': 'sulfur/Assimilatory sulfate reduction.fasta',
                    'SO': 'sulfur/Sulfide oxidation.fasta',
                    'SOX': 'sulfur/Sulfur oxidation, SOX system.fasta',
                    'S4I': 'sulfur/Sulfur oxidation, tetrathionate intermediate (S4I) pathway.fasta',
                    'SR': 'sulfur/Sulfur reduction.fasta',
                    'DSR': 'sulfur/Dissimilatory sulfate reduction.fasta'
                },
                'nitrogen_pathways': {
                    'ANR': 'nitrogen/Assimilatory nitrate reduction.fasta',
                    'DEN': 'nitrogen/Denitrification.fasta',
                    'DNR': 'nitrogen/Dissimilatory nitrate reduction.fasta',
                    'NIT': 'nitrogen/Nitrification.fasta'
                },
                'salt_tolerance': 'salt/salt_tolerance_data.fasta',
                'cultivation': 'cultivation/newprotein.fasta'
            },
            'pathway_names': {
                # Methane pathways
                'CO2-CH4': 'Methanogenesis from CO2',
                'METHYLAMINE-CH4': 'Methanogenesis from methylamine',
                'METHANOL-CH4': 'Methanogenesis from methanol',
                'METHANETHIOL-CH4': 'Methanogenesis from methanethiol',
                'ACETIC_ACID-CH4': 'Methanogenesis from acetic acid',
                'C16-CH4': 'Methanogenesis from C16 fatty acid',
                'CO-CH4': 'Methanogenesis from carbon monoxide',
                'FORMIC_ACID-CH4': 'Methanogenesis from formic acid',
                'METHOXY-CH4': 'Methanogenesis from methoxy compounds',
                'FATTY_ACID-CH4': 'Methanogenesis from fatty acids',
                'DIMETHYLAMINE-CH4': 'Methanogenesis from dimethylamine',
                'TRIMETHYLAMINE-CH4': 'Methanogenesis from trimethylamine',
                # Sulfur pathways
                'ASR': 'Assimilatory sulfate reduction',
                'SO': 'Sulfide oxidation',
                'SOX': 'Sulfur oxidation, SOX system',
                'S4I': 'Sulfur oxidation, tetrathionate intermediate (S4I) pathway',
                'SR': 'Sulfur reduction',
                'DSR': 'Dissimilatory sulfate reduction',
                # Nitrogen pathways
                'ANR': 'Assimilatory nitrate reduction',
                'DEN': 'Denitrification',
                'DNR': 'Dissimilatory nitrate reduction',
                'NIT': 'Nitrification',
                # Others
                'SALT_TOLERANCE': 'Salt tolerance',
                'CULTIVATION': 'Cultivation potential'
            },
            'reference_sequence_counts': {
                # Methane pathways
                'CO2-CH4': 12,
                'METHYLAMINE-CH4': 10,
                'METHANOL-CH4': 7,
                'METHANETHIOL-CH4': 6,
                'ACETIC_ACID-CH4': 14,
                'C16-CH4': 15,
                'CO-CH4': 11,
                'FORMIC_ACID-CH4': 15,
                'METHOXY-CH4': 9,
                'FATTY_ACID-CH4': 12,
                'DIMETHYLAMINE-CH4': 10,
                'TRIMETHYLAMINE-CH4': 10,
                # Sulfur pathways
                'ASR': 9,
                'SO': 2,
                'SOX': 7,
                'S4I': 5,
                'SR': 3,
                'DSR': 8,
                # Nitrogen pathways
                'ANR': 9,
                'DEN': 10,
                'DNR': 9,
                'NIT': 4,
                # Others
                'SALT_TOLERANCE': 15,
                'CULTIVATION': 14192
            },
            'output': {
                'base_dir': 'results',
                'formats': ['csv', 'json', 'html']
            },
            'logging': {
                'level': 'INFO',
                'format': '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                'file': 'metharct.log'
            }
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
        elif db_type == 'salt':
            rel_path = self.get('databases.salt_tolerance')
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
        
        # Salt tolerance
        salt_path = self.get('databases.salt_tolerance')
        if salt_path:
            paths['SALT_TOLERANCE'] = os.path.join(base_dir, salt_path)
        
        # Cultivation
        cult_path = self.get('databases.cultivation')
        if cult_path:
            paths['CULTIVATION'] = os.path.join(base_dir, cult_path)
        
        return paths
    
    def save_config(self, output_file: str):
        """Save current configuration to file"""
        with open(output_file, 'w', encoding='utf-8') as f:
            if output_file.endswith('.yaml') or output_file.endswith('.yml'):
                yaml.dump(self.config, f, default_flow_style=False, allow_unicode=True)
            elif output_file.endswith('.json'):
                json.dump(self.config, f, indent=2, ensure_ascii=False)
            else:
                raise ValueError(f"Unsupported output format: {output_file}")
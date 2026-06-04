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
                    'threads': 4,
                    'evalue': 1e-5,
                    'max_target_seqs': 1,
                    'identity_threshold': 30.0
                },
                'tome': {
                    'path': 'tome',
                    'threads': 4
                },
                'checkm2': {
                    'path': 'checkm2',
                    'threads': 4
                }
            },
            'databases': {
                'base_dir': 'data/databases',
                'methane_pathways': {
                    'CO2-CH4': 'methane/processed_CO2-CH4.fasta',
                    'JIAAN-CH4': 'methane/processed_JIAAN-CH4.fasta',
                    'JIACHUN-CH4': 'methane/processed_JIACHUN-CH4.fasta',
                    'JIALIUCHUN-CH4': 'methane/processed_JIALIUCHUN-CH4.fasta',
                    'YISUAN-CH4': 'methane/processed_YISUAN-CH4.fasta',
                    'C16-CH4': 'methane/processed_C16-CH4.fasta',
                    'CO-CH4': 'methane/processed_CO-CH4.fasta',
                    'JIASUAN-CH4': 'methane/processed_JIASUAN-CH4.fasta',
                    'JIAYANGJI-CH4': 'methane/processed_JIAYANGJI-CH4.fasta',
                    'ZHIFANGSUAN-CH4': 'methane/processed_ZHIFANGSUAN-CH4.fasta',
                    '2JIAAN-CH4': 'methane/processed_2JIAAN-CH4.fasta',
                    '3JIAAN-CH4': 'methane/processed_3JIAAN-CH4.fasta',
                    'Glycine betaine methanogenesis': 'methane/processed_甘氨酸甜菜碱产甲烷.fasta',
                    'Methylthiopropionate methanogenesis': 'methane/processed_硫代丙酸甲酯产甲烷.fasta',
                    'Tetramethylammonium methanogenesis': 'methane/四甲基铵产甲烷.fasta',
                    'Methanol dismutation methanogenesis': 'methane/processed_甲醇歧化产甲烷.fasta'
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
                'salt_tolerance': 'salt/naiyan_data.fasta',
                'cultivation_pathways': {
                    'amino_acid': 'cultivation/path.aa',
                    'vitamin': 'cultivation/path.vc',
                    'nucleotide': 'cultivation/path_hesuan',
                    'atp': 'cultivation/path_atp'
                }
            },
            'pathway_names': {
                # Methane pathways
                'CO2-CH4': 'CO2 reduction methanogenesis pathway',
                'JIAAN-CH4': 'Methylamine methanogenesis pathway',
                'JIACHUN-CH4': 'Methanol methanogenesis pathway',
                'JIALIUCHUN-CH4': 'Methanethiol methanogenesis pathway',
                'YISUAN-CH4': 'Acetate methanogenesis pathway',
                'C16-CH4': 'C16 methanogenesis pathway',
                'CO-CH4': 'CO methanogenesis pathway',
                'JIASUAN-CH4': 'Formate methanogenesis pathway',
                'JIAYANGJI-CH4': 'Methoxy methanogenesis pathway',
                'ZHIFANGSUAN-CH4': 'Fatty acid methanogenesis pathway',
                '2JIAAN-CH4': 'Dimethylamine methanogenesis pathway',
                '3JIAAN-CH4': 'Trimethylamine methanogenesis pathway',
                'Glycine betaine methanogenesis': 'Glycine betaine methanogenesis pathway',
                'Methylthiopropionate methanogenesis': 'Methylthiopropionate methanogenesis pathway',
                'Tetramethylammonium methanogenesis': 'Tetramethylammonium methanogenesis pathway',
                'Methanol dismutation methanogenesis': 'Methanol dismutation methanogenesis pathway',
                # Sulfur pathways
                'ASR': 'Assimilatory sulfate reduction pathway',
                'SO': 'Sulfide oxidation pathway',
                'SOX': 'Sulfur oxidation SOX system pathway',
                'S4I': 'Sulfur oxidation tetrathionate intermediate pathway',
                'SR': 'Sulfur reduction pathway',
                'DSR': 'Dissimilatory sulfate reduction pathway',
                # Nitrogen pathways
                'ANR': 'Assimilatory nitrate reduction pathway',
                'DEN': 'Denitrification pathway',
                'DNR': 'Dissimilatory nitrate reduction pathway',
                'NIT': 'Nitrification pathway',
                # Others
                'NAIYAN': 'Salt tolerance characteristics',
                'CULTIVATION': 'Cultivation assessment'
            },
            'reference_sequence_counts': {
                # Methane pathways
                'CO2-CH4': 15,
                'JIAAN-CH4': 10,
                'JIACHUN-CH4': 7,
                'JIALIUCHUN-CH4': 6,
                'YISUAN-CH4': 14,
                'Glycine betaine methanogenesis': 3,
                'Methylthiopropionate methanogenesis': 4,
                'Tetramethylammonium methanogenesis': 3,
                'Methanol dismutation methanogenesis': 15,
                'C16-CH4': 15,
                'CO-CH4': 11,
                'JIASUAN-CH4': 18,
                'JIAYANGJI-CH4': 9,
                'ZHIFANGSUAN-CH4': 12,
                '2JIAAN-CH4': 10,
                '3JIAAN-CH4': 10,
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
                'NAIYAN': 15,
                'CULTIVATION': 14192
            },
            'cultivation': {
                'diamond_path': 'diamond',
                'threads': 4,
                'evalue': 1e-5,
                'identity_threshold': 30.0,
                'pathway_weights': {
                    'amino_acid': 0.4,
                    'vitamin': 0.3,
                    'nucleotide': 0.2,
                    'atp': 0.1
                },
                'cultivability_thresholds': {
                    'high': 0.8,
                    'medium': 0.6,
                    'low': 0.4
                }
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
            paths['NAIYAN'] = os.path.join(base_dir, salt_path)
        
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
"""预测模块 - 加载模型并对特征进行温度预测"""

import os
import pickle
import warnings
import numpy as np

# 抑制 sklearn 版本不一致的警告（模型跨版本兼容）
warnings.filterwarnings('ignore', message='Trying to unpickle estimator')

# 模型文件路径（包内相对路径）
_MODEL_DIR = os.path.join(os.path.dirname(__file__), 'models')

# 模型文件名映射
_MODEL_FILES = {
    'T_min': 'T_min_BayesianRidge.pkl',
    'T_opt': 'T_opt_BayesianRidge.pkl',
    'T_max': 'T_max_Ridge(a=100).pkl',
}

# 全局模型缓存
_models = {}


def _load_model(target: str):
    """加载指定目标变量的模型（带缓存）"""
    if target not in _models:
        model_path = os.path.join(_MODEL_DIR, _MODEL_FILES[target])
        if not os.path.exists(model_path):
            raise FileNotFoundError(f"模型文件未找到: {model_path}")
        with open(model_path, 'rb') as f:
            _models[target] = pickle.load(f)
    return _models[target]


def predict(features: np.ndarray) -> dict:
    """对输入特征进行温度预测

    Parameters
    ----------
    features : np.ndarray
        形状为 (1, 420) 的特征矩阵

    Returns
    -------
    dict
        包含 T_min, T_opt, T_max 预测值的字典，单位：°C
    """
    results = {}
    for target in ['T_min', 'T_opt', 'T_max']:
        model = _load_model(target)
        pred = model.predict(features)
        results[target] = float(pred[0])
    return results


def predict_from_fasta(fasta_path: str) -> dict:
    """从 FASTA 文件直接预测生长温度范围

    Parameters
    ----------
    fasta_path : str
        FASTA 蛋白质序列文件路径

    Returns
    -------
    dict
        包含 T_min, T_opt, T_max 预测值的字典，单位：°C
    """
    from .features import extract_features_from_fasta
    features = extract_features_from_fasta(fasta_path)
    return predict(features)

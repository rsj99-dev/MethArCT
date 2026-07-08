"""特征提取模块 - 从蛋白质序列中提取 AAC 和 Dipeptide Composition 特征"""

import numpy as np

# 标准 20 种氨基酸（字母顺序）
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
AA_TO_IDX = {aa: i for i, aa in enumerate(AMINO_ACIDS)}
N_AA = len(AMINO_ACIDS)   # 20
N_FEATURES = N_AA + N_AA * N_AA  # 420


def _count_aa(sequence: str) -> np.ndarray:
    """统计一条蛋白质序列中各氨基酸的出现次数"""
    counts = np.zeros(N_AA, dtype=np.float64)
    for aa in sequence:
        if aa in AA_TO_IDX:
            counts[AA_TO_IDX[aa]] += 1
    return counts


def _count_dipeptide(sequence: str) -> np.ndarray:
    """统计一条蛋白质序列中二肽的出现次数"""
    counts = np.zeros((N_AA, N_AA), dtype=np.float64)
    for i in range(len(sequence) - 1):
        aa1 = sequence[i]
        aa2 = sequence[i + 1]
        if aa1 in AA_TO_IDX and aa2 in AA_TO_IDX:
            counts[AA_TO_IDX[aa1], AA_TO_IDX[aa2]] += 1
    return counts


def extract_features_from_sequences(sequences: list) -> np.ndarray:
    """从蛋白质序列列表中提取 420 维特征向量

    特征顺序: AAC (20) + DC (400)
    AAC: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
    DC:  AA, AC, AD, ..., AY, CA, CC, ..., YY

    Parameters
    ----------
    sequences : list of str
        蛋白质序列列表（氨基酸字符串）

    Returns
    -------
    np.ndarray
        形状为 (1, 420) 的特征矩阵（对整个蛋白质组求平均）
    """
    if not sequences:
        raise ValueError("输入的序列列表为空")

    n_seq = len(sequences)
    aac_accum = np.zeros(N_AA, dtype=np.float64)
    dc_accum = np.zeros((N_AA, N_AA), dtype=np.float64)
    total_aa = 0

    for seq in sequences:
        seq_len = len(seq)
        if seq_len == 0:
            continue
        # AAC: 按序列长度归一化
        aa_counts = _count_aa(seq)
        total_aa += np.sum(aa_counts)
        aac_accum += aa_counts

        # DC: 按 (seq_len - 1) 归一化后累加，最后对所有序列取平均
        dp_counts = _count_dipeptide(seq)
        total_dp = np.sum(dp_counts)
        if total_dp > 0:
            dc_accum += dp_counts / total_dp

    # AAC: 整组归一化为频率
    if total_aa > 0:
        aac = aac_accum / total_aa
    else:
        aac = np.zeros(N_AA)

    # DC: 对序列数取平均
    if n_seq > 0:
        dc = dc_accum / n_seq
    else:
        dc = np.zeros((N_AA, N_AA))

    # 拼接为 420 维向量
    dc_flat = dc.flatten()  # 行优先展开: AA, AC, AD, ..., AY, CA, CC, ..., YY
    features = np.concatenate([aac, dc_flat])
    return features.reshape(1, -1)


def extract_features_from_fasta(fasta_path: str) -> np.ndarray:
    """从 FASTA 文件中读取所有蛋白质序列并提取特征

    Parameters
    ----------
    fasta_path : str
        FASTA 文件路径 (.faa 或 .fasta)

    Returns
    -------
    np.ndarray
        形状为 (1, 420) 的特征矩阵
    """
    sequences = []
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_seq:
                    sequences.append(''.join(current_seq))
                    current_seq = []
            else:
                current_seq.append(line)

        # 最后一条序列
        if current_seq:
            sequences.append(''.join(current_seq))

    if not sequences:
        raise ValueError(f"FASTA 文件中未找到任何蛋白质序列: {fasta_path}")

    return extract_features_from_sequences(sequences)

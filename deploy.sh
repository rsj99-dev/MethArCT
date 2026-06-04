#!/bin/bash
# ============================================================================
# MethArCT 龙芯 (LoongArch64) 一键部署脚本
# ============================================================================
# 本脚本用于在龙芯架构上完整部署 MethArCT 分析工具
# 包括：Python 环境、依赖包、DIAMOND 编译、Tome 安装
# ============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="${SCRIPT_DIR}"
INSTALL_PREFIX="${PROJECT_DIR}/loongarch_install"

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log_info()  { echo -e "${BLUE}[INFO]${NC} $1"; }
log_ok()    { echo -e "${GREEN}[OK]${NC} $1"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }

echo "============================================================"
echo "  MethArCT 龙芯 (LoongArch64) 一键部署脚本"
echo "============================================================"
echo ""

# ------------------------------------------------------------------
# Step 0: 架构检查
# ------------------------------------------------------------------
ARCH=$(uname -m)
log_info "检测到系统架构: ${ARCH}"

if [[ "$ARCH" != "loongarch64" ]]; then
    log_warn "当前架构为 ${ARCH}，非 loongarch64"
    read -p "是否继续？(y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "已取消"
        exit 0
    fi
fi

# 检查操作系统
OS_INFO=""
if [ -f /etc/os-release ]; then
    OS_INFO=$(cat /etc/os-release | grep "^PRETTY_NAME" | cut -d'"' -f2)
    log_info "操作系统: ${OS_INFO}"
fi

# ------------------------------------------------------------------
# Step 1: 检查系统依赖
# ------------------------------------------------------------------
echo ""
echo "============================================================"
log_info "Step 1/6: 检查系统编译工具..."
echo "============================================================"

MISSING_TOOLS=()

check_tool() {
    if command -v "$1" &> /dev/null; then
        log_ok "$1: $(command -v $1)"
        return 0
    else
        log_error "$1 未安装"
        MISSING_TOOLS+=("$1")
        return 1
    fi
}

check_tool "gcc"
check_tool "g++"
check_tool "cmake"
check_tool "make"
check_tool "python3"
check_tool "pip3" || check_tool "pip"

if [ ${#MISSING_TOOLS[@]} -gt 0 ]; then
    echo ""
    log_error "以下工具缺失，请先安装:"
    for tool in "${MISSING_TOOLS[@]}"; do
        echo "  - $tool"
    done
    echo ""
    echo "安装方法:"
    echo "  Ubuntu/Loongnix: sudo apt-get install build-essential cmake python3 python3-pip"
    echo "  CentOS/RHEL:     sudo yum install gcc gcc-c++ cmake make python3 python3-pip"
    echo "  Conda:           conda install gcc_linux-loongarch64 gxx_linux-loongarch64 cmake make python"
    exit 1
fi

# ------------------------------------------------------------------
# Step 2: 创建 Python 虚拟环境
# ------------------------------------------------------------------
echo ""
echo "============================================================"
log_info "Step 2/6: 创建 Python 虚拟环境..."
echo "============================================================"

VENV_DIR="${INSTALL_PREFIX}/venv"

if [ -d "${VENV_DIR}" ]; then
    log_warn "虚拟环境已存在: ${VENV_DIR}"
    read -p "是否重新创建？(y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        rm -rf "${VENV_DIR}"
    fi
fi

if [ ! -d "${VENV_DIR}" ]; then
    python3 -m venv "${VENV_DIR}"
    log_ok "虚拟环境已创建: ${VENV_DIR}"
fi

# 激活虚拟环境
source "${VENV_DIR}/bin/activate"
log_info "Python: $(which python)"
log_info "Python 版本: $(python --version)"

# ------------------------------------------------------------------
# Step 3: 安装 Python 依赖
# ------------------------------------------------------------------
echo ""
echo "============================================================"
log_info "Step 3/6: 安装 Python 依赖包..."
echo "============================================================"

# 升级 pip
pip install --upgrade pip setuptools wheel

# 安装核心依赖 (龙芯架构下，conda-forge 提供的预编译包更可靠)
# 如果 pip 安装失败（无预编译 wheel），尝试从源码编译
install_package() {
    local pkg=$1
    log_info "安装 ${pkg}..."
    if pip install "${pkg}" 2>/dev/null; then
        log_ok "${pkg} 安装成功"
    else
        log_warn "${pkg} 预编译包不可用，尝试源码编译..."
        pip install --no-binary=:all: "${pkg}" || {
            log_error "${pkg} 安装失败，请手动安装"
            return 1
        }
    fi
}

# 科学计算基础包 (这些在龙芯上可能需要源码编译)
install_package "numpy>=1.21.0"
install_package "scipy>=1.7.0"
install_package "pandas>=1.3.0"
install_package "scikit-learn>=1.0.0"

# 纯 Python 包 (直接安装)
pip install \
    "biopython>=1.79" \
    "pyyaml>=5.4" \
    "tqdm>=4.60.0" \
    "requests>=2.25.0" \
    "joblib>=1.0.0" \
    "loguru>=0.6.0"

log_ok "Python 依赖安装完成"

# ------------------------------------------------------------------
# Step 4: 编译 DIAMOND
# ------------------------------------------------------------------
echo ""
echo "============================================================"
log_info "Step 4/6: 编译 DIAMOND..."
echo "============================================================"

DIAMOND_BIN="${INSTALL_PREFIX}/bin/diamond"

if [ -f "${DIAMOND_BIN}" ]; then
    log_ok "DIAMOND 已存在: ${DIAMOND_BIN}"
    "${DIAMOND_BIN}" version
else
    log_info "开始编译 DIAMOND..."
    cd "${PROJECT_DIR}"
    bash build_diamond.sh

    if [ -f "${PROJECT_DIR}/bin/diamond" ]; then
        mkdir -p "${INSTALL_PREFIX}/bin"
        cp "${PROJECT_DIR}/bin/diamond" "${INSTALL_PREFIX}/bin/diamond"
        log_ok "DIAMOND 编译安装成功"
    else
        log_error "DIAMOND 编译失败，请查看 build_diamond.sh 输出"
        echo "  你也可以手动编译后放到: ${INSTALL_PREFIX}/bin/diamond"
    fi
fi

# ------------------------------------------------------------------
# Step 5: 安装 Tome
# ------------------------------------------------------------------
echo ""
echo "============================================================"
log_info "Step 5/6: 安装 Tome (OGT 预测工具)..."
echo "============================================================"

TOME_DIR="${PROJECT_DIR}/Tome-1.1.0"

if [ -d "${TOME_DIR}" ]; then
    cd "${TOME_DIR}"
    pip install -e . 2>/dev/null || pip install .
    log_ok "Tome 安装完成"
else
    log_warn "Tome-1.1.0 目录不存在，跳过 Tome 安装"
    echo "  请将 Tome-1.1.0 放到项目根目录后重新运行"
fi

# ------------------------------------------------------------------
# Step 6: 安装 MethArCT
# ------------------------------------------------------------------
echo ""
echo "============================================================"
log_info "Step 6/6: 安装 MethArCT..."
echo "============================================================"

cd "${PROJECT_DIR}"
pip install -e . 2>/dev/null || pip install .

log_ok "MethArCT 安装完成"

# ------------------------------------------------------------------
# 生成龙芯专用配置文件
# ------------------------------------------------------------------
echo ""
echo "============================================================"
log_info "生成龙芯专用配置文件..."
echo "============================================================"

CONFIG_FILE="${PROJECT_DIR}/metharct_config.yaml"

cat > "${CONFIG_FILE}" << EOF
# MethArCT 龙芯 (LoongArch64) 专用配置
# 自动生成于部署脚本

tools:
  diamond:
    path: "${INSTALL_PREFIX}/bin/diamond"
    use_wsl: false
    threads: 4
    evalue: 1e-5
    max_target_seqs: 1
    identity_threshold: 30.0

  tome:
    path: "python_module"
    use_wsl: false
    threads: 4

  checkm2:
    path: "checkm2"
    use_wsl: false
    threads: 4

databases:
  base_dir: "data/databases"
EOF

log_ok "配置文件已生成: ${CONFIG_FILE}"

# ------------------------------------------------------------------
# 验证安装
# ------------------------------------------------------------------
echo ""
echo "============================================================"
echo "  安装验证"
echo "============================================================"

echo ""
log_info "验证 MethArCT..."
if python -c "import metharct; print(f'  版本: {metharct.__version__}')" 2>/dev/null; then
    log_ok "MethArCT 模块加载成功"
else
    log_warn "MethArCT 模块加载失败，请检查安装"
fi

echo ""
log_info "验证 DIAMOND..."
if [ -f "${INSTALL_PREFIX}/bin/diamond" ]; then
    "${INSTALL_PREFIX}/bin/diamond" version
    log_ok "DIAMOND 可用"
else
    log_warn "DIAMOND 未找到"
fi

echo ""
log_info "验证 Tome..."
if python -c "from tome.tome import load_model; print('  Tome 模块可用')" 2>/dev/null; then
    log_ok "Tome 可用"
else
    log_warn "Tome 不可用（OGT 预测功能将受限）"
fi

# ------------------------------------------------------------------
# 最终输出
# ------------------------------------------------------------------
echo ""
echo "============================================================"
echo "  MethArCT 龙芯部署完成！"
echo "============================================================"
echo ""
echo "激活环境:"
echo "  source ${VENV_DIR}/bin/activate"
echo ""
echo "运行分析:"
echo "  metharct analyze <input.fasta> --config ${CONFIG_FILE}"
echo ""
echo "或使用 Python API:"
echo "  python -c \"from metharct import MethArCTAnalyzer; ...\""
echo ""
echo "配置文件: ${CONFIG_FILE}"
echo "安装目录: ${INSTALL_PREFIX}"
echo ""

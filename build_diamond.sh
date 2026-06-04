#!/bin/bash
# ============================================================================
# DIAMOND 龙芯 (LoongArch64) 编译脚本
# ============================================================================
# DIAMOND 是 MethArCT 的核心依赖，需要从源码为龙芯架构编译
# 项目地址: https://github.com/bbuchfink/diamond
#
# 龙芯说明：DIAMOND 默认构建 x86 SIMD 优化版本 (SSE4.1/SSE4.2/AVX/AVX2)，
# 这些指令集在龙芯上不存在。本脚本会自动修补 CMakeLists.txt，
# 只编译 generic（通用标量）版本。
# ============================================================================

set -e

DIAMOND_VERSION="2.1.9"
BUILD_DIR="$(pwd)/diamond_build"
INSTALL_DIR="$(pwd)/bin"

echo "============================================"
echo "  DIAMOND ${DIAMOND_VERSION} 龙芯编译脚本"
echo "============================================"

# 检查架构
ARCH=$(uname -m)
echo "[INFO] 检测到架构: ${ARCH}"

if [[ "$ARCH" != "loongarch64" ]]; then
    echo "[WARN] 当前架构为 ${ARCH}，此脚本专为 loongarch64 设计"
fi

# 检查编译依赖
check_dependency() {
    if ! command -v "$1" &> /dev/null; then
        echo "[ERROR] 未找到 $1，请先安装"
        return 1
    fi
    echo "[OK] $1 已安装: $(command -v $1)"
    return 0
}

echo ""
echo "[1/6] 检查编译依赖..."
echo "--------------------------------------------"

check_dependency "cmake" || exit 1
check_dependency "gcc" || exit 1
check_dependency "g++" || exit 1
check_dependency "make" || exit 1
check_dependency "python3" || exit 1

echo ""
echo "[2/6] 下载 DIAMOND ${DIAMOND_VERSION} 源码..."
echo "--------------------------------------------"

mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

DIAMOND_TARBALL="v${DIAMOND_VERSION}.tar.gz"
DIAMOND_URL="https://github.com/bbuchfink/diamond/archive/refs/tags/${DIAMOND_TARBALL}"

if [ ! -f "${DIAMOND_TARBALL}" ]; then
    echo "从 ${DIAMOND_URL} 下载..."
    if command -v wget &> /dev/null; then
        wget -O "${DIAMOND_TARBALL}" "${DIAMOND_URL}"
    elif command -v curl &> /dev/null; then
        curl -L -o "${DIAMOND_TARBALL}" "${DIAMOND_URL}"
    else
        echo "[ERROR] 需要 wget 或 curl"
        exit 1
    fi
else
    echo "源码包已存在，跳过下载"
fi

echo ""
echo "[3/6] 解压源码..."
echo "--------------------------------------------"

if [ -d "diamond-${DIAMOND_VERSION}" ]; then
    rm -rf "diamond-${DIAMOND_VERSION}"
fi
tar xzf "${DIAMOND_TARBALL}"
cd "diamond-${DIAMOND_VERSION}"

echo ""
echo "[4/6] 配置 DIAMOND（禁用 x86 SIMD，使用 generic 架构）..."
echo "--------------------------------------------"

# DIAMOND 的 CMakeLists.txt 中 option(X86 "X86" ON) 默认开启 x86 SIMD 编译，
# 龙芯不支持这些指令集。通过 -DX86=OFF 告诉 CMake 跳过所有 x86 SIMD 目标，
# 只编译 generic（通用标量）版本。
echo "  [INFO] 将使用 -DX86=OFF 禁用 x86 SIMD 优化"

echo ""
echo "[5/6] 编译 DIAMOND（仅 generic 架构）..."
echo "--------------------------------------------"

# 获取 CPU 核心数
NPROC=$(nproc 2>/dev/null || echo 4)

# 创建构建目录
rm -rf build
mkdir -p build && cd build

# 龙芯编译选项
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}" \
    -DCMAKE_CXX_FLAGS="-O3" \
    -DX86=OFF

echo "使用 ${NPROC} 个核心编译..."
make -j${NPROC}

echo ""
echo "[6/6] 安装..."
echo "--------------------------------------------"

make install 2>/dev/null || true

# 验证安装
echo ""
echo "============================================"
echo "  编译结果验证"
echo "============================================"

DIAMOND_BIN=""
if [ -f "${INSTALL_DIR}/bin/diamond" ]; then
    DIAMOND_BIN="${INSTALL_DIR}/bin/diamond"
elif [ -f "diamond" ]; then
    mkdir -p "${INSTALL_DIR}/bin"
    cp diamond "${INSTALL_DIR}/bin/diamond"
    chmod +x "${INSTALL_DIR}/bin/diamond"
    DIAMOND_BIN="${INSTALL_DIR}/bin/diamond"
elif [ -f "../diamond" ]; then
    mkdir -p "${INSTALL_DIR}/bin"
    cp ../diamond "${INSTALL_DIR}/bin/diamond"
    chmod +x "${INSTALL_DIR}/bin/diamond"
    DIAMOND_BIN="${INSTALL_DIR}/bin/diamond"
fi

if [ -n "${DIAMOND_BIN}" ] && [ -f "${DIAMOND_BIN}" ]; then
    echo "[OK] DIAMOND 已安装到: ${DIAMOND_BIN}"
    "${DIAMOND_BIN}" version
    echo ""
    echo "使用方法:"
    echo "  1. 添加到 PATH: export PATH=\"${INSTALL_DIR}/bin:\$PATH\""
    echo "  2. 或在 metharct_config.yaml 中设置:"
    echo "     tools:"
    echo "       diamond:"
    echo "         path: \"${DIAMOND_BIN}\""
    echo ""
    echo "编译成功！"
else
    echo "[ERROR] DIAMOND 编译失败，请检查上方错误信息"
    echo ""
    echo "如需手动调试："
    echo "  cd ${BUILD_DIR}/diamond-${DIAMOND_VERSION}/build"
    echo "  cmake .. -DCMAKE_BUILD_TYPE=Release"
    echo "  make arch_generic -j${NPROC}"
    exit 1
fi

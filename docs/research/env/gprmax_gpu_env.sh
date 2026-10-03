# Environment for gprMax V4 GPU solves on this machine.
#
# Source this before running a GPU batch:
#   source artifacts/local_checks/cuda12_toolkit/env.sh
#
# Why it is needed. Two independent problems blocked the GPU path:
#
#  1. The system CUDA toolkit is 13.3, whose nvcc lists only
#     compute_100/103/110/120/121 (Blackwell). This machine's GPU is an
#     RTX 4090 Laptop with compute capability 8.9 (Ada), so nvcc 13.3 rejects
#     "-arch sm_89" and PyCUDA raises CompileError before any field array is
#     allocated. A project-local CUDA 12.4 nvcc (which does list sm_89) is
#     installed under cuda12_toolkit/cuda12.
#
#  2. The original error was simpler: nvcc needs the MSVC host compiler.
#     Without cl.exe it stops with "nvcc fatal : Cannot find compiler
#     'cl.exe' in PATH", and even with cl.exe it needs INCLUDE to point at the
#     MSVC and Windows SDK headers or crtdefs.h fails on a missing corecrt.h.
#     MSVC BuildTools 2022 and SDK 10.0.26100.0 are installed, but vcvars64.bat
#     cannot be invoked from this shell, so PATH, INCLUDE and LIB are set
#     explicitly here.
#
# PyCUDA's compile() defaults to nvcc="nvcc" and resolves it from PATH, so
# PATH is what actually redirects the compiler; CUDA_PATH alone is not used
# for nvcc lookup. INCLUDE and LIB use native Windows paths with semicolon
# separators because nvcc is a Windows program and does not understand POSIX
# paths. The system CUDA installation is NOT modified and nothing set here
# leaks outside this shell process.

# nvcc needs Windows-style paths, so build them with single quotes to keep the
# backslashes literal.
MSVC_INC='C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Tools\MSVC\14.44.35207\include'
MSVC_LIB='C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Tools\MSVC\14.44.35207\lib\x64'
SDK_INC_UCRT='C:\Program Files (x86)\Windows Kits\10\Include\10.0.26100.0\ucrt'
SDK_INC_UM='C:\Program Files (x86)\Windows Kits\10\Include\10.0.26100.0\um'
SDK_INC_SHARED='C:\Program Files (x86)\Windows Kits\10\Include\10.0.26100.0\shared'
SDK_LIB_UCRT='C:\Program Files (x86)\Windows Kits\10\Lib\10.0.26100.0\ucrt\x64'
SDK_LIB_UM='C:\Program Files (x86)\Windows Kits\10\Lib\10.0.26100.0\um\x64'

# POSIX form on purpose: REPO is interpolated into PATH, and a Windows-style
# "E:/..." value would be split at the colon by the path search.
export REPO="${REPO:-/e/automation_djh/automation_repo}"
export CUDA12="$REPO/artifacts/local_checks/cuda12_toolkit/cuda12"
MSVC_BIN_POSIX='/c/Program Files (x86)/Microsoft Visual Studio/2022/BuildTools/VC/Tools/MSVC/14.44.35207/bin/Hostx64/x64'

# Fail loudly rather than silently falling back to the system nvcc, which
# produced a confusing "sm_89 not supported" error instead of "no nvcc".
if [ ! -x "$CUDA12/bin/nvcc.exe" ]; then
  echo "ERROR: project CUDA 12.4 nvcc missing at $CUDA12/bin/nvcc.exe" >&2
  return 1 2>/dev/null || exit 1
fi
if [ ! -f "$MSVC_BIN_POSIX/cl.exe" ]; then
  echo "ERROR: MSVC cl.exe missing under $MSVC_BIN_POSIX" >&2
  return 1 2>/dev/null || exit 1
fi

# CUDA 12 first so plain "nvcc" resolves to the Ada-capable compiler.
export PATH="$CUDA12/bin:$MSVC_BIN_POSIX:$PATH"
export CUDA_PATH="$CUDA12"
export INCLUDE="$MSVC_INC;$SDK_INC_UCRT;$SDK_INC_UM;$SDK_INC_SHARED"
export LIB="$MSVC_LIB;$SDK_LIB_UCRT;$SDK_LIB_UM"

# gprMax reads this for its in-process CUBIN cache; default to this run dir.
export HS4_CUDA_CACHE_LOG="${HS4_CUDA_CACHE_LOG:-$PWD/cuda_cache.jsonl}"

# gprMax V4 GPU environment, project-local venv.
export HS4_PYTHON="$REPO/artifacts/local_checks/gprmax_v4_gpu_env/Scripts/python.exe"
export HS4_ENTRY="$REPO/scripts/gprmax_cached_cuda_entry.py"

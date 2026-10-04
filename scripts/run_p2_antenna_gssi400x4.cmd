@echo off
setlocal
if not defined HS4_VCVARS64 set "HS4_VCVARS64=E:\msvc2022bt\VC\Auxiliary\Build\vcvars64.bat"
if not defined HS4_CUDA_BIN set "HS4_CUDA_BIN=C:\cuda133\bin"
if not defined HS4_V4_PYTHON set "HS4_V4_PYTHON=E:\automation_djh\artifacts\local_checks\gprmax_v4_gpu_env\Scripts\python.exe"
call "%HS4_VCVARS64%"
if errorlevel 1 exit /b 1
set "PATH=%HS4_CUDA_BIN%;%PATH%"
"%HS4_V4_PYTHON%" "%~dp0p2_antenna_gssi400x4_bscan.py" %*
exit /b %errorlevel%

@echo off
setlocal
if not defined HS4_VCVARS64 set "HS4_VCVARS64=E:\sisual stdio 2022\VC\Auxiliary\Build\vcvars64.bat"
if not defined HS4_CUDA_BIN set "HS4_CUDA_BIN=C:\cuda118\bin"
if not defined HS4_V4_PYTHON set "HS4_V4_PYTHON=D:\gprmax_v4_gpu_env\Scripts\python.exe"
call "%HS4_VCVARS64%"
if errorlevel 1 exit /b 1
set "PATH=%HS4_CUDA_BIN%;%PATH%"
"%HS4_V4_PYTHON%" "%~dp0hs4_large_domain_dense.py" %*
exit /b %errorlevel%

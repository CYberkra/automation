@echo off
setlocal
set "PYTHONPATH=E:\djh\artifacts\local_checks\gprmax_v401_gpu_env\Lib\site-packages;E:\djh\artifacts\local_checks\gprmax_v4_gpu_env2\Lib\site-packages"
set "GPRMAX_PYTHON=E:\djh\artifacts\local_checks\py312conda\python.exe"
set "GPRMAX_VCVARS=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
set "GPRMAX_CUDA_BIN=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.3\bin"
cd /d E:\djh\line9_paper_mudstone_20261011_r2
call scripts\run_gprmax_v401.cmd scripts\line9_paper_mudstone_window.py %* --out execution
exit /b %errorlevel%

@echo off
setlocal
call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 exit /b 1
set "CUDA_BIN=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.3\bin"
if exist "%CUDA_BIN%\x64\curand64_10.dll" set "PATH=%CUDA_BIN%\x64;%PATH%"
if exist "%CUDA_BIN%\nvcc.exe" set "PATH=%CUDA_BIN%;%PATH%"
cd /d "E:\automation_djh\verify_runs\2026-10-02_halfspace_standard"
"E:\automation_djh\automation_repo\artifacts\local_checks\gprmax_v4_gpu_env\Scripts\python.exe" -u -c "import runpy,sys; sys.argv=['gprMax','%~1.in','--geometry-only']; runpy.run_module('gprMax',run_name='__main__')"
exit /b %errorlevel%

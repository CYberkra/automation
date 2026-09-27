@echo off
setlocal
set "GPRMAX_VS=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools"
if not exist "%GPRMAX_VS%\VC\Auxiliary\Build\vcvars64.bat" (
  echo "MSVC x64 Build Tools not found at %GPRMAX_VS%. No CPU fallback."
  exit /b 1
)
call "%GPRMAX_VS%\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 exit /b 1
set "CUDA_BIN=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.3\bin"
if exist "%CUDA_BIN%\x64\curand64_10.dll" set "PATH=%CUDA_BIN%\x64;%PATH%"
if exist "%CUDA_BIN%\nvcc.exe" set "PATH=%CUDA_BIN%;%PATH%"
"%~dp0..\artifacts\local_checks\gprmax_v4_gpu_env\Scripts\python.exe" "%~dp0run_approved_batch2d_slope_t3_co.py" %*
exit /b %errorlevel%

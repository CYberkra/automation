@echo off
setlocal
if not defined GPRMAX_PYTHON (
  echo Set GPRMAX_PYTHON to the verified V4 python.exe.
  exit /b 2
)
if not exist "%GPRMAX_PYTHON%" exit /b 2
if not defined GPRMAX_VCVARS (
  echo Set GPRMAX_VCVARS to the target-machine vcvars64.bat.
  exit /b 2
)
if not exist "%GPRMAX_VCVARS%" exit /b 2
if not defined GPRMAX_CUDA_BIN (
  echo Set GPRMAX_CUDA_BIN to the target-machine CUDA compiler bin directory.
  exit /b 2
)
if not exist "%GPRMAX_CUDA_BIN%" exit /b 2
call "%GPRMAX_VCVARS%"
if errorlevel 1 exit /b 1
set "PATH=%GPRMAX_CUDA_BIN%;%PATH%"
"%GPRMAX_PYTHON%" %*
exit /b %errorlevel%

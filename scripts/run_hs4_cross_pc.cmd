@echo off
setlocal
if not defined HS4_PYTHON (
  echo Set HS4_PYTHON to this machine's gprMax V4 Python executable.
  exit /b 2
)
if not exist "%HS4_PYTHON%" exit /b 2
if defined HS4_VCVARS (
  if not exist "%HS4_VCVARS%" exit /b 2
  call "%HS4_VCVARS%" >nul
  if errorlevel 1 exit /b 2
)
if defined HS4_CUDA_BIN set "PATH=%HS4_CUDA_BIN%;%PATH%"
if "%~1"=="run" (
  where cl >nul 2>nul
  if errorlevel 1 (
    echo Run in x64 Native Tools or set HS4_VCVARS.
    exit /b 2
  )
  where nvcc >nul 2>nul
  if errorlevel 1 exit /b 2
)
"%HS4_PYTHON%" "%~dp0hs4_cross_pc.py" %*
exit /b %errorlevel%

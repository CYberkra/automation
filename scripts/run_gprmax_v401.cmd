@echo off
setlocal
if not defined GPRMAX_PYTHON set "GPRMAX_PYTHON=%~dp0..\artifacts\local_checks\gprmax_v401_gpu_env\Scripts\python.exe"
if not exist "%GPRMAX_PYTHON%" (
  echo Set GPRMAX_PYTHON to the verified 4.0.1 environment.
  exit /b 2
)
if not defined GPRMAX_VCVARS (
  if exist "E:\sisual stdio 2022\VC\Auxiliary\Build\vcvars64.bat" (
    set "GPRMAX_VCVARS=E:\sisual stdio 2022\VC\Auxiliary\Build\vcvars64.bat"
  ) else (
    set "GPRMAX_VCVARS=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
  )
)
if not defined GPRMAX_CUDA_BIN (
  if exist "C:\cuda118\bin\nvcc.exe" (
    set "GPRMAX_CUDA_BIN=C:\cuda118\bin"
  ) else (
    set "GPRMAX_CUDA_BIN=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.3\bin"
  )
)
if not exist "%GPRMAX_VCVARS%" exit /b 2
if not exist "%GPRMAX_CUDA_BIN%\nvcc.exe" exit /b 2
call "%GPRMAX_VCVARS%" >nul
if errorlevel 1 exit /b 1
set "PATH=%GPRMAX_CUDA_BIN%;%PATH%"
"%GPRMAX_PYTHON%" -c "import gprMax; assert gprMax.__version__ == '4.0.1', gprMax.__version__"
if errorlevel 1 exit /b 1
"%GPRMAX_PYTHON%" %*
exit /b %errorlevel%

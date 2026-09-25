@echo off
setlocal
set "GPRMAX_VS="
for /f "usebackq tokens=*" %%i in (`"%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "GPRMAX_VS=%%i"
if not defined GPRMAX_VS (
  echo MSVC x64 Build Tools not found.
  exit /b 1
)
call "%GPRMAX_VS%\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 exit /b 1
rem ASCII aliases: nvcc/cl.exe cannot parse the non-ASCII or space-containing include paths.
set "PATH=C:\cuda118\bin;%PATH%"
set "PY=D:\gprmax_v4_gpu_env\Scripts\python.exe"
"%PY%" "%~dp0inspect_v4_runtime.py" --output "%~dp0..\artifacts\research_checks\2026-09-25_local_gpu_setup\runtime_identity.json" || exit /b 1
"%PY%" "%~dp0inspect_cuda_runtime.py" --output "%~dp0..\artifacts\research_checks\2026-09-25_local_gpu_setup\cuda_identity.json" || exit /b 1
exit /b 0

@echo off
setlocal
set "GPRMAX_VS="
for /f "usebackq tokens=*" %%i in (`"%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "GPRMAX_VS=%%i"
if not defined GPRMAX_VS (
  echo MSVC x64 Build Tools not found. No CPU fallback.
  exit /b 1
)
call "%GPRMAX_VS%\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 exit /b 1
"%~dp0..\artifacts\local_checks\gprmax_v4_gpu_env\Scripts\python.exe" "%~dp0run_approved_batch2d.py" %*
exit /b %errorlevel%

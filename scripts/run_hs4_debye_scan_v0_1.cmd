@echo off
setlocal
call "E:\msvc2022bt\VC\Auxiliary\Build\vcvars64.bat"
if errorlevel 1 exit /b 1
set "PATH=C:\cuda133\bin;C:\cuda133\bin\x64;%PATH%"
set "TEMP=E:\automation_djh\artifacts\local_checks\tmp"
set "TMP=E:\automation_djh\artifacts\local_checks\tmp"
set "MPLCONFIGDIR=E:\automation_djh\artifacts\local_checks\tmp\mpl"
"E:\automation_djh\artifacts\local_checks\gprmax_v4_gpu_env\Scripts\python.exe" "%~dp0hs4_debye_scan_v0_1.py" %*
exit /b %errorlevel%

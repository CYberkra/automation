@echo off
setlocal
call "E:\sisual stdio 2022\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 exit /b 1
rem Reuse existing ASCII aliases required by nvcc/MSVC preprocessing.
set "PATH=C:\cuda118\bin;%PATH%"
"D:\gprmax_v4_gpu_env\Scripts\python.exe" "%~dp0run_hs4t2d_relief_series.py" %*
exit /b %errorlevel%

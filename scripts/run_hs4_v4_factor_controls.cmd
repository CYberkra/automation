@echo off
setlocal
call "E:\sisual stdio 2022\VC\Auxiliary\Build\vcvars64.bat"
if errorlevel 1 exit /b 1
set "PATH=C:\cuda118\bin;%PATH%"
"D:\gprmax_v4_gpu_env\Scripts\python.exe" "%~dp0hs4_v4_factor_controls.py" %*
exit /b %errorlevel%

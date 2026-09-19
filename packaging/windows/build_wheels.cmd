@echo off
setlocal
rem Run from a conda-enabled Windows Command Prompt. Find VS automatically if
rem this is not already an x64 developer prompt.
where cl >nul 2>nul
if errorlevel 1 (
    call :setup_vs
    if errorlevel 1 exit /b 1
)
python "%~dp0build_wheels.py" %*
exit /b %errorlevel%

:setup_vs
set "MSANI_VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
if not exist "%MSANI_VSWHERE%" (
    echo ERROR: VS Build Tools not found. Use an x64 Native Tools Command Prompt.
    exit /b 1
)
for /f "usebackq tokens=*" %%i in (`"%MSANI_VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "MSANI_VSROOT=%%i"
if not defined MSANI_VSROOT (
    echo ERROR: Install the Visual Studio C++ build tools workload.
    exit /b 1
)
call "%MSANI_VSROOT%\Common7\Tools\VsDevCmd.bat" -arch=x64 -host_arch=x64
exit /b %errorlevel%

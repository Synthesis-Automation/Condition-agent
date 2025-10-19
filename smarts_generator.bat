@echo off
REM Reaction SMARTS Generator - Quick Launch Script
REM This script provides easy access to the SMARTS pattern generator CLI

setlocal

set "SCRIPT_DIR=%~dp0"
set "PROJECT_ROOT=%SCRIPT_DIR%"

REM Check if virtual environment exists
if exist "%PROJECT_ROOT%.venv\Scripts\python.exe" (
    set "PYTHON=%PROJECT_ROOT%.venv\Scripts\python.exe"
) else (
    set "PYTHON=python"
)

REM Show usage if no arguments
if "%~1"=="" (
    echo.
    echo 🧪 Reaction SMARTS Generator
    echo =============================
    echo.
    echo Usage:
    echo   smarts_generator.bat --interactive          Interactive mode
    echo   smarts_generator.bat --reaction "SMILES"    Single reaction
    echo   smarts_generator.bat --batch input.txt      Batch processing
    echo   smarts_generator.bat --check-rdkit          Check RDKit installation
    echo.
    echo Examples:
    echo   smarts_generator.bat --interactive
    echo   smarts_generator.bat --reaction "CCCCI.B2pin2>>CCCBpin" --output out.json
    echo   smarts_generator.bat --batch reactions.txt
    echo.
    goto :end
)

REM Run the CLI tool
"%PYTHON%" -m chemtools.protocol.smarts_generator_cli %*

:end
endlocal

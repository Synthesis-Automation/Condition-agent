@echo off
REM Windows batch script to run the protocol CLI tester
REM Sets up the Python path and runs the test script

REM Set UTF-8 encoding for emojis
chcp 65001 >nul

cd /d %~dp0
set PYTHONPATH=%~dp0
python tests\test_protocol_cli.py %*

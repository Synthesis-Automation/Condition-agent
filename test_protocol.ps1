# PowerShell script to run the protocol CLI tester
# Sets up the Python path and runs the test script

# Set UTF-8 encoding for emojis
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
$env:PYTHONIOENCODING = "utf-8"

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $scriptDir
$env:PYTHONPATH = $scriptDir

python tests\test_protocol_cli.py @args

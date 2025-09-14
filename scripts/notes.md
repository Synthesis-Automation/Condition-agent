Find PID: netstat -ano -p tcp | Select-String ':7860'
Inspect: Get-Process -Id <PID> | Select-Object Id, ProcessName, Path
Kill: Stop-Process -Id <PID> -Force


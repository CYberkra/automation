[Console]::OutputEncoding = [Text.Encoding]::UTF8
Write-Host "=== A. TrustedInstaller / TiWorker now ==="
Get-Process TrustedInstaller, TiWorker -ErrorAction SilentlyContinue |
  Select-Object Id, ProcessName, StartTime, CPU | Format-Table -AutoSize
Get-Service TrustedInstaller | Format-Table Status, Name -AutoSize

Write-Host "=== B. CBS boot-time phase (01:22-01:30) ==="
$lines = Get-Content C:\Windows\Logs\CBS\CBS.log -ErrorAction SilentlyContinue |
  Where-Object { $_ -match '2026-10-02 01:(2[2-9]|30)' }
Write-Host ("lines in boot window: " + ($lines | Measure-Object).Count)
$lines | Select-String -Pattern 'Startup|Failed|Error|0x8|pending' | Select-Object -First 25 | ForEach-Object { $_.Line }

Write-Host "=== C. CBS last 15 lines ==="
Get-Content C:\Windows\Logs\CBS\CBS.log -Tail 15 -ErrorAction SilentlyContinue

Write-Host "=== D. Setup event log (last 20) ==="
Get-WinEvent -FilterHashtable @{LogName='Setup'} -MaxEvents 20 -ErrorAction SilentlyContinue |
  Select-Object TimeCreated, Id, LevelDisplayName,
    @{N='Msg';E={ ($_.Message -split "`n")[0].Substring(0, [Math]::Min(140, ($_.Message -split "`n")[0].Length)) }} |
  Format-Table -AutoSize -Wrap

Write-Host "=== E. any TiWorker dumps ==="
Get-ChildItem C:\Windows\Logs\CBS\ -Filter *.dmp -ErrorAction SilentlyContinue | Select-Object Name, Length, LastWriteTime
Get-ChildItem C:\Windows\Logs\MoSetup -ErrorAction SilentlyContinue | Select-Object Name, LastWriteTime

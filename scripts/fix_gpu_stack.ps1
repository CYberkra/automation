$log = 'E:\automation_djh\fix_gpu_stack.log'
function W($m) { $line = ("[{0}] {1}" -f (Get-Date -Format 'yyyy-MM-dd HH:mm:ss'), $m); Add-Content -Path $log -Value $line; Write-Host $line }

W '=== fix_gpu_stack start ==='
W '--- step 1: sfc /scannow ---'
$p = Start-Process sfc.exe -ArgumentList '/scannow' -Wait -NoNewWindow -PassThru -RedirectStandardOutput "$log.sfc.out" -RedirectStandardError "$log.sfc.err"
W ("sfc exit code: " + $p.ExitCode)
Get-Content "$log.sfc.out" -Tail 8 | ForEach-Object { W $_ }

W '--- step 2: DISM RestoreHealth ---'
$p2 = Start-Process DISM.exe -ArgumentList '/Online /Cleanup-Image /RestoreHealth' -Wait -NoNewWindow -PassThru
W ("DISM exit code: " + $p2.ExitCode)

W '--- step 3: post-checks ---'
W ("pending.xml exists: " + (Test-Path C:\Windows\WinSxS\pending.xml))
W ("CBS RebootPending: " + (Test-Path 'HKLM:\SOFTWARE\Microsoft\Windows\CurrentVersion\Component Based Servicing\RebootPending'))
W '=== fix_gpu_stack end ==='

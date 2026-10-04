$p = Start-Process -FilePath 'E:\automation_djh\automation_repo\scripts\run_benchmark3d_r2_gpu.cmd' -ArgumentList '--execute' -WorkingDirectory 'E:\automation_djh\automation_repo' -WindowStyle Hidden -RedirectStandardOutput 'E:\automation_djh\batch_2d_stdout.log' -RedirectStandardError 'E:\automation_djh\batch_2d_stderr.log' -PassThru
Write-Host "2d batch launched pid $($p.Id)"
Start-Sleep -Seconds 25
Get-Content 'E:\automation_djh\batch_2d_stdout.log' -ErrorAction SilentlyContinue | Select-Object -First 3 | ForEach-Object { $_.Substring(0, [Math]::Min(90, $_.Length)) }

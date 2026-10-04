$p = Start-Process -FilePath 'E:\automation_djh\automation_repo\scripts\restore_t11_benchmark3d_gpu.cmd' -WorkingDirectory 'E:\automation_djh\automation_repo' -WindowStyle Hidden -RedirectStandardOutput 'E:\automation_djh\restore_t11_stdout.log' -RedirectStandardError 'E:\automation_djh\restore_t11_stderr.log' -PassThru
Write-Host "restore launched pid $($p.Id)"
Start-Sleep -Seconds 30
if (Test-Path 'E:\automation_djh\restore_t11_stdout.log') {
  Get-Content 'E:\automation_djh\restore_t11_stdout.log' | ForEach-Object { $_.Substring(0, [Math]::Min(110, $_.Length)) }
}

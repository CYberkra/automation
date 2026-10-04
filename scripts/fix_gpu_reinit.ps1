$log = 'E:\automation_djh\fix_gpu_stack.log'
function W($m) { $line = ("[{0}] {1}" -f (Get-Date -Format 'yyyy-MM-dd HH:mm:ss'), $m); Add-Content -Path $log -Value $line }

W '=== attempt 2: gpu device reinit ==='
$dev = Get-PnpDevice -Class Display | Where-Object { $_.FriendlyName -match 'RTX 4090' }
if (-not $dev) { W 'RTX 4090 device not found!'; exit 1 }
W ("device: " + $dev.InstanceId)

Disable-PnpDevice -InstanceId $dev.InstanceId -Confirm:$false -ErrorAction Continue
W ('disable returned: ' + $LASTEXITCODE)
Start-Sleep -Seconds 5
Enable-PnpDevice -InstanceId $dev.InstanceId -Confirm:$false -ErrorAction Continue
W ('enable returned: ' + $LASTEXITCODE)
Start-Sleep -Seconds 10

W '--- nvidia-smi test after reinit ---'
& nvidia-smi --query-gpu=name,driver_version,memory.used,memory.total --format=csv | Out-File -Append -FilePath $log -Encoding utf8
W ('nvidia-smi exit: ' + $LASTEXITCODE)
W '=== attempt 2 end ==='

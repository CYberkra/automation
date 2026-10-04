$log = 'E:\automation_djh\fix_gpu_stack.log'
$settings = 'HKLM:\SOFTWARE\Microsoft\WindowsUpdate\UX\Settings'
$start = (Get-Date).ToString('yyyy-MM-ddTHH:mm:ssZ')
$expiry = (Get-Date).AddDays(7).ToString('yyyy-MM-ddTHH:mm:ssZ')
Set-ItemProperty -Path $settings -Name 'PauseFeatureUpdatesStartTime' -Value $start
Set-ItemProperty -Path $settings -Name 'PauseQualityUpdatesStartTime' -Value $start
Set-ItemProperty -Path $settings -Name 'PauseFeatureUpdatesEndTime' -Value $expiry
Set-ItemProperty -Path $settings -Name 'PauseQualityUpdatesEndTime' -Value $expiry
Set-ItemProperty -Path $settings -Name 'PauseUpdatesExpiryTime' -Value $expiry
Add-Content $log "[$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')] WU paused until $expiry"

[Console]::OutputEncoding = [Text.Encoding]::UTF8
Write-Host "=== 1. NVIDIA driver version ==="
Get-CimInstance Win32_PnPSignedDriver -Filter "DeviceName like '%RTX 4090%'" |
  Select-Object DeviceName, DriverVersion, DriverDate | Format-List

Write-Host "=== 2. nvlddmkm service ==="
Get-Service nvlddmkm -ErrorAction SilentlyContinue | Format-Table Status, Name, StartType -AutoSize

Write-Host "=== 3. NVIDIA services ==="
Get-Service | Where-Object { $_.Name -match 'nvidia|NvContainer|NVDisplay' } |
  Select-Object Status, Name, DisplayName | Format-Table -AutoSize

Write-Host "=== 4. Recent nvlddmkm / Display errors ==="
Get-WinEvent -FilterHashtable @{LogName='System'; StartTime=(Get-Date).AddHours(-6)} -ErrorAction SilentlyContinue |
  Where-Object { $_.ProviderName -match 'nvlddmkm|Display|nvidia' -or $_.Message -match 'nvlddmkm|NVIDIA' } |
  Select-Object -First 15 TimeCreated, Id, LevelDisplayName, ProviderName,
    @{N='Msg';E={ ($_.Message -split "`n")[0].Substring(0, [Math]::Min(120, ($_.Message -split "`n")[0].Length)) }} |
  Format-Table -AutoSize -Wrap

Write-Host "=== 5. CBS.log recent errors ==="
$tail = Get-Content C:\Windows\Logs\CBS\CBS.log -Tail 400 -ErrorAction SilentlyContinue
$tail | Select-String -Pattern 'Error|Failed|0x8' | Select-Object -First 15 | ForEach-Object { $_.Line }

Write-Host "=== 6. RebootRequired flags ==="
$rr = 'HKLM:\SOFTWARE\Microsoft\Windows\CurrentVersion\WindowsUpdate\Auto Update\RebootRequired'
if (Test-Path $rr) { Get-ItemProperty $rr | Format-List } else { Write-Host "no WU RebootRequired key" }
$rr2 = 'HKLM:\SOFTWARE\Microsoft\Windows\CurrentVersion\Component Based Servicing\RebootPending'
if (Test-Path $rr2) { Write-Host "CBS RebootPending EXISTS" } else { Write-Host "no CBS RebootPending" }

Write-Host "=== 7. Installed LCU version ==="
Get-ItemProperty 'HKLM:\SOFTWARE\Microsoft\Windows NT\CurrentVersion' |
  Select-Object DisplayVersion, UBR, BuildLabEx | Format-List

Write-Host "=== 8. pending.xml transaction count ==="
[xml]$px = Get-Content C:\Windows\WinSxS\pending.xml
Write-Host ("transactions: " + $px.PendingTransaction.Transactions.Transaction.Count)

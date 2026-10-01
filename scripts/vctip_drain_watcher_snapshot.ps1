# vctip drain watcher (25 min) for snapshot_c1mx_v1 (2026-10-01)
# Same operational precedent as the b2 pilot watcher: kills MSVC telemetry
# orphans only. Never touches solver/nvcc/cl.
$log = 'E:\automation_djh\automation_repo\artifacts\research_checks\2026-10-01_vctip_drain_watcher.log'
$deadline = (Get-Date).AddMinutes(25)
while ((Get-Date) -lt $deadline) {
    $procs = Get-Process vctip -ErrorAction SilentlyContinue
    foreach ($p in $procs) {
        try {
            Stop-Process -Id $p.Id -Force -ErrorAction Stop
            Add-Content $log ("{0} killed vctip pid {1} (started {2}) [snapshot_c1mx_v1]" -f (Get-Date -Format 'HH:mm:ss'), $p.Id, $p.StartTime)
        } catch {}
    }
    Start-Sleep -Seconds 12
}
Add-Content $log ("{0} watcher (snapshot_c1mx_v1) deadline reached, exiting" -f (Get-Date -Format 'HH:mm:ss'))

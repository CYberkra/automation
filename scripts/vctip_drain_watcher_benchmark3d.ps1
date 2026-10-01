# vctip drain watcher for benchmark3d_r2_co (2026-10-01)
# Operational intervention only (registry precedent): kills MSVC telemetry orphans
# that keep supervision jobs non-empty after solver exit. Never touches solver/nvcc/cl.
$log = 'E:\automation_djh\automation_repo\artifacts\research_checks\2026-10-01_vctip_drain_watcher_benchmark3d.log'
$deadline = (Get-Date).AddHours(6)
while ((Get-Date) -lt $deadline) {
    $procs = Get-Process vctip -ErrorAction SilentlyContinue
    foreach ($p in $procs) {
        try {
            Stop-Process -Id $p.Id -Force -ErrorAction Stop
            Add-Content $log ("{0} killed vctip pid {1} (started {2})" -f (Get-Date -Format 'HH:mm:ss'), $p.Id, $p.StartTime)
        } catch {}
    }
    Start-Sleep -Seconds 12
}
Add-Content $log ("{0} watcher deadline reached, exiting" -f (Get-Date -Format 'HH:mm:ss'))

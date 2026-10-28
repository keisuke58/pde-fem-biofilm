# start_remote_control.ps1 -- IKMHIWI03: keep `claude remote-control` running for the Claude Code app (21 Oct 2026).
#
#   .\scripts\start_remote_control.ps1            # start it now in its own window, restarted whenever it exits
#   .\scripts\start_remote_control.ps1 -Install   # also start it at every logon (scheduled task, current user, no admin)
#   .\scripts\start_remote_control.ps1 -Uninstall # remove that scheduled task
#
# Its own window survives a closed VS Code tunnel tab; the loop restarts it after a crash or a lost connection
# (wait 30 s, log line in F:\remote_control.log); the scheduled task brings it back after a reboot (Windows updates,
# 14 Oct). ANSYS chains run detached (run_chain.ps1) either way. The PC must not sleep: check with `powercfg /a`.
param([string]$Repo = "C:\Users\nishioka\git\pde-fem-biofilm",
      [string]$Log = "F:\remote_control.log",
      [switch]$Install, [switch]$Uninstall, [switch]$Loop)

$task = "ClaudeRemoteControl"
if ($Uninstall) {
    Unregister-ScheduledTask -TaskName $task -Confirm:$false -ErrorAction SilentlyContinue
    Write-Host "removed scheduled task $task"; return
}
if ($Loop) {                                   # the restart loop (runs inside the window started below)
    Set-Location $Repo
    while ($true) {
        Add-Content $Log "$(Get-Date -Format s) start claude remote-control"
        claude remote-control
        Add-Content $Log "$(Get-Date -Format s) exited with code $LASTEXITCODE, restart in 30 s"
        Start-Sleep -Seconds 30
    }
}
$self = $MyInvocation.MyCommand.Path
if ($Install) {
    $act = New-ScheduledTaskAction -Execute "powershell.exe" `
        -Argument "-NoExit -ExecutionPolicy Bypass -File `"$self`" -Loop -Repo `"$Repo`" -Log `"$Log`""
    $trg = New-ScheduledTaskTrigger -AtLogOn -User $env:USERNAME
    $set = New-ScheduledTaskSettingsSet -AllowStartIfOnBatteries -DontStopIfGoingOnBatteries `
        -ExecutionTimeLimit ([TimeSpan]::Zero) -RestartCount 3 -RestartInterval (New-TimeSpan -Minutes 5)
    Register-ScheduledTask -TaskName $task -Action $act -Trigger $trg -Settings $set -Force | Out-Null
    Write-Host "scheduled task ${task}: starts at every logon of $env:USERNAME"
}
$running = Get-CimInstance Win32_Process -Filter "Name='powershell.exe'" |
    Where-Object { $_.CommandLine -like "*start_remote_control.ps1*-Loop*" }
if ($running) { Write-Host "already running (pid $($running.ProcessId -join ', '))"; return }
Start-Process powershell -ArgumentList '-NoExit', '-ExecutionPolicy', 'Bypass', '-File', "`"$self`"", '-Loop',
    '-Repo', "`"$Repo`"", '-Log', "`"$Log`""
Write-Host "started (log $Log)"

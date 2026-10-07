<#
wp3_auto.ps1 -- the unattended WP3 sequence on IKMHIWI03 (RUN_WP3_IKMHIWI03.md),
in a process detached from the calling shell (as run_chain.ps1):

  .\ansys_usermat\apdl\wp3_auto.ps1 [-AcceptLog F:\biofilm_upf_front\_chain_front_accept.log]
      [-ALog F:\biofilm_upf_wired\_chain_wp3_A.log]

1. waits until the acceptance chain of the front build has ended;
2. runs judge_front_accept.py, writes its report to
   ansys_usermat/apdl/results/2026-10-wp3/front_accept.txt, commits and pushes it;
3. if the gate passes: starts the block B chain (F:\biofilm_upf_front,
   _wp3_B_runs.txt from make_wp3_decks.py) behind the block A chain (-WaitFor
   its log), with JSON export and a push every 3 runs. If it fails: nothing
   more is started; the report says why.
Block A is started separately (run_chain.ps1 -Name wp3_A ..., it does not depend
on the gate). Log: F:\biofilm_upf_front\_wp3_auto.log.
#>
param(
    [string]$AcceptLog = 'F:\biofilm_upf_front\_chain_front_accept.log',
    [string]$ALog = 'F:\biofilm_upf_wired\_chain_wp3_A.log',
    [string]$Export = 'ansys_usermat\apdl\results\2026-10-wp3',
    [switch]$Worker
)
$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$front = 'F:\biofilm_upf_front'
$log = Join-Path $front '_wp3_auto.log'

if (-not $Worker) {
    $cmd = "powershell.exe -NoProfile -ExecutionPolicy Bypass -WindowStyle Hidden -File `"$PSCommandPath`" -Worker " +
           "-AcceptLog `"$AcceptLog`" -ALog `"$ALog`" -Export `"$Export`""
    $r = Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd; CurrentDirectory = $repo }
    if ($r.ReturnValue -ne 0) { throw "could not start (Win32_Process.Create returned $($r.ReturnValue))" }
    "wp3_auto started, PID $($r.ProcessId); log $log"
    exit 0
}

# PS 5.1: with Stop, git's progress lines on stderr (2>$null) would end the worker
$ErrorActionPreference = 'Continue'
function L($m) {
    $line = "$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') $m"
    for ($i = 0; $i -lt 5; $i++) { try { [IO.File]::AppendAllText($log, $line + "`r`n"); break } catch { Start-Sleep -Milliseconds 300 } }
}
try {
    Set-Location $repo
    $py = 'C:\Users\nishioka\AppData\Local\Programs\Python\Python312\python.exe'
    $env:Path = "C:\Users\nishioka\git\cmd;C:\Users\nishioka\git\mingw64\bin;" + $env:Path
    L "waiting for $AcceptLog"
    while (-not ((Test-Path $AcceptLog) -and ((Get-Content $AcceptLog -Raw) -match '(?m)^(\S+ )?\S+ \w+ end\s*$|EXCEPTION'))) { Start-Sleep 60 }
    L 'acceptance chain ended'

    $rep = Join-Path $repo "$Export\front_accept.txt"
    $o = & $py (Join-Path $repo 'ansys_usermat\apdl\judge_front_accept.py') --workdir $front --report $rep 2>&1
    $rc = $LASTEXITCODE
    $o | ForEach-Object { L "  judge: $_" }
    $verdict = if ($rc -eq 0) { 'gate passed' } else { 'gate failed' }

    $rel = ($Export -replace '\\', '/') + '/front_accept.txt'
    $msg = "Front-term build: acceptance on ds_fig7h ($verdict, automatic)`n`n" + (($o | ForEach-Object { "$_" }) -join "`n")
    $c = & (Join-Path $repo 'commit.ps1') -Files $rel -Message $msg 2>&1; $c | ForEach-Object { L "  commit: $_" }
    $br = (& git rev-parse --abbrev-ref HEAD).Trim()
    & git fetch origin $br 2>$null
    & git rebase --autostash FETCH_HEAD *> $null
    $p = & (Join-Path $repo 'push.ps1') -Branch $br 2>&1; $p | ForEach-Object { L "  push: $_" }

    if ($rc -ne 0) { L 'gate failed: block B not started'; L 'wp3_auto end'; exit 0 }

    $runs = Get-Content (Join-Path $front '_wp3_B_runs.txt') | Where-Object { $_.Trim() }
    $o = & powershell.exe -NoProfile -ExecutionPolicy Bypass -File (Join-Path $repo 'ansys_usermat\apdl\run_chain.ps1') `
        -Name wp3_B -Runs ($runs -join ',') -WorkDir $front -WaitFor $ALog -Export $Export -Push -PushEvery 3 2>&1
    $o | ForEach-Object { L "  run_chain: $_" }
    L 'wp3_auto end'
} catch { L "EXCEPTION: $_" }

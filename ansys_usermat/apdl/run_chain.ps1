<#
run_chain.ps1 -- run a list of wired decks one after another in a process that is
detached from the calling shell (a Claude session restart or a closed terminal
does not stop it), with a log, and optionally export the results as JSON and
push them when the chain ends.

  .\ansys_usermat\apdl\run_chain.ps1 -Name 1012 -Runs ds_a:2sp_case6, ds16_b:2sp_case3, ds24_c::300 `
      [-WaitFor F:\biofilm_upf_wired\_chain_other.log] [-Export results\2026-10-12_ansys] [-Push]

-Runs      deck names without .dat (in -WorkDir), each optionally ":case" (material
           server case for run_wired.ps1 -Case) and ":minutes" (timeout, default 150)
-WaitFor   a log of another chain or script: start only after it contains
           "<word> end" (chain end, after18 end, ...) or "EXCEPTION", and no ANSYS
           process is left (one ANSYS run at a time: one material-server port)
-Export    after the last run, export_runs_json.py writes <run>.json into this
           repo-relative folder; with -Push the JSON files are committed
           (commit.ps1, only those files) and the current branch is pushed (push.ps1)
-PushEvery n   with -Export -Push: export, commit and push after every n finished runs
           as well, so results show up on GitHub while a long chain is still running
-SkipDone  skip a run whose all_stress_<run>.csv is newer than its deck (to restart a
           chain after a reboot without repeating finished runs)

A run that fails within 3 minutes without writing any result (licence server
not reachable, network gone) is retried up to 3 times, 30 minutes apart; a run
that fails after solving, or stops with a FATAL message that is not about the
licence (e.g. out of memory), is not retried.

The caller returns at once and prints the PID and the log path. Progress:
<WorkDir>\_chain_<Name>.log (START / DONE rc= per run, run_wired's summary);
<WorkDir>\_chain_<Name>.transcript.txt shows the cause if the chain itself stops.
Each run goes through run_wired.ps1 in its own powershell process, so its
exit cannot end the chain. Found the hard way on 5 Oct: a chain started as a
background shell of the session died with the session; one started with
Start-Process but calling run_wired.ps1 in-process also stopped mid-way.
#>
param(
    [Parameter(Mandatory)] [string]$Name,
    [Parameter(Mandatory)] [string[]]$Runs,
    [string]$WaitFor = '',
    [string]$Export = '',
    [switch]$Push,
    [int]$PushEvery = 0,
    [switch]$SkipDone,
    [string]$WorkDir = 'F:\biofilm_upf_wired',
    [switch]$DryRun,          # log START/DONE without running ANSYS (to test the chain itself)
    [switch]$Worker,
    [string]$Branch = ''
)
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$log = Join-Path $WorkDir "_chain_$Name.log"

if (-not $Worker) {
    foreach ($r in $Runs) {
        $d = ($r -split ':')[0]
        if (-not (Test-Path (Join-Path $WorkDir "$d.dat"))) { throw "deck not found: $d.dat" }
        if (Test-Path (Join-Path $WorkDir "$d.lock")) { throw "$d.lock exists in $WorkDir (an aborted run?)" }
    }
    $env:Path = "C:\Users\nishioka\git\cmd;C:\Users\nishioka\git\mingw64\bin;" + $env:Path
    $br = (& git -C $repo rev-parse --abbrev-ref HEAD).Trim()
    $args_ = @('-NoProfile', '-ExecutionPolicy', 'Bypass', '-WindowStyle', 'Hidden', '-File', "`"$PSCommandPath`"",
               '-Worker', '-Name', $Name, '-Runs', ($Runs -join ','), '-WorkDir', "`"$WorkDir`"", '-Branch', $br)
    if ($WaitFor) { $args_ += @('-WaitFor', "`"$WaitFor`"") }
    if ($Export) { $args_ += @('-Export', "`"$Export`"") }
    if ($Push) { $args_ += '-Push' }
    if ($PushEvery) { $args_ += @('-PushEvery', "$PushEvery") }
    if ($SkipDone) { $args_ += '-SkipDone' }
    if ($DryRun) { $args_ += '-DryRun' }
    $r = Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{
        CommandLine = "powershell.exe $($args_ -join ' ')"; CurrentDirectory = $repo }
    if ($r.ReturnValue -ne 0) { throw "could not start the chain (Win32_Process.Create returned $($r.ReturnValue))" }
    "chain $Name started, PID $($r.ProcessId), $($Runs.Count) runs; log $log"
    exit 0
}

# ---------------- worker ----------------
Start-Transcript -Path (Join-Path $WorkDir "_chain_$Name.transcript.txt") -Append | Out-Null
function L($m) {
    $line = "$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') $m"
    for ($i = 0; $i -lt 5; $i++) { try { [IO.File]::AppendAllText($log, $line + "`r`n"); break } catch { Start-Sleep -Milliseconds 300 } }
}
$Runs = @($Runs -join ',' -split ',' | Where-Object { $_ })        # one comma-joined string from the launcher
function ExportPush($list, $tag) {
    $list = @($list)                       # one run arrives as a bare string, and @string splats per character
    if (-not $Export -or -not $list) { return }
    . (Join-Path $repo 'dev-env.ps1') | Out-Null
    $o = & python (Join-Path $repo 'ansys_usermat\apdl\export_runs_json.py') (Join-Path $repo $Export) @list 2>&1
    $o | ForEach-Object { L "  export: $_" }
    if (-not $Push) { return }
    $env:Path = "C:\Users\nishioka\git\cmd;C:\Users\nishioka\git\mingw64\bin;" + $env:Path
    $rel = $Export -replace '\\', '/'
    $files = $list | Where-Object { Test-Path (Join-Path $repo "$Export\$_.json") } | ForEach-Object { "$rel/$_.json" }
    if (-not $files) { L '  push: nothing to commit'; return }
    $summary = (Get-Content $log | Where-Object { $_ -match 'DONE|ERROR   MESSAGES|FATAL' } | Select-Object -Last 40) -join "`n"
    $msg = "ANSYS run chain $Name as JSON ($tag, automatic export)`n`nNo analysis yet. Runs: $($list -join ', ').`n`nChain log (last lines):`n$summary"
    $o = & (Join-Path $repo 'commit.ps1') -Files $files -Message $msg 2>&1; $o | ForEach-Object { L "  commit: $_" }
    & git fetch origin $Branch 2>$null
    & git rebase --autostash FETCH_HEAD *> $null
    $o = & (Join-Path $repo 'push.ps1') -Branch $Branch 2>&1; $o | ForEach-Object { L "  push: $_" }
}
try {
    Set-Location $repo
    L "chain start: $($Runs -join ', ')"
    if ($WaitFor) {
        L "waiting for $WaitFor"
        # "HH:mm:ss after18 end" (scripts) and "yyyy-MM-dd HH:mm:ss chain end" (chain logs)
        while (-not ((Test-Path $WaitFor) -and ((Get-Content $WaitFor -Raw) -match '(?m)^(\S+ )?\S+ \w+ end\s*$|EXCEPTION'))) { Start-Sleep 30 }
        while (Get-Process ANSYS -ErrorAction SilentlyContinue) { Start-Sleep 30 }
        Start-Sleep 10
    }
    $done = @(); $batch = @(); $k = 0
    foreach ($r in $Runs) {
        $k++
        $p = $r -split ':'
        $d = $p[0]; $case = if ($p.Count -gt 1) { $p[1] } else { '' }
        $tmo = if ($p.Count -gt 2 -and $p[2]) { $p[2] } else { '150' }
        $res = Join-Path $WorkDir "all_stress_$d.csv"
        if ($SkipDone -and (Test-Path $res) -and
            (Get-Item $res).LastWriteTime -gt (Get-Item (Join-Path $WorkDir "$d.dat")).LastWriteTime) {
            L "SKIP $d (finished earlier)"; continue
        }
        if ($DryRun) { L "START $d ($k/$($Runs.Count))"; Start-Sleep 2; L "DONE $d rc=0 (dry run)" }
        for ($try = 1; $try -le 4 -and -not $DryRun; $try++) {
            L "START $d ($k/$($Runs.Count))$(if ($case) { " case $case" }) timeout $tmo min$(if ($try -gt 1) { " retry $($try - 1)" })"
            $rc = -1; $t0 = Get-Date
            $lock = Join-Path $WorkDir "$d.lock"
            if ((Test-Path $lock) -and -not (Get-Process ANSYS -ErrorAction SilentlyContinue)) { Remove-Item $lock -Force; L "  removed stale $d.lock" }
            try {
                $o = & powershell.exe -NoProfile -ExecutionPolicy Bypass -File (Join-Path $repo 'ansys_usermat\apdl\run_wired.ps1') @(
                    '-Deck', "$d.dat", '-TimeoutMin', $tmo, '-WorkDir', $WorkDir) @(if ($case) { '-Case', $case }) 2>&1
                $rc = $LASTEXITCODE
                $o | ForEach-Object { L "  $_" }
            } catch { L "  exception: $_" }
            L "DONE $d rc=$rc"
            $quick = ((Get-Date) - $t0).TotalMinutes -lt 3
            $wrote = (Test-Path $res) -and (Get-Item $res).LastWriteTime -gt $t0
            $outf = Join-Path $WorkDir "out_$d.txt"
            $fm = if (Test-Path $outf) { @(Select-String $outf -Pattern '\*\*\* FATAL \*\*\*' -Context 0, 3) } else { @() }
            # the banner of every output says "LICENSORS": look at the FATAL message itself only
            $fatal = $fm.Count -gt 0 -and -not ($fm | Where-Object { ($_.Line + ' ' + ($_.Context.PostContext -join ' ')) -match 'licen[cs]' })
            if ($rc -eq 0 -or $wrote -or -not $quick -or $fatal -or $try -eq 4) { break }
            L "  failed within 3 min without results: waiting 30 min before retrying (licence/network?)"
            Start-Sleep 1800
        }
        $done += $d; $batch += $d
        if ($PushEvery -gt 0 -and $batch.Count -ge $PushEvery) { ExportPush $batch "$($done.Count) of $($Runs.Count)"; $batch = @() }
    }
    ExportPush $(if ($PushEvery -gt 0) { $batch } else { $done }) 'end'
    L 'chain end'
} catch {
    L "CHAIN EXCEPTION: $_"
} finally {
    Stop-Transcript | Out-Null
}

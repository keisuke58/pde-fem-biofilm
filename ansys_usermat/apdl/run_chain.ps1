<#
run_chain.ps1 -- run a list of wired decks one after another in a process that is
detached from the calling shell (a Claude session restart or a closed terminal
does not stop it), with a log, and optionally export the results as JSON and
push them when the chain ends.

  .\ansys_usermat\apdl\run_chain.ps1 -Name 1012 -Runs ds_a:2sp_case6, ds16_b:2sp_case3, ds24_c::300 `
      [-WaitFor F:\biofilm_upf_wired\_chain_other.log] [-Export results\2026-10-12_ansys] [-Push]

-Runs      deck names without .dat (in -WorkDir), each optionally ":case" (material
           server case for run_wired.ps1 -Case) and ":minutes" (timeout, default 150)
-WaitFor   a log of another chain: start only after it contains "chain end" or
           "CHAIN EXCEPTION" (one ANSYS run at a time: one material-server port)
-Export    after the last run, export_runs_json.py writes <run>.json into this
           repo-relative folder; with -Push the JSON files are committed
           (commit.ps1, only those files) and the current branch is pushed (push.ps1)

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
try {
    Set-Location $repo
    L "chain start: $($Runs -join ', ')"
    if ($WaitFor) {
        L "waiting for $WaitFor"
        while (-not ((Test-Path $WaitFor) -and ((Get-Content $WaitFor -Raw) -match 'chain end|CHAIN EXCEPTION'))) { Start-Sleep 30 }
        Start-Sleep 10
    }
    $done = @()
    foreach ($r in $Runs) {
        $p = $r -split ':'
        $d = $p[0]; $case = if ($p.Count -gt 1) { $p[1] } else { '' }
        $tmo = if ($p.Count -gt 2 -and $p[2]) { $p[2] } else { '150' }
        L "START $d$(if ($case) { " case $case" }) timeout $tmo min"
        $rc = -1
        if ($DryRun) { Start-Sleep 2; L "DONE $d rc=0 (dry run)"; $done += $d; continue }
        try {
            $o = & powershell.exe -NoProfile -ExecutionPolicy Bypass -File (Join-Path $repo 'ansys_usermat\apdl\run_wired.ps1') @(
                '-Deck', "$d.dat", '-TimeoutMin', $tmo, '-WorkDir', $WorkDir) @(if ($case) { '-Case', $case }) 2>&1
            $rc = $LASTEXITCODE
            $o | ForEach-Object { L "  $_" }
        } catch { L "  exception: $_" }
        L "DONE $d rc=$rc"
        $done += $d
    }
    if ($Export) {
        . (Join-Path $repo 'dev-env.ps1') | Out-Null
        $o = & python (Join-Path $repo 'ansys_usermat\apdl\export_runs_json.py') (Join-Path $repo $Export) @done 2>&1
        $o | ForEach-Object { L "  export: $_" }
        if ($Push) {
            $env:Path = "C:\Users\nishioka\git\cmd;C:\Users\nishioka\git\mingw64\bin;" + $env:Path
            $rel = $Export -replace '\\', '/'
            $files = $done | Where-Object { Test-Path (Join-Path $repo "$Export\$_.json") } | ForEach-Object { "$rel/$_.json" }
            $summary = (Get-Content $log | Where-Object { $_ -match 'DONE|ERROR   MESSAGES|FATAL' }) -join "`n"
            $msg = "ANSYS run chain $Name as JSON (automatic export)`n`nNo analysis yet. Runs: $($done -join ', ').`n`nChain log:`n$summary"
            $o = & (Join-Path $repo 'commit.ps1') -Files $files -Message $msg 2>&1; $o | ForEach-Object { L "  commit: $_" }
            & git fetch origin $Branch 2>$null
            & git rebase --autostash FETCH_HEAD *> $null
            $o = & (Join-Path $repo 'push.ps1') -Branch $Branch 2>&1; $o | ForEach-Object { L "  push: $_" }
        }
    }
    L 'chain end'
} catch {
    L "CHAIN EXCEPTION: $_"
} finally {
    Stop-Transcript | Out-Null
}

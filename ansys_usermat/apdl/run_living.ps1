<#
run_living.ps1 -- the coupling with the living share (A-19, 20 Oct 2026) in ANSYS: make the decks
(make_living_deck.py), run them one after another with the build in F:\biofilm_upf_living (usermat prop(28) = 11,
USSFin with LIVING_CONS / FRONT_GRAD / NUT2_TSTART), collect the fields and the point-model trace
(fields_living.py -> ansys_usermat/apdl/results/2026-10-living/) and, with -Push, commit and push the JSON.

  .\ansys_usermat\apdl\run_living.ps1 -Name chk8 -Runs "lv8_sep:8:--layout sep --front 0.1", ... [-Push]

-Runs: job:mesh:make_living_deck options. Detached from the session (WMI Win32_Process.Create), log
<WorkDir>\_living_<Name>.log; -Worker is the detached part. One ANSYS at a time, -np 1 (the partner's element hangs
under DMP; the point model reads its file at the first call). -After <log>: wait until that log ends with "end".
#>
param(
    [string]$Name = 'living',
    [string[]]$Runs = @(),
    [string]$WorkDir = 'F:\biofilm_upf_living',
    [string]$After = '',
    [switch]$Push,
    [switch]$Worker,
    [string]$Branch = '',
    [int]$Np = 1          # > 1: shared-memory parallel (-smp -np N), A-44; DMP hangs with the partner's element
)
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$env:Path = "C:\Users\nishioka\git\cmd;C:\Users\nishioka\git\mingw64\bin;" + $env:Path
$Runs = @($Runs | ForEach-Object { $_ -split ';' } | Where-Object { $_ })
if (-not $Worker) {
    foreach ($f in 'ANSYS.exe', 'base_w8_c6_g1_s015.dat', 'post_fields_living.mac', 'eco_fritsch_living.txt') {
        if (-not (Test-Path (Join-Path $WorkDir $f))) { throw "missing in ${WorkDir}: $f" }
    }
    $br = (& git -C $repo rev-parse --abbrev-ref HEAD).Trim()
    $runArg = ($Runs -join ';')
    $a = @('-NoProfile', '-ExecutionPolicy', 'Bypass', '-WindowStyle', 'Hidden', '-File', "`"$PSCommandPath`"", '-Worker',
           '-Name', $Name, '-WorkDir', "`"$WorkDir`"", '-Runs', "`"$runArg`"", '-Branch', $br)
    if ($After) { $a += @('-After', "`"$After`"") }
    if ($Push) { $a += '-Push' }
    if ($Np -gt 1) { $a += @('-Np', $Np) }
    $r = Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{
        CommandLine = "powershell.exe $($a -join ' ')"; CurrentDirectory = $repo }
    if ($r.ReturnValue -ne 0) { throw "could not start (Win32_Process.Create returned $($r.ReturnValue))" }
    Write-Output "run_living started, PID $($r.ProcessId); log $WorkDir\_living_$Name.log"
    return
}

$ErrorActionPreference = 'Continue'
Set-Location $repo
. (Join-Path $repo 'dev-env.ps1') | Out-Null
$py = (Get-Command python).Source
$log = Join-Path $WorkDir "_living_$Name.log"
function L($m) { [IO.File]::AppendAllText($log, "$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') $m`r`n") }
if ($After) {
    L "waiting for $After"
    while (-not ((Test-Path $After) -and ((Get-Content $After -Tail 1) -match ' end$'))) { Start-Sleep 60 }
}
L "start: $($Runs -join ' | ')"
$done = @()
$exe = "$env:AWP_ROOT222\ANSYS\bin\winx64\ANSYS222.exe"
foreach ($spec in $Runs) {
    $parts = $spec -split ':', 3
    $job, $mesh, $opts = $parts[0], $parts[1], $parts[2]
    # --eco FILE (a file name in WorkDir): this run's point-model file instead of eco_fritsch_living.txt
    $eco = 'eco_fritsch_living.txt'
    if ($opts -match '--eco\s+(\S+)') { $eco = $Matches[1]; $opts = ($opts -replace '--eco\s+\S+', '').Trim() }
    # --trace-elems a,b (build living8t and later, A-53): the traces only for these elements
    $env:BIOFILM_TRACE_ELEMS = ''
    if ($opts -match '--trace-elems\s+(\S+)') { $env:BIOFILM_TRACE_ELEMS = $Matches[1]; $opts = ($opts -replace '--trace-elems\s+\S+', '').Trim() }
    $base = Join-Path $WorkDir "base_w${mesh}_c6_g1_s015.dat"
    $deck = Join-Path $WorkDir "$job.dat"
    $o = & $py ansys_usermat\apdl\make_living_deck.py $base $deck @($opts -split ' ' | Where-Object { $_ }) 2>&1
    $o | ForEach-Object { L "  deck: $_" }
    if ($LASTEXITCODE -ne 0) { L "  ${job}: deck FAILED"; continue }
    Push-Location $WorkDir
    Remove-Item comp_trace.csv, fields_*.csv, all_stress.csv, nut_field.csv, "$job.lock" -ErrorAction SilentlyContinue
    Get-ChildItem -Filter "$job.*" | Where-Object { $_.Extension -notin '.dat' } | Remove-Item -ErrorAction SilentlyContinue
    $env:BIOFILM_ECO_CASE = Join-Path $WorkDir $eco
    $t0 = Get-Date
    L "START $job ($mesh^3: $opts; point model $eco)"
    # as run_calibrated.ps1: a hidden window and a poll (-NoNewWindow -Wait exits at once in the detached process)
    $par = if ($Np -gt 1) { @('-smp', '-np', "$Np") } else { @('-np', '1') }
    $p = Start-Process $exe -WorkingDirectory $WorkDir -ArgumentList (@('-b') + $par + @('-j', $job, '-custom',
        '.\ANSYS.exe', '-i', "$job.dat", '-o', "out_$job.txt")) -PassThru -WindowStyle Hidden
    $null = $p.Handle
    $limit = @{ '8' = 60; '16' = 600; '24' = 2400; '32' = 6000 }[($mesh -replace '\D', '')]   # 16q: quarter base deck (A-55)
    if ($opts -match '--deltim\s+([0-9.eE+-]+)') { $limit *= [Math]::Max(1.0, 0.00125 / [double]$Matches[1]) }   # smaller increments run longer
    while (-not $p.HasExited) {
        Start-Sleep 20
        if (((Get-Date) - $t0).TotalMinutes -gt $limit) {
            Get-CimInstance Win32_Process -Filter "ParentProcessId=$($p.Id)" | ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }
            Stop-Process -Id $p.Id -Force -ErrorAction SilentlyContinue; L "  killed after $limit min"; break
        }
    }
    $err = (Select-String "out_$job.txt" -Pattern 'NUMBER OF ERROR\s+MESSAGES ENCOUNTERED=\s*(\d+)' |
            Select-Object -Last 1 | ForEach-Object { $_.Matches[0].Groups[1].Value })
    L ("DONE $job exit $($p.ExitCode), errors $err, {0:N1} min" -f ((Get-Date) - $t0).TotalMinutes)
    $keep = Join-Path $WorkDir "csv_$job"
    New-Item -ItemType Directory -Force $keep | Out-Null
    Get-ChildItem -Filter '*.csv' | Copy-Item -Destination $keep -Force
    Remove-Item "$job.lock" -ErrorAction SilentlyContinue
    Pop-Location
    if (-not (Test-Path (Join-Path $keep 'fields_4.csv')) -or $err -ne '0') { L "  $job FAILED, not exported"; continue }
    $o = & $py ansys_usermat\apdl\fields_living.py $job $keep $deck 2>&1
    $o | ForEach-Object { L "  $_" }
    $done += "ansys_usermat/apdl/results/2026-10-living/fields_$job.json"
    if (Test-Path "ansys_usermat\apdl\results\2026-10-living\${job}_pm.json") {
        $done += "ansys_usermat/apdl/results/2026-10-living/${job}_pm.json"
    }
}
if ($Push -and $done.Count) {
    $c = & (Join-Path $repo 'commit.ps1') -Files $done -Message "Coupling with the living share in ANSYS (A-19): $Name (automatic)" 2>&1
    $c | ForEach-Object { L "  commit: $_" }
    for ($i = 1; $i -le 5; $i++) {
        & git fetch origin $Branch 2>$null
        & git merge --no-edit -m "Merge $Branch" FETCH_HEAD *> $null
        $q = & powershell.exe -NoProfile -ExecutionPolicy Bypass -File (Join-Path $repo 'push.ps1') -Branch $Branch 2>&1
        L "  push try ${i}: $(($q | Select-Object -Last 1))"
        if (($q | ForEach-Object { "$_" }) -match '^OK:') { break }
        Start-Sleep 120
    }
}
L "end"

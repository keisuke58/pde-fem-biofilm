<#
run_calibrated.ps1 -- research idea 7: the calibrated five-species point model
(TMCMC final MAP, gate off, c* = 25) in the partner element, 8^3, one run per
condition. Start it by hand when a condition's MAP file is on GitHub (no polling).

  .\ansys_usermat\apdl\run_calibrated.ps1 -Conditions CS [-Push]

Runs detached (WMI Win32_Process.Create, like run_chain.ps1) and returns at once.
Log: <WorkDir>\_cal5.log. For each condition given, once, it looks for
  keisuke58/Tmcmc202601, branch claude/gate-off-map-check,
  docs/revision/generated/final_theta_MAP/<C>.json
and, when the file is there (pinned to the last commit that touched it):
  1. write_eco_cfg.py --theta-json (kept as a record of the constants);
  2. make_wired_deck.py from base_w8_c6_g1_s015.dat (the two-species composition
     deck: Klempt 2024 stiffness, consumption 1, beta 0.02, s 0.15, phi_cap 0.9,
     T* = 1.0) with prop(8:27) = theta and prop(37) = 5 -> w8_<C>_g1_s015.dat;
  3. ANSYS with the native exe (surface fix) in WorkDir, -np 1, WITHOUT
     BIOFILM_ECO_CASE: the native build compares the deck's theta with the
     file's bit for bit, and APDL's reading of a 16-digit decimal differs from
     Fortran's in the last bit (10 Oct: every call failed, keycut). Without the
     file it takes theta from the deck and its defaults n = 5, c* = 25,
     alpha* = 0, eta = 1, which are the file's values (ecology_constants.py);
  4. export_runs_json.py (stress, nutrient) and export_cal5_json.py (point-model
     trace, all digits) into -Out; with -Push commit and push to the current
     branch.
A condition whose file is not there yet is skipped and logged.
The work dir must hold the native exe with its DLLs and base_w8_c6_g1_s015.dat
(F:\biofilm_upf_ch5, set up 10 Oct from F:\biofilm_upf_nativefix).
#>
param(
    [string[]]$Conditions = @('CH', 'CS', 'DS', 'DH'),
    [string]$WorkDir = 'F:\biofilm_upf_ch5',
    [string]$Out = 'ansys_usermat\apdl\results\2026-10-cal5',
    [switch]$Push,
    [switch]$Worker,
    [string]$Branch = ''
)
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$env:Path = "C:\Users\nishioka\git\cmd;C:\Users\nishioka\git\mingw64\bin;" + $env:Path
if (-not $Worker) {
    foreach ($f in 'ANSYS.exe', 'base_w8_c6_g1_s015.dat', 'post_all_stress.mac') {
        if (-not (Test-Path (Join-Path $WorkDir $f))) { throw "missing in ${WorkDir}: $f" }
    }
    $br = (& git -C $repo rev-parse --abbrev-ref HEAD).Trim()
    $a = @('-NoProfile', '-ExecutionPolicy', 'Bypass', '-WindowStyle', 'Hidden', '-File', "`"$PSCommandPath`"", '-Worker',
           '-Conditions', ($Conditions -join ','), '-WorkDir', "`"$WorkDir`"", '-Out', "`"$Out`"",
           '-Branch', $br)
    if ($Push) { $a += '-Push' }
    $r = Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{
        CommandLine = "powershell.exe $($a -join ' ')"; CurrentDirectory = $repo }
    if ($r.ReturnValue -ne 0) { throw "could not start (Win32_Process.Create returned $($r.ReturnValue))" }
    Write-Output "run_calibrated started, PID $($r.ProcessId); log $WorkDir\_cal5.log"
    return
}

$ErrorActionPreference = 'Continue'
Set-Location $repo
. (Join-Path $repo 'dev-env.ps1') | Out-Null
$py = (Get-Command python).Source
$log = Join-Path $WorkDir '_cal5.log'
function L($m) { [IO.File]::AppendAllText($log, "$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') $m`r`n") }
$api = 'https://api.github.com/repos/keisuke58/Tmcmc202601'
$mapBranch = 'claude/gate-off-map-check'
$mapDir = 'docs/revision/generated/final_theta_MAP'
$outs = 'all_stress.csv', 'elem_stress.csv', 'nut_field.csv', 'comp_trace.csv', 'phi_trace.csv', 'age_trace.csv', 'pm_trace.csv'
$Conditions = @($Conditions -split ',' | Where-Object { $_ })

function Find-Map($c) {
    # the last commit on the branch that touched <C>.json, or $null
    try {
        $cm = Invoke-RestMethod "$api/commits?sha=$mapBranch&path=$mapDir/$c.json&per_page=1" -TimeoutSec 60
        if (-not $cm -or @($cm).Count -eq 0) { return $null }
        return @($cm)[0].sha
    } catch { L "  $c lookup failed: $($_.Exception.Message)"; return $null }
}

function Push-Files($files, $msg) {
    $c = & (Join-Path $repo 'commit.ps1') -Files $files -Message $msg 2>&1
    $c | ForEach-Object { L "  commit: $_" }
    for ($i = 1; $i -le 5; $i++) {
        & git fetch origin $Branch 2>$null
        & git merge --no-edit -m "Merge $Branch" FETCH_HEAD *> $null
        $q = & powershell.exe -NoProfile -ExecutionPolicy Bypass -File (Join-Path $repo 'push.ps1') -Branch $Branch 2>&1
        L "  push try ${i}: $(($q | Select-Object -Last 1))"
        if (($q | ForEach-Object { "$_" }) -match '^OK:') { return }
        Start-Sleep 120
    }
}

function Run-Condition($c, $sha) {
    $job = "w8_${c}_g1_s015"
    $tj = Join-Path $WorkDir "${c}_theta_MAP_$($sha.Substring(0,7)).json"
    Invoke-WebRequest "https://raw.githubusercontent.com/keisuke58/Tmcmc202601/$sha/$mapDir/$c.json" -OutFile $tj -UseBasicParsing
    $o = & $py abaqus_composition\write_eco_cfg.py --theta-json $tj (Join-Path $WorkDir "eco_$c.txt") 2>&1
    L "  $o"
    $th = (Get-Content (Join-Path $WorkDir "eco_$c.txt"))[3].Trim().Split(' ')
    if ($th.Count -ne 20) { throw "eco_$c.txt: $($th.Count) theta values" }
    $props = ((0..19 | ForEach-Object { "$($_ + 8)=$($th[$_])" }) -join ',') + ',37=5'
    $o = & $py ansys_usermat\apdl\make_wired_deck.py (Join-Path $WorkDir 'base_w8_c6_g1_s015.dat') (Join-Path $WorkDir "$job.dat") `
        --props $props --post both --post-elem 220 2>&1
    if ($LASTEXITCODE -ne 0) { throw "make_wired_deck: $o" }
    foreach ($f in $outs) { $p = Join-Path $WorkDir $f; if (Test-Path $p) { Remove-Item $p -Force } }
    Remove-Item Env:\BIOFILM_ECO_CASE -ErrorAction SilentlyContinue
    L "START $job (theta $c.json @ $($sha.Substring(0,7)))"
    $t0 = Get-Date
    $p = Start-Process "$env:AWP_ROOT222\ANSYS\bin\winx64\ANSYS222.exe" -WorkingDirectory $WorkDir `
        -ArgumentList '-b', '-np', '1', '-j', $job, '-custom', '.\ANSYS.exe', '-i', "$job.dat", '-o', "out_$job.txt" -PassThru -WindowStyle Hidden
    $null = $p.Handle
    while (-not $p.HasExited) {
        Start-Sleep 20
        if (((Get-Date) - $t0).TotalMinutes -gt 60) {
            Get-CimInstance Win32_Process -Filter "ParentProcessId=$($p.Id)" | ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }
            Stop-Process -Id $p.Id -Force -ErrorAction SilentlyContinue; L '  killed after 60 min'; break
        }
    }
    $outTxt = Join-Path $WorkDir "out_$job.txt"
    $err = (Select-String $outTxt -Pattern 'ERROR\s+MESSAGES ENCOUNTERED=\s*(\d+)' | ForEach-Object { $_.Matches[0].Groups[1].Value }) -join ''
    $cut = (Select-String (Join-Path $WorkDir "$job.err") -Pattern 'bisection key' -ErrorAction SilentlyContinue | Measure-Object).Count
    foreach ($f in $outs) { $q = Join-Path $WorkDir $f; if (Test-Path $q) { Copy-Item $q (Join-Path $WorkDir ($f -replace '\.csv$', "_$job.csv")) -Force } }
    L "DONE $job exit $($p.ExitCode), errors $err, bisection warnings $cut, $([math]::Round(((Get-Date) - $t0).TotalMinutes, 1)) min"
    if ($err -ne '0' -or $cut -gt 0 -or -not (Test-Path (Join-Path $WorkDir "comp_trace_$job.csv"))) { L "  $job FAILED, not exported"; return $false }
    $env:BIOFILM_WORKDIR = $WorkDir
    $o = & $py ansys_usermat\apdl\export_runs_json.py $Out $job 2>&1; $o | ForEach-Object { L "  export: $_" }
    $o = & $py ansys_usermat\apdl\export_cal5_json.py $Out $job --theta-json $tj --source "Tmcmc202601@${sha}:$mapDir/$c.json" 2>&1
    $o | ForEach-Object { L "  export: $_" }
    if ($Push) {
        $files = @("$Out\$job.json", "$Out\${job}_pm.json") -replace '\\', '/'
        Push-Files $files "Calibrated five-species $c in the element, 8^3 (automatic)"
    }
    return $true
}

L "cal5 start: $($Conditions -join ','), branch $Branch"
try {
    foreach ($c in $Conditions) {
        $sha = Find-Map $c
        if (-not $sha) { L "  $c.json not on GitHub yet: skipped"; continue }
        try { $null = Run-Condition $c $sha } catch { L "  $c EXCEPTION: $_" }
    }
    L 'cal5 end'
} catch { L "EXCEPTION: $_" }

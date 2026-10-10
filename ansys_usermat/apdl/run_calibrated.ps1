<#
run_calibrated.ps1 -- research idea 7: the calibrated five-species point model
(TMCMC final MAP, gate off, c* = 25) in the partner element, 8^3, one run per
condition. Start it by hand when a condition's MAP file is on GitHub (no polling).

  .\ansys_usermat\apdl\run_calibrated.ps1 -Conditions CS [-Weights 1,1,1,1.5,2 -Tag fnpg] [-Push]
  .\ansys_usermat\apdl\run_calibrated.ps1 -Conditions CS -Growth [-Kappa 8.822731e-3] -Tag k0882 -Push

-Growth (10 Oct, P3 stage 2): mode 9 of the call-site fragment, alpha from the
point model's live biomass, the point model free at each Gauss point from the
condition's Day-1 composition (the calibration's config.json, taken from the
earlier _posterior run of the same condition on master), s = 0.2375 (T* = 1 is
Day 21), deltim 1/95, local nutrient on, BIOFILM_ECO_CASE set (newton 12 1e-20,
theta_tol 1e-12). Job w8_<C>_cal[_Tag]. Without -Growth: mode 7 as before
(composition only, Eq. 36 growth, optional -Weights), job w8_<C>_g1_s015[_Tag].
The work dir F:\biofilm_upf_cal5 holds the executable with mode 9 (built
10 Oct from F:\biofilm_upf_nativefix + the repo's fragments and ecology_native.f).

Runs detached (WMI Win32_Process.Create, like run_chain.ps1) and returns at once.
Log: <WorkDir>\_cal5.log. For each condition given, once, it looks for
  keisuke58/Tmcmc202601, branch claude/gate-off-map-check,
  docs/revision/generated/final_theta_MAP/<C>.json
and, when the file is there (pinned to the last commit that touched it):
  1. write_eco_cfg.py --theta-json (kept as a record of the constants);
  2. make_wired_deck.py from base_w8_c6_g1_s015.dat (the two-species composition
     deck: Klempt 2024 stiffness, consumption 1, beta 0.02, s 0.15, phi_cap 0.9,
     T* = 1.0) with prop(8:27) = theta and prop(37) = 5 -> w8_<C>_g1_s015.dat;
     -Weights f1..f5 (So, An, Vd, Fn, Pg; an assumption, not from a paper) adds the
     species-weighted growth law prop(36) = 1, prop(38:42) = f: alpha_dot =
     k_alpha phi sum_i f_i phi_i/sum(phi) -> w8_<C>_g1_s015_<Tag>.dat;
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
    [string]$WorkDir = 'F:\biofilm_upf_cal5',
    [string]$Out = 'ansys_usermat\apdl\results\2026-10-cal5',
    [string]$Weights = '',    # f1,f2,f3,f4,f5 for prop(38:42); empty = Eq. 36 unchanged
    [switch]$Growth,          # mode 9: alpha from the point model's live biomass (P3 stage 2)
    [double]$Kappa = 8.822731e-3,   # prop(7) of -Growth: CH at T* = 1 as Eq. 36 (10 Oct)
    [double]$YoungBio = 0,    # E_bio in MPa (deck YOUNG_BIO); 0 = the base deck's 1e-5 (10 Pa, Klempt 2024)
    [string]$Tag = '',
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
           '-Conditions', ($Conditions -join ','), '-WorkDir', "`"$WorkDir`"", '-Out', "`"$Out`"", '-Weights', "`"$Weights`"", '-Tag', "`"$Tag`"",
           '-Branch', $br)
    if ($Push) { $a += '-Push' }
    if ($Growth) { $a += @('-Growth', '-Kappa', $Kappa.ToString('R')) }
    if ($YoungBio -gt 0) { $a += @('-YoungBio', $YoungBio.ToString('R')) }
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

# the calibration's config.json (metadata.phi_init_exp, the Day-1 composition;
# the same loader as the ult runs, whose config.json is only on the GPU host)
$cfgPath = @{
    CS = 'data_5species/_runs/commensal_static_posterior/config.json'
    CH = 'data_5species/_runs/commensal_hobic_posterior/config.json'
    DS = 'data_5species/_runs/dysbiotic_static_posterior/config.json'
    DH = 'data_5species/_runs/Dysbiotic_HOBIC_20260226_041232/config.json'
}

function Run-Condition($c, $sha) {
    $job = if ($Growth) { "w8_${c}_cal" } else { "w8_${c}_g1_s015" }
    if ($Tag) { $job += "_$Tag" }
    $tj = Join-Path $WorkDir "${c}_theta_MAP_$($sha.Substring(0,7)).json"
    Invoke-WebRequest "https://raw.githubusercontent.com/keisuke58/Tmcmc202601/$sha/$mapDir/$c.json" -OutFile $tj -UseBasicParsing
    $eco = Join-Path $WorkDir "eco_$c.txt"
    $mk = @('--post', 'both', '--post-elem', '220')
    # stiffness sensitivity (P3_LITERATURE.ja.md: Pattem 2018, 0.5-46 kPa by AFM):
    # E(phi) = (phi^2 + f) E_bio, only E_bio changes
    if ($YoungBio -gt 0) { $mk += @('--set', "YOUNG_BIO=$($YoungBio.ToString('R'))") }
    if ($Growth) {
        # P3 stage 2 (10 Oct): the point model free at each point from the
        # Day-1 state, s = 0.2375 (T* = 1 is Day 21), deltim = 1/95 so one
        # coupling step is 25 x 1e-4 of point-model time, the paper's Newton
        # control, the deck-vs-file theta check with a tolerance
        $cfg = Join-Path $WorkDir "${c}_config.json"
        Invoke-WebRequest "https://raw.githubusercontent.com/keisuke58/Tmcmc202601/master/$($cfgPath[$c])" -OutFile $cfg -UseBasicParsing
        $eco = Join-Path $WorkDir "eco_${c}_cal.txt"
        $o = & $py abaqus_composition\write_eco_cfg.py --theta-json $tj --phi-init-config $cfg `
            --newton 12,1e-20 --theta-tol 1e-12 $eco 2>&1
        $mk += @('--deltim', (1.0 / 95.0).ToString('R'), '--time', '1.0')
    } else {
        $o = & $py abaqus_composition\write_eco_cfg.py --theta-json $tj $eco 2>&1
    }
    L "  $o"
    $th = (Get-Content $eco)[3].Trim().Split(' ')
    if ($th.Count -ne 20) { throw "$($eco): $($th.Count) theta values" }
    $props = ((0..19 | ForEach-Object { "$($_ + 8)=$($th[$_])" }) -join ',') + ',37=5'
    if ($Growth) { $props += ",7=$($Kappa.ToString('R')),28=9,31=0.2375,33=1.0" }
    if ($Weights) {
        if ($Growth) { throw "-Weights is the species-weighted Eq. 36, not for -Growth" }
        $w = @($Weights -split '[ ,]+' | Where-Object { $_ })
        if ($w.Count -ne 5) { throw "-Weights needs 5 values, got $($w.Count)" }
        $props += ',36=1,' + ((0..4 | ForEach-Object { "$($_ + 38)=$($w[$_])" }) -join ',')
    }
    $o = & $py ansys_usermat\apdl\make_wired_deck.py (Join-Path $WorkDir 'base_w8_c6_g1_s015.dat') (Join-Path $WorkDir "$job.dat") `
        --props $props @mk 2>&1
    if ($LASTEXITCODE -ne 0) { throw "make_wired_deck: $o" }
    foreach ($f in $outs) { $p = Join-Path $WorkDir $f; if (Test-Path $p) { Remove-Item $p -Force } }
    # ANSYS appends to <job>.err: an earlier run's warnings would count for this one
    Remove-Item (Join-Path $WorkDir "$job.err") -Force -ErrorAction SilentlyContinue
    if ($Growth) { $env:BIOFILM_ECO_CASE = $eco } else { Remove-Item Env:\BIOFILM_ECO_CASE -ErrorAction SilentlyContinue }
    L "START $job (theta $c.json @ $($sha.Substring(0,7)), weights '$Weights', growth $Growth kappa $Kappa)"
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
        $what = if ($Growth) { "growth from the calibrated point model, $c, 8^3, kappa $Kappa" } else { "Calibrated five-species $c in the element, 8^3" }
        Push-Files $files "$what (automatic)"
    }
    return $true
}

L "cal5 start: $($Conditions -join ','), weights '$Weights' growth $Growth kappa $Kappa tag '$Tag', branch $Branch"
try {
    foreach ($c in $Conditions) {
        $sha = Find-Map $c
        if (-not $sha) { L "  $c.json not on GitHub yet: skipped"; continue }
        try { $null = Run-Condition $c $sha } catch { L "  $c EXCEPTION: $_" }
    }
    L 'cal5 end'
} catch { L "EXCEPTION: $_" }

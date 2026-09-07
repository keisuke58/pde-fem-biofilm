<#
.SYNOPSIS
    Recompile Oliver's usercm.inc-dependent UPF pool and run a species-trial
    deck against it, always the safe way -- automates the two mechanical
    mistakes that caused both real bugs in the n=3/4/5 growth-term trials
    (N3_GROWTH_TRIAL.md): a partial usercm.inc recompile, and running a
    freshly-rebuilt Ussfin under the default DMP launch instead of
    -smp -np 1.

.DESCRIPTION
    Every edit to usercm.inc (adding a new species' scalars/offsets) is only
    safe if ALL FIVE files that #include it are recompiled together:
    USolBeg, Ussfin, Usermat_P21-V21_v222, NEM_UserData_P21_V05,
    userdata_P21-V21_Conection_Test. Forgetting one leaves it linked against
    a stale, layout-inconsistent copy of the common block -- this produced a
    SIG$SEGV once and a NEM D-Matrix NaN once, both silent about the real
    cause. This script always recompiles the fixed five, in the same
    link_v222.ps1 -Sources call, so that mistake is no longer possible.

    Separately: ANY freshly-recompiled Ussfin.obj hangs silently and
    indefinitely under ANSYS's default DMP (distributed-memory-parallel /
    MPI) launch mode, reproducible even on a completely unmodified deck.
    The fix is -smp -np 1 (serial, no MPI). This script always launches
    that way.

    This script does NOT generate the Fortran edits for a new species --
    only steps 5/6 (build, run) of N3_GROWTH_TRIAL.md's per-species
    checklist are automated here. Steps 1-4 (usercm.inc/USolBeg/Ussfin
    edits, deck parameters) still need to be written by hand for a new
    species, the same way n=3/4/5 were.

.PARAMETER WorkDir
    Directory holding Oliver's UPF source, .lrf/.def/.manifest scaffolding,
    and the deck to run. Not committed to this repo -- see RUNBOOK.md.

.PARAMETER DeckPath
    Path to the .dat deck to run (absolute, or relative to WorkDir).

.PARAMETER JobName
    ANSYS -j jobname. Also used to name the output log (<JobName>_out.txt)
    and to detect/clear a stale <JobName>.lock before running.

.PARAMETER SkipCompile
    Skip the recompile step and just run the deck against whatever
    ANSYS.exe is already in WorkDir. Use only when you are certain no
    usercm.inc-dependent file has changed since the last run through this
    script -- otherwise you are back to the exact hazard this script exists
    to prevent.

.EXAMPLE
    .\run_species_trial.ps1 -WorkDir F:\biofilm_upf_wired `
        -DeckPath ds_oliver_wired_n5trial.dat -JobName n5run3

.EXAMPLE
    # Re-run only (no recompile) -- e.g. re-checking determinism
    .\run_species_trial.ps1 -WorkDir F:\biofilm_upf_wired `
        -DeckPath ds_oliver_wired_n5trial.dat -JobName n5run4 -SkipCompile
#>
param(
    [Parameter(Mandatory=$true)][string]$WorkDir,
    [Parameter(Mandatory=$true)][string]$DeckPath,
    [Parameter(Mandatory=$true)][string]$JobName,
    [switch]$SkipCompile,
    [string]$AnsysRoot = $env:AWP_ROOT222
)

$ErrorActionPreference = "Stop"

# The fixed list of files that #include usercm.inc. This must stay in sync
# with N3_GROWTH_TRIAL.md's "recompile all 5 together" rule -- if Oliver's
# pool ever gains or loses a file that includes usercm.inc, update both
# this list and that doc in the same commit.
$UsercmDependentSources = @(
    "USolBeg_P21-V21_Conection_Test.F",
    "Ussfin_P21-V21_Conection_Test.F",
    "Usermat_P21-V21_v222.F",
    "NEM_UserData_P21_V05.F",
    "userdata_P21-V21_Conection_Test.f"
)

if (-not (Test-Path $WorkDir)) { throw "WorkDir does not exist: $WorkDir" }
if (-not $AnsysRoot) { throw "AWP_ROOT222 is not set and -AnsysRoot was not given." }

$deckFull = if ([System.IO.Path]::IsPathRooted($DeckPath)) { $DeckPath } else { Join-Path $WorkDir $DeckPath }
if (-not (Test-Path $deckFull)) { throw "Deck not found: $deckFull" }

if (-not $SkipCompile) {
    Write-Output "=== Recompiling all 5 usercm.inc-dependent files together ==="
    $linkScript = Join-Path $PSScriptRoot "link_v222.ps1"
    & $linkScript -WorkDir $WorkDir -Sources $UsercmDependentSources -AnsysRoot $AnsysRoot
    if ($LASTEXITCODE -ne 0) { throw "link_v222.ps1 failed -- see log above. Not attempting a run against a bad/missing ANSYS.exe." }
} else {
    Write-Output "=== -SkipCompile: reusing existing ANSYS.exe in $WorkDir (not recompiled) ==="
}

# A stale lock from a prior run of the same jobname blocks a fresh one.
$lockFile = Join-Path $WorkDir "$JobName.lock"
if (Test-Path $lockFile) {
    Write-Output "Removing stale lock: $lockFile"
    Remove-Item $lockFile -Force
}

$outFile = Join-Path $WorkDir "${JobName}_out.txt"
$ansysExe = Join-Path $AnsysRoot "ANSYS\bin\winx64\ANSYS222.exe"
$customExe = Join-Path $WorkDir "ANSYS.exe"
if (-not (Test-Path $customExe)) { throw "No ANSYS.exe in $WorkDir -- run without -SkipCompile first." }

Write-Output ""
Write-Output "=== Running $DeckPath as job '$JobName' (-smp -np 1) ==="
Push-Location $WorkDir
try {
    & $ansysExe -smp -np 1 -b -custom .\ANSYS.exe -j $JobName -i $deckFull -o $outFile
    $exitCode = $LASTEXITCODE
} finally {
    Pop-Location
}

if (-not (Test-Path $outFile)) {
    throw "ANSYS produced no output file at $outFile (exit code $exitCode) -- check for a launch-level failure."
}

$log = Get-Content $outFile -Raw
$errorLine = ($log -split "`r?`n" | Select-String "NUMBER OF ERROR" | Select-Object -Last 1)
$completedLine = ($log -split "`r?`n" | Select-String "RUN COMPLETED" | Select-Object -Last 1)

Write-Output ""
Write-Output "=== Result ==="
if ($errorLine) { Write-Output $errorLine.Line.Trim() } else { Write-Output "(no 'NUMBER OF ERROR' line found -- run may have aborted early)" }
if ($completedLine) { Write-Output $completedLine.Line.Trim() } else { Write-Output "(no 'RUN COMPLETED' banner found)" }

$errorCountMatch = [regex]::Match($log, "NUMBER OF ERROR\s+MESSAGES ENCOUNTERED=\s*(\d+)")
if ($errorCountMatch.Success -and $errorCountMatch.Groups[1].Value -eq "0" -and $completedLine) {
    Write-Output ""
    Write-Output "PASS: 0 errors, run completed. Log: $outFile"
    exit 0
} else {
    Write-Output ""
    Write-Output "FAIL (or inconclusive) -- inspect $outFile directly."
    exit 1
}

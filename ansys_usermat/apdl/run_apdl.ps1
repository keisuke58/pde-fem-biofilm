<#
.SYNOPSIS
    Run one APDL deck through the custom-built ANSYS.exe, with the disk-hygiene
    steps that were being done by hand all session: clean stale scratch files
    before AND after, check free space before running, and summarize
    errors/warnings from the log instead of requiring a manual grep.

.DESCRIPTION
    This machine (IKMHIWI03) hit 0 bytes free mid-session more than once from
    ANSYS scratch files (file*.esav/.full/.db/.rdb/.rst/.err) left behind by a
    crashed or completed run. This script is the checklist that grew out of
    that, made repeatable instead of re-typed by hand each time:
      1. clean scratch from any previous run
      2. report free space; ABORT if it is below -MinFreeGB (default 0.3 GB)
         rather than let a run crash mid-solve from disk-full
      3. copy the deck from the repo into the ANSYS working directory
      4. run it through the custom-built ANSYS.exe
      5. summarize NUMBER OF ERROR/WARNING MESSAGES from the log
      6. clean scratch again, report free space after

.PARAMETER Deck
    Deck filename, relative to this script's own directory
    (ansys_usermat/apdl/), e.g. t_growth_free.dat

.PARAMETER WorkDir
    The writable ANSYS working directory holding the custom-built ANSYS.exe
    and its runtime DLLs. Defaults to this machine's actual build location
    (see RUNBOOK.md) but can be overridden for a different machine/user.

.PARAMETER MinFreeGB
    Refuse to run if free space on the WorkDir's drive is below this many GB.
    Default 0.3 -- below that, a real solve has previously died mid-write.

.EXAMPLE
    .\run_apdl.ps1 -Deck t_growth_free.dat
    .\run_apdl.ps1 -Deck t_growth_cylinder_shell.dat -MinFreeGB 0.5
#>
param(
    [Parameter(Mandatory=$true)][string]$Deck,
    [string]$WorkDir = "C:\Users\nishioka\work\biofilm_upf",
    [double]$MinFreeGB = 0.3
)

$ErrorActionPreference = "Stop"
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$deckPath = Join-Path $scriptDir $Deck
if (-not (Test-Path $deckPath)) {
    throw "Deck not found: $deckPath"
}
$ansysExe = Join-Path $WorkDir "ANSYS.exe"
if (-not (Test-Path $ansysExe)) {
    throw "Custom ANSYS.exe not found in $WorkDir -- build it first (see RUNBOOK.md)."
}
if (-not $env:AWP_ROOT222) {
    throw "AWP_ROOT222 is not set -- ANSYS 2022 R2 environment not initialised."
}

function Get-FreeGB {
    [math]::Round((Get-PSDrive ($WorkDir.Substring(0,1))).Free / 1GB, 2)
}

function Clear-Scratch {
    Remove-Item (Join-Path $WorkDir "file*.esav"), (Join-Path $WorkDir "file*.full"),
                (Join-Path $WorkDir "file.db"), (Join-Path $WorkDir "file.rdb"),
                (Join-Path $WorkDir "file*.rst"), (Join-Path $WorkDir "file*.err") `
                -Force -ErrorAction SilentlyContinue
}

Write-Output "== run_apdl.ps1: $Deck =="

Clear-Scratch
$freeBefore = Get-FreeGB
Write-Output "Free space before run: $freeBefore GB"
if ($freeBefore -lt $MinFreeGB) {
    throw "Only $freeBefore GB free (< -MinFreeGB $MinFreeGB) -- refusing to run. " +
          "Clean up more (this machine has a known VSS shadow-copy issue where " +
          "deletions don't reliably free space -- see CLAUDE.md / draft email to Timo) " +
          "or lower -MinFreeGB explicitly if you are sure."
}

$outLog = [IO.Path]::ChangeExtension($Deck, $null).TrimEnd('.') + "_out.txt"
Copy-Item $deckPath (Join-Path $WorkDir $Deck) -Force

Push-Location $WorkDir
try {
    & "$env:AWP_ROOT222\ANSYS\bin\winx64\ANSYS222.exe" -b -custom .\ANSYS.exe -i $Deck -o $outLog
    $exitCode = $LASTEXITCODE
} finally {
    Pop-Location
}

$logPath = Join-Path $WorkDir $outLog
Write-Output "Exit code: $exitCode"
if (Test-Path $logPath) {
    $errLine = Select-String -Path $logPath -Pattern "NUMBER OF ERROR" | Select-Object -Last 1
    $warnLine = Select-String -Path $logPath -Pattern "NUMBER OF WARNING" | Select-Object -Last 1
    if ($errLine) { Write-Output $errLine.Line.Trim() }
    if ($warnLine) { Write-Output $warnLine.Line.Trim() }
    if (-not $errLine) {
        Write-Output "(no 'NUMBER OF ERROR' line found in the log -- check $outLog directly, e.g. a crash before that point)"
    }
} else {
    Write-Output "WARNING: expected log $logPath was not created."
}

Clear-Scratch
$freeAfter = Get-FreeGB
Write-Output "Free space after cleanup: $freeAfter GB"
Write-Output "Log: $logPath"

<#
.SYNOPSIS
    Run one APDL deck through the custom-built ANSYS.exe, with the disk-hygiene
    steps that were being done by hand all session: clean stale scratch files
    before AND after, check free space before running, and summarize
    errors/warnings from the log instead of requiring a manual grep.

.DESCRIPTION
    This machine (IKMHIWI03) hit 0 bytes free mid-session more than once from
    ANSYS scratch files (file*.esav/.full/.db/.rdb/.rst/.err) left behind by a
    crashed or completed run -- all on C:, which turned out to be chronically
    near-full and, per a same-day investigation, not reliably recoverable by
    deleting things (looks like VSS retaining the blocks; see the disk-space
    memory / draft email to Timo). Fixed properly on 2026-08-20 by moving the
    whole build (ANSYS.exe + runtime DLLs) to F:\biofilm_upf -- an otherwise
    untouched local data drive with ~3.7 TB free, confirmed working (same
    t_growth_free.dat result as the C: copy: 0 errors, 2 benign warnings).
    The disk-hygiene steps below are kept as defense in depth, not because
    F: is expected to fill up:
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
    and its runtime DLLs. Defaults to F:\biofilm_upf_kusepy, the 2026-09-29
    build (usermat_biofilm.f + ecology bridge), which is safe under -smp.
    F:\biofilm_upf holds the 8/19 build: give -Np 1 with it -- it solves
    non-deterministically under -smp -np > 1 (CLAUDE.md).

.PARAMETER Np
    Threads (shared-memory parallel, -smp -np N). Default 4. 1 = plain run.

.PARAMETER Ecology
    The deck calls the ecology bridge: start the material server
    (ansys_usermat/coupling/material_server.py, 127.0.0.1:8765) if nothing is
    listening there, and stop it again afterwards if this script started it.
    Without the server the bridge silently degrades to field mode.

.PARAMETER MinFreeGB
    Refuse to run if free space on the WorkDir's drive is below this many GB.
    Default 0.3 -- below that, a real solve has previously died mid-write.

.EXAMPLE
    .\run_apdl.ps1 -Deck t_growth_free.dat
    .\run_apdl.ps1 -Deck t_growth_cylinder_shell.dat -MinFreeGB 0.5
#>
param(
    [Parameter(Mandatory=$true)][string]$Deck,
    [string]$WorkDir = "F:\biofilm_upf_kusepy",
    [int]$Np = 4,
    [switch]$Ecology,
    [double]$MinFreeGB = 0.3
)
if ($Np -gt 1 -and $WorkDir.TrimEnd('\') -ieq "F:\biofilm_upf") {
    throw "F:\biofilm_upf's ANSYS.exe is not thread-safe -- use -Np 1 with it."
}

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
    # file.* and, under distributed runs, file0.*, file1.* ...
    Get-ChildItem $WorkDir -File -ErrorAction SilentlyContinue |
        Where-Object { $_.Name -match '^file\d*\.(esav|full|db|rdb|rst|err|r\d{3}|ldhi|stat|mntr|page|osav|emat|log|DSP|BCS|PCS)$' } |
        Remove-Item -Force -ErrorAction SilentlyContinue
}

function Test-Server {
    $c = New-Object Net.Sockets.TcpClient
    try { $c.Connect("127.0.0.1", 8765); $true } catch { $false } finally { $c.Close() }
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

$server = $null
if ($Ecology -and -not (Test-Server)) {
    $repoRoot = Split-Path -Parent (Split-Path -Parent $scriptDir)
    . (Join-Path $repoRoot "dev-env.ps1") | Out-Null
    $server = Start-Process python -ArgumentList "`"$(Join-Path $repoRoot 'ansys_usermat\coupling\material_server.py')`"" `
        -WorkingDirectory $repoRoot -WindowStyle Hidden -PassThru
    for ($i = 0; $i -lt 60 -and -not (Test-Server); $i++) { Start-Sleep -Milliseconds 500 }
    if (-not (Test-Server)) { $server | Stop-Process -Force; throw "material server did not come up on 8765" }
    Write-Output "Material server started (pid $($server.Id))"
} elseif ($Ecology) {
    Write-Output "Material server already listening on 8765"
}

$par = if ($Np -gt 1) { @("-smp", "-np", "$Np") } else { @() }
Push-Location $WorkDir
try {
    & "$env:AWP_ROOT222\ANSYS\bin\winx64\ANSYS222.exe" -b @par -custom .\ANSYS.exe -i $Deck -o $outLog
    $exitCode = $LASTEXITCODE
} finally {
    Pop-Location
    if ($server) { $server | Stop-Process -Force -ErrorAction SilentlyContinue; Write-Output "Material server stopped" }
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

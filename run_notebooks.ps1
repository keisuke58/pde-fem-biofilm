<#
.SYNOPSIS
    Re-executes every committed Jupyter notebook in the repo and reports
    which ones still run clean.

.DESCRIPTION
    This repo's verification notebooks (ansys_usermat/apdl,
    ansys_usermat/crosscheck, ansys_usermat/coupling,
    umat_tangent_test/abaqus_1elem) carry baked-in outputs -- their value is
    that the numbers shown are real, not stale. This script finds every
    *.ipynb under the repo (skipping .ipynb_checkpoints), re-executes each
    with nbconvert --execute --inplace, and summarizes pass/fail so drift
    (code or data changing under a notebook) gets caught instead of silently
    leaving old output in place.

    Some notebooks depend on machine-specific state that may not exist on
    every clone (an F:-based ANSYS/Abaqus work directory, gfortran on PATH,
    a licensed abaqus.bat) -- those are expected to fail cleanly elsewhere
    and are noted as such, not silently skipped.

.PARAMETER Path
    Root to search for notebooks. Defaults to the repo root (this script's
    own directory).

.EXAMPLE
    .\run_notebooks.ps1
#>
param(
    [string]$Path = $PSScriptRoot
)

. "$PSScriptRoot\dev-env.ps1"

$notebooks = Get-ChildItem -Path $Path -Recurse -Filter "*.ipynb" |
    Where-Object { $_.FullName -notmatch '\.ipynb_checkpoints' } |
    Sort-Object FullName

if (-not $notebooks) {
    Write-Output "No notebooks found under $Path"
    exit 0
}

Write-Output "== run_notebooks.ps1: $($notebooks.Count) notebook(s) found =="

$results = @()
$ErrorActionPreference = "Continue"
foreach ($nb in $notebooks) {
    $rel = Resolve-Path -Relative $nb.FullName
    Write-Output ""
    Write-Output "-- $rel --"
    python -m nbconvert --to notebook --execute --inplace $nb.FullName 2>&1 | Out-Null
    $ok = ($LASTEXITCODE -eq 0)
    $results += [PSCustomObject]@{ Notebook = $rel; OK = $ok }
    Write-Output $(if ($ok) { "  OK" } else { "  FAILED (exit $LASTEXITCODE)" })
}
$ErrorActionPreference = "Stop"

Write-Output ""
Write-Output "== summary =="
$results | ForEach-Object {
    $mark = if ($_.OK) { "PASS" } else { "FAIL" }
    Write-Output ("  {0,-4} {1}" -f $mark, $_.Notebook)
}

$failed = $results | Where-Object { -not $_.OK }
if ($failed) {
    Write-Output ""
    Write-Output "$($failed.Count) of $($results.Count) notebook(s) failed to re-execute."
    exit 1
}
Write-Output ""
Write-Output "All $($results.Count) notebooks re-executed clean."
exit 0

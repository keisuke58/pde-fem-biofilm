<#
verify_all.ps1 -- the verification chain of the living-share coupling in one command (A-66, 24 Oct 2026), everything
that runs without a licence: Python manufactured solutions (8/16/32^3, order of phi and c), the Python mesh series of
the separate example with antibiotic (r h/4 8..64, TVD 16..64, LW + sensor 16..64), the GCI collation of all series
(gci_a52.py, which also picks up the Abaqus and ANSYS results that exist) and the convergence figure.

  .\verify_all.ps1 [-Quick]      # -Quick: 8/16^3 MMS and meshes up to 32^3 (about 15 min); full about 2 h
Results: ansys_usermat/apdl/results/2026-10-living/gci_a52.{md,json}, papers/coupling_paper/fig_convergence.pdf,
log verify_all.log in the repository root.
#>
param([switch]$Quick)
$ErrorActionPreference = 'Continue'
. (Join-Path $PSScriptRoot 'dev-env.ps1') | Out-Null
Set-Location $PSScriptRoot
$log = Join-Path $PSScriptRoot 'verify_all.log'
function L($m) { $s = "$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') $m"; Write-Output $s; [IO.File]::AppendAllText($log, "$s`r`n") }
$cl = 'ansys_usermat\apdl\coupled_living_consumption.py'
$base = @('--front', '1', '--phi-min', '0.001', '--v2', '--ansys-bc', '--start-position', '--alpha0', '100', '--gab', '1',
          '--gc', '1', '--tab', '0.333333', '--fortran', '--fast', '--wall-grad', 'mirror', '--only', 'separate_ab')
L "start (quick = $Quick)"
# 1. manufactured solutions (Python): order of phi and c
$mms = if ($Quick) { 8, 16 } else { 8, 16, 32 }
foreach ($n in $mms) {
    $nf = if ($n -eq 32) { @() } else { @('--nf', "$n") }
    $o = & python $cl @nf --mms --artdiff --dtc 0.005 --tend 0.05 --only mixed --wall-grad mirror 2>&1
    L ("  MMS {0}^3: {1}" -f $n, (($o | Select-String 'L2 error') -join ' '))
}
# 2. mesh series
$meshes = if ($Quick) { @(@('artdiff', 8), @('artdiff', 16), @('artdiff', 32), @('tvd', 16), @('tvd', 32), @('lw', 16), @('lw', 32)) }
          else { @(@('artdiff', 8), @('artdiff', 16), @('artdiff', 24), @('artdiff', 32), @('artdiff', 48), @('artdiff', 64),
                   @('tvd', 16), @('tvd', 32), @('tvd', 64), @('lw', 16), @('lw', 24), @('lw', 32), @('lw', 48), @('lw', 64)) }
foreach ($m in $meshes) {
    $scheme, $n = $m
    $nf = if ($n -eq 32) { @() } else { @('--nf', "$n") }
    $opt = switch ($scheme) { 'artdiff' { @('--artdiff', '--dtc', '0.00125') } 'tvd' { @('--tvd') } 'lw' { @('--lw', '--sensor', '--sensor-pow', '2', '--dtc', '0.00125') } }
    $t0 = Get-Date
    $o = & python $cl @nf @base @opt 2>&1
    L ("  {0} {1}^3: {2} ({3:N1} min)" -f $scheme, $n, ($o | Select-Object -Last 1), ((Get-Date) - $t0).TotalMinutes)
}
# 3. GCI of all series and the figure
$o = & python ansys_usermat\apdl\gci_a52.py 2>&1
$o | Where-Object { "$_" -match '^\| (Python|Abaqus|ANSYS)' } | ForEach-Object { L "  $_" }
$o = & python papers\coupling_paper\fig_convergence.py 2>&1
L "  $($o | Select-Object -Last 1)"
L "end"

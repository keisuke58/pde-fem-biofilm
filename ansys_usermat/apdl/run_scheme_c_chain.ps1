<#
run_scheme_c_chain.ps1 -- a list of run_calibrated.ps1 runs (scheme (c) and its
reference (a), CLOUD_TO_IKMHIWI03_2026-10-14_C.ja.md) one after another, in a
process detached from the session (WMI Win32_Process.Create, as run_chain.ps1),
so a session restart does not stop it. Each run pushes its JSON when it ends.

  .\ansys_usermat\apdl\run_scheme_c_chain.ps1 -Set B [-DryRun]
  .\ansys_usermat\apdl\run_scheme_c_chain.ps1 -Set A6

Sets (the order is the instruction's):
  A6  antibiotic from T* = 0.5 (prop(54)), b = 1 with alpha* 0.32 and 0.65,
      So x10, anaerobes x10, (c) and (a); the pressure feedback test
      (prop(52:53)); all CH 8^3. Done 14 Oct; kept to rerun.
  B   the long list: CS 16^3 (a), (c); DS, DH 8^3 and 16^3 (skipped and logged
      while their MAP is missing); CH 24^3 (a), (c); CH antibiotic 16^3 (a), (c);
      kappa 2.8 on 16^3 (stopped at T* 0.35 on 8^3: checks whether the mesh
      changes that); DH anaerobes x10 8^3 when the DH MAP is there.
Minutes per run (IKMHIWI03, -np 1): 8^3 1-2, 16^3 (a) 10 / (c) 23, 24^3 about
40 / 90 (estimates, written to the log before each run).
Exe: F:\biofilm_upf_cal5sc (the cal5 build plus alpha* per point, prop(52:54),
off unless set; it reproduces F:\biofilm_upf_cal5 to the last digit).
Log: <WorkDir>\_cal5.log (lines "chain <set> ...").
#>
param(
    [Parameter(Mandatory = $true)][ValidateSet('A6', 'B')][string]$Set,
    [string]$WorkDir = 'F:\biofilm_upf_cal5sc',
    [switch]$DryRun,
    [switch]$Worker
)
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$rc = Join-Path $PSScriptRoot 'run_calibrated.ps1'
$log = Join-Path $WorkDir '_cal5.log'
function L($m) { [IO.File]::AppendAllText($log, "$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') $m`r`n") }

$b1 = @('-B', '1,1,1,1,1')
$sets = @{
    A6 = @(
        @{ n = 'pressure feedback (c)'; m = 2; a = @('-Conditions', 'CH', '-Growth') + $b1 + @('-ExtraProps', '52=0.32,53=2.4e-10,54=0', '-Tag', 'k0882_sc') },
        @{ n = 'A-6 (c) b=1 0.32'; m = 2; a = @('-Conditions', 'CH', '-Growth', '-AlphaStar', '0.32') + $b1 + @('-ExtraProps', '54=0.5', '-Tag', 'k0882_ab05') },
        @{ n = 'A-6 (a) b=1 0.32'; m = 2; a = @('-Conditions', 'CH', '-AlphaStar', '0.32') + $b1 + @('-ExtraProps', '54=0.5', '-Tag', 'ab05') },
        @{ n = 'A-6 (c) b=1 0.65'; m = 2; a = @('-Conditions', 'CH', '-Growth', '-AlphaStar', '0.65') + $b1 + @('-ExtraProps', '54=0.5', '-Tag', 'k0882_ab05_a065') },
        @{ n = 'A-6 (a) b=1 0.65'; m = 2; a = @('-Conditions', 'CH', '-AlphaStar', '0.65') + $b1 + @('-ExtraProps', '54=0.5', '-Tag', 'ab05_a065') },
        @{ n = 'A-6 (c) So x10'; m = 2; a = @('-Conditions', 'CH', '-Growth', '-AlphaStar', '0.32', '-B', '10,1,1,1,1', '-ExtraProps', '54=0.5', '-Tag', 'k0882_ab05_So') },
        @{ n = 'A-6 (a) So x10'; m = 2; a = @('-Conditions', 'CH', '-AlphaStar', '0.32', '-B', '10,1,1,1,1', '-ExtraProps', '54=0.5', '-Tag', 'ab05_So') },
        @{ n = 'A-6 (c) anaerobes x10'; m = 2; a = @('-Conditions', 'CH', '-Growth', '-AlphaStar', '0.32', '-B', '1,1,10,10,10', '-ExtraProps', '54=0.5', '-Tag', 'k0882_ab05_ana') },
        @{ n = 'A-6 (a) anaerobes x10'; m = 2; a = @('-Conditions', 'CH', '-AlphaStar', '0.32', '-B', '1,1,10,10,10', '-ExtraProps', '54=0.5', '-Tag', 'ab05_ana') }
    )
    B = @(
        @{ n = 'CS 16^3 (a)'; m = 10; a = @('-Conditions', 'CS', '-Mesh', '16') },
        @{ n = 'CS 16^3 (c)'; m = 25; a = @('-Conditions', 'CS', '-Mesh', '16', '-Growth', '-Tag', 'k0882') },
        @{ n = 'DS, DH 8^3 (a) (skipped without MAP)'; m = 4; a = @('-Conditions', 'DS,DH') },
        @{ n = 'DS, DH 8^3 (c) (skipped without MAP)'; m = 4; a = @('-Conditions', 'DS,DH', '-Growth', '-Tag', 'k0882') },
        @{ n = 'DS, DH 16^3 (a) (skipped without MAP)'; m = 20; a = @('-Conditions', 'DS,DH', '-Mesh', '16') },
        @{ n = 'DS, DH 16^3 (c) (skipped without MAP)'; m = 50; a = @('-Conditions', 'DS,DH', '-Mesh', '16', '-Growth', '-Tag', 'k0882') },
        @{ n = 'DH 8^3 (c) anaerobes x10 from T* 0.5 (skipped without MAP)'; m = 2; a = @('-Conditions', 'DH', '-Growth', '-AlphaStar', '0.32', '-B', '1,1,10,10,10', '-ExtraProps', '54=0.5', '-Tag', 'k0882_ab05_ana') },
        @{ n = 'CH 24^3 (a)'; m = 40; a = @('-Conditions', 'CH', '-Mesh', '24') },
        @{ n = 'CH 24^3 (c)'; m = 90; a = @('-Conditions', 'CH', '-Mesh', '24', '-Growth', '-Tag', 'k0882') },
        @{ n = 'CH 16^3 (c) antibiotic from T* 0'; m = 25; a = @('-Conditions', 'CH', '-Mesh', '16', '-Growth', '-AlphaStar', '0.32') + $b1 + @('-Tag', 'k0882_ab') },
        @{ n = 'CH 16^3 (a) antibiotic from T* 0'; m = 10; a = @('-Conditions', 'CH', '-Mesh', '16', '-AlphaStar', '0.32') + $b1 + @('-Tag', 'ab') },
        @{ n = 'CH 16^3 (c) kappa 2.8'; m = 25; a = @('-Conditions', 'CH', '-Mesh', '16', '-Growth', '-Kappa', '2.8', '-Tag', 'k2800') }
    )
}
$runs = $sets[$Set]

if (-not $Worker) {
    $tot = ($runs | ForEach-Object { $_.m } | Measure-Object -Sum).Sum
    Write-Output "set $Set`: $($runs.Count) runs, about $tot min if every MAP is there"
    $runs | ForEach-Object { Write-Output ("  {0,-55} {1,4} min  {2}" -f $_.n, $_.m, ($_.a -join ' ')) }
    if ($DryRun) { return }
    $cmd = "powershell.exe -NoProfile -ExecutionPolicy Bypass -WindowStyle Hidden -File `"$PSCommandPath`" -Worker -Set $Set -WorkDir `"$WorkDir`""
    $r = Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd; CurrentDirectory = $repo }
    if ($r.ReturnValue -ne 0) { throw "could not start (Win32_Process.Create returned $($r.ReturnValue))" }
    Write-Output "started, PID $($r.ProcessId); log $log"
    return
}

$env:Path = "C:\Users\nishioka\git\cmd;C:\Users\nishioka\git\mingw64\bin;" + $env:Path
$br = (& git -C $repo rev-parse --abbrev-ref HEAD).Trim()
L "chain $Set start: $($runs.Count) runs, branch $br, exe $WorkDir"
foreach ($r in $runs) {
    L "chain $Set next: $($r.n), estimate $($r.m) min"
    & powershell.exe -NoProfile -ExecutionPolicy Bypass -File $rc -Worker -Branch $br -WorkDir $WorkDir -Push @($r.a)
}
L "chain $Set end"

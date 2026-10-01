<#
run_wired.ps1 -- run one deck through the partner-element ANSYS.exe in
F:\biofilm_upf_wired, with the material server around it.

  .\ansys_usermat\apdl\run_wired.ps1 -Deck ds_mode7_case3_e220.dat [-Case 2sp_case3]
      [-ActiveSpecies N] [-Job name] [-Judge]

- always -np 1 (the default 4-rank DMP hangs with this element);
- its own job name (default: the deck name), so a stale file.lock from an
  aborted run under another name does not block it;
- moves comp_trace.csv / pm_trace.csv / phi_trace.csv aside first
  (<name>_prev_<time>.csv), so the trace holds this run only;
- starts material_server.py (--case or --active-species), stops it afterwards
  even on failure;
- prints the error/warning counts; -Judge runs judge_comp_trace.py on the
  composition trace and copies it to comp_trace_<job>.csv.
#>
param(
    [Parameter(Mandatory)] [string]$Deck,
    [string]$Case = '',
    [int]$ActiveSpecies = 0,
    [string]$Job = '',
    [string]$WorkDir = 'F:\biofilm_upf_wired',
    [switch]$Judge
)
$ErrorActionPreference = 'Stop'
$repo = Resolve-Path (Join-Path $PSScriptRoot '..\..')
. (Join-Path $repo 'dev-env.ps1') | Out-Null
if (-not $Job) { $Job = [IO.Path]::GetFileNameWithoutExtension($Deck) }
if (-not (Test-Path (Join-Path $WorkDir $Deck))) { throw "deck not found: $Deck" }

$ts = Get-Date -Format 'HHmmssfff'
foreach ($f in 'comp_trace.csv', 'pm_trace.csv', 'phi_trace.csv') {
    $p = Join-Path $WorkDir $f
    if (Test-Path $p) { Move-Item $p ($p -replace '\.csv$', "_prev_$ts.csv") }
}

$srvArgs = @((Join-Path $repo 'ansys_usermat\coupling\material_server.py'))
if ($Case) { $srvArgs += @('--case', $Case) }
elseif ($ActiveSpecies) { $srvArgs += @('--active-species', "$ActiveSpecies") }
$py = (Get-Command python).Source
$srv = Start-Process $py -ArgumentList $srvArgs -PassThru -WindowStyle Hidden `
    -RedirectStandardOutput "$env:TEMP\ms_out.txt" -RedirectStandardError "$env:TEMP\ms_err.txt"
try {
    $up = $false
    for ($i = 0; $i -lt 30 -and -not $up; $i++) {
        Start-Sleep 1
        if ($srv.HasExited) { throw "material server exited: $(Get-Content "$env:TEMP\ms_err.txt" -Raw)" }
        $c = New-Object Net.Sockets.TcpClient
        try { $c.Connect('127.0.0.1', 8765); $up = $true } catch {} finally { $c.Close() }
    }
    if (-not $up) { throw 'material server did not open port 8765' }

    Push-Location $WorkDir
    try {
        $out = "out_$Job.txt"
        $t0 = Get-Date
        & "$env:AWP_ROOT222\ANSYS\bin\winx64\ANSYS222.exe" -b -np 1 -j $Job -custom .\ANSYS.exe -i $Deck -o $out
        $rc = $LASTEXITCODE
        "ANSYS exit $rc, $([int]((Get-Date) - $t0).TotalSeconds) s, output $out"
        Select-String $out -Pattern 'NUMBER OF (ERROR|WARNING)' | ForEach-Object { $_.Line.Trim() }
    } finally { Pop-Location }
} finally {
    if (-not $srv.HasExited) { Stop-Process -Id $srv.Id -Confirm:$false }
    'material server stopped'
}

if ($Judge) {
    $trace = Join-Path $WorkDir 'comp_trace.csv'
    if (-not (Test-Path $trace)) { throw 'no comp_trace.csv written' }
    Copy-Item $trace (Join-Path $WorkDir "comp_trace_$Job.csv")
    $jargs = @((Join-Path $PSScriptRoot 'judge_comp_trace.py'), $trace)
    if ($Case) { $jargs += @('--case', $Case) }
    & $py @jargs
    exit $LASTEXITCODE
}
exit $rc

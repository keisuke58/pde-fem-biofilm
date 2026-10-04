<#
run_wired.ps1 -- run one deck through the partner-element ANSYS.exe in
F:\biofilm_upf_wired, with the material server around it.

  .\ansys_usermat\apdl\run_wired.ps1 -Deck ds_mode7_case3_e220.dat [-Case 2sp_case3]
      [-ActiveSpecies N] [-Job name] [-TimeoutMin 60] [-Judge]

- always -np 1 (the default 4-rank DMP hangs with this element);
- its own job name (default: the deck name), so a stale file.lock from an
  aborted run under another name does not block it; a stale <job>.lock stops
  the script (it is never deleted here: pick another -Job);
- refuses to start when port 8765 is taken (a leftover server would answer
  with its own case / species count);
- copies the post-processing macros (callsite/post_*.mac) into the work dir;
- moves comp/pm/phi traces and all_stress.csv / elem_stress.csv aside first
  (<name>_prev_<time>.csv), so they hold this run only, and afterwards copies
  each one written to <name>_<job>.csv;
- starts material_server.py (--case or --active-species), stops it afterwards
  even on failure;
- kills ANSYS (and its child processes) after -TimeoutMin minutes;
- exit code non-zero when ANSYS fails, reports errors, or -Judge fails.
#>
param(
    [Parameter(Mandatory)] [string]$Deck,
    [string]$Case = '',
    [int]$ActiveSpecies = 0,
    [string]$Job = '',
    [string]$WorkDir = 'F:\biofilm_upf_wired',
    [double]$TimeoutMin = 60,
    [switch]$Judge
)
$ErrorActionPreference = 'Stop'
$repo = Resolve-Path (Join-Path $PSScriptRoot '..\..')
. (Join-Path $repo 'dev-env.ps1') | Out-Null
if (-not $Job) { $Job = [IO.Path]::GetFileNameWithoutExtension($Deck) }
if (-not (Test-Path (Join-Path $WorkDir $Deck))) { throw "deck not found: $Deck" }
if (Test-Path (Join-Path $WorkDir "$Job.lock")) {
    throw "$Job.lock exists in $WorkDir (an aborted run?): use another -Job"
}

function Test-Port8765 {
    $c = New-Object Net.Sockets.TcpClient
    try { $c.Connect('127.0.0.1', 8765); return $true } catch { return $false } finally { $c.Close() }
}
if (Test-Port8765) { throw 'port 8765 already in use: another material server is running; stop it first' }

function Stop-Tree([int]$id) {
    Get-CimInstance Win32_Process -Filter "ParentProcessId=$id" | ForEach-Object { Stop-Tree $_.ProcessId }
    Stop-Process -Id $id -Force -Confirm:$false -ErrorAction SilentlyContinue
}

Copy-Item (Join-Path $PSScriptRoot 'callsite\post_*.mac') $WorkDir -Force
$outputs = 'comp_trace.csv', 'age_trace.csv', 'pm_trace.csv', 'phi_trace.csv', 'all_stress.csv', 'elem_stress.csv'
$ts = Get-Date -Format 'HHmmssfff'
foreach ($f in $outputs) {
    $p = Join-Path $WorkDir $f
    if (Test-Path $p) { Move-Item $p ($p -replace '\.csv$', "_prev_$ts.csv") }
}

$srvArgs = @((Join-Path $repo 'ansys_usermat\coupling\material_server.py'))
if ($Case) { $srvArgs += @('--case', $Case) }
elseif ($ActiveSpecies) { $srvArgs += @('--active-species', "$ActiveSpecies") }
$py = (Get-Command python).Source
$srv = Start-Process $py -ArgumentList $srvArgs -PassThru -WindowStyle Hidden `
    -RedirectStandardOutput "$env:TEMP\ms_out.txt" -RedirectStandardError "$env:TEMP\ms_err.txt"
$rc = 1
try {
    $up = $false
    for ($i = 0; $i -lt 30 -and -not $up; $i++) {
        Start-Sleep 1
        if ($srv.HasExited) { throw "material server exited: $(Get-Content "$env:TEMP\ms_err.txt" -Raw)" }
        $up = Test-Port8765
    }
    if (-not $up) { throw 'material server did not open port 8765' }

    $out = "out_$Job.txt"
    $t0 = Get-Date
    $ans = Start-Process "$env:AWP_ROOT222\ANSYS\bin\winx64\ANSYS222.exe" -WorkingDirectory $WorkDir `
        -ArgumentList '-b', '-np', '1', '-j', $Job, '-custom', '.\ANSYS.exe', '-i', $Deck, '-o', $out `
        -PassThru -NoNewWindow
    $null = $ans.Handle                                  # keeps ExitCode readable after exit (PS 5.1)
    if (-not $ans.WaitForExit([int]($TimeoutMin * 60000))) {
        Stop-Tree $ans.Id
        throw "ANSYS did not finish within $TimeoutMin min: killed (see $out)"
    }
    $rc = $ans.ExitCode
    "ANSYS exit $rc, $([int]((Get-Date) - $t0).TotalSeconds) s, output $out"
    $counts = Select-String (Join-Path $WorkDir $out) -Pattern 'NUMBER OF (ERROR|WARNING)'
    $counts | ForEach-Object { $_.Line.Trim() }
    $err = $counts | Where-Object { $_.Line -match 'ERROR\s+MESSAGES ENCOUNTERED=\s*(\d+)' } |
        ForEach-Object { [int]$Matches[1] }
    if ($rc -eq 0 -and $err -gt 0) { $rc = 2 }
} finally {
    if (-not $srv.HasExited) { Stop-Process -Id $srv.Id -Confirm:$false }
    'material server stopped'
}

foreach ($f in $outputs) {
    $p = Join-Path $WorkDir $f
    if (Test-Path $p) {
        $dst = $p -replace '\.csv$', "_$Job.csv"
        Copy-Item $p $dst
        "wrote $(Split-Path $dst -Leaf)"
    }
}

if ($Judge) {
    $trace = Join-Path $WorkDir 'comp_trace.csv'
    if (-not (Test-Path $trace)) { throw 'no comp_trace.csv written' }
    $jargs = @((Join-Path $PSScriptRoot 'judge_comp_trace.py'), $trace)
    if ($Case) { $jargs += @('--case', $Case) }
    & $py @jargs
    if ($LASTEXITCODE -ne 0 -and $rc -eq 0) { $rc = $LASTEXITCODE }
}
exit $rc

<#
run_comp.ps1 -- build the composition UMAT (make_umat.py) and run one Abaqus job
with the material server, on IKMHIWI03 (Abaqus 2024, ifort 2025.3 + VS 18).

  .\abaqus_composition\run_comp.ps1 -Inp abaqus_composition\one_elem.inp -Case 2sp_case6 [-Port 8766]

- work dir F:\abaqus_work\comp_<job> (never C:);
- compiles the C shim (ansys_usermat\coupling\biofilm_py_eval.c) with cl and
  links it through a job-local abaqus_v6.env (link_sl += the object);
- starts material_server.py on -Port (default 8766: the ANSYS chains use 8765)
  with --case, sets BIOFILM_PY_PORT for Abaqus, stops the server afterwards;
- the fragments write comp_trace.csv into the work dir, so
  judge_comp_trace.py can check the run exactly as it checks an ANSYS run.
#>
param(
    [Parameter(Mandatory)] [string]$Inp,
    [string]$Case = '2sp_case6',
    [int]$Port = 8766,
    [string]$WorkRoot = 'F:\abaqus_work'
)
$ErrorActionPreference = 'Stop'
$repo = Resolve-Path (Join-Path $PSScriptRoot '..')
# no dev-env.ps1 here: its C:\msys64\usr\bin (GNU find, link) breaks vcvars64 and the MSVC link
$py = 'C:\Users\nishioka\AppData\Local\Programs\Python\Python312\python.exe'
$env:Path = ($env:Path -split ';' | Where-Object { $_ -notmatch 'msys64' }) -join ';'   # also when the caller sourced it
$job = [IO.Path]::GetFileNameWithoutExtension($Inp)
$wd = Join-Path $WorkRoot "comp_$job"
if (Test-Path $wd) { Remove-Item $wd -Recurse -Force }
New-Item -ItemType Directory -Force $wd | Out-Null
Copy-Item (Join-Path $repo $Inp) $wd

if (-not (Get-Command ifort -ErrorAction SilentlyContinue)) {
    # as link_v222.ps1: vcvars64.bat needs vswhere.exe on PATH, or it leaves LIB unset
    $env:Path = 'C:\Program Files (x86)\Microsoft Visual Studio\Installer;' + $env:Path
    $vcvars = 'C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat'
    $ifortEnv = 'C:\Program Files (x86)\Intel\oneAPI\compiler\2025.3\env\vars.bat'
    $dump = cmd /c "call `"$vcvars`" & call `"$ifortEnv`" & set" 2>&1
    foreach ($l in $dump) { if ("$l" -match '^([^=]+)=(.*)$') { [Environment]::SetEnvironmentVariable($Matches[1], $Matches[2], 'Process') } }
    if (-not (Get-Command ifort -ErrorAction SilentlyContinue)) {
        throw "ifort still not on PATH after vcvars64 + vars.bat ($($dump.Count) lines from cmd; first: $($dump | Select-Object -First 3))"
    }
}
& $py (Join-Path $repo 'abaqus_composition\make_umat.py') (Join-Path $wd 'umat_comp.for')
Push-Location $wd
try {
    $o = & cl /nologo /c /O2 /MD /GS- /Zl (Join-Path $repo 'ansys_usermat\coupling\biofilm_py_eval.c') /Fo:biofilm_py_eval.obj 2>&1
    if ($LASTEXITCODE -ne 0) { $o; throw 'cl failed on the C shim' }
    $obj = (Join-Path $wd 'biofilm_py_eval.obj') -replace '\\', '/'
    # the fragments' trace OPEN/WRITE pull in the static Intel runtime (libifcoremt)
    # unless the object asks for the DLL one, which Abaqus links (LIBIFCOREMD)
    @("compile_fortran += ['/libs:dll', '/threads']", "link_sl += ['$obj']") | Set-Content abaqus_v6.env -Encoding ascii

    $srv = Start-Process $py -ArgumentList @((Join-Path $repo 'ansys_usermat\coupling\material_server.py'),
        '--port', "$Port", '--case', $Case) -PassThru -WindowStyle Hidden `
        -RedirectStandardOutput (Join-Path $wd 'ms_out.txt') -RedirectStandardError (Join-Path $wd 'ms_err.txt')
    try {
        $up = $false
        for ($i = 0; $i -lt 60 -and -not $up; $i++) {
            Start-Sleep 1
            if ($srv.HasExited) { throw "material server exited: $(Get-Content (Join-Path $wd 'ms_err.txt') -Raw)" }
            $c = New-Object Net.Sockets.TcpClient
            try { $c.Connect('127.0.0.1', $Port); $up = $true } catch {} finally { $c.Close() }
        }
        if (-not $up) { throw "material server did not open port $Port" }
        $env:BIOFILM_PY_PORT = "$Port"
        $ErrorActionPreference = 'Continue'
        & abaqus job=$job input="$job.inp" user=umat_comp.for cpus=1 interactive ask_delete=OFF 2>&1 | Tee-Object -FilePath (Join-Path $wd 'abaqus_out.txt') | Select-Object -Last 15
        $ErrorActionPreference = 'Stop'
    } finally {
        if (-not $srv.HasExited) { Stop-Process -Id $srv.Id -Force }
    }
    $sta = Join-Path $wd "$job.sta"
    if ((Test-Path $sta) -and (Select-String $sta -Pattern 'COMPLETED SUCCESSFULLY' -Quiet)) { "PASS: $job completed" }
    else { "NOT COMPLETE: see $wd ($job.msg / $job.log)" }
    "work dir: $wd"
} finally { Pop-Location }

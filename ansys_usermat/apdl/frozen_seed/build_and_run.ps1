<#
build_and_run.ps1 -- request A-16 on IKMHIWI03: build the test exe
(usermat_frozen_wrap.f + the unchanged usermat_biofilm.f, its subroutine
usermat renamed usermat_core in the build directory only), run
frozen_seed.dat and write results/2026-10-uel/seed_frozen_ansys_<N>.json.

    .\ansys_usermat\apdl\frozen_seed\build_and_run.ps1 [-N 8] [-NoBuild]

Build directory F:\biofilm_upf_frozen (scaffolding copied from
F:\biofilm_upf_kusepy: DLLs, ansys.lrf, the C shim object of the py hook).
#>
param([int]$N = 8, [switch]$NoBuild, [string]$WorkDir = 'F:\biofilm_upf_frozen')
$ErrorActionPreference = 'Stop'
$repo = Resolve-Path (Join-Path $PSScriptRoot '..\..\..')
$here = $PSScriptRoot
$py = 'C:\Users\nishioka\AppData\Local\Programs\Python\Python312\python.exe'
if (-not $NoBuild) {
    if (-not (Test-Path $WorkDir)) {
        New-Item -ItemType Directory $WorkDir | Out-Null
        Get-ChildItem F:\biofilm_upf_kusepy -File | Where-Object {
            $_.Extension -in '.dll', '.lrf', '.def', '.manifest' } |
            Copy-Item -Destination $WorkDir
    }
    Get-ChildItem $WorkDir -Filter *.obj | Remove-Item
    Copy-Item (Join-Path $here 'py_bridge_stub.f') $WorkDir
    $core = (Get-Content (Join-Path $repo 'ansys_usermat\usermat_biofilm.f') -Raw) -replace '(?m)^(\s+subroutine\s+)usermat\(', '${1}usermat_core('
    if ($core -notmatch 'subroutine usermat_core\(') { throw 'rename of usermat failed' }
    Set-Content (Join-Path $WorkDir 'usermat_biofilm_core.f') $core -NoNewline -Encoding ascii
    Copy-Item (Join-Path $here 'usermat_frozen_wrap.f') $WorkDir
    # the stub module first: usermat_biofilm_core.f uses it
    & (Join-Path $repo 'ansys_usermat\apdl\link_v222.ps1') -WorkDir $WorkDir `
        -Sources py_bridge_stub.f, usermat_biofilm_core.f, usermat_frozen_wrap.f | Select-Object -Last 3
    if (-not (Test-Path (Join-Path $WorkDir 'ANSYS.exe'))) { throw "no ANSYS.exe in $WorkDir (see _link_v222_tmp.log)" }
}
$job = "frozen_seed_$N"
Get-ChildItem $WorkDir -Filter "$job.*" | Remove-Item -Force
(Get-Content (Join-Path $here 'frozen_seed.dat')) -replace '^N = 8$', "N = $N" | Set-Content (Join-Path $WorkDir "$job.inp") -Encoding ascii
Push-Location $WorkDir
try {
    Remove-Item frozen_elems.csv, frozen_gp.csv -ErrorAction SilentlyContinue
    $t0 = Get-Date
    & "$env:AWP_ROOT222\ANSYS\bin\winx64\ANSYS222.exe" -b -smp -np 1 -custom .\ANSYS.exe -j $job -i "$job.inp" -o "$job.out" | Out-Null
    "{0}: exit {1}, {2:N1} min, errors: {3}" -f $job, $LASTEXITCODE, ((Get-Date) - $t0).TotalMinutes,
        ((Select-String "$job.out" -Pattern '\*\*\* ERROR' | Measure-Object).Count)
} finally { Pop-Location }
& $py (Join-Path $here 'frozen_seed_json.py') $WorkDir $N

<#
post_living_dense.ps1 -- the 40 fields of a finished run of the coupling with the living share (addendum 5, the
animation): post_fields_living_dense.mac on <WorkDir>\<Job>.rst in <WorkDir>\post_dense_<Job> (so a chain running in
WorkDir is not disturbed), then fields_living.py <Job> <dir> <deck> 40 -> results/2026-10-living/fields_<Job>_dense.json.

  .\ansys_usermat\apdl\post_living_dense.ps1 -Job v2s_16_sep_ab [-WorkDir F:\biofilm_upf_living3]
#>
param(
    [Parameter(Mandatory = $true)][string]$Job,
    [string]$WorkDir = 'F:\biofilm_upf_living3'
)
$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$rst = Join-Path $WorkDir "$Job.rst"
if (-not (Test-Path $rst)) { throw "no $rst" }
$d = Join-Path $WorkDir "post_dense_$Job"
New-Item -ItemType Directory -Force $d | Out-Null
Copy-Item (Join-Path $PSScriptRoot 'post_fields_living_dense.mac') $d -Force
Set-Content -Encoding ascii (Join-Path $d 'post_dense.dat') "/POST1`r`nFILE,$Job,rst,$WorkDir\`r`n/INPUT,post_fields_living_dense,mac`r`nFINISH`r`n"
$exe = "$env:AWP_ROOT222\ANSYS\bin\winx64\ANSYS222.exe"
Push-Location $d
& $exe -b -np 1 -j postd -i post_dense.dat -o out_post.txt
$code = $LASTEXITCODE
Pop-Location
if (-not (Test-Path (Join-Path $d 'fields_40.csv'))) { throw "post failed (exit $code), see $d\out_post.txt" }
Set-Location $repo
. (Join-Path $repo 'dev-env.ps1') | Out-Null
& python ansys_usermat\apdl\fields_living.py $Job $d (Join-Path $WorkDir "$Job.dat") 40

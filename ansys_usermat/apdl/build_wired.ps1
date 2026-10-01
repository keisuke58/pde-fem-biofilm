<#
build_wired.ps1 -- rebuild the partner-element ANSYS.exe in F:\biofilm_upf_wired
after paste_fragments.py, with the module files compiled first:
  usermat_py_hook.f  (module biofilm_py_bridge)
  split_rates.f      (module biofilm_split)
  everything else    (*.f / *.F, minus the unused Conection_Test usermat)
Then checks that the new ANSYS.exe exists and is newer than the start.

  .\ansys_usermat\apdl\build_wired.ps1 [-WorkDir F:\biofilm_upf_wired]
#>
param([string]$WorkDir = 'F:\biofilm_upf_wired')
$ErrorActionPreference = 'Stop'
$first = @('usermat_py_hook.f', 'split_rates.f')
$skip  = @('Usermat_P21-V21_Conection_Test.F') + $first
foreach ($f in $first) {
    if (-not (Test-Path (Join-Path $WorkDir $f))) { throw "missing $f in $WorkDir" }
}
$rest = Get-ChildItem $WorkDir -File |
    Where-Object { $_.Extension -in '.f', '.F' -and $skip -notcontains $_.Name } |
    ForEach-Object Name
$t0 = Get-Date
# retries: the link intermittently ends with no new ANSYS.exe and no error
# message (1 Oct: once right after copying the 390 MB exe, once twice in a
# row; the next attempt passed) -- up to 3 attempts, 30 s apart
for ($try = 1; $try -le 3; $try++) {
    & (Join-Path $PSScriptRoot 'link_v222.ps1') -WorkDir $WorkDir -Sources ($first + $rest)
    $exe = Get-Item (Join-Path $WorkDir 'ANSYS.exe') -ErrorAction SilentlyContinue
    if ($exe -and $exe.LastWriteTime -ge $t0) { break }
    if ($try -lt 3) { "link produced no new ANSYS.exe (attempt $try); retrying in 30 s"; Start-Sleep 30 }
}
if (-not $exe -or $exe.LastWriteTime -lt $t0) { throw 'ANSYS.exe was not rebuilt' }
"OK: $($exe.FullName)  $($exe.Length) bytes  $($exe.LastWriteTime)"

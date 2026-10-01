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
# one retry: right after the 390 MB exe has been copied (e.g. as a backup) the
# link failed once with no error message and passed on the retry (1 Oct)
for ($try = 1; $try -le 2; $try++) {
    & (Join-Path $PSScriptRoot 'link_v222.ps1') -WorkDir $WorkDir -Sources ($first + $rest)
    $exe = Get-Item (Join-Path $WorkDir 'ANSYS.exe') -ErrorAction SilentlyContinue
    if ($exe -and $exe.LastWriteTime -ge $t0) { break }
    if ($try -eq 1) { 'link produced no new ANSYS.exe; retrying once in 10 s'; Start-Sleep 10 }
}
if (-not $exe -or $exe.LastWriteTime -lt $t0) { throw 'ANSYS.exe was not rebuilt' }
"OK: $($exe.FullName)  $($exe.Length) bytes  $($exe.LastWriteTime)"

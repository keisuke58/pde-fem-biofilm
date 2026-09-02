<#
.SYNOPSIS
    Non-interactive compile + link of a v222 custom UPF ANSYS.exe, bypassing
    ANSCUST.BAT's interactive prompts entirely.

.DESCRIPTION
    ANSCUST.BAT (the ANSYS-supplied UPF build script) is genuinely
    interactive -- it uses a bundled ASK.EXE that reads the console
    directly, not stdin, so it cannot be scripted or piped (confirmed in
    RUNBOOK.md). This script reproduces what ANSCUST.BAT does for the
    compile and link steps only, using the exact macro/flag set read out of
    ANSCUST.BAT itself (CUSTMACROS/FMACS/FSWITCH), so no human needs to be
    at the console for a rebuild.

    Two environment-setup bugs had to be worked around to get here (see
    V222_PORT_INSTRUCTIONS.md 1.6 and the 2026-09-02 link section for the
    full story):

    1. `setvars.bat` and even a bare `vcvars64.bat` call shell out to
       vswhere.exe, which is not itself on PATH
       (C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe).
       Without it, vcvars64.bat prints "[vcvarsall.bat] Environment
       initialized" and looks like it succeeded, but LIB/INCLUDE are left
       unset -- this script puts vswhere.exe's directory on PATH first.
    2. Chaining `call vcvars64.bat && set LIB=...;%LIB%&& ...` in one
       cmd /c line expands %LIB% at PARSE time (before vcvars64.bat has run),
       silently discarding everything vcvars64.bat set. This script writes
       a small temp .bat file instead, so each line's variable expansion
       happens after the previous line actually executed.
    3. A stale ANSYS.exe/.lib/.exp/.map from a prior link attempt gets
       swept up by ansys.lrf's own `*.lib` wildcard and collides with the
       new output (LNK1149). This script deletes all four before every
       link, since ANSCUST.BAT itself only auto-deletes ANSYS.exe.

.PARAMETER WorkDir
    Directory containing the .obj files to link (and, for a fresh build,
    the source files to compile first). Must already contain a copy of a
    working custom-UPF template directory's DLLs, ansys.lrf, ansysex.def,
    and app.manifest (e.g. copied from a prior successful ANSCUST.BAT run --
    see RUNBOOK.md step 0c/1). This script does not create that scaffolding.

.PARAMETER Sources
    Optional list of .f/.F source files (in WorkDir) to compile before
    linking. Omit to link whatever .obj files are already in WorkDir.

.PARAMETER MklInclude
    Include path for MKL's .fi Fortran interface files, only needed if a
    source being compiled uses MKL (e.g. PARDISO). Defaults to the payload
    extracted per V222_PORT_INSTRUCTIONS.md's MKL section
    (F:\mkl_payload\_installdir\mkl\2026.1\include) if present.

.PARAMETER MklLib
    Matching .lib directory for the link step. Same default pattern as
    -MklInclude.

.EXAMPLE
    # Compile everything in a working dir and link a fresh ANSYS.exe
    .\link_v222.ps1 -WorkDir F:\biofilm_upf_oliver -Sources (Get-ChildItem F:\biofilm_upf_oliver -Filter *.f,*.F)

.EXAMPLE
    # Just relink from existing .obj files (e.g. after copying new ones in)
    .\link_v222.ps1 -WorkDir F:\biofilm_upf_link
#>
param(
    [Parameter(Mandatory=$true)][string]$WorkDir,
    [string[]]$Sources = @(),
    [string]$MklInclude = "F:\mkl_payload\_installdir\mkl\2026.1\include",
    [string]$MklLib = "F:\mkl_payload\_installdir\mkl\2026.1\lib",
    [string]$AnsysRoot = $env:AWP_ROOT222,
    [string]$VsVcvars64 = "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat",
    [string]$IntelCompilerVars = "C:\Program Files (x86)\Intel\oneAPI\compiler\2025.3\env\vars.bat",
    [string]$VsWhereDir = "C:\Program Files (x86)\Microsoft Visual Studio\Installer"
)

$ErrorActionPreference = "Stop"

if (-not $AnsysRoot) { throw "AWP_ROOT222 is not set and -AnsysRoot was not given." }
if (-not (Test-Path $WorkDir)) { throw "WorkDir does not exist: $WorkDir" }

# A stale ANSYS.lib/.exp/.map from a prior link (successful or not) gets
# swept up by ansys.lrf's own `*.lib` wildcard and collides with the new
# ANSYS.exe/.lib output (LNK1149: output filename identical to input) --
# see RUNBOOK.md's documented fix. ANSCUST.BAT only auto-deletes ANSYS.exe,
# not the others, so this must be done explicitly on every relink.
foreach ($ext in "exe","lib","exp","map") {
    $stale = Join-Path $WorkDir "ANSYS.$ext"
    if (Test-Path $stale) { Remove-Item $stale -Force }
}

$incAnsys = Join-Path $AnsysRoot "ansys\customize\include"
$incMpi = Join-Path $AnsysRoot "commonfiles\MPI\Intel\2021.6.0\winx64\include"
$libAnsysCustom = Join-Path $AnsysRoot "ansys\Custom\Lib\winx64"

$macros = "/DNOSTDCALL /DARGTRAIL /DPCWIN64_SYS /DPCWINX64_SYS /DPCWINNT_SYS /DCADOE_ANSYS"
$fmacs = "/D__EFL /DFORTRAN"
$switch = "/O2 /fpp /4Yportlib /auto /c /Fo.\ /MD /watch:source"

$batPath = Join-Path $WorkDir "_link_v222_tmp.bat"
$logPath = Join-Path $WorkDir "_link_v222_tmp.log"

$lines = @(
    "@echo off",
    "set PATH=$VsWhereDir;%PATH%",
    "call `"$VsVcvars64`" >nul 2>&1",
    "call `"$IntelCompilerVars`" >nul 2>&1",
    "set LIB=$libAnsysCustom;$MklLib;%LIB%",
    "cd /d `"$WorkDir`"",
    "echo === LIB === > `"$logPath`"",
    "echo %LIB%>> `"$logPath`""
)

if ($Sources.Count -gt 0) {
    $srcList = ($Sources -join " ")
    $lines += "echo === COMPILE ===>> `"$logPath`""
    $lines += "ifort /nologo $macros $fmacs $switch /I`"$incAnsys`" /I`"$incMpi`" /I`"$MklInclude`" $srcList >> `"$logPath`" 2>&1"
}

$lines += "echo === LINK ===>> `"$logPath`""
$lines += "link @ansys.lrf >> `"$logPath`" 2>&1"

$lines -join "`r`n" | Set-Content -Path $batPath -Encoding ASCII

& cmd /c "`"$batPath`""
$log = Get-Content $logPath -Raw
Write-Output $log

$exePath = Join-Path $WorkDir "ANSYS.exe"
if (Test-Path $exePath) {
    $item = Get-Item $exePath
    Write-Output ""
    Write-Output "ANSYS.exe produced: $($item.Length) bytes, $($item.LastWriteTime)"
    exit 0
} else {
    Write-Output ""
    Write-Output "No ANSYS.exe produced -- check the log above for the fatal error."
    exit 1
}

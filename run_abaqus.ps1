<#
.SYNOPSIS
    Run one Abaqus job (optionally with a user subroutine) in an F:-based
    working directory, never on C:.

.DESCRIPTION
    Mirrors ansys_usermat/apdl/run_apdl.ps1's reasoning: C: on IKMHIWI03 is
    chronically near-full (VSS, admin needed to actually fix -- see the
    disk-space memory), while F:\ is a separate local data drive with
    ~3.7 TB free. Abaqus scratch (.odb/.sta/.msg/.dat/.prt/.com/.SMABulk/...)
    and the user-subroutine compile step (Fortran + MSVC, its own build
    directory) both land wherever the job is launched from -- so this script
    always launches from an F:-based per-job directory instead of wherever
    the .inp happens to live in the repo (which is on C:).

    Only the .inp (and the user Fortran source, if given) are copied to F:;
    nothing is copied back automatically, since Abaqus results (.odb
    especially) can be large and are usually inspected in place. The
    -WorkDir path is printed at the end.

.PARAMETER InpFile
    Path to the .inp file, relative to the repo root (this script's own
    directory), e.g. umat_tangent_test/abaqus_1elem/umat_visco_1elem.inp

.PARAMETER UserFile
    Path to a Fortran user-subroutine source, relative to the repo root,
    e.g. umat_biofilm_visco.f. Optional -- omit for a job with no UMAT.

.PARAMETER WorkDir
    F:-based root for per-job working directories. Each job gets its own
    subdirectory named after the .inp's basename, so repeat runs don't
    collide. Default F:\abaqus_work.

.PARAMETER Cpus
    Passed through to Abaqus's cpus= argument. Default 1 (these are meant
    to be small single-element comparison jobs, not production runs --
    see the Abaqus-compute-location memory: heavy jobs run at Keio, not
    here).

.EXAMPLE
    .\run_abaqus.ps1 -InpFile umat_tangent_test/abaqus_1elem/umat_visco_1elem.inp -UserFile umat_biofilm_visco.f
#>
param(
    [Parameter(Mandatory=$true)][string]$InpFile,
    [string]$UserFile,
    [string]$WorkDir = "F:\abaqus_work",
    [int]$Cpus = 1
)

$repoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$inpPath = Join-Path $repoRoot $InpFile
if (-not (Test-Path $inpPath)) {
    throw "InpFile not found: $inpPath"
}
if (-not (Get-Command abaqus -ErrorAction SilentlyContinue)) {
    throw "abaqus not found on PATH -- expected C:\SIMULIA\Commands on this machine."
}

if ($UserFile -and -not (Get-Command ifort -ErrorAction SilentlyContinue)) {
    # A user subroutine needs ifort on PATH to compile. Plain vcvars64.bat +
    # Intel's env\vars.bat via setvars.bat does NOT reliably work on this
    # machine (see RUNBOOK.md Step 0b -- setvars.bat shells out to a bare
    # vswhere.exe that isn't itself on PATH). The sequence that does work,
    # verified when building the ANSYS USERMAT, calls the two env scripts
    # directly. Batch files set env vars in a cmd.exe child process that
    # doesn't propagate to PowerShell on its own, so run them in cmd, dump
    # the resulting environment with `set`, and import each VAR=VALUE line.
    Write-Output "ifort not on PATH -- initialising the Intel/VS toolchain (see RUNBOOK.md Step 0b) ..."
    $vcvars = "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat"
    $ifortEnv = "C:\Program Files (x86)\Intel\oneAPI\compiler\2025.3\env\vars.bat"
    if (-not (Test-Path $vcvars) -or -not (Test-Path $ifortEnv)) {
        throw "Expected toolchain paths not found ($vcvars / $ifortEnv) -- re-check RUNBOOK.md Step 0b, versions may have changed."
    }
    $envDump = cmd /c "call `"$vcvars`" && call `"$ifortEnv`" && set" 2>&1
    foreach ($line in $envDump) {
        if ($line -match '^([^=]+)=(.*)$') {
            [System.Environment]::SetEnvironmentVariable($matches[1], $matches[2], "Process")
        }
    }
    if (-not (Get-Command ifort -ErrorAction SilentlyContinue)) {
        throw "ifort still not found after toolchain init -- something changed vs. RUNBOOK.md Step 0b."
    }
    Write-Output "ifort on PATH now: $((Get-Command ifort).Source)"
}

$jobName = [IO.Path]::GetFileNameWithoutExtension($InpFile)
$jobDir = Join-Path $WorkDir $jobName
New-Item -ItemType Directory -Path $jobDir -Force | Out-Null
Copy-Item $inpPath $jobDir -Force

$userArg = @()
if ($UserFile) {
    $userPath = Join-Path $repoRoot $UserFile
    if (-not (Test-Path $userPath)) {
        throw "UserFile not found: $userPath"
    }
    # abaqus user= takes a path; pointing it straight at the repo copy is
    # fine (Abaqus only reads it to compile), no need to copy the source in.
    $userArg = @("user=$userPath")
}

Write-Output "== abaqus job=$jobName (in $jobDir) =="
Push-Location $jobDir
try {
    $prevEAP = $ErrorActionPreference
    $ErrorActionPreference = "Continue"
    & abaqus job=$jobName input="$jobName.inp" @userArg cpus=$Cpus interactive ask_delete=OFF 2>&1 |
        ForEach-Object { Write-Output $_ }
    $abqExit = $LASTEXITCODE
    $ErrorActionPreference = $prevEAP
} finally {
    Pop-Location
}

$staPath = Join-Path $jobDir "$jobName.sta"
if (Test-Path $staPath) {
    if (Select-String -Path $staPath -Pattern "COMPLETED SUCCESSFULLY" -Quiet) {
        Write-Output "PASS: $jobName COMPLETED SUCCESSFULLY"
    } else {
        Write-Output "NOT COMPLETE: check $staPath"
    }
} else {
    Write-Output "No .sta produced -- job likely failed before starting the analysis. abaqus exit code: $abqExit"
}
Write-Output "Job directory: $jobDir"

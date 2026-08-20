<#
.SYNOPSIS
    Run the Python test suite with this machine's toolchain on PATH and the
    known-environment-limited tests skipped, instead of retyping the same
    flags by hand each time.

.DESCRIPTION
    Two tests are excluded here, both confirmed environment limitations of
    THIS machine, not code bugs (see CLAUDE.md / memory for the
    investigation): tests/test_coupling_shim.py needs POSIX socket headers
    (arpa/inet.h) that this machine's native-Windows mingw toolchain does
    not ship, and tests/test_viscoelastic*.py need scipy, which is not
    installed here (installing it was judged too disk-risky on a machine
    that has hit 0 bytes free -- see the disk-space memory). Both run fine
    in CI (ubuntu-latest, full requirements.txt).

    -All also runs the dual-solver equivalence suite
    (ansys_usermat/crosscheck/), which is otherwise easy to forget since
    it lives outside tests/.

.PARAMETER All
    Also run ansys_usermat/crosscheck/crosscheck.py and adversarial.py.

.EXAMPLE
    .\run_tests.ps1
    .\run_tests.ps1 -All
#>
param(
    [switch]$All
)

$ErrorActionPreference = "Continue"
$repoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
. (Join-Path $repoRoot "dev-env.ps1")

Push-Location $repoRoot
try {
    Write-Output "== pytest tests/ (skipping scipy- and POSIX-socket-limited tests) =="
    python -m pytest tests/ -q --tb=short --continue-on-collection-errors `
        --ignore=tests/test_coupling_shim.py `
        --deselect=tests/test_viscoelastic.py::test_2d_viscoelastic_relaxation `
        --deselect=tests/test_viscoelastic.py::test_2d_viscoelastic_time_history `
        --ignore=tests/test_viscoelastic_2d.py
    $testExit = $LASTEXITCODE

    if ($All) {
        Write-Output ""
        Write-Output "== ansys_usermat/crosscheck (dual-solver equivalence) =="
        python -m pytest ansys_usermat/crosscheck/crosscheck.py ansys_usermat/crosscheck/adversarial.py -q --tb=short
        $crossExit = $LASTEXITCODE
    } else {
        $crossExit = 0
    }
} finally {
    Pop-Location
}

Write-Output ""
if ($testExit -eq 0 -and $crossExit -eq 0) {
    Write-Output "run_tests.ps1: all green (within this machine's known limitations -- see script header)."
} else {
    Write-Output "run_tests.ps1: FAILURES -- tests exit=$testExit, crosscheck exit=$crossExit (see output above)."
}
exit ([int]($testExit -ne 0 -or $crossExit -ne 0))

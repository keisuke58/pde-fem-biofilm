<#
.SYNOPSIS
    Dot-source this to put git (MSYS64), the portable Python, and the
    portable gfortran on PATH for the current PowerShell command.

.DESCRIPTION
    This harness's shell state (env vars, PATH) does not persist between
    tool calls -- only the working directory does -- so every command that
    needs git/python/gfortran has been re-prepending PATH by hand, all
    session, three separate paths retyped each time. This is that one line,
    reusable.

    Must be DOT-SOURCED, not just run, so the PATH change applies to the
    calling shell rather than a child scope that exits immediately:

        . .\dev-env.ps1
        python -m pytest tests/ -q

    or inline in one command:

        . .\dev-env.ps1; python -m pytest tests/ -q

.NOTES
    Paths are this machine's actual install locations (see CLAUDE.md /
    memory for how they got there): MSYS64 git (no git on PATH otherwise),
    a per-user Python 3.12.7 install, and a portable WinLibs gfortran.
    Silently skips any of the three that isn't present, rather than
    erroring, so this stays usable if one toolchain is missing.
#>

$paths = @(
    "C:\msys64\usr\bin",
    "C:\Users\nishioka\AppData\Local\Programs\Python\Python312",
    "C:\Users\nishioka\AppData\Local\Programs\Python\Python312\Scripts",
    "C:\Users\nishioka\tools\mingw64\bin"
) | Where-Object { Test-Path $_ }

$env:Path = ($paths -join ";") + ";" + $env:Path
Write-Output ("dev-env: added " + $paths.Count + " path(s) -- " + ($paths -join ", "))

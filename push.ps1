<#
.SYNOPSIS
    Push the current branch to origin using the GitHub PAT in .env, without
    ever printing the token.

.DESCRIPTION
    This machine has no credential helper configured (see CLAUDE.md), so a
    plain `git push` prompts for username/password on the console and hangs
    in a non-interactive session. The working pattern all session has been:
    read GITHUB_PAT from .env via raw file I/O (never the Read tool, so it
    never lands in a transcript), build an inline
    https://x-access-token:<token>@github.com/... push URL, push, then
    redact the token from anything printed. This script is that pattern,
    written once instead of retyped by hand for every push.

    Also prepends C:\msys64\usr\bin to PATH for the call, since there is no
    git on PATH otherwise on this machine (again, see CLAUDE.md) -- without
    it git's https/credential-helper subprocesses fail with "shared
    libraries" errors.

.PARAMETER Branch
    Branch to push. Defaults to the current branch.

.PARAMETER Remote
    Repo URL host path, e.g. "keisuke58/pde-fem-biofilm". Defaults to what
    this repo's origin already points at.

.EXAMPLE
    .\push.ps1
    .\push.ps1 -Branch feature/foo
#>
param(
    [string]$Branch,
    [string]$Remote = "keisuke58/pde-fem-biofilm"
)

$ErrorActionPreference = "Stop"
$repoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$envPath = Join-Path $repoRoot ".env"
if (-not (Test-Path $envPath)) {
    throw ".env not found at $envPath -- expected a GITHUB_PAT=... line there."
}

$env:Path = "C:\msys64\usr\bin;" + $env:Path

if (-not $Branch) {
    $Branch = (git -C $repoRoot rev-parse --abbrev-ref HEAD).Trim()
}

$envContent = [System.IO.File]::ReadAllText($envPath)
$tokenLine = ($envContent -split "`n" | Where-Object { $_ -match '^GITHUB_PAT=' })
if (-not $tokenLine) {
    throw "No GITHUB_PAT= line found in .env"
}
$token = ($tokenLine -replace '^GITHUB_PAT=', '').Trim()
if (-not $token) {
    throw "GITHUB_PAT in .env is empty"
}

$pushUrl = "https://x-access-token:$token@github.com/$Remote.git"

Write-Output "Pushing branch '$Branch' to $Remote ..."
$result = git -C $repoRoot push $pushUrl $Branch 2>&1 | Out-String
Write-Output ($result -replace [regex]::Escape($token), '***')

# Keep the local origin/<branch> tracking ref in sync so `git status` doesn't
# report stale "ahead by N commits" after a push done via this inline URL
# instead of the configured 'origin' remote.
git -C $repoRoot fetch origin $Branch 2>&1 | Out-Null

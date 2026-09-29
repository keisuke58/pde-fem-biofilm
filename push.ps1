<#
.SYNOPSIS
    Push HEAD to origin/master with the PAT from .env, never printing the
    token, then confirm against the GitHub API that the remote really moved.

.DESCRIPTION
    No credential helper is configured on this machine, and the local
    tracking ref refs/remotes/origin/master has a persistent rename-lock
    (CLAUDE.md), so `git status` cannot be trusted to say whether a push
    landed. This script:
      1. reads GITHUB_PAT from .env (repo root)
      2. pushes HEAD to master over https with the token in the URL
      3. prints git's output with the token redacted
      4. asks the GitHub API for the remote master sha and compares it
         with the local HEAD
    git's normal progress output goes to stderr, which PowerShell 5.1 turns
    into error records; the push is judged by git's exit code and by the
    API comparison, not by stderr.

.PARAMETER Branch
    Remote branch to push to. Default master.

.EXAMPLE
    .\push.ps1
#>
param([string]$Branch = "master")

$repoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$env:Path = "C:\Users\nishioka\git\cmd;C:\Users\nishioka\git\mingw64\bin;$env:Path"
Push-Location $repoRoot
try {
    $envFile = Join-Path $repoRoot ".env"
    if (-not (Test-Path $envFile)) { throw ".env not found at $envFile" }
    $line = Get-Content $envFile | Where-Object { $_ -match '^\s*GITHUB_PAT\s*=' } | Select-Object -First 1
    if (-not $line) { throw "GITHUB_PAT not set in .env" }
    $tok = ($line -replace '^\s*GITHUB_PAT\s*=\s*', '').Trim().Trim('"').Trim("'")

    $head = (& git rev-parse HEAD).Trim()
    Write-Output "== push $($head.Substring(0,7)) -> origin/$Branch =="
    $url = "https://x-access-token:$tok@github.com/keisuke58/pde-fem-biofilm.git"
    $prev = $ErrorActionPreference
    $ErrorActionPreference = "Continue"
    $out = & git push $url "HEAD:$Branch" 2>&1 | Out-String
    $code = $LASTEXITCODE
    $ErrorActionPreference = $prev
    ($out -replace [regex]::Escape($tok), '***').Trim() -split "`n" |
        ForEach-Object { $_ -replace '^\s*git(\.exe)?\s*:\s*', '' } |
        Where-Object { $_ -notmatch '^\s*(\+ |In Zeile|In line|In .*(Zeichen|char)|CategoryInfo|FullyQualifiedErrorId|~)' -and $_.Trim() } |
        ForEach-Object { Write-Output "  $($_.Trim())" }
    if ($code -ne 0) { throw "git push exited $code" }

    $remote = (Invoke-RestMethod "https://api.github.com/repos/keisuke58/pde-fem-biofilm/commits/$Branch").sha
    if ($remote -eq $head) {
        Write-Output "OK: GitHub $Branch = $($remote.Substring(0,7)) (matches local HEAD)"
    } else {
        Write-Output "MISMATCH: GitHub $Branch = $($remote.Substring(0,7)), local HEAD = $($head.Substring(0,7))"
    }
} finally {
    Remove-Variable tok, url -ErrorAction SilentlyContinue
    Pop-Location
}

<#
.SYNOPSIS
    Stage exactly the named files and commit them, with the message written
    as UTF-8 without a BOM.

.DESCRIPTION
    The two rules this repo's working tree forces on every commit
    (CLAUDE.md): never `git add -A` / `commit -a` -- ~480 files show pure
    CRLF/LF churn -- and write the message through a file, because
    PowerShell 5.1's Set-Content/Out-File put a BOM in front of it that
    then shows up in the commit subject. This script does both, refuses to
    commit if anything other than the named files is already staged, and
    refuses an AI-attribution trailer.

.PARAMETER Files
    Paths to stage (repo-relative). Deletions are staged too.

.PARAMETER Message
    Full commit message: subject line, blank line, body.

.PARAMETER Push
    Also run push.ps1 afterwards.

.EXAMPLE
    .\commit.ps1 -Files thesis_ch5/ch5_ansys_contribution.tex -Message "Ch. 5: fix a typo"
    .\commit.ps1 -Files a.py,b.png -Message $msg -Push
#>
param(
    [Parameter(Mandatory=$true)][string[]]$Files,
    [Parameter(Mandatory=$true)][string]$Message,
    [switch]$Push
)

$repoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$env:Path = "C:\Users\nishioka\git\cmd;C:\Users\nishioka\git\mingw64\bin;$env:Path"
Push-Location $repoRoot
try {
    if ($Message -match '(?im)^\s*Co-Authored-By:.*(claude|anthropic)') {
        throw "Message carries an AI co-author trailer -- not allowed in this repo (CLAUDE.md)."
    }
    $pre = @(& git diff --cached --name-only)
    $want = $Files | ForEach-Object { $_.Replace('\', '/') }
    $stray = $pre | Where-Object { $_ -and ($want -notcontains $_) }
    if ($stray) {
        throw "Already staged but not named: $($stray -join ', ') -- unstage them or name them."
    }

    $prev = $ErrorActionPreference
    $ErrorActionPreference = "Continue"
    & git add -A -- $Files 2>&1 | Where-Object { $_ -notmatch 'LF will be replaced|CRLF will be replaced' } |
        ForEach-Object { Write-Output "  $_" }
    $ErrorActionPreference = $prev
    if ($LASTEXITCODE -ne 0) { throw "git add failed" }

    $staged = @(& git diff --cached --name-only)
    if (-not $staged) { Write-Output "Nothing to commit for the named files."; return }
    $tmp = Join-Path $repoRoot ".git\COMMIT_EDITMSG_TMP"
    $text = $Message.Replace("`r`n", "`n")
    if (-not $text.EndsWith("`n")) { $text += "`n" }
    [IO.File]::WriteAllText($tmp, $text, (New-Object Text.UTF8Encoding($false)))
    & git commit -q -F $tmp
    if ($LASTEXITCODE -ne 0) { throw "git commit failed" }
    Remove-Item $tmp -ErrorAction SilentlyContinue
    Write-Output ("committed " + (& git log --oneline -1) + "  [" + ($staged -join ', ') + "]")
} finally {
    Pop-Location
}
if ($Push) { & (Join-Path $repoRoot "push.ps1") }

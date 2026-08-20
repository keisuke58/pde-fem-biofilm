<#
.SYNOPSIS
    Build slides_1005.tex, check it against the 20-page cap, and clean up
    the LaTeX intermediate files -- the sequence that's been done by hand
    (pdflatex, grep the page count, rm the aux/log/nav/out/snm/toc files)
    every time a frame was added this session.

.PARAMETER Tex
    The .tex file to build. Defaults to slides_1005.tex.

.PARAMETER MaxPages
    Page-count cap to warn against. Defaults to 20 (see slides_1005.tex's
    own header comment -- keep this in sync if that constraint ever changes).

.PARAMETER KeepPdf
    If set, also copies the built PDF next to the script with a
    timestamped name, for quick before/after comparison. Off by default --
    the PDF stays wherever pdflatex put it either way.

.EXAMPLE
    .\build_slides.ps1
    .\build_slides.ps1 -Tex other_deck.tex -MaxPages 15
#>
param(
    [string]$Tex = "slides_1005.tex",
    [int]$MaxPages = 20,
    [switch]$KeepPdf
)

$ErrorActionPreference = "Stop"
$repoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$texPath = Join-Path $repoRoot $Tex
if (-not (Test-Path $texPath)) {
    throw "Not found: $texPath"
}
$base = [IO.Path]::GetFileNameWithoutExtension($Tex)

Push-Location $repoRoot
try {
    Write-Output "== pdflatex $Tex =="
    $prevEAP = $ErrorActionPreference
    $ErrorActionPreference = "Continue"
    $out = & pdflatex -interaction=nonstopmode $Tex 2>&1 | Out-String
    $texExit = $LASTEXITCODE
    $ErrorActionPreference = $prevEAP

    $errLines = $out -split "`n" | Where-Object { $_ -match '^!' }
    if ($errLines) {
        Write-Output "LaTeX errors:"
        $errLines | ForEach-Object { Write-Output "  $_" }
    }
    if ($texExit -ne 0 -and -not $errLines) {
        Write-Output "pdflatex exited $texExit with no '!' error lines -- check the log if the PDF looks wrong."
    }

    $pdfPath = Join-Path $repoRoot "$base.pdf"
    if (Test-Path $pdfPath) {
        $pageMatch = Select-String -InputObject $out -Pattern "Output written on .* \((\d+) pages?"
        if ($pageMatch) {
            $pages = [int]$pageMatch.Matches[0].Groups[1].Value
            Write-Output "Pages: $pages (cap $MaxPages)"
            if ($pages -gt $MaxPages) {
                Write-Output "OVER THE CAP by $($pages - $MaxPages) page(s) -- trim a frame before adding more."
            }
        } else {
            Write-Output "Built, but couldn't parse the page count from pdflatex output -- check manually."
        }
        if ($KeepPdf) {
            $stamp = Get-Date -Format "yyyyMMdd_HHmmss"
            Copy-Item $pdfPath (Join-Path $repoRoot "${base}_$stamp.pdf")
            Write-Output "Also saved a timestamped copy: ${base}_$stamp.pdf"
        }
    } else {
        Write-Output "No PDF was produced -- build failed. See errors above."
    }
} finally {
    Remove-Item (Join-Path $repoRoot "$base.aux"), (Join-Path $repoRoot "$base.log"),
                (Join-Path $repoRoot "$base.nav"), (Join-Path $repoRoot "$base.out"),
                (Join-Path $repoRoot "$base.snm"), (Join-Path $repoRoot "$base.toc") `
                -Force -ErrorAction SilentlyContinue
    Pop-Location
}

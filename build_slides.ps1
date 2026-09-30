<#
.SYNOPSIS
    Build a LaTeX deck (or the chapter-5 check document), report errors,
    frames that run off the page and the page count, and clean up the
    intermediate files.

.DESCRIPTION
    Runs the engine twice (references), then reads the .log -- not the
    console -- for:
      * '!' error lines (a build can still write a PDF past an error;
        that is how a broken ch. 5 went unnoticed once)
      * 'Overfull \vbox' -- a beamer frame whose content runs off the page
      * undefined references (second pass only)
      * the page count, against -MaxPages
    The engine is lualatex for files that load luatexja (the Japanese
    deck), pdflatex otherwise. The build runs in the .tex file's own
    directory, so thesis_ch5/_build_check.tex resolves its relative paths.

.PARAMETER Tex
    .tex file(s), repo-relative. Default slides_1005.tex.

.PARAMETER All
    Build the three meeting decks: slides_1001_oliver.tex, slides_1005.tex,
    slides_muramatsu_ja.tex.

.PARAMETER Chapter
    Build thesis_ch5/_build_check.tex (chapter 5 alone) and delete its PDF
    afterwards -- a check, not a deliverable.

.PARAMETER MaxPages
    Page cap to warn against. Default 20 (slides_1005.tex's own limit).

.PARAMETER KeepPdf
    Also save a timestamped copy of each PDF for before/after comparison.

.EXAMPLE
    .\build_slides.ps1 -All
    .\build_slides.ps1 -Chapter
    .\build_slides.ps1 -Tex slides_1001_oliver.tex -MaxPages 12
#>
param(
    [string[]]$Tex = @("slides_1005.tex"),
    [switch]$All,
    [switch]$Chapter,
    [int]$MaxPages = 20,
    [switch]$KeepPdf
)

$repoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
if ($All) { $Tex = @("slides_1001_oliver.tex", "slides_1005.tex", "slides_muramatsu_ja.tex") }
if ($Chapter) { $Tex = @("thesis_ch5/_build_check.tex") }

$bad = 0
foreach ($t in $Tex) {
    $texPath = Join-Path $repoRoot $t
    if (-not (Test-Path $texPath)) { Write-Output "Not found: $t"; $bad++; continue }
    $dir = Split-Path -Parent $texPath
    $name = Split-Path -Leaf $texPath
    $base = [IO.Path]::GetFileNameWithoutExtension($name)
    $engine = if (Select-String -Path $texPath -Pattern 'luatexja' -Quiet) { "lualatex" } else { "pdflatex" }

    Write-Output "== $engine $t =="
    Push-Location $dir
    try {
        $prev = $ErrorActionPreference
        $ErrorActionPreference = "Continue"
        & $engine -interaction=nonstopmode $name 2>&1 | Out-Null
        & $engine -interaction=nonstopmode $name 2>&1 | Out-Null
        $ErrorActionPreference = $prev

        $log = "$base.log"
        if (-not (Test-Path $log)) { Write-Output "  no log written -- engine did not run"; $bad++; continue }
        $errs = Select-String -Path $log -Pattern '^!' | ForEach-Object { $_.Line.Trim() }
        $over = Select-String -Path $log -Pattern 'Overfull \\vbox \(([\d.]+)pt too high\) detected at line (\d+)' |
            ForEach-Object { "{0}pt at line {1}" -f $_.Matches[0].Groups[1].Value, $_.Matches[0].Groups[2].Value }
        $undef = Select-String -Path $log -Pattern "(Reference|Citation) [``']([^']+)' .*undefined" |
            ForEach-Object { $_.Matches[0].Groups[2].Value } | Sort-Object -Unique
        $pages = Select-String -Path $log -Pattern 'Output written on .* \((\d+) pages?' |
            ForEach-Object { [int]$_.Matches[0].Groups[1].Value } | Select-Object -Last 1

        if ($errs) { $bad++; Write-Output "  ERRORS:"; $errs | ForEach-Object { Write-Output "    $_" } }
        if ($over) { $bad++; Write-Output "  RUNS OFF THE PAGE (frame ending at):"; $over | ForEach-Object { Write-Output "    $_" } }
        if ($undef) { Write-Output "  undefined (expected in the chapter check for other chapters' labels): $($undef -join ', ')" }
        if ($pages) {
            $cap = if ($Chapter) { "" } else { " (cap $MaxPages)" }
            Write-Output "  pages: $pages$cap"
            if (-not $Chapter -and $pages -gt $MaxPages) { $bad++; Write-Output "  OVER THE CAP by $($pages - $MaxPages)" }
        } else { $bad++; Write-Output "  no PDF produced" }
        if (-not $errs -and -not $over -and $pages) { Write-Output "  OK" }

        if ($KeepPdf -and (Test-Path "$base.pdf")) {
            Copy-Item "$base.pdf" ("{0}_{1}.pdf" -f $base, (Get-Date -Format "yyyyMMdd_HHmmss"))
        }
    } finally {
        Remove-Item "$base.aux", "$base.log", "$base.nav", "$base.out", "$base.snm", "$base.toc" `
            -Force -ErrorAction SilentlyContinue
        if ($Chapter) { Remove-Item "$base.pdf" -Force -ErrorAction SilentlyContinue }
        Pop-Location
    }
}
if ($bad) { Write-Output "== $bad problem(s) =="; exit 1 } else { Write-Output "== all clean ==" }

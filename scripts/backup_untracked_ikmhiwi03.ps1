<#
backup_untracked_ikmhiwi03.ps1 -- zip what IKMHIWI03 holds outside git before the machine goes at the end of December 2026
(the user, 22 Oct 2026). Plan and what goes where: BACKUP_PLAN_2026-12.ja.md.

  .\scripts\backup_untracked_ikmhiwi03.ps1                 # dry run: lists every zip, its files and sizes
  .\scripts\backup_untracked_ikmhiwi03.ps1 -Go [-Rst] [-Figdata]   # writes the zips to F:\backup_2026-12
                                                           # (-Rst: also the u8_ .rst; -Figdata: the Python fields figdata\*.npz)

Zips (Zip64 through System.IO.Compression, so files over 2 GB are fine):
  private_felix.zip     F:\felix_private as it is                       -- private: never to a public place
  private_element.zip   the partner element's sources of every build (.f .F .f90 .inc .h) and the patch folders
                                                                        -- private: never to a public place
  work_text.zip         per work folder: decks (.dat), scripts, macros, point-model files, lap_cons_*.txt, run and chain
                        logs, automatic chains (auto_*), small JSON; Abaqus: living_gp.txt, .sta, .msg, .dat; the
                        20 Oct backup; the Python runs' scripts and logs in F:\auto_*
  python_fields.zip     (-Figdata) the Python 3D fields figdata\*.npz (not figdata\ckpt)
  rst_main.zip         (-Rst) the .rst of the main runs (u8_*)
Left out on purpose (can be recomputed or are scratch): .rst/.odb except -Rst, comp_trace.csv, phi_trace.csv, front_dump,
nem_dump, out_*.txt over 20 MB, solver scratch (.esav .full .emat .r00* .nr* .page .osav .stt .mdl .prt .SMABulk ...).
#>
param([switch]$Go, [switch]$Rst, [switch]$Figdata, [string]$Out = 'F:\backup_2026-12')
Add-Type -AssemblyName System.IO.Compression, System.IO.Compression.FileSystem
$work = Get-ChildItem F:\ -Directory | Where-Object { $_.Name -like 'biofilm_upf*' }
$srcExt = '.f', '.F', '.f90', '.inc', '.h'
$textExt = '.dat', '.ps1', '.py', '.mac', '.txt', '.md', '.json', '.log', '.bat', '.patch', '.lrf'
$skip = '^(comp_trace|phi_trace|front_dump|nem_dump|ansys_L_trace)'
$sets = [ordered]@{}
$sets['private_felix.zip'] = @(Get-ChildItem F:\felix_private -Recurse -File -ErrorAction SilentlyContinue)
$sets['private_element.zip'] = @($work | ForEach-Object { Get-ChildItem $_.FullName -File | Where-Object { $_.Extension -cin $srcExt } }) +
    @(Get-ChildItem F:\oliver_element_patches* -Recurse -File -ErrorAction SilentlyContinue)
$sets['work_text.zip'] = @($work | ForEach-Object { Get-ChildItem $_.FullName -Recurse -File -ErrorAction SilentlyContinue |
        Where-Object { $_.Extension -in $textExt -and $_.Name -notmatch $skip -and -not ($_.Name -like 'out_*' -and $_.Length -gt 20MB) -and $_.Length -lt 500MB } |
        # lap_cons: only where it was built (the *_8n/_16n/_24n folders, nem*), not the copies in the run folders, not the bad 24^3 one
        Where-Object { $_.Name -notlike 'lap_cons*' -or (($_.Directory.Name -match '_\d+n$|^nem') -and $_.Name -notmatch 'BAD') } }) +
    @(Get-ChildItem F:\abaqus_work -Recurse -File -ErrorAction SilentlyContinue | Where-Object { $_.Name -eq 'living_gp.txt' -or $_.Extension -in '.sta', '.msg', '.dat', '.inp' } |
        Where-Object { $_.Length -lt 200MB }) +
    @(Get-ChildItem F:\backup_2026-10-20 -Recurse -File -ErrorAction SilentlyContinue) +
    # the Python runs' unattended scripts and logs (F:\auto_night, F:\auto_a33, ...), not only the auto_* inside the work folders
    @(Get-ChildItem F:\ -Directory -Filter 'auto_*' | ForEach-Object { Get-ChildItem $_.FullName -Recurse -File -ErrorAction SilentlyContinue |
        Where-Object { $_.Extension -in $textExt -and $_.Length -lt 500MB } })
# (-Figdata) the Python 3D fields figdata\*.npz (git ignores *.npz; figures and videos are drawn from them, 192^3 takes hours
# to recompute); not the checkpoints in figdata\ckpt. The zip keeps the path below the repository.
$repo = Split-Path $PSScriptRoot -Parent
if ($Figdata) { $sets['python_fields.zip'] = @(Get-ChildItem (Join-Path $repo 'figdata') -File -Filter '*.npz') }
# (-Rst) only the .rst behind the thesis numbers (user, 28 Oct; all u8_ would be 169 GB): the upwind front on 8^3/16^3/24^3
# (seed stress 2.14 / 4.56 / 5.15e-4 Pa, extrapolated 5.8e-4) and the constant added diffusion r h/4 on 16^3 (about a quarter)
$rstThesis = 'biofilm_upf_living11_8\up8_8_sep_ab.rst', 'biofilm_upf_living11_16\up16_16_sep_ab.rst',
             'biofilm_upf_living11_24\up24_24_sep_ab.rst', 'biofilm_upf_living8_16a\u8_16_sep_ab.rst'
if ($Rst) { $sets['rst_main.zip'] = @($rstThesis | ForEach-Object { Get-Item (Join-Path 'F:\' $_) -ErrorAction SilentlyContinue }) }
foreach ($k in $sets.Keys) {
    $fs = $sets[$k]
    "{0,-22} {1,6} files {2,9:N1} MB" -f $k, $fs.Count, (($fs | Measure-Object Length -Sum).Sum / 1MB)
}
if (-not $Go) { 'dry run; -Go writes the zips to ' + $Out; return }
New-Item -ItemType Directory -Force $Out | Out-Null
foreach ($k in $sets.Keys) {
    $zp = Join-Path $Out $k
    if (Test-Path $zp) { Remove-Item $zp }
    $z = [System.IO.Compression.ZipFile]::Open($zp, 'Create')
    try {
        foreach ($f in $sets[$k]) {
            $rel = if ($f.FullName.StartsWith($repo)) { $f.FullName.Substring($repo.Length + 1) } else { $f.FullName.Substring(3) }   # path below the repository or F:\
            [void][System.IO.Compression.ZipFileExtensions]::CreateEntryFromFile($z, $f.FullName, $rel, 'Optimal')
        }
    } finally { $z.Dispose() }
    "{0}: {1:N1} MB" -f $zp, ((Get-Item $zp).Length / 1MB)
}

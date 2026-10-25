# Download EQ3/6 v8.0a (software, data0 databases, test library) and the CHNOSZ OBIGT
# data files into a folder OUTSIDE the repository and OUTSIDE OneDrive.
#
#   powershell -ExecutionPolicy Bypass -File scripts\fetch_references.ps1
#
# Sources and licenses:
#   EQ3/6   https://github.com/LLNL/EQ3_6      software BSD-3-Clause; data0 file license unconfirmed
#   OBIGT   https://github.com/jedick/CHNOSZ   GPL-3.0
# Uses curl.exe because some Python installs here have a stale CA bundle.
param([string]$Dest = "$env:LOCALAPPDATA\Reaction\external")

$ErrorActionPreference = "Stop"
$zips = @{
    "EQ3_6v80a Windows PC.zip"      = "eq36_win"
    "EQ3_6v8.0a TestLibrary PC.zip" = "eq36_tests"
}
foreach ($name in $zips.Keys) {
    $dir = Join-Path $Dest $zips[$name]
    New-Item -ItemType Directory -Force -Path $dir | Out-Null
    $zip = Join-Path $dir "download.zip"
    $url = "https://raw.githubusercontent.com/LLNL/EQ3_6/main/" + [uri]::EscapeDataString($name)
    curl.exe -sSL -o $zip $url
    Expand-Archive -Path $zip -DestinationPath $dir -Force
    Remove-Item $zip
}

$obigt = Join-Path $Dest "obigt"
New-Item -ItemType Directory -Force -Path $obigt | Out-Null
foreach ($f in "AD","Berman_cr","DEW","H2O_aq","SLOP98-a","SLOP98-b","SUPCRT92",
               "inorganic_aq","inorganic_cr","inorganic_gas","inorganic_liq","organic_aq","organic_gas") {
    curl.exe -sSL -o (Join-Path $obigt "$f.csv") `
        "https://raw.githubusercontent.com/jedick/CHNOSZ/main/inst/extdata/OBIGT/$f.csv"
}
# Parameters of the Berman (1988) minerals listed in Berman_cr.csv
$berman = Join-Path $obigt "Berman"
New-Item -ItemType Directory -Force -Path $berman | Out-Null
foreach ($f in "BDat17_2017","Ber88_1988","Ber90_1990","DS10_2010","FDM+14_2014","Got04_2004",
               "JUN92_1992","SHD91_1991","VGT92_1992","VPT01_2001","VPV05_2005","ZS92_1992") {
    curl.exe -sSL -o (Join-Path $berman "$f.csv") `
        ("https://raw.githubusercontent.com/jedick/CHNOSZ/main/inst/extdata/Berman/" + [uri]::EscapeDataString("$f.csv"))
}

Write-Host "Done. To run the validation tests:"
Write-Host "  `$env:REACTION_EQ36_DIR   = '$Dest\eq36_win\Eq3_6v8.0a'"
Write-Host "  `$env:REACTION_EQ36_TESTS = '$Dest\eq36_tests'"
Write-Host "  `$env:REACTION_OBIGT_DIR  = '$obigt'"
Write-Host "  python -m pytest"

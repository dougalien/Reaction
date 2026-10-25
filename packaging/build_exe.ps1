# Build dist\Reaction.exe on Windows.  Run from the repository root:  powershell -File packaging\build_exe.ps1
$ErrorActionPreference = "Stop"
python -m pip install --upgrade pip
python -m pip install numpy scipy pyinstaller
python -m PyInstaller --noconfirm --clean packaging\reaction.spec
Write-Host "Built dist\Reaction.exe"

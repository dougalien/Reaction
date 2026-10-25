# PyInstaller recipe for Reaction.exe: one file holding Python, numpy, scipy, the app and its own data
# tables.  The thermodynamic databases are not inside it; they stay in %LOCALAPPDATA%\Reaction\external.
# Build with packaging/build_exe.ps1 (or: pyinstaller packaging/reaction.spec from the repository root).
from pathlib import Path

ROOT = Path(SPECPATH).parent
SRC = ROOT / "src" / "reaction"

a = Analysis(
    [str(ROOT / "packaging" / "reaction_exe.py")],
    pathex=[str(ROOT / "src")],
    datas=[(str(SRC / "data" / "*.tsv"), "reaction/data"),
           (str(SRC / "app" / "static"), "reaction/app/static")],
    excludes=["matplotlib", "pandas", "tkinter", "pytest", "IPython"],
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz, a.scripts, a.binaries, a.datas,
    name="Reaction",
    console=True,          # the window shows where the app runs; closing it stops the app
    upx=False,
)

# PyInstaller recipe for Reaction.exe: one file holding Python, numpy, scipy, the app and its own data
# tables, including the OBIGT and EQ3/6 reference data in reaction/data/reference (see its NOTICE.md).
# Build with packaging/build_exe.ps1 (or: pyinstaller packaging/reaction.spec from the repository root).
from pathlib import Path

ROOT = Path(SPECPATH).parent
SRC = ROOT / "src" / "reaction"

a = Analysis(
    [str(ROOT / "packaging" / "reaction_exe.py")],
    pathex=[str(ROOT / "src")],
    datas=[(str(SRC / "data" / "*.tsv"), "reaction/data"),
           (str(SRC / "data" / "*.dat"), "reaction/data"),
           (str(SRC / "data" / "reference"), "reaction/data/reference"),
           (str(SRC / "app" / "static"), "reaction/app/static")],
    excludes=["matplotlib", "pandas", "tkinter", "pytest", "IPython"],
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz, a.scripts, a.binaries, a.datas,
    name="Reaction",
    console=True,          # the window shows where the app runs; closing it stops the app
    upx=False,
    icon=str(ROOT / "packaging" / "reaction.ico"),           # the Reaction R in the Dougalien Systems style (packaging/icon/make_icons.py)
    version=str(ROOT / "packaging" / "version_info.txt"),  # company and product in Properties > Details
)

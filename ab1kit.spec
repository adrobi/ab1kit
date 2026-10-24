# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec for the ab1kit one-file Windows executable.

Build:  pyinstaller ab1kit.spec --noconfirm
"""

a = Analysis(
    ["chrom_cli.py"],
    pathex=[],
    binaries=[],
    datas=[],
    hiddenimports=["chrom_core", "ui", "i18n", "translations_en"],
    hookspath=[],
    excludes=["tkinter", "matplotlib", "numpy.tests"],
    noarchive=False,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name="ab1kit",
    debug=False,
    strip=False,
    upx=False,
    console=False,
)

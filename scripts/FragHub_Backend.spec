# -*- mode: python ; coding: utf-8 -*-
import os
import sys

# Icône adaptée à la plateforme : .icns sur macOS, .ico sous Windows et Linux.
# Les icônes sont dans GUI/app/assets (chemin relatif au dossier de ce .spec).
IS_MAC = sys.platform == "darwin"
ICON = os.path.join(SPECPATH, '..', 'GUI', 'app', 'assets', 'FragHub_icon.' + ('icns' if IS_MAC else 'ico'))


a = Analysis(
    ['FragHub.py'],
    pathex=['.'],
    binaries=[],
    datas=[('../datas', 'datas'), ('GUI/assets', 'GUI/assets')],
    hiddenimports=['uvicorn.logging', 'uvicorn.loops', 'uvicorn.loops.auto', 'uvicorn.protocols.http.auto', 'uvicorn.protocols.websockets.auto', 'mzspeclib_converter', 'rdkit_worker'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='FragHub_Backend',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=[ICON],
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='FragHub_Backend',
)
# Le bundle .app n'a de sens que sur macOS
if IS_MAC:
    app = BUNDLE(
        coll,
        name='FragHub_Backend.app',
        icon=ICON,
        bundle_identifier=None,
    )

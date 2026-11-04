"""Reconstruit tout FragHub : back-end (Rust + PyInstaller), front-end (Nuxt), puis l'installeur Electron.

Usage :  python build_all.py                (tout)
         python build_all.py --skip-rust    (ne recompile pas le Rust)
         python build_all.py --skip-backend (ni Rust ni PyInstaller : front + installeur)
         python build_all.py --skip-frontend
         python build_all.py --no-installer (s'arrête avant electron-builder)
         python build_all.py --dry-run      (affiche les commandes sans rien lancer)

Résultat : GUI/dist_electron/ (.dmg sur macOS, .exe sur Windows)
"""
import sys

from build_backend import build_pyinstaller, build_rust
from build_frontend import build_frontend, npm
from build_utils import (BuildError, DRY_RUN, GUI_DIR, IS_WINDOWS, SCRIPTS_DIR,
                         main_wrapper, run, step)


def build_installer():
    step("Installeur Electron (electron-builder)")
    backend = SCRIPTS_DIR / "dist" / "FragHub_Backend" / ("FragHub_Backend.exe" if IS_WINDOWS else "FragHub_Backend")
    if not backend.exists() and not DRY_RUN:
        raise BuildError(f"Back-end compilé introuvable ({backend}). Lancez d'abord python build_backend.py.")
    if not (GUI_DIR / ".output" / "public" / "index.html").exists() and not DRY_RUN:
        raise BuildError("Front-end compilé introuvable (GUI/.output/public). Lancez d'abord python build_frontend.py.")
    run([npm(), "run", "build:electron"], cwd=GUI_DIR, title="Création de l'installeur dans GUI/dist_electron")


def build_all():
    """Reconstruit tout FragHub : back-end, front-end, installeur Electron.

Usage :  python build_all.py [--skip-rust] [--skip-backend] [--skip-frontend] [--no-installer] [--dry-run]"""
    args = sys.argv
    if "--skip-backend" not in args:
        if "--skip-rust" not in args:
            build_rust()
        build_pyinstaller()
    if "--skip-frontend" not in args:
        build_frontend()
    if "--no-installer" not in args:
        build_installer()


if __name__ == "__main__":
    sys.exit(main_wrapper(build_all))

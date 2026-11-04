"""Recompile le back-end FragHub : moteur Rust (maturin) puis exécutable Python (PyInstaller).

Usage :  python build_backend.py            (Rust + PyInstaller)
         python build_backend.py --rust-only (juste le module Rust, pour tester en dev)
         python build_backend.py --dry-run   (affiche les commandes sans rien lancer)

Résultat : scripts/dist/FragHub_Backend/
"""
import sys

from build_utils import (BuildError, RUST_DIR, SCRIPTS_DIR, VENV_DIR, ensure_venv,
                         main_wrapper, require_tool, run, step, venv_bin)


def build_rust():
    step("Back-end 1/2 : moteur Rust (maturin develop --release)")
    require_tool("cargo", "Installez Rust : https://rustup.rs")
    ensure_venv()
    maturin = venv_bin("maturin")
    cmd = [maturin, "develop", "--release"] if maturin.exists() else [venv_bin("python"), "-m", "maturin", "develop", "--release"]
    # VIRTUAL_ENV : indique à maturin d'installer le module dans le .venv du projet
    run(cmd, cwd=RUST_DIR, title="Compilation et installation de fraghub_rust dans .venv",
        env_extra={"VIRTUAL_ENV": str(VENV_DIR)})


def build_pyinstaller():
    step("Back-end 2/2 : exécutable Python (PyInstaller)")
    ensure_venv()
    spec = SCRIPTS_DIR / "FragHub_Backend.spec"
    if not spec.exists():
        raise BuildError(f"Fichier introuvable : {spec}")
    run([venv_bin("python"), "-m", "PyInstaller", "-y", spec.name], cwd=SCRIPTS_DIR,
        title="Création de scripts/dist/FragHub_Backend")


def build_backend():
    """Recompile le back-end FragHub : Rust (maturin) puis PyInstaller.

Usage :  python build_backend.py [--rust-only] [--dry-run]"""
    build_rust()
    if "--rust-only" not in sys.argv:
        build_pyinstaller()


if __name__ == "__main__":
    sys.exit(main_wrapper(build_backend))

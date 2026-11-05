"""Recompile le back-end FragHub : moteur Rust (maturin) puis exécutable Python (PyInstaller).

Usage :  python build_backend.py            (Rust + PyInstaller)
         python build_backend.py --rust-only (juste le module Rust, pour tester en dev)
         python build_backend.py --dry-run   (affiche les commandes sans rien lancer)

Résultat : scripts/dist/FragHub_Backend/
"""
import sys

from build_utils import (BuildError, DRY_RUN, IS_WINDOWS, RDKIT_PREFIX, RUST_DIR, SCRIPTS_DIR, VENV_DIR, ensure_venv,
                         main_wrapper, require_tool, run, step, venv_bin)

# Commande d'installation de RDKit C++, rappelée si on ne le trouve pas
_RDKIT_INSTALL = ("mamba create -p \"{0}\" -c conda-forge librdkit-dev=2026.09.1 libboost-headers=1.90.0"
                  .format(RDKIT_PREFIX))


def rdkit_lib_dirs():
    """Dossiers de RDKit C++ : (en-têtes, DLL ou .so). Erreur claire s'il n'est pas installé."""
    base = RDKIT_PREFIX / "Library" if IS_WINDOWS else RDKIT_PREFIX
    include = base / "include" / "rdkit"
    libs = base / ("bin" if IS_WINDOWS else "lib")
    if not include.is_dir():
        raise BuildError(f"RDKit C++ introuvable dans {RDKIT_PREFIX}.\n"
                         f"Installez-le avec :\n  {_RDKIT_INSTALL}\n"
                         "ou indiquez son emplacement avec la variable d'environnement RDKIT_PREFIX.")
    return include, libs


def build_rust():
    step("Back-end 1/2 : moteur Rust (maturin develop --release)")
    require_tool("cargo", "Installez Rust : https://rustup.rs")
    rdkit_lib_dirs()
    ensure_venv()
    maturin = venv_bin("maturin")
    cmd = [maturin, "develop", "--release"] if maturin.exists() else [venv_bin("python"), "-m", "maturin", "develop", "--release"]
    # VIRTUAL_ENV : indique à maturin d'installer le module dans le .venv du projet
    run(cmd, cwd=RUST_DIR, title="Compilation et installation de fraghub_rust dans .venv",
        env_extra={"VIRTUAL_ENV": str(VENV_DIR), "RDKIT_PREFIX": str(RDKIT_PREFIX)})


def build_pyinstaller():
    step("Back-end 2/2 : exécutable Python (PyInstaller)")
    ensure_venv()
    spec = SCRIPTS_DIR / "FragHub_Backend.spec"
    if not spec.exists():
        raise BuildError(f"Fichier introuvable : {spec}")
    _, rdkit_libs = rdkit_lib_dirs()
    # Le .spec ajoute les DLL de RDKit C++ au PATH pour que PyInstaller les copie dans l'exécutable
    run([venv_bin("python"), "-m", "PyInstaller", "-y", spec.name], cwd=SCRIPTS_DIR,
        title="Création de scripts/dist/FragHub_Backend",
        env_extra={"RDKIT_PREFIX": str(RDKIT_PREFIX)})
    # Vérification : sans les DLL de RDKit, l'exécutable ne démarrerait pas sur une machine sans mamba
    internal = SCRIPTS_DIR / "dist" / "FragHub_Backend" / "_internal"
    if IS_WINDOWS and not DRY_RUN and not any(internal.glob("RDKitGraphMol*.dll")):
        raise BuildError(f"Les DLL de RDKit C++ n'ont pas été copiées dans {internal} (cherchées dans {rdkit_libs}).")


def build_backend():
    """Recompile le back-end FragHub : Rust (maturin) puis PyInstaller.

Usage :  python build_backend.py [--rust-only] [--dry-run]"""
    build_rust()
    if "--rust-only" not in sys.argv:
        build_pyinstaller()


if __name__ == "__main__":
    sys.exit(main_wrapper(build_backend))

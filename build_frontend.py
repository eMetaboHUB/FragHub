"""Recompile l'interface FragHub (Nuxt) dans GUI/.output/public.

Usage :  python build_frontend.py           (npm install si besoin, puis npm run generate)
         python build_frontend.py --dry-run (affiche les commandes sans rien lancer)
"""
import sys

from build_utils import GUI_DIR, main_wrapper, require_tool, run, step


def npm():
    return require_tool("npm", "Installez Node.js : https://nodejs.org")


def install_deps_if_needed():
    if not (GUI_DIR / "node_modules").exists():
        run([npm(), "install"], cwd=GUI_DIR, title="Installation des dépendances (npm install)")


def build_frontend():
    """Recompile l'interface FragHub (Nuxt) dans GUI/.output/public.

Usage :  python build_frontend.py [--dry-run]"""
    step("Front-end : interface Nuxt (npm run generate)")
    install_deps_if_needed()
    run([npm(), "run", "generate"], cwd=GUI_DIR, title="Génération de GUI/.output/public")


if __name__ == "__main__":
    sys.exit(main_wrapper(build_frontend))

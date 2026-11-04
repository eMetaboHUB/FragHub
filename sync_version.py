"""Synchronise le numéro de version de FragHub dans tout le projet.

SOURCE UNIQUE : pyproject.toml (racine), ligne `version = "..."` de la section [project].

Usage :
    python sync_version.py            recopie la version de pyproject.toml partout
    python sync_version.py 2.1.0      change la version dans pyproject.toml, puis la recopie partout
    python sync_version.py --check    vérifie seulement (code de sortie 1 si un fichier n'est pas à jour)

Les scripts build_backend.py, build_frontend.py et build_all.py l'appellent automatiquement.

Volontairement NON modifiés : CHANGELOG.md (historique), la phrase "FragHub 2.0 introduces…" du README,
et les liens de téléchargement Zenodo du README (ils pointent vers des fichiers déjà publiés : un nouveau
numéro sans nouveau dépôt Zenodo donnerait des liens morts).
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "pyproject.toml"
VERSION_RE = r"\d+\.\d+\.\d+(?:[-.+][0-9A-Za-z.]+)?"

# (fichier, motif, nombre d'occurrences attendu)
# Chaque motif capture : (avant)(version)(après). Seul le groupe "version" est remplacé.
TARGETS = [
    ("GUI/package.json", r'(\A\{\s*"name":\s*"FragHub",\s*"version":\s*")(' + VERSION_RE + r')(")', 1),
    ("GUI/package-lock.json", r'(\A\{\s*"name":\s*"FragHub",\s*"version":\s*")(' + VERSION_RE + r')(")', 1),
    ("GUI/package-lock.json", r'("":\s*\{\s*"name":\s*"FragHub",\s*"version":\s*")(' + VERSION_RE + r')(")', 1),
    ("scripts/fraghub_rust/Cargo.toml", r'(\[package\]\s*name\s*=\s*"fraghub_rust"\s*version\s*=\s*")(' + VERSION_RE + r')(")', 1),
    ("scripts/fraghub_rust/Cargo.lock", r'(name = "fraghub_rust"\r?\nversion = ")(' + VERSION_RE + r')(")', 1),
    ("scripts/fraghub_rust/pyproject.toml", r'(\[project\]\s*name\s*=\s*"fraghub-rust"\s*version\s*=\s*")(' + VERSION_RE + r')(")', 1),
    ("uv.lock", r'(name = "fraghub"\r?\nversion = ")(' + VERSION_RE + r')(")', 1),
    ("uv.lock", r'(name = "fraghub-rust"\r?\nversion = ")(' + VERSION_RE + r')(")', 1),
    ("README.md", r'(alt="FragHub Icon" /> )(' + VERSION_RE + r')()', 1),
    ("README.md", r'(^# FragHub\s+\()(' + VERSION_RE + r')(\))', 1),
    ("documentation/cmd_pyinstaller.txt", r'(FragHub_(?:Windows|macOS|Linux)_)(' + VERSION_RE + r')(_)', 5),
]


def _read(path):
    # newline="" : conserve les fins de ligne d'origine (LF ou CRLF), pour ne pas réécrire tout le fichier sous Windows
    with open(path, encoding="utf-8", newline="") as f:
        return f.read()


def _write(path, text):
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write(text)


def read_source_version():
    text = _read(SOURCE)
    m = re.search(r'^\[project\][^\[]*?^version\s*=\s*"(' + VERSION_RE + r')"', text, re.M | re.S)
    if not m:
        sys.exit(f"✘ Impossible de lire la version dans {SOURCE} (section [project], ligne version = \"x.y.z\").")
    return m.group(1)


def write_source_version(new_version):
    text = _read(SOURCE)
    new_text, n = re.subn(r'(^\[project\][^\[]*?^version\s*=\s*")(' + VERSION_RE + r')(")',
                          lambda m: m.group(1) + new_version + m.group(3), text, count=1, flags=re.M | re.S)
    if n != 1:
        sys.exit(f"✘ Impossible d'écrire la version dans {SOURCE}.")
    _write(SOURCE, new_text)


def sync(version, check_only=False, quiet=False):
    """Recopie `version` dans tous les fichiers cibles. Retourne la liste des fichiers pas à jour / modifiés."""
    changed, problems = [], []
    contents = {}
    for rel, pattern, expected in TARGETS:
        path = ROOT / rel
        if not path.exists():
            problems.append(f"{rel} : fichier introuvable")
            continue
        text = contents.get(rel) or _read(path)
        regex = re.compile(pattern, re.M)
        found = regex.findall(text)
        if len(found) != expected:
            problems.append(f"{rel} : {len(found)} occurrence(s) trouvée(s) au lieu de {expected} (motif à revoir)")
            contents[rel] = text
            continue
        if any(v != version for _, v, _ in found):
            text = regex.sub(lambda m: m.group(1) + version + m.group(3), text)
            if rel not in changed:
                changed.append(rel)
        contents[rel] = text

    if not check_only:
        for rel in changed:
            _write(ROOT / rel, contents[rel])

    if not quiet or changed or problems:
        verb = "pas à jour" if check_only else "mis à jour"
        if changed:
            print(f"Version {version} — fichiers {verb} :")
            for rel in changed:
                print(f"  • {rel}")
        elif not quiet:
            print(f"Version {version} — tous les fichiers sont déjà à jour.")
        for p in problems:
            print(f"  ⚠ {p}")
    return changed, problems


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    check_only = "--check" in sys.argv
    if "-h" in sys.argv or "--help" in sys.argv:
        print(__doc__.strip())
        return 0
    if args:
        new_version = args[0]
        if not re.fullmatch(VERSION_RE, new_version):
            sys.exit(f"✘ Version invalide : {new_version!r} (format attendu : 2.1.0)")
        if check_only:
            sys.exit("✘ --check ne s'utilise pas avec un numéro de version.")
        write_source_version(new_version)
        print(f"pyproject.toml → version {new_version}")
    version = read_source_version()
    changed, problems = sync(version, check_only=check_only)
    if problems:
        return 1
    if check_only and changed:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

# Documentation Développeur : FragHub

Ce guide détaille les étapes nécessaires pour installer l'environnement de développement, compiler le backend Python, et générer l'installateur de l'application Electron sous macOS.

---

## 1. Prérequis et Configuration Python

Avant de commencer, assurez-vous que les outils suivants sont installés sur votre machine :
* **Node.js** : [Télécharger Node.js](https://nodejs.org/) (Version LTS recommandée).
* **Python 3.12** : [Télécharger Python 3.12](https://www.python.org/downloads/release/python-3128/).

### Initialisation de l'environnement Python
Pour isoler les dépendances du projet et les gérer de manière optimale, nous utilisons `uv`. Ouvrez votre terminal dans le dossier contenant le code Python et exécutez les commandes suivantes :

```bash
# 1. Installation de uv (si ce n'est pas déjà fait)
curl -LsSf https://astral.sh/uv/install.sh | sh

# 2. Synchronisation de l'environnement et installation des dépendances
# Cela va créer automatiquement le dossier .venv et récupérer les dépendances via pyproject.toml / uv.lock
uv sync

# 3. Activation de l'environnement virtuel
source .venv/bin/activate
```

---

## 2. Installation du Projet

Commencez par récupérer le code source, puis installez les dépendances nécessaires au Front-End.

```bash
# 1. Cloner le dépôt localement
git clone <URL_DU_DEPOT_GIT>
cd FragHub

# 2. Se rendre dans le dossier du Front-End (Interface Graphique)
cd GUI

# 3. Installer les paquets Node.js
npm install
```

> ⚠️ **NOTE IMPORTANTE POUR MACOS : BUG D'INSTALLATION ELECTRON**
> Il a été constaté que la commande classique `npm install` bug actuellement sous macOS pour l'installation d'Electron. Pour contourner ce problème, vous devez installer Electron manuellement dans vos `node_modules` directement depuis le dépôt Git officiel.
> Vous pouvez utiliser la commande suivante (toujours dans le dossier `GUI`) :
> ```bash
> npm install git+[https://github.com/electron/electron.git](https://github.com/electron/electron.git)
> ```

---

## 3. Compilation du Backend (Python)

Le cœur logique de l'application est développé en Python. Pour qu'Electron puisse l'exécuter de manière autonome, il doit être empaqueté avec PyInstaller.

1. Placez-vous dans le dossier `scripts`.
2. Exécutez la commande suivante (notez l'utilisation des deux-points `:` comme séparateur pour macOS) :

```bash
uv run pyinstaller --noconfirm --onedir --noconsole --paths . --icon="GUI/assets/FragHub_icon.icns" --name="FragHub_Backend" --add-data="../datas:datas" --add-data="GUI/assets:GUI/assets" --hidden-import=uvicorn.logging --hidden-import=uvicorn.loops --hidden-import=uvicorn.loops.auto --hidden-import=uvicorn.protocols.http.auto --hidden-import=uvicorn.protocols.websockets.auto --hidden-import=mzspeclib_converter --hidden-import=rdkit_worker FragHub.py
```

> **📌 Notes importantes concernant le Backend :**
> * **Versionnage :** Gardez `--name=FragHub_Backend` tel quel : Electron cherche un exécutable portant exactement ce nom. Pour changer le numéro de version de FragHub, voir la section 6.
> * **Icône :** Sous macOS, le format standard pour les icônes est `.icns`.
> * **Fichiers générés :** PyInstaller va créer le backend compilé sous forme de dossier ou de bundle `.app` (incluant l'exécutable et le dossier `_internal`) dans `scripts/dist/FragHub_Backend`.
> * ⚠️ **Règle d'or :** Cette commande PyInstaller doit être obligatoirement relancée **à chaque fois que vous modifiez un fichier du code source Python ou des données externes tel que les bases de données**.

---

## 4. Mode Développement (Dev)

Pour travailler sur l'interface graphique en temps réel (avec rechargement à chaud) et tester l'application :

1. Placez-vous dans le dossier `GUI`.
2. Dans un premier terminal, lancez le serveur de développement Nuxt :
```bash
npm run dev
```
3. Dans un second terminal (toujours dans `GUI`), lancez l'application Electron en mode dev :
```bash
npm run electron:dev
```
*(Vous pouvez également utiliser `npm run electron` pour lancer l'enveloppe sans les outils de développement étendus).*

---

## 5. Build et Création de l'Installeur (Production)

Pour générer la version distribuable de l'application (le fichier `.dmg` ou `.app` d'installation pour l'utilisateur final) :

1. Placez-vous dans le dossier `GUI`.
2. Générez les fichiers statiques du Front-End :
```bash
npm run generate
```
> * **Fichiers générés :** Le code compilé du Front-End sera placé dans le dossier caché `.output/public/`.
> * ⚠️ **Règle d'or :** La commande `npm run generate` doit être refaite **à chaque fois que vous modifiez le code du Front-End** (fichiers `.vue`, `main.js`, etc.) avant de lancer la création de l'installeur.

3. Lancez la construction de l'installeur Electron final pour macOS :
```bash
npm run build:electron -- --mac
```
> * **Fichiers générés :** Une fois le processus terminé, l'installeur macOS (fichier `.dmg`) prêt à être distribué se trouvera dans le dossier `dist_electron/`, situé à la racine de votre dossier `GUI`.

---

## 6. Changer le numéro de version

Le numéro de version de FragHub est défini à **un seul endroit** : la ligne `version` de la section `[project]` du fichier `pyproject.toml`, à la racine du projet.

```toml
[project]
name = "fraghub"
version = "2.0.0"
```

**Ne modifiez pas la version à la main ailleurs** : le script `sync_version.py` la recopie dans tous les autres fichiers.

### Procédure
Depuis la racine du projet :
```bash
python sync_version.py 2.1.0
```
Cette commande écrit `2.1.0` dans `pyproject.toml`, puis la recopie partout. Si vous préférez modifier `pyproject.toml` vous-même, lancez ensuite `python sync_version.py` (sans numéro).

Pour vérifier que tous les fichiers sont cohérents, sans rien modifier :
```bash
python sync_version.py --check
```

> * **Automatique :** `build_all.py`, `build_backend.py` et `build_frontend.py` lancent cette synchronisation au démarrage. Un build utilise donc toujours la version de `pyproject.toml`.
> * **Fichiers mis à jour :** `GUI/package.json` (qui donne aussi son nom à l'installeur, ex. `FragHub-2.1.0-arm64.dmg`), `GUI/package-lock.json`, `scripts/fraghub_rust/Cargo.toml`, `scripts/fraghub_rust/Cargo.lock`, `scripts/fraghub_rust/pyproject.toml`, `uv.lock`, le titre du `README.md` et `documentation/cmd_pyinstaller.txt`.
> * **Volontairement non modifiés :** `CHANGELOG.md` (ajoutez vous-même la nouvelle entrée), la phrase « FragHub 2.0 introduces… » du README, et les liens de téléchargement Zenodo du README (à mettre à jour à la main une fois la nouvelle version publiée sur Zenodo, car l'identifiant du dépôt change).
> * ⚠️ **Si le script affiche un avertissement ⚠**, un fichier cible a changé de structure et la version n'y a pas été trouvée : vérifiez le fichier et, si besoin, ajustez son motif dans la liste `TARGETS` en haut de `sync_version.py`.

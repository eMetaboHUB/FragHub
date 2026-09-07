# Developer Documentation: FragHub

This guide details the steps required to set up the development environment, compile the Python backend, and generate the Electron application installer on Linux.

---

## 1. Prerequisites and Python Configuration

Before you begin, make sure the following tools are installed on your machine:
* **Node.js**: [Download Node.js](https://nodejs.org/) (LTS version recommended).
* **Python 3.12**: [Download Python 3.12](https://www.python.org/downloads/release/python-3128/).

### Python Environment Initialization
To isolate project dependencies and manage them efficiently, we use `uv`. Open your terminal in the folder containing the Python code and run the following commands:

```bash
# 1. Install uv (if not already installed)
curl -LsSf https://astral.sh/uv/install.sh | sh

# 2. Sync the environment and install dependencies
# This automatically creates the .venv and resolves the dependencies from pyproject.toml / uv.lock
uv sync

# 3. Activate the virtual environment
source .venv/bin/activate
```

---

## 2. Project Installation

Start by retrieving the source code, then install the Front-End dependencies.

```bash
# 1. Clone the repository locally
git clone [https://github.com/eMetaboHUB/FragHub.git](https://github.com/eMetaboHUB/FragHub.git)
cd FragHub

# 2. Navigate to the Front-End folder (Graphical Interface)
cd GUI

# 3. Install Node.js packages
npm install
```

---

## 3. Backend Compilation (Python)

The core logic of the application is developed in Python. For Electron to run it autonomously, it must be packaged with PyInstaller.

1. Navigate to the `scripts` folder.
2. Run the following command (note the use of `:` as separator for Linux):

```bash
uv run pyinstaller --noconfirm --onedir --noconsole --paths . --icon="GUI/assets/FragHub_icon.ico" --name="FragHub_Backend" --add-data="../datas:datas" --add-data="GUI/assets:GUI/assets" --hidden-import=uvicorn.logging --hidden-import=uvicorn.loops --hidden-import=uvicorn.loops.auto --hidden-import=uvicorn.protocols.http.auto --hidden-import=uvicorn.protocols.websockets.auto --hidden-import=mzspeclib_converter --hidden-import=rdkit_worker FragHub.py
```

> **📌 Important notes regarding the Backend:**
> * **Versioning:** Remember to adapt the `--name=FragHub_Backend` argument to match the current version number if the application evolves.
> * **Generated files:** PyInstaller will create the compiled backend as a folder (including the executable and the `_internal` folder) in `scripts/dist/FragHub_Backend`.
> * ⚠️ **Golden rule:** This PyInstaller command **must be re-run every time you modify a Python source file or external data such as databases**.

---

## 4. Development Mode (Dev)

To work on the graphical interface in real time (with hot reloading) and test the application:

1. Navigate to the `GUI` folder.
2. In a first terminal, start the Nuxt development server:
```bash
npm run dev
```
3. In a second terminal (still in `GUI`), launch the Electron application in dev mode:
```bash
npm run electron:dev
```
*(You can also use `npm run electron` to launch the wrapper without the extended development tools).*

---

## 5. Build and Installer Creation (Production)

To generate the distributable version of the application (the `.dmg` or `.app` installation file for the end user):

1. Navigate to the `GUI` folder.
2. Generate the static Front-End files:
```bash
npm run generate
```
> * **Generated files:** The compiled Front-End code will be placed in the hidden `.output/public/` folder.
> * ⚠️ **Golden rule:** The `npm run generate` command **must be re-run every time you modify Front-End code** (`.vue` files, `main.js`, etc.) before launching the installer creation.

3. Launch the final Electron installer build for Linux:
```bash
npm run build:electron -- --linux
```
> * **Generated files:** Once the process is complete, the Linux installer (`.AppImage` or `.deb` file) ready for distribution will be located in the `dist_electron/` folder, at the root of your `GUI` folder.
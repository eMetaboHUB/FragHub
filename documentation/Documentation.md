# Developer Documentation: FragHub

This guide details the steps required to set up the development environment, compile the Python backend, and generate the Electron application installer.

---

## 1. Prerequisites and Python Configuration

Before you begin, make sure the following tools are installed on your machine:
*   **Node.js**: [Download Node.js](https://nodejs.org/) (LTS version recommended).
*   **Python 3.12**: [Download Python 3.12](https://www.python.org/downloads/release/python-3128/).

### Python Environment Initialization
To isolate project dependencies and manage them efficiently, we use `uv`. Open your terminal in the folder containing the Python code and run the following commands:

```bash
# 1. Install uv (if not already installed)
curl -LsSf https://astral.sh/uv/install.sh | sh
# On Windows, you can use: powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"

# 2. Sync the environment and install dependencies
# This automatically creates the .venv and resolves the dependencies from pyproject.toml / uv.lock
uv sync

# 3. Activate the virtual environment
# On Windows:
.venv\Scripts\activate
# On macOS / Linux:
source .venv/bin/activate
```

---

## 2. Project Installation

Start by retrieving the source code, then install the Front-End dependencies.

```bash
# 1. Clone the repository locally
git clone https://github.com/eMetaboHUB/FragHub.git
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
2. Run the following command:

```bash
uv run pyinstaller --noconfirm --onedir --noconsole --paths . --icon="GUI/assets/FragHub_icon.ico" --name="FragHub_Backend" --add-data="../datas;datas" --add-data="GUI/assets;GUI/assets" --hidden-import=uvicorn.logging --hidden-import=uvicorn.loops --hidden-import=uvicorn.loops.auto --hidden-import=uvicorn.protocols.http.auto --hidden-import=uvicorn.protocols.websockets.auto --hidden-import=mzspeclib_converter --hidden-import=rdkit_worker FragHub.py
```

> **📌 Important notes regarding the Backend:**
> *   `--paths .` and the two `--hidden-import` flags for `mzspeclib_converter` and `rdkit_worker` are required: these two files live next to `FragHub.py` in `scripts/` and PyInstaller's static analysis does not reliably auto-detect sibling local modules without an explicit search path. Omitting them produces a backend that crashes at startup with `ModuleNotFoundError`.
> *   **Generated files:** PyInstaller will create the compiled backend as a folder (including the executable and the `_internal` folder) in `scripts/dist/FragHub_Backend`.
> *   ⚠️ **Golden rule:** This PyInstaller command **must be re-run every time you modify a Python source file or external data such as databases**.

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

To generate the distributable version of the application (the installation `.exe` for the end user):

1. Navigate to the `GUI` folder.
2. Generate the static Front-End files:
```bash
npm run generate
```
> *   **Generated files:** The compiled Front-End code will be placed in the hidden `.output/public/` folder.
> *   ⚠️ **Golden rule:** The `npm run generate` command **must be re-run every time you modify Front-End code** (`.vue` files, `main.js`, etc.) before launching the installer creation.

3. Launch the final Electron installer build:
```bash
npm run build:electron -- --win
```
> *   **Generated files:** Once the process is complete, the Windows installer ready for distribution will be located in the `dist_electron/` folder, at the root of your `GUI` folder.

---

## 6. Changing the Version Number

The FragHub version number is defined in **one single place**: the `version` line of the `[project]` section in the `pyproject.toml` file at the root of the project.

```toml
[project]
name = "fraghub"
version = "2.0.0"
```

**Do not edit the version by hand anywhere else**: the `sync_version.py` script copies it into all the other files.

### Procedure
From the root of the project:
```bash
python sync_version.py 2.1.0
```
This command writes `2.1.0` into `pyproject.toml`, then copies it everywhere else. If you prefer to edit `pyproject.toml` yourself, run `python sync_version.py` afterwards (without a number).

To check that every file is consistent, without modifying anything:
```bash
python sync_version.py --check
```

> * **Automatic:** `build_all.py`, `build_backend.py` and `build_frontend.py` run this synchronization at startup, so a build always uses the version from `pyproject.toml`.
> * **Updated files:** `GUI/package.json` (which also names the installer, e.g. `FragHub-2.1.0-arm64.dmg`), `GUI/package-lock.json`, `scripts/fraghub_rust/Cargo.toml`, `scripts/fraghub_rust/Cargo.lock`, `scripts/fraghub_rust/pyproject.toml`, `uv.lock`, the title of `README.md` and `documentation/cmd_pyinstaller.txt`.
> * **Intentionally left unchanged:** `CHANGELOG.md` (add the new entry yourself), the "FragHub 2.0 introduces…" sentence in the README, and the Zenodo download links in the README (update them by hand once the new release is published on Zenodo, since the record ID changes).
> * ⚠️ **If the script shows a ⚠ warning**, a target file has changed its structure and the version could not be found in it: check the file and, if needed, adjust its pattern in the `TARGETS` list at the top of `sync_version.py`.

---

## 7. RDKit C++ (native calculation engine)

The Rust engine calls RDKit in C++ (`scripts/fraghub_rust/build.rs`, `src/rdkit_native/`). RDKit C++ is not installed by `uv`: it comes from conda-forge.

### Installation (once)
With `mamba` (or `micromamba`, same syntax):
```powershell
mamba create -p C:\Users\<name>\rdkit_cpp -c conda-forge librdkit-dev=2026.09.1 libboost-headers=1.90.0
```
* **Location:** `build.rs` and the build scripts look for RDKit in the `RDKIT_PREFIX` environment variable, otherwise in `%USERPROFILE%\rdkit_cpp` (headers in `Library\include`, DLLs in `Library\bin`).
* **In development**, `FragHub.py` adds `%USERPROFILE%\rdkit_cpp\Library\bin` to the DLL search path.
* **Compiler:** Visual Studio Build Tools (MSVC, C++20).

### Who needs it?
* **End users of the installer / executable: nothing to install.** PyInstaller copies the RDKit C++ libraries into the executable (`scripts/dist/FragHub_Backend/_internal`), and `build_backend.py` stops with an error if they are missing.
* **Anyone who compiles FragHub, or runs it from source** (`python scripts/FragHub.py --cli …`): RDKit C++ must be installed. The Rust module is linked to it and does not load without it, even with `FRAGHUB_RDKIT_BACKEND=python` (that variable only selects the calculation engine).

### Calculation engine
By default, RDKit calculations are done in C++, directly from Rust. The former Python engine (`rdkit_worker.py`) is still available and gives identical results:
```bash
FRAGHUB_RDKIT_BACKEND=python python scripts/FragHub.py --cli ...
```

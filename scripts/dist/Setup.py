from PyQt6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QPushButton, QCheckBox, QProgressBar,
    QFileDialog, QLabel, QMessageBox
)
from PyQt6.QtGui import QIcon
from PyQt6.QtCore import Qt, QSize, QThread, pyqtSignal
import os
import sys
import zipfile
import subprocess
from pathlib import Path
import ctypes  # For setting AppUserModelID (Windows Taskbar Icon)
import platform
import shutil

if platform.system() == "Windows":
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("FragHub.Installer")

# Si le fichier est exécuté comme un exécutable PyInstaller
if getattr(sys, 'frozen', False):
    BASE_DIR = sys._MEIPASS
else:
    # Si le fichier est exécuté comme un script Python
    BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__)))


class InstallerThread(QThread):
    """
    A QThread worker that extracts a ZIP archive to a specified directory.
    It uses Python's zipfile module on Windows and the 'unzip' command
    on other operating systems (macOS, Linux).
    """
    progress_changed = pyqtSignal(int) # Signal emitted during progress (currently only 100% at end)
    installation_complete = pyqtSignal() # Signal emitted on successful completion
    error_occurred = pyqtSignal(str)   # Signal emitted if an error occurs

    def __init__(self, zip_path: Path, install_dir: Path):
        """
        Initializes the thread.

        Args:
            zip_path: Path object pointing to the source ZIP file.
            install_dir: Path object for the destination directory where files will be extracted.
        """
        super().__init__()
        self.zip_path = zip_path
        self.install_dir = install_dir

    def run(self):
        """
        Executes the extraction logic when the thread starts.
        Determines the OS and calls the appropriate extraction method.
        """
        try:
            # Ensure the installation directory exists, create it if necessary.
            # parents=True creates any necessary parent directories.
            # exist_ok=True prevents an error if the directory already exists.
            self.install_dir.mkdir(parents=True, exist_ok=True)

            if platform.system() == "Windows":
                # --- Logic for Windows using the zipfile module ---
                try:
                    with zipfile.ZipFile(self.zip_path, 'r') as zip_ref:
                        # extractall extracts all members from the archive into install_dir
                        zip_ref.extractall(self.install_dir)
                    # If extraction completes without exceptions:
                    self.progress_changed.emit(100) # Emit 100% upon completion
                    self.installation_complete.emit()
                except zipfile.BadZipFile:
                    self.error_occurred.emit(f"Error: The file '{self.zip_path.name}' is not a valid ZIP archive.")
                except Exception as e:
                    # Catch other potential zipfile errors (permissions, disk full, etc.)
                    self.error_occurred.emit(f"Error during extraction (zipfile): {e}")

            else:
                # --- Logic for macOS/Linux using the unzip command ---
                try:
                    command = ['unzip', '-o', str(self.zip_path), '-d', str(self.install_dir)]
                    # Use subprocess.run for simpler handling of completed processes
                    # capture_output=True gets stdout and stderr
                    # text=True decodes stdout/stderr as text (requires Python 3.7+)
                    # check=False prevents raising an exception on failure; we check returncode manually
                    process = subprocess.run(command, capture_output=True, text=True, check=False)

                    if process.returncode == 0:
                        self.progress_changed.emit(100)
                        self.installation_complete.emit()
                    else:
                        # Use stderr if it contains anything, otherwise a generic message
                        error_message = process.stderr or f"The unzip command failed with return code {process.returncode}"
                        self.error_occurred.emit(f"Error during extraction (unzip): {error_message.strip()}")

                except FileNotFoundError:
                    # Occurs if the 'unzip' command is not found in the system's PATH
                    self.error_occurred.emit("Error: The 'unzip' command was not found. Please ensure it is installed and in the system PATH.")
                except Exception as e:
                    # Catch other potential subprocess-related errors
                    self.error_occurred.emit(f"Error executing the unzip command: {e}")

        except Exception as e:
            # Catch potential errors even before starting extraction
            # e.g., permission error when creating install_dir
            self.error_occurred.emit(f"An unexpected error occurred: {e}")


class InstallerApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("FragHub Installer")
        self.resize(400, 300)

        # Set the window icon
        window_icon_path = os.path.join(BASE_DIR, r"setup_gui/assets/FragHub_Python_icon.ico")
        if os.path.exists(window_icon_path):
            self.setWindowIcon(QIcon(window_icon_path))

        # Main layout
        layout = QVBoxLayout()

        # Add title label
        self.info_label = QLabel("Choose Install Directory:")
        self.info_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.info_label.setStyleSheet("font-size: 14px; font-weight: bold; margin-bottom: 10px;")
        layout.addWidget(self.info_label)

        # Button for directory selection
        self.select_dir_button = QPushButton()
        self.configure_directory_button(self.select_dir_button,
                                        os.path.join(BASE_DIR, "setup_gui/assets/directory.png"), QSize(64, 64))
        self.select_dir_button.clicked.connect(self.select_directory)
        layout.addWidget(self.select_dir_button, alignment=Qt.AlignmentFlag.AlignCenter)

        # Directory label
        self.selected_dir_label = QLabel("No directory selected")
        self.selected_dir_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.selected_dir_label)

        # Add the "Create Shortcut" option for Windows only
        if platform.system() == "Windows":
            self.create_shortcut_checkbox = QCheckBox("Create Desktop Shortcut")
            self.create_shortcut_checkbox.setChecked(True)  # Checked by default
            layout.addWidget(self.create_shortcut_checkbox)

        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setValue(0)
        self.progress_bar.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.progress_bar.setTextVisible(False)  # Hide text inside the progress bar
        layout.addWidget(self.progress_bar)

        # Label for percentage (under the progress bar)
        self.progress_label = QLabel("0%")
        self.progress_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.progress_label.setStyleSheet("font-size: 12px; color: #0078d7;")  # Add some styling for visibility
        layout.addWidget(self.progress_label)

        # Install button (changed to Exit after install)
        self.install_button = QPushButton("Install")
        layout.addWidget(self.install_button)
        self.install_button.clicked.connect(self.install)

        # Set main layout
        self.setLayout(layout)

        # Variables
        self.selected_directory = None
        self.installer_thread = None  # Thread pour l'installation

    def configure_directory_button(self, button: QPushButton, image_path: str, size: QSize):
        """Configures the 'Choose install directory' button with an icon."""
        if os.path.exists(image_path):
            button.setIcon(QIcon(image_path))
            button.setIconSize(size)
            button.setFixedSize(size.width() + 20, size.height() + 20)
            button.setStyleSheet(
                """
                QPushButton {
                    border: 2px solid #0078d7;
                    border-radius: 8px;
                    background-color: transparent;
                }
                QPushButton:hover {
                    border-color: #005fb8;
                }
                QPushButton:pressed {
                    border-color: #004c94;
                }
                """
            )

    def select_directory(self):
        """Opens a dialog box to select a directory."""
        directory = QFileDialog.getExistingDirectory(self, "Select Install Directory")
        if directory:
            self.selected_directory = directory
            self.selected_dir_label.setText(f"Selected Directory: {directory}")
        else:
            self.selected_directory = None
            self.selected_dir_label.setText("No directory selected")

    def install(self):
        """Démarre l'installation en supprimant uniquement les fichiers liés, si nécessaire."""
        if not self.selected_directory:
            self.selected_dir_label.setText("Please select a directory first!")
            return

        zip_path = None

        # Recherche du fichier ZIP contenant "FragHub"
        for file in Path(__file__).parent.glob("*.zip"):
            if "FragHub" in file.name:
                zip_path = file
                break

        # Vérifier si le fichier a été trouvé
        if not zip_path or not zip_path.exists():
            self.selected_dir_label.setText("FragHub_1.3.0.zip not found!")
            return

        install_dir = Path(self.selected_directory)

        # Vérifier l'existence d'une installation précédente
        if self.previous_installation_exists(install_dir, zip_path):
            response = QMessageBox.question(
                self,
                "Existing Installation Detected",
                "A previous installation of FragHub was detected.\n"
                "Do you want to remove it before proceeding?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            )

            if response == QMessageBox.StandardButton.Yes:
                # Supprimer les fichiers liés à l'application uniquement
                try:
                    self.remove_previous_installation(install_dir, zip_path)
                    self.selected_dir_label.setText("Previous installation removed successfully.")
                except Exception as e:
                    self.selected_dir_label.setText(f"Failed to remove previous installation: {e}")
                    return

            else:
                self.selected_dir_label.setText("Installation cancelled.")
                return

        # Lancer l'installation dans un thread séparé
        self.installer_thread = InstallerThread(zip_path, install_dir)
        self.installer_thread.progress_changed.connect(self.update_progress)
        self.installer_thread.installation_complete.connect(self.installation_complete)
        self.installer_thread.error_occurred.connect(self.installation_failed)
        self.installer_thread.start()
        self.selected_dir_label.setText("Installing...")

    def previous_installation_exists(self, install_dir: Path, zip_path: Path) -> bool:
        """Vérifie si les fichiers de l'installation précédente existent."""
        if not install_dir.exists():
            return False

        # Liste des fichiers contenus dans le ZIP
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_files = zip_ref.namelist()

        # Vérifier si au moins un fichier du ZIP existe dans le dossier d'installation
        for file in zip_files:
            if (install_dir / file).exists():
                return True

        return False

    def remove_previous_installation(self, install_dir: Path, zip_path: Path):
        """Supprime uniquement les fichiers liés à l'application."""
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_files = zip_ref.namelist()

        for file in zip_files:
            file_path = install_dir / file
            try:
                if file_path.is_file() or file_path.is_symlink():
                    file_path.unlink()  # Supprime le fichier ou le lien symbolique
                elif file_path.is_dir() and file_path.exists():
                    shutil.rmtree(file_path)  # Supprime le dossier
            except Exception as e:
                print(f"Failed to delete {file_path}: {e}")

    def update_progress(self, progress):
        """Updates the progress bar and percentage label."""
        self.progress_bar.setValue(progress)
        self.progress_label.setText(f"{progress}%")

    def installation_complete(self):
        """Handles completion of installation."""
        self.selected_dir_label.setText("Installation complete!")
        self.progress_bar.setValue(100)
        self.progress_label.setText("100%")

        # Create desktop shortcut if requested and on Windows
        if platform.system() == "Windows" and hasattr(self, "create_shortcut_checkbox"):
            if self.create_shortcut_checkbox.isChecked():
                self.create_shortcut()

        # Change the button to Exit
        self.install_button.setText("Exit")
        self.install_button.clicked.disconnect()  # Remove the old signal
        self.install_button.clicked.connect(self.exit_app)

    def installation_failed(self, error_message):
        """Handles errors during installation."""
        self.selected_dir_label.setText(f"Error: {error_message}")

    def create_shortcut(self):
        """Creates a desktop shortcut depending on the OS."""
        if platform.system() != "Windows":
            self.selected_dir_label.setText("Shortcut creation is only supported on Windows.")
            return

        target = None
        for item in Path(self.selected_directory).iterdir():
            if "FragHub" in item.name and item.is_dir():
                target = item / f"{item.name}.exe"
                break

        desktop = Path(os.path.join(os.path.expanduser("~"), "Desktop"))

        if not target or not target.exists():
            self.selected_dir_label.setText("Target for shortcut not found!")
            return

        self.create_windows_shortcut(desktop / "FragHub_1.3.0.lnk", target)

    def create_windows_shortcut(self, shortcut_path, target):
        """Creates a Windows shortcut."""
        import win32com.client  # Windows-specific
        shell = win32com.client.Dispatch("WScript.Shell")
        shortcut = shell.CreateShortcut(str(shortcut_path))
        shortcut.TargetPath = str(target)
        shortcut.WorkingDirectory = str(target.parent)
        shortcut.IconLocation = str(target)
        shortcut.save()

    def create_linux_macos_shortcut(self, shortcut_path, target):
        """Creates a Linux/macOS shortcut."""
        content = f"""
        [Desktop Entry]
        Version=1.0
        Name=FragHub
        Comment=Launcher for FragHub
        Exec={target}
        Path={target.parent}
        Terminal=false
        Type=Application
        Icon=utilities-terminal
        """
        with open(shortcut_path, "w") as shortcut_file:
            shortcut_file.write(content.strip())
        os.chmod(shortcut_path, 0o755)  # Make it executable

    def exit_app(self):
        """Quits the application."""
        QApplication.quit()


if __name__ == "__main__":
    # Define the AppUserModelID for taskbar icon visibility

    app = QApplication([])

    # Set the application icon
    if platform.system() == "Darwin":
        app.setWindowIcon(QIcon(os.path.join(BASE_DIR, r"setup_gui/assets/FragHub_Python_icon.icns")))
    else:
        app.setWindowIcon(QIcon(os.path.join(BASE_DIR, r"setup_gui/assets/FragHub_Python_icon.ico")))


    # Show the installer GUI
    window = InstallerApp()
    window.show()

    # Run the application
    app.exec()

import os
import zipfile
from pathlib import Path
import platform
from PyQt6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QPushButton, QCheckBox,
    QProgressBar, QFileDialog, QLabel, QMessageBox
)
from PyQt6.QtGui import QIcon
from PyQt6.QtCore import Qt, QSize

# Ajout pour définir l'AppUserModelID
import ctypes

# Import Windows-specific modules conditionally
if platform.system() == "Windows":
    try:
        import win32com.client  # PyWin32 for creating Windows shortcuts
    except ImportError:
        win32com = None  # If PyWin32 is not installed, disable shortcut creation functionality


class InstallerApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("FragHub Installer")
        self.resize(400, 300)

        # Set the window and taskbar icon
        window_icon_path = r"D:\Axel\PYTHON\FragHub\scripts\dist\setup_gui\assets\FragHub_Python_icon.ico"
        if os.path.exists(window_icon_path):
            self.setWindowIcon(QIcon(window_icon_path))  # Set icon for both window and taskbar
        else:
            QMessageBox.warning(self, "Missing Icon", f"Taskbar and window icon not found at: {window_icon_path}")

        # Main layout
        layout = QVBoxLayout()

        # Label for the directory selection
        self.info_label = QLabel("Choose Install Directory:")
        self.info_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.info_label.setStyleSheet("font-size: 14px; font-weight: bold; margin-bottom: 10px;")
        layout.addWidget(self.info_label)

        # Button for selecting directory (visible, no white background)
        self.select_dir_button = QPushButton()
        self.configure_directory_button(self.select_dir_button, "setup_gui/assets/directory.png", QSize(64, 64))
        self.select_dir_button.clicked.connect(self.select_directory)
        layout.addWidget(self.select_dir_button, alignment=Qt.AlignmentFlag.AlignCenter)

        # Label to display the selected directory
        self.selected_dir_label = QLabel("No directory selected")
        self.selected_dir_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.selected_dir_label.setWordWrap(True)
        layout.addWidget(self.selected_dir_label)

        # Checkbox to create a shortcut
        self.create_shortcut_checkbox = QCheckBox("Create Desktop Shortcut")
        layout.addWidget(self.create_shortcut_checkbox)

        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setValue(0)
        self.progress_bar.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.progress_bar)

        # Install button
        self.install_button = QPushButton("Install")
        self.install_button.clicked.connect(self.install)
        layout.addWidget(self.install_button)

        # Set main layout
        self.setLayout(layout)

        # Variables
        self.selected_directory = None

    def configure_directory_button(self, button: QPushButton, image_path: str, size: QSize):
        """
        Configures the 'Choose install directory' button as a visible button without a white background.
        """
        if os.path.exists(image_path):
            button.setIcon(QIcon(image_path))
            button.setIconSize(size)
            button.setFixedSize(size.width() + 20, size.height() + 20)  # Bigger button with padding
            button.setStyleSheet(
                """
                QPushButton {
                    border: 2px solid #0078d7; /* Blue border */
                    border-radius: 8px;       /* Rounded corners */
                    background-color: transparent; /* No background */
                }
                QPushButton:hover {
                    border-color: #005fb8; /* Darker blue on hover */
                }
                QPushButton:pressed {
                    border-color: #004c94; /* Even darker blue when pressed */
                }
                """
            )
        else:
            QMessageBox.warning(self, "Missing Icon", f"Could not find the icon at '{image_path}'.")

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
        """Performs the installation: extracts files and creates a shortcut."""
        if not self.selected_directory:
            QMessageBox.warning(self, "Error", "Please select a directory!")
            return

        # Path to the ZIP file
        zip_path = Path(__file__).with_name("FragHub.zip")
        if not zip_path.exists():
            QMessageBox.critical(self, "Error", f"File {zip_path.name} not found.")
            return

        # Extract the ZIP file
        try:
            self.progress_bar.setValue(0)
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                files = zip_ref.namelist()
                total_files = len(files)
                for i, file in enumerate(files):
                    zip_ref.extract(file, self.selected_directory)
                    self.progress_bar.setValue(int((i + 1) / total_files * 100))

            QMessageBox.information(self, "Success", "Extraction complete!")

            # Shortcut creation if requested
            if self.create_shortcut_checkbox.isChecked():
                self.create_shortcut()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred: {e}")

    def create_shortcut(self):
        """Creates a desktop shortcut depending on the OS."""
        target = Path(self.selected_directory) / "FragHub" / "FragHub.exe"
        desktop = Path(os.path.join(os.path.expanduser("~"), "Desktop"))

        if not target.exists():
            QMessageBox.critical(self, "Error", f"Target file for shortcut not found: {target}")
            return

        if platform.system() == "Windows":
            self.create_windows_shortcut(desktop / "FragHub.lnk", target)
        elif platform.system() in ["Linux", "Darwin"]:  # Linux or macOS
            self.create_linux_macos_shortcut(desktop / "FragHub.desktop", target)
        else:
            QMessageBox.warning(self, "Unsupported OS", "Shortcut creation is not supported on this operating system.")

    def create_windows_shortcut(self, shortcut_path, target):
        """Creates a Windows shortcut."""
        if win32com:
            try:
                shell = win32com.client.Dispatch("WScript.Shell")
                shortcut = shell.CreateShortcut(str(shortcut_path))
                shortcut.TargetPath = str(target)
                shortcut.WorkingDirectory = str(target.parent)
                shortcut.IconLocation = str(target)
                shortcut.save()
                QMessageBox.information(self, "Success", "Shortcut created on the desktop.")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Couldn't create shortcut: {e}")
        else:
            QMessageBox.critical(self, "Error", "PyWin32 is not installed. Cannot create shortcut on Windows.")

    def create_linux_macos_shortcut(self, shortcut_path, target):
        """Creates a shortcut on Linux or macOS (as a .desktop file)."""
        try:
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
            os.chmod(shortcut_path, 0o755)  # Make the file executable
            QMessageBox.information(self, "Success", "Shortcut created on the Desktop.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Couldn't create shortcut: {e}")


if __name__ == "__main__":
    # Définir l'AppUserModelID (obligatoire pour l'icône dans la barre des tâches)
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("FragHub.Installer")

    app = QApplication([])

    # Définir l'icône globale pour l'application
    app.setWindowIcon(QIcon(r"D:\Axel\PYTHON\FragHub\scripts\dist\setup_gui\assets\FragHub_Python_icon.ico"))

    # Créer la fenêtre principale
    window = InstallerApp()
    window.show()

    # Lancer l'application
    app.exec()

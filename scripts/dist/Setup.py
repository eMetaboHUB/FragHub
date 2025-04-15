from PyQt6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QPushButton, QCheckBox, QProgressBar,
    QFileDialog, QLabel
)
from PyQt6.QtGui import QIcon
from PyQt6.QtCore import Qt, QSize, QThread, pyqtSignal
import os
import zipfile
from pathlib import Path
import platform
import ctypes  # For setting AppUserModelID (Windows Taskbar Icon)


class InstallerThread(QThread):
    """Thread responsable de l'extraction des fichiers et de l'avancement de la barre de progression."""
    progress_changed = pyqtSignal(int)  # Signal pour mettre à jour la progression
    installation_complete = pyqtSignal()  # Signal pour indiquer que l'installation est terminée
    error_occurred = pyqtSignal(str)  # Signal pour gérer les erreurs

    def __init__(self, zip_path: Path, install_dir: Path):
        super().__init__()
        self.zip_path = zip_path
        self.install_dir = install_dir

    def run(self):
        """Extrait les fichiers du ZIP et met à jour la progression."""
        try:
            with zipfile.ZipFile(self.zip_path, 'r') as zip_ref:
                files = zip_ref.namelist()
                total_files = len(files)

                for i, file in enumerate(files):
                    zip_ref.extract(file, self.install_dir)
                    progress = int((i + 1) / total_files * 100)
                    self.progress_changed.emit(progress)  # Émet le signal pour mettre à jour la barre

            self.installation_complete.emit()  # Émet un signal quand c'est terminé

        except Exception as e:
            self.error_occurred.emit(str(e))  # Émet un signal en cas d'erreur


class InstallerApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("FragHub Installer")
        self.resize(400, 300)

        # Set the window icon
        window_icon_path = r"D:\Axel\PYTHON\FragHub\scripts\dist\setup_gui\assets\FragHub_Python_icon.ico"
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
        self.configure_directory_button(self.select_dir_button, "setup_gui/assets/directory.png", QSize(64, 64))
        self.select_dir_button.clicked.connect(self.select_directory)
        layout.addWidget(self.select_dir_button, alignment=Qt.AlignmentFlag.AlignCenter)

        # Directory label
        self.selected_dir_label = QLabel("No directory selected")
        self.selected_dir_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.selected_dir_label)

        # Checkbox for shortcut creation
        self.create_shortcut_checkbox = QCheckBox("Create Desktop Shortcut")
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
        """Starts the installation in a separate thread."""
        if not self.selected_directory:
            self.selected_dir_label.setText("Please select a directory first!")
            return

        zip_path = Path(__file__).with_name("FragHub.zip")
        if not zip_path.exists():
            self.selected_dir_label.setText("FragHub.zip not found!")
            return

        # Start installation in a separate thread
        self.installer_thread = InstallerThread(zip_path, Path(self.selected_directory))
        self.installer_thread.progress_changed.connect(self.update_progress)
        self.installer_thread.installation_complete.connect(self.installation_complete)
        self.installer_thread.error_occurred.connect(self.installation_failed)
        self.installer_thread.start()  # Lance le thread d'installation

        self.selected_dir_label.setText("Installing...")

    def update_progress(self, progress):
        """Updates the progress bar and percentage label."""
        self.progress_bar.setValue(progress)
        self.progress_label.setText(f"{progress}%")

    def installation_complete(self):
        """Handles completion of installation."""
        self.selected_dir_label.setText("Installation complete!")
        self.progress_bar.setValue(100)
        self.progress_label.setText("100%")

        # Create desktop shortcut if requested
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
        target = Path(self.selected_directory) / "FragHub" / "FragHub.exe"
        desktop = Path(os.path.join(os.path.expanduser("~"), "Desktop"))

        if not target.exists():
            self.selected_dir_label.setText("Target for shortcut not found!")
            return

        if platform.system() == "Windows":
            self.create_windows_shortcut(desktop / "FragHub.lnk", target)
        elif platform.system() in ["Linux", "Darwin"]:  # Linux or macOS
            self.create_linux_macos_shortcut(desktop / "FragHub.desktop", target)

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
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("FragHub.Installer")

    app = QApplication([])

    # Set the application icon
    app.setWindowIcon(QIcon(r"D:\Axel\PYTHON\FragHub\scripts\dist\setup_gui\assets\FragHub_Python_icon.ico"))

    # Show the installer GUI
    window = InstallerApp()
    window.show()

    # Run the application
    app.exec()

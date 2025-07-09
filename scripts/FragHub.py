import sys
import os
import ctypes
import time
import platform
import traceback

from PyQt6.QtWidgets import QApplication, QMessageBox
from PyQt6.QtGui import QPixmap, QIcon
from PyQt6.QtCore import pyqtSignal, QObject, QThread

# Logique de démarrage et gestion des erreurs depuis leur nouvel emplacement
from scripts.GUI.error_handler import exception_hook
from scripts.GUI.splash_screen import LoadingSplashScreen
# Import de la fenêtre principale depuis son nouvel emplacement
from scripts.GUI.main_GUI import MainWindow

if getattr(sys, 'frozen', False):
    BASE_DIR = sys._MEIPASS
else:
    BASE_DIR = os.path.abspath(os.path.dirname(__file__))

if platform.system() == "Windows":
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("FragHub")


class StartupWorker(QObject):
    finished = pyqtSignal(object)
    error = pyqtSignal(str)
    update_splash_message = pyqtSignal(str, int)

    def __init__(self, base_dir):
        super().__init__()
        self._base_dir = base_dir

    def run_startup_tasks(self):
        try:
            self.update_splash_message.emit("Loading FragHub, please wait...", 20)
            time.sleep(1)
            # L'import de MAIN reste correct car on part de la racine
            from scripts.MAIN import MAIN as imported_main
            self.update_splash_message.emit("Initializing main window", 20)
            time.sleep(1)
            self.finished.emit(imported_main)
        except ImportError as e:
            self.error.emit(f"Failed to import 'scripts.MAIN': {e}\n{traceback.format_exc()}")
        except Exception as e:
            self.error.emit(f"Unexpected error during startup tasks: {e}\n{traceback.format_exc()}")


def run_GUI():
    # Ajoute le dossier 'scripts' au path pour assurer que les imports relatifs fonctionnent
    scripts_path = os.path.join(BASE_DIR, 'scripts')
    if scripts_path not in sys.path:
        sys.path.insert(0, scripts_path)

    sys.excepthook = exception_hook
    app = QApplication(sys.argv if hasattr(sys, 'argv') else [''])

    # CORRECTION: Chemin des icônes mis à jour
    scripts_gui_assets = os.path.join(BASE_DIR, "GUI", "assets")
    if platform.system() == "Darwin":
        app_icon_path = os.path.join(BASE_DIR, scripts_gui_assets, "FragHub_Python_icon.icns")
    else:
        app_icon_path = os.path.join(BASE_DIR, scripts_gui_assets, "FragHub_Python_icon.ico")
    if not os.path.exists(app_icon_path):
        app_icon_path = os.path.join(BASE_DIR, scripts_gui_assets, "FragHub_icon.png")
    if os.path.exists(app_icon_path):
        app.setWindowIcon(QIcon(app_icon_path))

    # CORRECTION: Chemin du splash screen mis à jour
    splash_pix_path = os.path.join(BASE_DIR, scripts_gui_assets, "FragHub_icon.png")
    splash_pixmap = QPixmap(splash_pix_path)
    splash = LoadingSplashScreen(splash_pixmap)
    splash.showMessage("Loading FragHub...")
    splash.show()

    startup_thread = QThread()
    startup_worker = StartupWorker(BASE_DIR)
    startup_worker.moveToThread(startup_thread)

    shared_state = {'main_window': None, 'main_function': None, 'splash_screen': splash}

    def on_startup_complete(imported_main_function):
        if not QApplication.instance():
            return
        shared_state['main_function'] = imported_main_function
        try:
            shared_state['main_window'] = MainWindow(main_function_ref=shared_state['main_function'])
            shared_state['main_window'].show()
        except Exception as e:
            on_startup_error(f"Failed to create main window: {e}\n{traceback.format_exc()}")
            return

        if shared_state['splash_screen']:
            shared_state['splash_screen'].close()
            shared_state['splash_screen'] = None

    def on_startup_error(error_message):
        if shared_state['splash_screen']:
            shared_state['splash_screen'].close()
            shared_state['splash_screen'] = None
        QMessageBox.critical(None, "Fatal Startup Error", f"Could not start FragHub:\n{error_message}")
        if QApplication.instance():
            QApplication.instance().quit()

    def update_splash(message, font_size=12):
        if shared_state['splash_screen']:
            shared_state['splash_screen'].showMessage(message, font_size=font_size)

    startup_thread.started.connect(startup_worker.run_startup_tasks)
    startup_worker.finished.connect(on_startup_complete)
    startup_worker.error.connect(on_startup_error)
    startup_worker.update_splash_message.connect(update_splash)

    startup_worker.finished.connect(startup_thread.quit)
    startup_worker.error.connect(startup_thread.quit)
    startup_worker.finished.connect(startup_worker.deleteLater)
    startup_worker.error.connect(startup_worker.deleteLater)
    startup_thread.finished.connect(startup_thread.deleteLater)

    startup_thread.start()
    exit_code = app.exec()
    sys.exit(exit_code)


if __name__ == "__main__":
    run_GUI()
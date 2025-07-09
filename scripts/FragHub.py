import sys
import os
import ctypes
import time
import platform
import traceback
from threading import Thread

from main_worker import run_main_in_worker
from error_handler import show_error_message, exception_hook
from splash_screen import LoadingSplashScreen

if getattr(sys, 'frozen', False):
    BASE_DIR = sys._MEIPASS
else:
    BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__)))

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QHBoxLayout,
    QWidget, QPushButton, QLabel, QTabWidget, QMessageBox
)
from PyQt6.QtGui import QFont, QPixmap, QIcon
from PyQt6.QtCore import Qt, pyqtSignal, QObject, QThread

try:
    from scripts.GUI.tabs.tab_input import InputTab
    from scripts.GUI.tabs.tab_output import OutputTab
    from scripts.GUI.tabs.tab_filters import FiltersTab
    from scripts.GUI.tabs.tab_output_settings import OutputSettingTab
    from scripts.GUI.tabs.tab_projects import ProjectsTab
    from scripts.progress_window import ProgressWindow
except ImportError as e:
    print(f"Erreur critique d'importation des composants UI: {e}")
    sys.exit(f"Erreur critique d'importation: {e}")

if platform.system() == "Windows":
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("FragHub")

class MainWindow(QMainWindow):
    update_progress_signal = pyqtSignal(int)
    error_occurred_signal = pyqtSignal(str)
    task_finished_signal = pyqtSignal()

    def __init__(self, main_function_ref):
        super().__init__()
        self.setWindowFlags(
            Qt.WindowType.Window | Qt.WindowType.WindowMinimizeButtonHint | Qt.WindowType.WindowCloseButtonHint
        )
        self.MAIN_function = main_function_ref
        self.setWindowTitle("FragHub 1.3.2")
        self.setGeometry(100, 100, 1280, 720)
        main_layout = QVBoxLayout()
        banner = QLabel()
        icon_path = os.path.join(BASE_DIR, "GUI", "assets", "FragHub_icon.png")
        pixmap = QPixmap(icon_path)
        if not pixmap.isNull():
            pixmap_scaled = pixmap.scaled(200, 200, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation)
            banner.setPixmap(pixmap_scaled)
        else:
            banner.setText("FragHub Logo")
        banner.setAlignment(Qt.AlignmentFlag.AlignCenter)
        main_layout.addWidget(banner)
        self.tabs = QTabWidget()
        self.input_tab = InputTab(); self.output_tab = OutputTab(); self.filters_tab = FiltersTab(); self.output_settings_tab = OutputSettingTab(); self.projects_tab = ProjectsTab()
        self.tabs.addTab(self.input_tab, "INPUT"); self.tabs.addTab(self.output_tab, "OUTPUT"); self.tabs.addTab(self.filters_tab, "Filters settings"); self.tabs.addTab(self.output_settings_tab, "Output settings"); self.tabs.addTab(self.projects_tab, "Projects settings")
        main_layout.addWidget(self.tabs)
        self.output_tab.output_directory_changed.connect(self.projects_tab.output_directory_changed_signal)
        self.start_button = QPushButton("START")
        self.start_button.setFixedSize(140, 60)
        self.start_button.setFont(QFont("Arial", 16, QFont.Weight.Bold))
        self.start_button.clicked.connect(self.open_progress_window)
        main_layout.addWidget(self.start_button, alignment=Qt.AlignmentFlag.AlignCenter)
        central_widget = QWidget(); central_widget.setLayout(main_layout); self.setCentralWidget(central_widget)
        self.running = False; self.progress_window = None; self.thread = None; self.stop_thread_flag = False
        self.error_occurred_signal.connect(self.handle_execution_error)
        self.task_finished_signal.connect(self.handle_task_finished)

    def open_progress_window(self):
        if not self.running:
            self.showMinimized()
            self.setEnabled(False)
            self.progress_window = ProgressWindow(parent=self)
            # AJOUT: Connexion du signal d'arrêt
            self.progress_window.stop_requested_signal.connect(self.handle_stop_request)
            self.progress_window.show()
            self.start_execution()

    # AJOUT: Méthode pour gérer la demande d'arrêt
    def handle_stop_request(self):
        if self.running:
            self.stop_thread_flag = True

    def start_execution(self):
        if self.running: print("Un traitement est déjà en cours."); return
        self.running = True; self.stop_thread_flag = False
        if self.progress_window: self.progress_window.progress_bar_widget.update_total_items(total=100, completed=0); self.progress_window.progress_bar_widget.update_progress_bar(0)
        if self.thread and self.thread.is_alive(): self.stop_thread_flag = True; self.thread.join(); self.thread = None
        callbacks = {'progress': self.progress_window.update_progress_signal.emit, 'total_items': self.progress_window.update_total_signal.emit, 'prefix': self.progress_window.update_prefix_signal.emit, 'item_type': self.progress_window.update_item_type_signal.emit, 'step': self.progress_window.update_step_signal.emit, 'completion': self.progress_window.completion_callback.emit, 'deletion': self.progress_window.deletion_callback.emit}
        signals = {'error': self.error_occurred_signal, 'finished': self.task_finished_signal}
        stop_flag_provider = lambda: self.stop_thread_flag
        self.thread = Thread(target=run_main_in_worker, args=(self.MAIN_function, callbacks, signals, stop_flag_provider), daemon=True, name="MainWorkerThread_FragHub")
        self.thread.start()

    def handle_task_finished(self):
        self.running = False
        # Ne pas réinitialiser le drapeau ici, laissez clean_exit le faire
        # self.stop_thread_flag = False
        self.thread = None
        if self.progress_window and not self.stop_thread_flag:
            # La tâche s'est terminée normalement
            pass
        elif self.progress_window and self.stop_thread_flag:
            # La tâche a été arrêtée par l'utilisateur
            self.progress_window.close()
        self.stop_thread_flag = False


    def handle_execution_error(self, traceback_str):
        print("handle_execution_error appelé.")
        if self.progress_window: self.progress_window.close()
        show_error_message(parent=self, title="Execution Error", message=traceback_str, on_close=lambda: self.clean_exit(force_quit_app=False))

    def clean_exit(self, force_quit_app=True):
        if self.running and self.thread and self.thread.is_alive(): print("Demande d'arrêt au thread worker..."); self.stop_thread_flag = True; self.thread.join(timeout=3.0)
        if self.progress_window: self.progress_window.close(); self.progress_window = None
        self.running = False
        if force_quit_app: QApplication.instance().quit()
        else: self.setEnabled(True); self.showNormal(); self.activateWindow()

    def closeEvent(self, event):
        if self.running:
            reply = QMessageBox.question(self, 'Quitter FragHub ?', "Un traitement est en cours. Êtes-vous sûr de vouloir quitter ?", QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No, QMessageBox.StandardButton.No)
            if reply == QMessageBox.StandardButton.Yes: self.clean_exit(force_quit_app=True); event.accept()
            else: event.ignore()
        else: event.accept()

class StartupWorker(QObject):
    finished = pyqtSignal(object); error = pyqtSignal(str); update_splash_message = pyqtSignal(str, int)
    def __init__(self, base_dir): super().__init__(); self._base_dir = base_dir
    def run_startup_tasks(self):
        try:
            self.update_splash_message.emit("Loading FragHub, please wait...", 20); time.sleep(1)
            from scripts.MAIN import MAIN as imported_main
            self.update_splash_message.emit("Initializing main window", 20); time.sleep(1)
            self.finished.emit(imported_main)
        except ImportError as e: self.error.emit(f"Failed to import 'scripts.MAIN': {e}\n{traceback.format_exc()}")
        except Exception as e: self.error.emit(f"Unexpected error during startup tasks: {e}\n{traceback.format_exc()}")

def run_GUI():
    sys.excepthook = exception_hook
    app = QApplication(sys.argv if hasattr(sys, 'argv') else [''])
    if platform.system() == "Darwin": app_icon_path = os.path.join(BASE_DIR, "GUI", "assets", "FragHub_Python_icon.icns")
    else: app_icon_path = os.path.join(BASE_DIR, "GUI", "assets", "FragHub_Python_icon.ico")
    if not os.path.exists(app_icon_path): app_icon_path = os.path.join(BASE_DIR, "GUI", "assets", "FragHub_icon.png")
    if os.path.exists(app_icon_path): app.setWindowIcon(QIcon(app_icon_path))
    splash_pix_path = os.path.join(BASE_DIR, "GUI", "assets", "FragHub_icon.png"); splash_pixmap = QPixmap(splash_pix_path); splash = LoadingSplashScreen(splash_pixmap); splash.showMessage("Loading FragHub..."); splash.show()
    startup_thread = QThread(); startup_worker = StartupWorker(BASE_DIR); startup_worker.moveToThread(startup_thread); shared_state = {'main_window': None, 'main_function': None, 'splash_screen': splash}
    def on_startup_complete(imported_main_function):
        if not QApplication.instance(): return
        shared_state['main_function'] = imported_main_function
        try:
            shared_state['main_window'] = MainWindow(main_function_ref=shared_state['main_function'])
            shared_state['main_window'].show()
        except Exception as e: on_startup_error(f"Failed to create main window: {e}\n{traceback.format_exc()}"); return
        if shared_state['splash_screen']: shared_state['splash_screen'].close(); shared_state['splash_screen'] = None
    def on_startup_error(error_message):
        if shared_state['splash_screen']: shared_state['splash_screen'].close(); shared_state['splash_screen'] = None
        QMessageBox.critical(None, "Fatal Startup Error", f"Could not start FragHub:\n{error_message}")
        if QApplication.instance(): QApplication.instance().quit()
    def update_splash(message, font_size=12):
        if shared_state['splash_screen']: shared_state['splash_screen'].showMessage(message, font_size=font_size)
    startup_thread.started.connect(startup_worker.run_startup_tasks); startup_worker.finished.connect(on_startup_complete); startup_worker.error.connect(on_startup_error); startup_worker.update_splash_message.connect(update_splash); startup_worker.finished.connect(startup_thread.quit); startup_worker.error.connect(startup_thread.quit); startup_worker.finished.connect(startup_worker.deleteLater); startup_worker.error.connect(startup_worker.deleteLater); startup_thread.finished.connect(startup_thread.deleteLater)
    startup_thread.start()
    exit_code = app.exec()
    sys.exit(exit_code)

if __name__ == "__main__":
    run_GUI()
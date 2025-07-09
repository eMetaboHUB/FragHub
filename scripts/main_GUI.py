import os
import sys
from threading import Thread

from PyQt6.QtWidgets import (
    QMainWindow, QVBoxLayout, QHBoxLayout, QWidget, QPushButton, QLabel, QTabWidget, QMessageBox, QApplication
)
from PyQt6.QtGui import QFont, QPixmap, QIcon
from PyQt6.QtCore import Qt, pyqtSignal

# Import des composants nécessaires
from main_worker import run_main_in_worker
from error_handler import show_error_message
from scripts.GUI.tabs.tab_input import InputTab
from scripts.GUI.tabs.tab_output import OutputTab
from scripts.GUI.tabs.tab_filters import FiltersTab
from scripts.GUI.tabs.tab_output_settings import OutputSettingTab
from scripts.GUI.tabs.tab_projects import ProjectsTab
from scripts.progress_window import ProgressWindow

if getattr(sys, 'frozen', False):
    BASE_DIR = sys._MEIPASS
else:
    # CORRECTION : Le '..' a été retiré pour pointer vers le bon dossier racine du projet
    BASE_DIR = os.path.abspath(os.path.dirname(__file__))


class MainWindow(QMainWindow):
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

        # Le chemin d'accès devrait maintenant être correct grâce à la correction de BASE_DIR
        icon_path = os.path.join(BASE_DIR, "GUI", "assets", "FragHub_icon.png")
        pixmap = QPixmap(icon_path)

        if not pixmap.isNull():
            pixmap_scaled = pixmap.scaled(200, 200, Qt.AspectRatioMode.KeepAspectRatio,
                                          Qt.TransformationMode.SmoothTransformation)
            banner.setPixmap(pixmap_scaled)
        # CORRECTION : Le bloc "else" avec le texte a été retiré comme demandé.

        banner.setAlignment(Qt.AlignmentFlag.AlignCenter)
        main_layout.addWidget(banner)

        self.tabs = QTabWidget()
        self.input_tab = InputTab()
        self.output_tab = OutputTab()
        self.filters_tab = FiltersTab()
        self.output_settings_tab = OutputSettingTab()
        self.projects_tab = ProjectsTab()

        self.tabs.addTab(self.input_tab, "INPUT")
        self.tabs.addTab(self.output_tab, "OUTPUT")
        self.tabs.addTab(self.filters_tab, "Filters settings")
        self.tabs.addTab(self.output_settings_tab, "Output settings")
        self.tabs.addTab(self.projects_tab, "Projects settings")
        main_layout.addWidget(self.tabs)

        self.output_tab.output_directory_changed.connect(self.projects_tab.output_directory_changed_signal)

        self.start_button = QPushButton("START")
        self.start_button.setFixedSize(140, 60)
        self.start_button.setFont(QFont("Arial", 16, QFont.Weight.Bold))
        self.start_button.clicked.connect(self.open_progress_window)
        main_layout.addWidget(self.start_button, alignment=Qt.AlignmentFlag.AlignCenter)

        central_widget = QWidget()
        central_widget.setLayout(main_layout)
        self.setCentralWidget(central_widget)

        self.running = False
        self.progress_window = None
        self.thread = None
        self.stop_thread_flag = False

        self.error_occurred_signal.connect(self.handle_execution_error)
        self.task_finished_signal.connect(self.handle_task_finished)

    def open_progress_window(self):
        if not self.running:
            self.showMinimized()
            self.setEnabled(False)
            self.progress_window = ProgressWindow(parent=self)
            self.progress_window.stop_requested_signal.connect(self.handle_stop_request)
            self.progress_window.show()
            self.start_execution()

    def handle_stop_request(self):
        if self.running:
            self.stop_thread_flag = True

    def start_execution(self):
        if self.running:
            return

        self.running = True
        self.stop_thread_flag = False

        if self.progress_window:
            self.progress_window.progress_bar_widget.update_total_items(total=100, completed=0)
            self.progress_window.progress_bar_widget.update_progress_bar(0)

        if self.thread and self.thread.is_alive():
            self.stop_thread_flag = True
            self.thread.join()
            self.thread = None

        callbacks = {
            'progress': self.progress_window.update_progress_signal.emit,
            'total_items': self.progress_window.update_total_signal.emit,
            'prefix': self.progress_window.update_prefix_signal.emit,
            'item_type': self.progress_window.update_item_type_signal.emit,
            'step': self.progress_window.update_step_signal.emit,
            'completion': self.progress_window.completion_callback.emit,
            'deletion': self.progress_window.deletion_callback.emit
        }
        signals = {'error': self.error_occurred_signal, 'finished': self.task_finished_signal}
        stop_flag_provider = lambda: self.stop_thread_flag

        self.thread = Thread(target=run_main_in_worker,
                             args=(self.MAIN_function, callbacks, signals, stop_flag_provider), daemon=True,
                             name="MainWorkerThread_FragHub")
        self.thread.start()

    def handle_task_finished(self):
        self.running = False
        self.thread = None
        if self.progress_window and self.stop_thread_flag:
            self.progress_window.close()
        self.stop_thread_flag = False

    def handle_execution_error(self, traceback_str):
        if self.progress_window:
            self.progress_window.close()
        show_error_message(parent=self, title="Execution Error", message=traceback_str,
                           on_close=lambda: self.clean_exit(force_quit_app=False))

    def clean_exit(self, force_quit_app=True):
        if self.running and self.thread and self.thread.is_alive():
            self.stop_thread_flag = True
            self.thread.join(timeout=3.0)
        if self.progress_window:
            self.progress_window.close()
            self.progress_window = None
        self.running = False
        if force_quit_app:
            QApplication.instance().quit()
        else:
            self.setEnabled(True)
            self.showNormal()
            self.activateWindow()

    def closeEvent(self, event):
        if self.running:
            reply = QMessageBox.question(self, 'Quitter FragHub ?',
                                         "Un traitement est en cours. Êtes-vous sûr de vouloir quitter ?",
                                         QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                                         QMessageBox.StandardButton.No)
            if reply == QMessageBox.StandardButton.Yes:
                self.clean_exit(force_quit_app=True)
                event.accept()
            else:
                event.ignore()
        else:
            event.accept()
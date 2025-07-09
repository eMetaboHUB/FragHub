import sys
import os
import traceback
import ctypes
import time
import platform
from threading import Thread

# Import de la fonction worker depuis le nouveau fichier
from main_worker import run_main_in_worker

# Configuration de BASE_DIR
if getattr(sys, 'frozen', False):
    BASE_DIR = sys._MEIPASS
else:
    BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__)))

# Imports PyQt6
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QHBoxLayout,
    QWidget, QPushButton, QLabel, QProgressBar, QTabWidget, QMessageBox, QSplashScreen,
    QSpacerItem, QSizePolicy, QTextEdit
)
from PyQt6.QtGui import QFont, QPixmap, QIcon, QPainter, QColor, QPen
from PyQt6.QtCore import (
    Qt, pyqtSignal, QSize, QTimer, QRectF, QObject, QThread
)

# Imports des composants UI
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


def show_error_message(parent, title, message, on_close=None):
    """Affiche une erreur dans une QMessageBox élargie avec un QTextEdit."""
    msg_box = QMessageBox(parent)
    msg_box.setIcon(QMessageBox.Icon.Critical)
    msg_box.setWindowTitle(title)
    msg_box.setText("An error occurred during execution.")
    msg_box.setStandardButtons(QMessageBox.StandardButton.Ok)

    text_area = QTextEdit()
    text_area.setText(message)
    text_area.setReadOnly(True)
    text_area.setStyleSheet("font-family: Consolas, monospace; font-size: 10pt;")
    text_area.setMinimumSize(700, 350)

    try:
        grid_layout = msg_box.layout()
        grid_layout.addWidget(text_area, grid_layout.rowCount(), 0, 1, grid_layout.columnCount())
        msg_box.setMinimumSize(750, 450)
    except Exception as layout_e:
        print(
            f"Warning: Could not add QTextEdit to QMessageBox layout ({layout_e}). Using setDetailedText as fallback.")
        msg_box.setDetailedText(message)

    clicked_button = msg_box.exec()

    if clicked_button == QMessageBox.StandardButton.Ok and on_close:
        try:
            on_close()
        except Exception as callback_e:
            print(f"ERREUR: Exception dans le callback on_close de show_error_message: {callback_e}")
            traceback.print_exc()


class StreamCapturer:
    """Capture le flux de sortie (stdout / stderr) et le stocke dans un buffer."""

    def __init__(self, original_stream, signal=None):
        self.original_stream = original_stream
        self.signal = signal
        self.buffer = []

    def write(self, message):
        self.buffer.append(message)
        if self.signal:
            self.signal.emit(message)
        self.original_stream.write(message)

    def flush(self):
        self.original_stream.flush()

    def get_captured_text(self):
        return ''.join(self.buffer)

    def clear(self):
        self.buffer = []


class InfiniteSpinnerWidget(QWidget):
    """Widget pour une animation de chargement infinie."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._angle = 0
        self._timer = QTimer(self)
        self._timer.timeout.connect(self._update_angle)
        self._timer.setInterval(25)

        self._spinner_color = QColor(Qt.GlobalColor.white)
        self._pen_width = 3.0
        self.setFixedSize(40, 40)

    def _update_angle(self):
        self._angle = (self._angle - 10) % 360
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        padding = self._pen_width / 2.0
        rect = QRectF(padding, padding, float(self.width()) - self._pen_width, float(self.height()) - self._pen_width)
        if rect.width() <= 0 or rect.height() <= 0:
            return

        pen = QPen(self._spinner_color, self._pen_width, Qt.PenStyle.SolidLine)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        painter.setPen(pen)
        painter.drawArc(rect, self._angle * 16, 90 * 16)

    def start_animation(self):
        if not self._timer.isActive():
            self._timer.start()

    def stop_animation(self):
        if self._timer.isActive():
            self._timer.stop()

    def set_color(self, color: QColor):
        self._spinner_color = color
        self.update()

    def hideEvent(self, event):
        self.stop_animation()
        super().hideEvent(event)

    def showEvent(self, event):
        self.start_animation()
        super().showEvent(event)


class LoadingSplashScreen(QWidget):
    """Écran de chargement personnalisé avec un spinner."""

    def __init__(self, icon_pixmap: QPixmap, parent=None):
        super().__init__(parent)
        self.setWindowFlags(
            Qt.WindowType.SplashScreen | Qt.WindowType.WindowStaysOnTopHint | Qt.WindowType.FramelessWindowHint)
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 20, 20, 20)
        self.icon_label = QLabel()
        if not icon_pixmap.isNull():
            self.icon_label.setPixmap(icon_pixmap)
        self.icon_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.icon_label)
        layout.addSpacerItem(QSpacerItem(20, 15, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed))
        self.spinner = InfiniteSpinnerWidget(self)
        spinner_layout = QHBoxLayout()
        spinner_layout.addStretch()
        spinner_layout.addWidget(self.spinner)
        spinner_layout.addStretch()
        layout.addLayout(spinner_layout)
        layout.addSpacerItem(QSpacerItem(20, 15, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed))
        self.message_label = QLabel("Loading FragHub, please wait...")
        self.message_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.message_label.setStyleSheet("QLabel { color : white; font-weight: bold; }")
        layout.addWidget(self.message_label)
        self.setLayout(layout)
        self.adjustSize()

        if QApplication.primaryScreen():
            center_point = QApplication.primaryScreen().availableGeometry().center()
            self.move(center_point.x() - self.width() // 2, center_point.y() - self.height() // 2)

    def showMessage(self, message, font_size=28, alignment=Qt.AlignmentFlag.AlignCenter, color=Qt.GlobalColor.white):
        self.message_label.setText(message)
        self.message_label.setAlignment(alignment)
        self.message_label.setStyleSheet(f"QLabel {{ color: white; font-weight: bold; font-size: {font_size}px; }}")

    def closeEvent(self, event):
        self.spinner.stop_animation()
        super().closeEvent(event)


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
            pixmap_scaled = pixmap.scaled(
                200, 200, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation
            )
            banner.setPixmap(pixmap_scaled)
        else:
            banner.setText("FragHub Logo")
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

        # --- Connexion des signaux aux slots ---
        self.error_occurred_signal.connect(self.handle_execution_error)
        self.task_finished_signal.connect(self.handle_task_finished)

    def open_progress_window(self):
        if not self.running:
            self.showMinimized()
            self.setEnabled(False)
            self.progress_window = ProgressWindow(parent=self)
            self.progress_window.setWindowFlags(
                Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint | Qt.WindowType.Tool
            )
            self.progress_window.show()
            self.start_execution()

    def start_execution(self):
        """Démarre le traitement en lançant le worker dans un nouveau thread."""
        if self.running:
            print("Un traitement est déjà en cours.")
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

        # Préparation des arguments pour la fonction du worker
        callbacks = {
            'progress': self.progress_window.update_progress_signal.emit,
            'total_items': self.progress_window.update_total_signal.emit,
            'prefix': self.progress_window.update_prefix_signal.emit,
            'item_type': self.progress_window.update_item_type_signal.emit,
            'step': self.progress_window.update_step_signal.emit,
            'completion': self.progress_window.completion_callback.emit,
            'deletion': self.progress_window.deletion_callback.emit
        }
        signals = {
            'error': self.error_occurred_signal,
            'finished': self.task_finished_signal
        }
        stop_flag_provider = lambda: self.stop_thread_flag

        # Lancement du nouveau thread avec la fonction worker comme cible
        self.thread = Thread(
            target=run_main_in_worker,
            args=(self.MAIN_function, callbacks, signals, stop_flag_provider),
            daemon=True,
            name="MainWorkerThread_FragHub"
        )
        self.thread.start()

    def handle_task_finished(self):
        """Slot pour nettoyer l'état après la fin (normale ou non) de la tâche."""
        print("Signal de fin de tâche reçu. Nettoyage en cours...")

        if self.progress_window and self.progress_window.stop_button.text() != 'FINISH':
            try:
                self.progress_window.close()
            except Exception as e:
                print(f"Erreur à la fermeture de progress_window : {e}")

        self.progress_window = None
        self.running = False
        self.thread = None
        self.stop_thread_flag = False

        self.setEnabled(True)
        self.showNormal()
        self.activateWindow()

    def handle_execution_error(self, traceback_str):
        """Slot pour afficher les erreurs venant du thread worker."""
        print("handle_execution_error appelé.")
        show_error_message(
            parent=self,
            title="Execution Error",
            message=traceback_str,
            on_close=self.clean_exit
        )

    def clean_exit(self, force_quit_app=True):
        """Tente d'arrêter le thread worker proprement et quitte l'application."""
        if self.running and self.thread and self.thread.is_alive():
            print("Demande d'arrêt au thread worker...")
            self.stop_thread_flag = True
            self.thread.join(timeout=3.0)

        if self.progress_window:
            self.progress_window.close()
            self.progress_window = None

        self.running = False
        if force_quit_app:
            QApplication.instance().quit()
        else:
            self.show()
            self.activateWindow()

    def closeEvent(self, event):
        """Gère la fermeture de la fenêtre principale."""
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


class StartupWorker(QObject):
    """Effectue les tâches de démarrage longues dans un thread séparé."""
    finished = pyqtSignal(object)
    error = pyqtSignal(str)
    update_splash_message = pyqtSignal(str, int)

    def __init__(self, base_dir):
        super().__init__()
        self._base_dir = base_dir

    def run_startup_tasks(self):
        """Exécute l'importation et autres tâches longues."""
        try:
            self.update_splash_message.emit("Loading FragHub, please wait...", 20)
            time.sleep(1)
            from scripts.MAIN import MAIN as imported_main
            self.update_splash_message.emit("Initializing main window", 20)
            time.sleep(1)
            self.finished.emit(imported_main)
        except ImportError as e:
            error_msg = f"Failed to import 'scripts.MAIN': {e}\n{traceback.format_exc()}"
            self.error.emit(error_msg)
        except Exception as e:
            error_msg = f"Unexpected error during startup tasks: {e}\n{traceback.format_exc()}"
            self.error.emit(error_msg)


def exception_hook(exctype, value, tb):
    """Hook global pour capturer les exceptions non gérées."""
    error_message = ''.join(traceback.format_exception(exctype, value, tb))
    sys.stderr.write(f"Uncaught exception:\n{error_message}")
    if QApplication.instance():
        message_box = QMessageBox()
        message_box.setIcon(QMessageBox.Icon.Critical)
        message_box.setWindowTitle("Unhandled Exception")
        message_box.setText("An critical unexpected error occurred!")
        message_box.setDetailedText(error_message)
        message_box.setMinimumSize(600, 300)
        try:
            text_edit = message_box.findChild(QTextEdit)
            if text_edit:
                text_edit.setMinimumSize(550, 200)
        except Exception:
            pass
        message_box.exec()
    sys.exit(1)


def run_GUI():
    """Fonction principale pour lancer l'application GUI."""
    sys.excepthook = exception_hook
    app = QApplication(sys.argv if hasattr(sys, 'argv') else [''])

    if platform.system() == "Darwin":
        app_icon_path = os.path.join(BASE_DIR, "GUI", "assets", "FragHub_Python_icon.icns")
    else:
        app_icon_path = os.path.join(BASE_DIR, "GUI", "assets", "FragHub_Python_icon.ico")
    if not os.path.exists(app_icon_path):
        app_icon_path = os.path.join(BASE_DIR, "GUI", "assets", "FragHub_icon.png")
    if os.path.exists(app_icon_path):
        app.setWindowIcon(QIcon(app_icon_path))

    splash_pix_path = os.path.join(BASE_DIR, "GUI", "assets", "FragHub_icon.png")
    splash_pixmap = QPixmap(splash_pix_path)
    splash = LoadingSplashScreen(splash_pixmap)
    splash.showMessage("Loading FragHub...")
    splash.show()

    startup_thread = QThread()
    startup_worker = StartupWorker(BASE_DIR)
    startup_worker.moveToThread(startup_thread)
    shared_state = {'main_window': None, 'main_function': None, 'splash_screen': splash}

    def on_startup_complete(imported_main_function):
        if not QApplication.instance(): return
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
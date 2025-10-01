from PyQt6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QProgressBar, QWidget, QSizePolicy,
    QPushButton, QMainWindow, QSplitter, QTabWidget, QScrollArea, QApplication
)
from PyQt6.QtGui import QFont, QPixmap, QIcon
from PyQt6.QtCore import Qt, pyqtSignal
import sys
import time
import ctypes
import os
import platform

if getattr(sys, 'frozen', False):
    BASE_DIR = sys._MEIPASS
else:
    # CORRECTION: Remonte de 2 niveaux (GUI -> scripts -> racine)
    BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

if platform.system() == "Windows":
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("FragHub")


def format_time(time_in_seconds):
    hours, remainder = divmod(time_in_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{int(hours):02}:{int(minutes):02}:{int(seconds):02}"


class ProgressBarWidget(QWidget):
    progress_completed_signal = pyqtSignal(str, str)

    def __init__(self, update_progress_signal, update_total_signal, update_prefix_signal, update_item_type_signal,
                 parent=None):
        super().__init__(parent)
        layout = QHBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(5)
        self.progress_prefix = QLabel("Début...")
        self.progress_prefix.setFont(QFont("Arial", 12))
        self.progress_prefix.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        layout.addWidget(self.progress_prefix)
        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(100)
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(False)
        self.progress_bar.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.progress_bar.setStyleSheet("""
            QProgressBar { height: 24px; border: 1px solid #000; border-radius: 4px; background: #e0e0e0; }
            QProgressBar::chunk { background-color: #3b8dff; border-radius: 4px; }
        """)
        layout.addWidget(self.progress_bar)
        self.progress_suffix = QLabel("0.00%")
        self.progress_suffix.setFont(QFont("Arial", 12))
        self.progress_suffix.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.progress_suffix.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        layout.addWidget(self.progress_suffix)
        self.setLayout(layout)
        self.total_items = 100
        self.start_time = time.time()
        self.item_type = "items"
        update_progress_signal.connect(self.update_progress_bar)
        update_total_signal.connect(self.update_total_items)
        update_prefix_signal.connect(self.update_progress_prefix)
        update_item_type_signal.connect(self.update_item_type)

    def update_item_type(self, item_type):
        self.item_type = item_type

    def update_progress_prefix(self, prefix_text):
        self.progress_prefix.setText(prefix_text)

    def update_progress_bar(self, progress):
        self.progress_bar.setValue(progress)
        progress_percent = (progress / self.total_items) * 100 if self.total_items > 0 else 0
        elapsed_time = time.time() - self.start_time
        items_per_second = progress / elapsed_time if elapsed_time > 0 else 0
        remaining_items = self.total_items - progress
        estimated_time_left = remaining_items / items_per_second if items_per_second > 0 else 0
        self.progress_suffix.setText(
            f"{progress_percent:.2f}% | {progress}/{self.total_items} {self.item_type} "
            f"[{format_time(elapsed_time)} < {format_time(estimated_time_left)}, {items_per_second:.2f} {self.item_type}/s]"
        )
        if progress >= self.total_items:
            self.progress_completed_signal.emit(self.progress_prefix.text(), self.progress_suffix.text())

    def update_total_items(self, total, completed=0):
        self.total_items = total
        self.progress_bar.setMaximum(total)
        self.start_time = time.time()
        self.update_progress_bar(completed)


class ProgressWindow(QMainWindow):
    update_progress_signal = pyqtSignal(int)
    update_total_signal = pyqtSignal(int, int)
    update_prefix_signal = pyqtSignal(str)
    update_item_type_signal = pyqtSignal(str)
    update_step_signal = pyqtSignal(str)
    completion_callback = pyqtSignal(str)
    deletion_callback = pyqtSignal(str)
    stop_requested_signal = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.main_window_ref = parent
        self.setWindowTitle("FragHub 1.4.1")
        # CORRECTION: Chemin de l'icône
        icon_path = os.path.join(BASE_DIR, "GUI", "assets", "FragHub_icon.png")
        self.setWindowIcon(QIcon(icon_path))
        self.setGeometry(100, 100, 1280, 720)
        self.setWindowFlags(
            Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint | Qt.WindowType.Tool | Qt.WindowType.WindowMinimizeButtonHint)
        self.setAttribute(Qt.WidgetAttribute.WA_ShowWithoutActivating, False)

        banner = QLabel()
        # CORRECTION: Chemin du pixmap
        pixmap = QPixmap(icon_path).scaled(200, 200,
                                           Qt.AspectRatioMode.KeepAspectRatio,
                                           Qt.TransformationMode.SmoothTransformation)
        banner.setPixmap(pixmap)
        banner.setAlignment(Qt.AlignmentFlag.AlignCenter)
        splitter = QSplitter(Qt.Orientation.Vertical)
        self.top_tab_widget = QTabWidget()
        self.report_tab = QWidget()
        self.report_layout = QVBoxLayout()
        self.report_widget = QWidget()
        self.report_content = QVBoxLayout()
        self.report_widget.setLayout(self.report_content)
        self.report_content.addStretch()
        self.report_widget.setMinimumWidth(1)
        self.report_scroll = QScrollArea()
        self.report_scroll.setWidgetResizable(True)
        self.report_scroll.setWidget(self.report_widget)
        self.report_layout.addWidget(self.report_scroll)
        self.report_tab.setLayout(self.report_layout)
        self.top_tab_widget.addTab(self.report_tab, "Report")
        splitter.addWidget(self.top_tab_widget)
        self.bottom_tab_widget = QTabWidget()
        self.progress_tab = QWidget()
        self.progress_layout = QVBoxLayout()
        self.progress_bar_widget = ProgressBarWidget(self.update_progress_signal, self.update_total_signal,
                                                     self.update_prefix_signal, self.update_item_type_signal)
        self.progress_bar_widget.progress_completed_signal.connect(self.add_to_report)
        self.progress_layout.addWidget(self.progress_bar_widget)
        self.progress_tab.setLayout(self.progress_layout)
        self.bottom_tab_widget.addTab(self.progress_tab, "Progress")
        splitter.addWidget(self.bottom_tab_widget)
        self.bottom_tab_widget.setMinimumHeight(60)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 0)

        main_layout = QVBoxLayout()
        main_layout.addWidget(banner)
        main_layout.addWidget(splitter)

        self.finish_button = QPushButton("FINISH")
        self.finish_button.setFixedSize(120, 40)
        self.finish_button.setStyleSheet(
            "background-color: green; color: white; font-weight: bold; font-size: 14px; padding: 10px; border-radius: 5px;")
        self.finish_button.clicked.connect(self.finish_button_clicked)
        self.finish_button.hide()

        self.stop_button = QPushButton("STOP")
        self.stop_button.setFixedSize(120, 40)
        self.stop_button.setStyleSheet(
            "background-color: red; color: white; font-weight: bold; font-size: 14px; padding: 10px; border-radius: 5px;")
        self.stop_button.clicked.connect(self.stop_button_clicked)

        button_layout = QHBoxLayout()
        button_layout.addStretch()
        button_layout.addWidget(self.stop_button)
        button_layout.addWidget(self.finish_button)
        button_layout.addStretch()
        main_layout.addLayout(button_layout)
        container = QWidget()
        container.setLayout(main_layout)
        self.setCentralWidget(container)

        self.update_step_signal.connect(self.add_step_to_report)
        self.completion_callback.connect(self.handle_completion)
        self.deletion_callback.connect(self.add_deletion_to_report)

    def stop_button_clicked(self):
        self.stop_requested_signal.emit()
        self.stop_button.setEnabled(False)
        self.stop_button.setText("STOPPING...")

    def finish_button_clicked(self):
        self.close()

    def closeEvent(self, event):
        if self.main_window_ref and self.main_window_ref.running:
            self.stop_requested_signal.emit()
        if self.main_window_ref:
            self.main_window_ref.setEnabled(True)
            self.main_window_ref.showNormal()
            self.main_window_ref.activateWindow()
            self.main_window_ref.progress_window = None
        super().closeEvent(event)

    def handle_completion(self, completion_message):
        while self.progress_layout.count() > 0:
            item = self.progress_layout.takeAt(0)
            if widget := item.widget():
                widget.deleteLater()
        message_label = QLabel(completion_message)
        message_label.setFont(QFont("Arial", 16, QFont.Weight.Bold))
        message_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.progress_layout.addWidget(message_label)
        self.stop_button.hide()
        self.finish_button.show()

    def add_to_report(self, prefix_text, suffix_text):
        report_layout = QHBoxLayout()
        report_layout.setContentsMargins(10, 5, 10, 5)
        report_layout.setSpacing(10)
        prefix_label = QLabel(prefix_text)
        prefix_label.setFont(QFont("Arial", 10, QFont.Weight.Bold))
        prefix_label.setAlignment(Qt.AlignmentFlag.AlignLeft)
        report_layout.addWidget(prefix_label)
        fake_progress_bar = QProgressBar()
        fake_progress_bar.setMinimum(0)
        fake_progress_bar.setMaximum(100)
        fake_progress_bar.setValue(100)
        fake_progress_bar.setTextVisible(False)
        fake_progress_bar.setStyleSheet(
            """QProgressBar { height: 18px; border: 1px solid #000; border-radius: 4px; background: #e0e0e0; } QProgressBar::chunk { background-color: #3b8dff; border-radius: 4px; }""")
        report_layout.addWidget(fake_progress_bar)
        suffix_label = QLabel(suffix_text)
        suffix_label.setFont(QFont("Arial", 10))
        suffix_label.setAlignment(Qt.AlignmentFlag.AlignRight)
        report_layout.addWidget(suffix_label)
        report_widget = QWidget()
        report_widget.setLayout(report_layout)
        self.report_content.insertWidget(self.report_content.count() - 1, report_widget)
        self.report_scroll.verticalScrollBar().setValue(self.report_scroll.verticalScrollBar().maximum())

    def add_step_to_report(self, step_message):
        new_step = QLabel(step_message)
        new_step.setFont(QFont("Arial", 12, QFont.Weight.Bold))
        new_step.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.report_content.insertWidget(self.report_content.count() - 1, new_step)
        self.report_scroll.verticalScrollBar().setValue(self.report_scroll.verticalScrollBar().maximum())

    def add_deletion_to_report(self, deletion_message):
        new_deletion = QLabel(deletion_message)
        new_deletion.setFont(QFont("Arial", 12, QFont.Weight.Normal))
        new_deletion.setAlignment(Qt.AlignmentFlag.AlignCenter)
        new_deletion.setStyleSheet("color: red;")
        self.report_content.insertWidget(self.report_content.count() - 1, new_deletion)
        self.report_scroll.verticalScrollBar().setValue(self.report_scroll.verticalScrollBar().maximum())
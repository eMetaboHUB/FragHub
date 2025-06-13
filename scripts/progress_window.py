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

# If the file is executed as a PyInstaller executable
if getattr(sys, 'frozen', False):
    BASE_DIR = sys._MEIPASS
else:
    # If the file is executed as a Python script
    BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__)))

if platform.system() == "Windows":
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("FragHub")

def format_time(time_in_seconds):
    """Format the time in seconds to HH:mm:ss"""
    hours, remainder = divmod(time_in_seconds, 3600)  # Convert to hours
    minutes, seconds = divmod(remainder, 60)  # Convert to minutes and seconds
    return f"{int(hours):02}:{int(minutes):02}:{int(seconds):02}"


class ProgressBarWidget(QWidget):
    """
    Custom widget for a progress bar dynamically updated through signals.
    """

    # Signal to indicate that progress has reached 100%
    progress_completed_signal = pyqtSignal(str, str)  # (prefix_text, suffix_text)

    def __init__(self, update_progress_signal, update_total_signal, update_prefix_signal, update_item_type_signal,
                 parent=None):
        super().__init__(parent)

        # Main layout to organize the components of the bar
        layout = QHBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)  # No inner margins
        layout.setSpacing(5)  # Spacing between elements

        # Progress bar prefix
        self.progress_prefix = QLabel("Start...")
        self.progress_prefix.setFont(QFont("Arial", 12))  # Text size set to 12
        self.progress_prefix.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        layout.addWidget(self.progress_prefix)

        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(100)
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(False)  # Disable text on the bar itself
        self.progress_bar.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

        # Stylesheet: make the bar thicker
        self.progress_bar.setStyleSheet("""
            QProgressBar {
                height: 24px;
                border: 1px solid #000;
                border-radius: 4px;
                background: #e0e0e0;
            }
            QProgressBar::chunk {
                background-color: #3b8dff;
                border-radius: 4px;
            }
        """)
        layout.addWidget(self.progress_bar)

        # Suffix containing statistics
        self.progress_suffix = QLabel("0.00%")
        self.progress_suffix.setFont(QFont("Arial", 12))  # Text size set to 12
        self.progress_suffix.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.progress_suffix.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        layout.addWidget(self.progress_suffix)

        # Apply the layout to the widget
        self.setLayout(layout)

        # Variables to manage the state
        self.total_items = 100  # Default value
        self.start_time = time.time()
        self.item_type = "items"

        # Signal connections
        update_progress_signal.connect(self.update_progress_bar)
        update_total_signal.connect(self.update_total_items)
        update_prefix_signal.connect(self.update_progress_prefix)
        update_item_type_signal.connect(self.update_item_type)

    def update_item_type(self, item_type):
        """Dynamically update the item unit for display."""
        self.item_type = item_type

    def update_progress_prefix(self, prefix_text):
        """Dynamically update the prefix text."""
        self.progress_prefix.setText(prefix_text)

    def update_progress_bar(self, progress):
        """Update the bar value and display stats."""
        self.progress_bar.setValue(progress)

        # Calculate progress (percentage, estimated time, speed)
        progress_percent = (progress / self.total_items) * 100 if self.total_items > 0 else 0
        elapsed_time = time.time() - self.start_time
        items_per_second = progress / elapsed_time if elapsed_time > 0 else 0
        remaining_items = self.total_items - progress
        estimated_time_left = remaining_items / items_per_second if items_per_second > 0 else 0

        # Update the progress suffix
        self.progress_suffix.setText(
            f"{progress_percent:.2f}% | {progress}/{self.total_items} {self.item_type} "
            f"[{format_time(elapsed_time)} < {format_time(estimated_time_left)}, {items_per_second:.2f} {self.item_type}/s]"
        )

        # Check if the bar reaches 100%
        if progress >= self.total_items:
            # Emit a signal to the ProgressWindow to add a report
            self.progress_completed_signal.emit(self.progress_prefix.text(), self.progress_suffix.text())

    def update_total_items(self, total, completed=0):
        """Update the total number of items and reset if necessary."""
        self.total_items = total
        self.progress_bar.setMaximum(total)

        # Reset the elapsed time
        self.start_time = time.time()  # Current time as a new reference

        # Update the bar with the new totals and current progress (completed)
        self.update_progress_bar(completed)

class ProgressWindow(QMainWindow):
    """
    Main window with a Progress tab (at the bottom) and a dynamic Report tab (at the top).
    """

    # PyQt signals for updating the progress bar
    update_progress_signal = pyqtSignal(int)  # Current progress
    update_total_signal = pyqtSignal(int, int)  # Total and current step
    update_prefix_signal = pyqtSignal(str)  # Prefix text
    update_item_type_signal = pyqtSignal(str)  # Type of items (files, steps, etc.)
    update_step_signal = pyqtSignal(str)  # New signal to display a step in the "Report" tab
    completion_callback = pyqtSignal(str)  # New signal to indicate completion
    deletion_callback = pyqtSignal(str)  # New signal for deletion actions

    def __init__(self, parent=None):
        super().__init__(parent)

        # Store a reference to the main window
        self.main_window_ref = parent  # Assumes 'parent' is an instance of MainWindow.

        # Window title and icon
        self.setWindowTitle("FragHub 1.3.2")
        self.setWindowIcon(QIcon(os.path.join(BASE_DIR, "GUI/assets/FragHub_icon.png")))
        self.setGeometry(100, 100, 1280, 720)

        # Configure the window to appear independently in the taskbar
        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowMinimizeButtonHint)
        self.setAttribute(Qt.WidgetAttribute.WA_ShowWithoutActivating, False)  # Display without stealing focus

        # Add a banner with the icon
        banner = QLabel()
        pixmap = QPixmap(os.path.join(BASE_DIR, "GUI/assets/FragHub_icon.png")).scaled(
            200, 200, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation
        )
        banner.setPixmap(pixmap)
        banner.setAlignment(Qt.AlignmentFlag.AlignCenter)

        # Splitter to divide the space between the two widgets
        splitter = QSplitter(Qt.Orientation.Vertical)

        # First QTabWidget (top) with a "Report" tab
        self.top_tab_widget = QTabWidget()

        # Create the "Report" tab
        self.report_tab = QWidget()
        self.report_layout = QVBoxLayout()
        self.report_widget = QWidget()
        self.report_content = QVBoxLayout()
        self.report_widget.setLayout(self.report_content)

        # Add a stretchable space to manage dynamic content
        self.report_content.addStretch()

        # Set a minimum size to enable scrolling
        self.report_widget.setMinimumWidth(1)

        # Add a scrollable area for the report
        self.report_scroll = QScrollArea()
        self.report_scroll.setWidgetResizable(True)
        self.report_scroll.setWidget(self.report_widget)
        self.report_layout.addWidget(self.report_scroll)

        self.report_tab.setLayout(self.report_layout)
        self.top_tab_widget.addTab(self.report_tab, "Report")

        splitter.addWidget(self.top_tab_widget)

        # "Progress" tab (at the bottom)
        self.bottom_tab_widget = QTabWidget()
        self.progress_tab = QWidget()
        self.progress_layout = QVBoxLayout()

        # Progress bar widget
        self.progress_bar_widget = ProgressBarWidget(
            self.update_progress_signal,
            self.update_total_signal,
            self.update_prefix_signal,
            self.update_item_type_signal
        )
        self.progress_bar_widget.progress_completed_signal.connect(self.add_to_report)  # Connect to report
        self.progress_layout.addWidget(self.progress_bar_widget)

        self.progress_tab.setLayout(self.progress_layout)
        self.bottom_tab_widget.addTab(self.progress_tab, "Progress")

        splitter.addWidget(self.bottom_tab_widget)

        # Configure the splitter's sizes
        self.bottom_tab_widget.setMinimumHeight(60)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 0)

        # Main layout
        main_layout = QVBoxLayout()
        main_layout.addWidget(banner)
        main_layout.addWidget(splitter)

        # STOP button
        self.stop_button = QPushButton("STOP")
        self.stop_button.setFixedSize(120, 40)
        self.stop_button.setStyleSheet("background-color: red; color: white; font-weight: bold; font-size: 12px;")
        self.stop_button.clicked.connect(self.finish_button_clicked)  # Connect to specific handling

        # Center the button
        button_layout = QHBoxLayout()
        button_layout.addStretch()
        button_layout.addWidget(self.stop_button)
        button_layout.addStretch()
        main_layout.addLayout(button_layout)

        # Apply the main layout
        container = QWidget()
        container.setLayout(main_layout)
        self.setCentralWidget(container)

        # Connect signals to actions
        self.update_step_signal.connect(self.add_step_to_report)
        self.completion_callback.connect(self.handle_completion)
        self.deletion_callback.connect(self.add_deletion_to_report)

    def finish_button_clicked(self):
        """Handle the FINISH button in the progress window."""
        if self.stop_button.text() == "FINISH":
            # Reset in the main_window to allow a new start
            if self.main_window_ref:
                self.main_window_ref.running = False
                self.main_window_ref.stop_thread_flag = False  # Reset to enable Start
                self.main_window_ref.thread = None  # Reset the previous thread

            # Close the progress window (return to the main window)
            self.close()
            if self.main_window_ref:
                self.main_window_ref.show()

    def add_to_report(self, prefix_text, suffix_text):
        """
        Add a new line to the report (prefix + static progress bar + suffix).
        """
        # Create a horizontal container for prefix, progress bar, and suffix
        report_layout = QHBoxLayout()
        report_layout.setContentsMargins(10, 5, 10, 5)
        report_layout.setSpacing(10)

        # Prefix
        prefix_label = QLabel(prefix_text)
        prefix_label.setFont(QFont("Arial", 10, QFont.Weight.Bold))
        prefix_label.setAlignment(Qt.AlignmentFlag.AlignLeft)
        report_layout.addWidget(prefix_label)

        # Fake progress bar
        fake_progress_bar = QProgressBar()
        fake_progress_bar.setMinimum(0)
        fake_progress_bar.setMaximum(100)
        fake_progress_bar.setValue(100)  # Set a fixed value (e.g., 50%)
        fake_progress_bar.setTextVisible(False)  # Disable displaying text

        # Apply the same style as your main progress bar
        fake_progress_bar.setStyleSheet("""
            QProgressBar {
                height: 18px;
                border: 1px solid #000;
                border-radius: 4px;
                background: #e0e0e0;
            }
            QProgressBar::chunk {
                background-color: #3b8dff;
                border-radius: 4px;
            }
        """)
        report_layout.addWidget(fake_progress_bar)

        # Suffix
        suffix_label = QLabel(suffix_text)
        suffix_label.setFont(QFont("Arial", 10))
        suffix_label.setAlignment(Qt.AlignmentFlag.AlignRight)
        report_layout.addWidget(suffix_label)

        # Create a container widget to encapsulate the layout
        report_widget = QWidget()
        report_widget.setLayout(report_layout)

        # Add the widget to the report content layout
        self.report_content.insertWidget(self.report_content.count() - 1, report_widget)

        # Force scrolling to the bottom to display the last entry
        self.report_scroll.verticalScrollBar().setValue(
            self.report_scroll.verticalScrollBar().maximum()
        )

    def add_step_to_report(self, step_message):
        """
        Add a new step to the "Report" tab.
        """
        new_step = QLabel(step_message)
        new_step.setFont(QFont("Arial", 12, QFont.Weight.Bold))
        new_step.setAlignment(Qt.AlignmentFlag.AlignCenter)  # Center both horizontally and vertically

        # Add the step before the stretchable elements
        self.report_content.insertWidget(self.report_content.count() - 1, new_step)

        # Force scrolling to the bottom
        self.report_scroll.verticalScrollBar().setValue(
            self.report_scroll.verticalScrollBar().maximum()
        )

    def handle_completion(self, completion_message):
        """
        Display a final message in the Progress tab, replace the progress bar,
        and change the STOP button to FINISH while keeping its behavior.
        """
        # Remove all existing widgets in the progress_layout
        while self.progress_layout.count() > 0:
            item = self.progress_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()  # Properly delete each widget

        # Create a QLabel to display the final message
        message_label = QLabel(completion_message)
        message_label.setFont(QFont("Arial", 16, QFont.Weight.Bold))
        message_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.progress_layout.addWidget(message_label)

        # Change the text and style of the "STOP" button
        self.stop_button.setText("FINISH")  # Change the button text
        self.stop_button.setStyleSheet(
            "background-color: green; color: white; font-weight: bold; font-size: 14px; padding: 10px; border-radius: 5px;"
        )  # Make the button green with white text

    def add_deletion_to_report(self, deletion_message):
        """
        Add a deletion message to the "Report" tab.
        """
        new_deletion = QLabel(deletion_message)
        new_deletion.setFont(QFont("Arial", 12, QFont.Weight.Normal))  # Normal text for deletions
        new_deletion.setAlignment(Qt.AlignmentFlag.AlignCenter)  # Center both horizontally and vertically
        new_deletion.setStyleSheet("color: red;")  # For example: red text to indicate a deletion action

        # Add the message before the stretchable elements
        self.report_content.insertWidget(self.report_content.count() - 1, new_deletion)

        # Force scrolling to the bottom
        self.report_scroll.verticalScrollBar().setValue(
            self.report_scroll.verticalScrollBar().maximum()
        )

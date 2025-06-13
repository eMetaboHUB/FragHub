import sys
import os
import traceback
import ctypes
import time
import platform
from threading import Thread

# Configuration of BASE_DIR
if getattr(sys, 'frozen', False):
    BASE_DIR = sys._MEIPASS
else:
    BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__)))

# PyQt6 Imports
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QHBoxLayout,
    QWidget, QPushButton, QLabel, QProgressBar, QTabWidget, QMessageBox, QSplashScreen,
    QSpacerItem, QSizePolicy, QTextEdit
)
from PyQt6.QtGui import QFont, QPixmap, QIcon, QPainter, QColor, QPen
from PyQt6.QtCore import (
    Qt, pyqtSignal, QSize, QTimer, QRectF, QObject, QThread
)

# Imports of UI components (unchanged)
try:
    from scripts.GUI.tabs.tab_input import InputTab
    from scripts.GUI.tabs.tab_output import OutputTab
    from scripts.GUI.tabs.tab_filters import FiltersTab
    from scripts.GUI.tabs.tab_output_settings import OutputSettingTab
    from scripts.GUI.tabs.tab_projects import ProjectsTab
    from scripts.progress_window import ProgressWindow
except ImportError as e:
    print(f"Critical error importing UI components: {e}")
    # In a real application, show a QMessageBox here before exiting
    sys.exit(f"Critical import error: {e}")


def show_error_message(parent, title, message, on_close=None):  # Add the parameter on_close=None
    """Displays an error in an expanded QMessageBox with a QTextEdit."""
    msg_box = QMessageBox(parent)
    msg_box.setIcon(QMessageBox.Icon.Critical)
    msg_box.setWindowTitle(title)
    msg_box.setText("An error occurred during execution.")
    msg_box.setStandardButtons(QMessageBox.StandardButton.Ok)

    # Option 2: Add a QTextEdit (more control over size/appearance)
    text_area = QTextEdit()
    text_area.setText(message)
    text_area.setReadOnly(True)
    text_area.setStyleSheet("font-family: Consolas, monospace; font-size: 10pt;")
    text_area.setMinimumSize(700, 350)  # Adjust the minimum size

    try:
        grid_layout = msg_box.layout()
        grid_layout.addWidget(text_area, grid_layout.rowCount(), 0, 1, grid_layout.columnCount())
        msg_box.setMinimumSize(750, 450)  # Adjust the dialog size
    except Exception as layout_e:
        print(f"Warning: Could not add QTextEdit to QMessageBox layout ({layout_e}). Using setDetailedText as fallback.")
        msg_box.setDetailedText(message)  # Fallback solution

    # Show the dialog and wait for the user to click a button
    clicked_button = msg_box.exec()

    # --- NEW: Execute the callback AFTER the dialog is closed, if OK was clicked ---
    if clicked_button == QMessageBox.StandardButton.Ok and on_close:
        try:
            on_close()  # Call the function passed as an argument
        except Exception as callback_e:
            # Log an error if the callback itself fails
            print(f"ERROR: Exception in on_close callback of show_error_message: {callback_e}")
            traceback.print_exc()


class StreamCapturer:
    """
    Captures the output stream (stdout / stderr) and stores it in a buffer.
    """

    def __init__(self, original_stream, signal=None):
        self.original_stream = original_stream  # Original stdout or stderr
        self.signal = signal  # Optional signal to transmit captured messages
        self.buffer = []  # Storage for captured messages

    def write(self, message):
        # Store the message in the buffer
        self.buffer.append(message)
        # If a signal is defined, emit it
        if self.signal:
            self.signal.emit(message)
        # Also send the message to the original stream for debugging via terminal
        self.original_stream.write(message)

    def flush(self):
        # Necessary to prevent errors during system calls via sys.stdout/sys.stderr
        self.original_stream.flush()

    def get_captured_text(self):
        # Returns all captured content
        return ''.join(self.buffer)

    def clear(self):
        # Resets the buffer
        self.buffer = []


# --- Custom Spinner Widget (Unchanged from your corrected version) ---
class InfiniteSpinnerWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._angle = 0
        self._timer = QTimer(self)
        self._timer.timeout.connect(self._update_angle)
        self._timer.setInterval(25)  # ~40 FPS

        self._spinner_color = QColor(Qt.GlobalColor.white)
        self._pen_width = 3.0

        self.setFixedSize(40, 40)

    def _update_angle(self):
        """Updates the angle to rotate the animation clockwise."""
        self._angle = (self._angle - 10) % 360  # Clockwise
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        padding = self._pen_width / 2.0
        rect_x = padding
        rect_y = padding
        rect_w = float(self.width()) - self._pen_width
        rect_h = float(self.height()) - self._pen_width

        if rect_w <= 0 or rect_h <= 0:
            return

        draw_rect = QRectF(rect_x, rect_y, rect_w, rect_h)
        pen = QPen(self._spinner_color, self._pen_width, Qt.PenStyle.SolidLine)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        painter.setPen(pen)

        start_angle = self._angle * 16
        span_angle = 90 * 16
        painter.drawArc(draw_rect, start_angle, span_angle)

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
# --- End of Custom Spinner Widget ---


# --- LoadingSplashScreen Class (Unchanged) ---
class LoadingSplashScreen(QWidget):
    #... (Identical to your version)...
    def __init__(self, icon_pixmap: QPixmap, parent=None):
        super().__init__(parent)
        self.setWindowFlags(Qt.WindowType.SplashScreen | Qt.WindowType.WindowStaysOnTopHint | Qt.WindowType.FramelessWindowHint)
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 20, 20, 20)

        self.icon_label = QLabel()
        if not icon_pixmap.isNull():
            self.icon_label.setPixmap(icon_pixmap)
        else:
            self.icon_label.setText("Icon Error")
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
        """Displays a message in the splash screen with a configurable font size."""
        self.message_label.setText(message)
        self.message_label.setAlignment(alignment)
        self.message_label.setStyleSheet(f"""
            QLabel {{
                color: white;
                font-weight: bold;
                font-size: {font_size}px;  /* Apply the specified size */
            }}
        """)

    def closeEvent(self, event):
        self.spinner.stop_animation()
        super().closeEvent(event)
# --- End of LoadingSplashScreen ---


# --- Definition of AppUserModelID (Unchanged) ---
if platform.system() == "Windows":
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("FragHub")


# --- MainWindow Class (Unchanged) ---
class MainWindow(QMainWindow):
    #... (Identical to your version)...
    update_progress_signal = pyqtSignal(int)

    # Signal emitted when an error occurs in the worker thread
    error_occurred_signal = pyqtSignal(str)
    # Signal emitted when the task (MAIN) is finished (success or error)
    task_finished_signal = pyqtSignal()

    def __init__(self, main_function_ref):
        super().__init__()
        self.MAIN_function = main_function_ref
        self.setWindowTitle("FragHub 1.3.2")
        self.setGeometry(100, 100, 1280, 720)

        main_layout = QVBoxLayout()
        banner = QLabel()
        # Using os.path.join for compatibility
        icon_path = os.path.join(BASE_DIR, "GUI", "assets", "FragHub_icon.png")
        pixmap = QPixmap(icon_path)
        if not pixmap.isNull():
            pixmap_scaled = pixmap.scaled(
                200, 200, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation
            )
            banner.setPixmap(pixmap_scaled)
        else:
            print(f"Warning: Could not load banner image at {icon_path}")
            banner.setText("FragHub Logo") # Fallback text
        banner.setAlignment(Qt.AlignmentFlag.AlignCenter)
        main_layout.addWidget(banner)

        self.tabs = QTabWidget()
        self.input_tab = InputTab() # Keep a reference if needed
        self.output_tab = OutputTab()
        self.filters_tab = FiltersTab() # Keep a reference if needed
        self.output_settings_tab = OutputSettingTab() # Keep a reference if needed
        self.projects_tab = ProjectsTab()

        self.tabs.addTab(self.input_tab, "INPUT")
        self.tabs.addTab(self.output_tab, "OUTPUT")
        self.tabs.addTab(self.filters_tab, "Filters settings")
        self.tabs.addTab(self.output_settings_tab, "Output settings")
        self.tabs.addTab(self.projects_tab, "Projects settings")
        main_layout.addWidget(self.tabs)

        # Connecting the signal for the output directory
        self.output_tab.output_directory_changed.connect(self.projects_tab.output_directory_changed_signal)

        self.start_button = QPushButton("START")
        self.start_button.setFixedSize(140, 60)
        self.start_button.setFont(QFont("Arial", 16, QFont.Weight.Bold))
        self.start_button.clicked.connect(self.open_progress_window)
        main_layout.addWidget(self.start_button, alignment=Qt.AlignmentFlag.AlignCenter)

        central_widget = QWidget()
        central_widget.setLayout(main_layout)
        self.setCentralWidget(central_widget)

        # --- State variables ---
        self.running = False
        self.progress_window = None
        self.thread = None  # Reference to the worker thread (type: threading.Thread)
        self.stop_thread_flag = False  # Flag to request the thread to stop

        # --- Connecting signals to slots ---
        self.error_occurred_signal.connect(self.handle_execution_error)


    def open_progress_window(self):
        if not self.running:
            # Hide the main window during execution
            self.hide()
            # Create and show the progress window
            # Make sure ProgressWindow is either a QWidget or QDialog
            self.progress_window = ProgressWindow(parent=self) # parent=None to make it independent
            self.progress_window.show()
            # Start execution in a separate thread
            self.start_execution()

    def start_execution(self):
        """Starts processing when the user presses START."""
        if self.running:
            print("A process is already running.")
            return

        # State reset for a new process
        self.running = True
        self.stop_thread_flag = False  # Reset the stop flag just in case
        if self.progress_window:
            self.progress_window.progress_bar_widget.update_total_items(total=100, completed=0)
            self.progress_window.progress_bar_widget.update_progress_bar(0)

        # Stop a previous thread if it was still active
        if self.thread and self.thread.is_alive():
            self.stop_thread_flag = True
            self.thread.join()  # Wait for the thread to cleanly finish
            self.thread = None

        # Launch a new thread
        self.thread = Thread(target=self.run_main_function, daemon=True, name="MainWorkerThread_FragHub")
        self.thread.start()

    def _cleanup_after_error_or_finish(self):
        """
        Cleans up state after an error (or could be used for a normal end if necessary).
        Closes the progress window, resets the state, and shows the main window.
        """
        print("Running cleanup (_cleanup_after_error_or_finish)...")
        # Close the progress window if it exists
        if self.progress_window:
            try:
                print("Closing progress window during cleanup.")
                self.progress_window.close()
                self.progress_window = None
            except Exception as e:
                print(f"Error while closing progress window during cleanup: {e}")

        # Reset the state (as in handle_task_finished)
        self.running = False
        self.thread = None
        self.stop_thread_flag = False
        print("Showing the main window after cleanup.")
        self.show()
        self.activateWindow()


    def run_main_function(self):
        """
        Wrapper executed in the worker thread. Calls MAIN and handles signals.
        DO NOT directly interact with the UI here (except via emitted signals).
        """
        try:
            if not self.progress_window:
                print("Warning: Progress window not initialized before calling MAIN.")
                # Could raise an exception or continue cautiously

            # Call the main function with callbacks and the stop flag
            # Ensure your MAIN function supports 'stop_flag'
            self.MAIN_function(
                progress_callback=self.progress_window.update_progress_signal.emit,
                total_items_callback=self.progress_window.update_total_signal.emit,
                prefix_callback=self.progress_window.update_prefix_signal.emit,
                item_type_callback=self.progress_window.update_item_type_signal.emit,
                step_callback=self.progress_window.update_step_signal.emit,
                completion_callback=self.progress_window.completion_callback.emit,
                deletion_callback=self.progress_window.deletion_callback.emit,
                # Pass a lambda function checking the stop flag
                stop_flag=lambda: self.stop_thread_flag
            )
        except InterruptedError as ie:
            print(f"Execution interrupted by user: {ie}")
            # Emit a specific error signal or just end cleanly
            # In this case, we don’t consider it an "error" to display to the user
            # The finally block will execute to clean up.
        except Exception:  # Capture ALL other exceptions
            full_traceback = traceback.format_exc()
            print(f"Error during MAIN execution in thread:\n{full_traceback}")
            # Emit the error signal to the GUI thread
            self.error_occurred_signal.emit(full_traceback)
        finally:
            # This block always executes (success, error, interruption)
            # Emit the finish signal so the GUI thread can clean up
            self.task_finished_signal.emit()

    def handle_execution_error(self, traceback_str):
        """
        Slot executed in the GUI thread to display errors from the worker.
        After clicking OK in the error message, the application will close cleanly.
        """
        print("handle_execution_error called.")
        # DO NOT close the progress window here immediately.

        # Show an error dialog.
        # Pass self.clean_exit as a callback to execute after clicking OK.
        show_error_message(
            parent=self,  # The parent can be self (MainWindow), even if it is hidden
            title="Execution Error",
            message=traceback_str,
            on_close=self.clean_exit  # Pass the reference to the clean_exit method
            # By default, clean_exit(force_quit_app=True) will quit the application.
        )

    # --- Methods for clean shutdown ---

    def clean_exit(self, force_quit_app=True):
        """
        Attempts to cleanly stop the worker thread and exit the application.
        :param force_quit_app: If True, quits QApplication at the end.
        """
        if self.running and self.thread and self.thread.is_alive():
            print("Requesting stop from the worker thread...")
            self.stop_thread_flag = True  # Notify the thread to stop

            # Wait a little for the thread to finish via the stop_flag
            wait_start_time = time.time()
            timeout_seconds = 3.0  # Wait max 3 seconds
            while self.thread.is_alive() and (time.time() - wait_start_time) < timeout_seconds:
                QApplication.processEvents()  # Allow the event loop to run
                time.sleep(0.05)  # Small pause

            if self.thread.is_alive():
                print(f"WARNING: Worker thread '{self.thread.name}' did not finish after {timeout_seconds}s.")
                # Here, we could decide to force stop or just proceed with UI cleanup

        # Close the progress window if it still exists
        if self.progress_window:
            try:
                self.progress_window.close()
                self.progress_window = None
            except Exception as e:
                print(f"Error while closing progress_window in clean_exit: {e}")

        self.running = False  # Ensure the state is correct

        if force_quit_app:
            QApplication.instance().quit()
        else:
            # If we don’t quit the app (e.g., just cancel the task),
            # ensure the main window is visible
            self.show()
            self.activateWindow()

    def closeEvent(self, event):
        """Handles the main window close event triggered by the user (clicking X)."""
        if self.running:
            reply = QMessageBox.question(self, 'Quit FragHub ?',
                                         "Are you sure you want to quit?",
                                         QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                                         QMessageBox.StandardButton.No)
            if reply == QMessageBox.StandardButton.Yes:
                self.clean_exit(force_quit_app=True)  # Attempt clean shutdown and exit
                event.accept()  # Accept the window close event

        else:
            # No need to call clean_exit here if we are just quitting the app
            # QApplication.quit() will be triggered by normal window close
            event.accept()  # Accept the close event


# --- End of MainWindow ---


# --- Worker for Startup Tasks (NEW) ---
class StartupWorker(QObject):
    """Performs long startup tasks in a separate thread."""
    # Signal emitted when tasks are finished, passes the imported MAIN function
    finished = pyqtSignal(object)
    # Signal emitted in case of an error
    error = pyqtSignal(str)
    # Signal to update the splash screen message
    update_splash_message = pyqtSignal(str, int)

    def __init__(self, base_dir):
        super().__init__()
        self._base_dir = base_dir

    def showMessage(self, message, font_size=28, alignment=Qt.AlignmentFlag.AlignCenter, color=Qt.GlobalColor.white):
        """Displays a message in the splash screen with a configurable font size."""
        self.message_label.setText(message)
        self.message_label.setAlignment(alignment)
        self.message_label.setStyleSheet(f"""
            QLabel {{
                color: white;
                font-weight: bold;
                font-size: {font_size}px;  /* Dynamic font size adjustment */
            }}
        """)

    def run_startup_tasks(self):
        """Executes import and other long tasks."""
        try:
            # --- Potentially Long Task 1: Import ---
            self.update_splash_message.emit("Loading FragHub, please wait...", 20)  # Message 1
            time.sleep(1)  # Simulates a slight delay to show the spinner

            # Import the main component
            from scripts.MAIN import MAIN as imported_main

            # --- Potentially Long Task 2: Preparation ---
            self.update_splash_message.emit("Initializing main window", 20)  # Message 2
            time.sleep(1)

            # Emit the finished signal on success
            self.finished.emit(imported_main)

        except ImportError as e:
            error_msg = f"Failed to import 'scripts.MAIN': {e}\n{traceback.format_exc()}"
            print(f"ERROR in worker thread: {error_msg}")
            self.error.emit(error_msg)  # Emit the error signal
        except Exception as e:
            error_msg = f"Unexpected error during startup tasks: {e}\n{traceback.format_exc()}"
            print(f"ERROR in worker thread: {error_msg}")
            self.error.emit(error_msg)  # Emit the error signal

# --- Fin de StartupWorker ---


# --- Global Exception Hook ---
def exception_hook(exctype, value, tb):
    error_message = ''.join(traceback.format_exception(exctype, value, tb))
    sys.stderr.write(f"Uncaught exception:\n{error_message}")
    # Attempt to display a QMessageBox even for unhandled errors
    try:
        # Check if QApplication exists before attempting to create a widget
        if QApplication.instance():
            message_box = QMessageBox()
            message_box.setIcon(QMessageBox.Icon.Critical)
            message_box.setWindowTitle("Unhandled Exception")
            message_box.setText("A critical unexpected error occurred!")
            message_box.setDetailedText(error_message)
            # Set a minimum size for the dialog box
            message_box.setMinimumSize(600, 300)
            # Attempt to add detailed text into a scrollable area (if setDetailedText is not enough)
            try:
                text_edit = message_box.findChild(QTextEdit)
                if text_edit:
                    text_edit.setMinimumSize(550, 200)
            except Exception:
                pass  # Ignore if we cannot find/modify the QTextEdit
            message_box.exec()
        else:
            print("QApplication not available, unable to display QMessageBox for the unhandled error.")
    except Exception as e:
        print(f"Unable to display QMessageBox for the unhandled error: {e}")
    sys.exit(1)


# --- run_GUI Function (MODIFIED to use threading) ---
def run_GUI():
    # Set the global exception hook *before* creating QApplication
    sys.excepthook = exception_hook

    app_args = sys.argv if hasattr(sys, 'argv') else ['']
    app = QApplication(app_args)

    # Application icon configuration (unchanged)
    if platform.system() == "Darwin":
        app_icon_path = os.path.join(BASE_DIR, "GUI", "assets", "FragHub_Python_icon.icns")
    else:
        app_icon_path = os.path.join(BASE_DIR, "GUI", "assets", "FragHub_Python_icon.ico")
    if not os.path.exists(app_icon_path):
        app_icon_path = os.path.join(BASE_DIR, "GUI", "assets", "FragHub_icon.png")
    if os.path.exists(app_icon_path):
        app.setWindowIcon(QIcon(app_icon_path))
    else:
        print(f"Warning: Application icon not found at {app_icon_path} or fallback path.")

    # --- Splash Screen Setup ---
    splash_pix_path = os.path.join(BASE_DIR, "GUI", "assets", "FragHub_icon.png")
    splash_pixmap = QPixmap(splash_pix_path)
    if splash_pixmap.isNull(): splash_pixmap = QPixmap(1, 1); splash_pixmap.fill(Qt.GlobalColor.transparent)

    splash = LoadingSplashScreen(splash_pixmap)  # Use a placeholder class if needed
    splash.showMessage("Loading FragHub...")  # Initial message
    splash.show()

    # --- Threading Setup for Startup Tasks ---
    startup_thread = QThread()
    startup_worker = StartupWorker(BASE_DIR)  # Use a placeholder class if needed
    startup_worker.moveToThread(startup_thread)
    shared_state = {'main_window': None, 'main_function': None, 'splash_screen': splash}

    # --- Local Slots ---
    def on_startup_complete(imported_main_function):
        if not QApplication.instance(): return  # Safeguard if the app was quit in the meantime
        shared_state['main_function'] = imported_main_function
        try:
            # Create the MainWindow after worker success
            shared_state['main_window'] = MainWindow(main_function_ref=shared_state['main_function'])
            shared_state['main_window'].show()
        except Exception as e:
            print(f"FATAL ERROR during MainWindow creation: {e}")
            on_startup_error(f"Failed to create main window: {e}\n{traceback.format_exc()}")
            return  # Important: Do not close the splash if the error occurs here

        if shared_state['splash_screen']:
            shared_state['splash_screen'].close()
            shared_state['splash_screen'] = None

    def on_startup_error(error_message):
        if not QApplication.instance(): return
        if shared_state['splash_screen']:
            shared_state['splash_screen'].close()
            shared_state['splash_screen'] = None
        QMessageBox.critical(None, "Fatal Startup Error", f"Could not start FragHub:\n{error_message}")
        app.quit()

    def update_splash(message, font_size=12):  # Default size adjusted
        if shared_state['splash_screen']:
            # Call the showMessage method of the actual splash object
            try:
                # Assume that LoadingSplashScreen has a showMessage(msg) method
                shared_state['splash_screen'].showMessage(message)  # Adjust if the method takes different params
            except Exception as e:
                print(f"Could not update splash message: {e}")


    # --- Connect Signals and Slots ---
    startup_thread.started.connect(startup_worker.run_startup_tasks)
    startup_worker.finished.connect(on_startup_complete)
    startup_worker.error.connect(on_startup_error)
    startup_worker.update_splash_message.connect(update_splash)  # Ensure StartupWorker emits this signal

    # Cleanup (unchanged)
    startup_worker.finished.connect(startup_thread.quit)
    startup_worker.error.connect(startup_thread.quit)  # Quit on error too
    startup_worker.finished.connect(startup_worker.deleteLater)
    startup_worker.error.connect(startup_worker.deleteLater)
    startup_thread.finished.connect(startup_thread.deleteLater)

    # --- Start the thread ---
    startup_thread.start()

    # --- Main Loop ---
    exit_code = app.exec()
    sys.exit(exit_code)


# --- Main Execution Block (Unchanged) ---
if __name__ == "__main__":
    # You can add logs or initializations here if necessary
    run_GUI()
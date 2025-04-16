import sys
import os

# Si le fichier est exécuté comme un exécutable PyInstaller
if getattr(sys, 'frozen', False):
    BASE_DIR = sys._MEIPASS
else:
    # Si le fichier est exécuté comme un script Python
    BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__)))


from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QHBoxLayout,
    QWidget, QPushButton, QLabel, QProgressBar, QTabWidget, QMessageBox, QSplashScreen
)
from PyQt6.QtGui import QFont, QPixmap, QIcon
from PyQt6.QtCore import Qt, pyqtSignal
from threading import Thread
import time
from scripts.GUI.tabs.tab_input import InputTab
from scripts.GUI.tabs.tab_output import OutputTab
from scripts.GUI.tabs.tab_filters import FiltersTab
from scripts.GUI.tabs.tab_output_settings import OutputSettingTab
from scripts.GUI.tabs.tab_projects import ProjectsTab
from scripts.progress_window import ProgressWindow


import traceback
import ctypes

ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("FragHub")


class MainWindow(QMainWindow):
    # Signal pour mettre à jour la barre de progression
    update_progress_signal = pyqtSignal(int)

    def __init__(self):
        super().__init__()
        self.setWindowTitle("FragHub 1.3.0")
        self.setWindowIcon(QIcon(os.path.join(BASE_DIR,"GUI/assets/FragHub_icon.png")))
        self.setGeometry(100, 100, 1280, 720)

        # Création du layout principal
        main_layout = QVBoxLayout()

        # Ajouter le logo FragHub en haut
        banner = QLabel()
        pixmap = QPixmap(os.path.join(BASE_DIR,"GUI/assets/FragHub_icon.png")).scaled(
            200, 200, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation
        )
        banner.setPixmap(pixmap)
        banner.setAlignment(Qt.AlignmentFlag.AlignCenter)
        main_layout.addWidget(banner)

        # Créez les onglets principaux
        self.tabs = QTabWidget()

        # Créer des instances explicites de OutputTab et ProjectsTab
        self.output_tab = OutputTab()
        self.projects_tab = ProjectsTab()

        # Ajouter les instances d'onglet à QTabWidget
        self.tabs.addTab(InputTab(), "INPUT")
        self.tabs.addTab(self.output_tab, "OUTPUT")
        self.tabs.addTab(FiltersTab(), "Filters settings")
        self.tabs.addTab(OutputSettingTab(), "Output settings")
        self.tabs.addTab(self.projects_tab, "Projects settings")
        main_layout.addWidget(self.tabs)

        # Connecter le signal (après avoir ajouté les onglets)
        self.output_tab.output_directory_changed.connect(self.projects_tab.output_directory_changed_signal)

        # Bouton START/STOP
        self.start_button = QPushButton("START")
        self.start_button.setFixedSize(140, 60)
        self.start_button.setFont(QFont("Arial", 16, QFont.Weight.Bold))
        self.start_button.clicked.connect(self.open_progress_window)
        main_layout.addWidget(self.start_button, alignment=Qt.AlignmentFlag.AlignCenter)

        # Ajouter le layout principal au widget central
        central_widget = QWidget()
        central_widget.setLayout(main_layout)
        self.setCentralWidget(central_widget)

        # Variables de contrôle
        self.running = False
        self.progress_window = None
        self.thread = None

    def open_progress_window(self):
        """
        Affiche la fenêtre de progression et masque la fenêtre principale (indépendante).
        """
        if not self.running:
            self.hide()  # Cache simplement la fenêtre principale.

            # Crée une nouvelle fenêtre de progression sans parent
            self.progress_window = ProgressWindow(parent=None)
            self.progress_window.show()  # Affiche la fenêtre de progression

            # Commencer le processus d'exécution
            self.start_execution()

    def start_execution(self):
        """
        Lance l'exécution de MAIN() dans un thread séparé.
        """
        self.running = True
        self.thread = Thread(target=self.run_main_function, daemon=True)
        self.thread.start()

    def run_main_function(self):
        """
        Exécution de la tâche principale avec gestion de la progression.
        """
        try:
            # Lancement du processus principal
            MAIN(
                progress_callback=self.progress_window.update_progress_signal.emit,
                total_items_callback=self.progress_window.update_total_signal.emit,
                prefix_callback=self.progress_window.update_prefix_signal.emit,
                item_type_callback=self.progress_window.update_item_type_signal.emit,
                step_callback=self.progress_window.update_step_signal.emit,
                completion_callback=self.progress_window.completion_callback.emit,
                deletion_callback=self.progress_window.deletion_callback.emit,
            )
        except Exception as e:
            traceback.print_exc()
        finally:
            # Restaurer la fenêtre principale et terminer la progression
            self.show()  # Réaffiche la fenêtre principale
            if self.progress_window:
                self.progress_window.close()  # Ferme la fenêtre de progression
                self.progress_window = None
            self.running = False


def run_GUI():
    # (5a) Initialize QApplication
    app = QApplication(sys.argv)

    # --- Early Application Setup ---
    # Set AppUserModelID for Windows Taskbar behavior
    try:
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("FragHub")
    except AttributeError:
        print("Note: Could not set AppUserModelID (not on Windows or ctypes issue).")
        pass # Ignore if not on Windows or ctypes fails

    # Set Application Icon (used by OS for window, taskbar)
    app_icon_path = os.path.join(BASE_DIR, "GUI/assets/FragHub_Python_icon.ico") # Or.png
    if os.path.exists(app_icon_path):
        app.setWindowIcon(QIcon(app_icon_path))
    else:
        print(f"Warning: Application icon not found at {app_icon_path}")

    # (5b) Create and configure the LoadingWindow (QSplashScreen)
    splash = None # Initialize splash to None
    splash_pix_path = os.path.join(BASE_DIR, "GUI/assets/FragHub_icon.png")
    try:
        splash_pixmap = QPixmap(splash_pix_path)
        if splash_pixmap.isNull():
            print(f"Warning: Could not load splash image at {splash_pix_path}. Using text-only splash.")
            splash = QSplashScreen() # Empty splash
            splash.setWindowFlags(Qt.WindowType.WindowStaysOnTopHint | Qt.WindowType.FramelessWindowHint)
            splash.setEnabled(False)
            splash.showMessage("Loading FragHub, please wait...",
                               Qt.AlignmentFlag.AlignCenter, Qt.GlobalColor.black)
        else:
            # Optional scaling if needed:
            # splash_pixmap = splash_pixmap.scaled(300, 300, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation)
            splash = QSplashScreen(splash_pixmap, Qt.WindowType.WindowStaysOnTopHint | Qt.WindowType.FramelessWindowHint)
            splash.setEnabled(False)
            splash.showMessage("Loading FragHub, please wait...",
                               Qt.AlignmentFlag.AlignBottom | Qt.AlignmentFlag.AlignCenter,
                               Qt.GlobalColor.white) # Adjust color for visibility

        splash.show() # Schedule the splash screen to be shown

    except Exception as e:
        print(f"Error creating splash screen: {e}")
        traceback.print_exc()
        # Continue without splash if creation fails

    # (5c) Force GUI event loop processing to paint the splash screen
    # This is ESSENTIAL before the blocking import. It processes pending events,
    # including painting the splash screen, ensuring it's visible to the user.
    app.processEvents()

    # (5d) Execute the delayed import within error handling
    MAIN_function = None # Variable to hold the imported function/class
    try:
        # *** The critical deferred import ***
        from scripts.MAIN import MAIN as imported_main
        MAIN_function = imported_main # Assign to variable

        if splash: # Optional: Update message after successful import
           splash.showMessage("Initializing main window...",
                              Qt.AlignmentFlag.AlignBottom | Qt.AlignmentFlag.AlignCenter,
                              Qt.GlobalColor.white)
           app.processEvents() # Process message update

    except ImportError as e:
        print(f"FATAL ERROR: Could not import required module 'scripts.MAIN': {e}")
        traceback.print_exc()
        # Show critical error message to the user
        QMessageBox.critical(None, "Fatal Error",
                             f"Failed to load a critical application component (scripts.MAIN).\nError: {e}\nThe application will now exit.")
        if splash:
            splash.close() # Close splash before exiting
        sys.exit(1) # Exit application

    except Exception as e: # Catch other potential errors during import/module execution
        print(f"FATAL ERROR: An unexpected error occurred during MAIN module import: {e}")
        traceback.print_exc()
        QMessageBox.critical(None, "Fatal Error",
                             f"An unexpected error occurred during application startup.\nError: {e}\nThe application will now exit.")
        if splash:
            splash.close()
        sys.exit(1)

    # (5e) Create the main application window
    # The MainWindow class definition needs access to the imported MAIN.
    # Since the import happened *before* instantiation, methods within MainWindow
    # that reference 'MAIN' (as imported originally) should now work correctly.
    # Ensure the MainWindow class definition uses the name 'MAIN' as expected.
    # If MainWindow needs the imported object passed explicitly, modify its __init__.
    # Assuming MainWindow finds 'MAIN' via the now-completed import:
    try:
        # Need to ensure MainWindow class definition is available here
        # It might need to be imported or defined before run_GUI
        from scripts.MAIN import MAIN # Example if MainWindow is in another file
        # Or if MainWindow class is defined in *this* file *before* run_GUI:
        # (Assuming MainWindow class definition is visible here)

        # Pass the imported function if MainWindow's __init__ requires it
        # window = MainWindow(main_function=MAIN_function)
        # Or if MainWindow accesses MAIN globally (as in original code):
        window = MainWindow() # Assuming MainWindow uses the name 'MAIN' internally

    except NameError:
         print("FATAL ERROR: MainWindow class definition not found.")
         traceback.print_exc()
         QMessageBox.critical(None, "Fatal Error", "Could not find the main application window definition.\nThe application will now exit.")
         if splash:
             splash.close()
         sys.exit(1)
    except Exception as e:
         print(f"FATAL ERROR: Failed to initialize MainWindow: {e}")
         traceback.print_exc()
         QMessageBox.critical(None, "Fatal Error", f"Failed to initialize the main application window.\nError: {e}\nThe application will now exit.")
         if splash:
             splash.close()
         sys.exit(1)


    # (5f) Close the LoadingWindow using finish() for smooth transition
    if splash:
        splash.finish(window) # Waits for 'window' to be shown before closing splash

    # (5g) Show the MainWindow
    window.show()

    # (5h) Start the main application event loop
    app.exec()

# --- Main execution block ---
if __name__ == "__main__":
    # Need to import the MainWindow class definition if it's in another file
    # Example: from main_window_module import MainWindow
    # Or ensure the class MainWindow(...): definition exists above run_GUI()

    # Ensure the MainWindow class definition exists before calling run_GUI
    # For this example, assume MainWindow is defined in this file *above* run_GUI
    class MainWindow(QMainWindow): # Placeholder for the actual MainWindow definition
        #... (Full MainWindow class code from user query, ensuring it uses 'MAIN')
        update_progress_signal = pyqtSignal(int)

        def __init__(self):
            super().__init__()
            self.setWindowTitle("FragHub 1.3.0")
            self.setWindowIcon(QIcon(os.path.join(BASE_DIR,"GUI/assets/FragHub_icon.png")))
            self.setGeometry(100, 100, 1280, 720)

            main_layout = QVBoxLayout()
            banner = QLabel()
            pixmap = QPixmap(os.path.join(BASE_DIR,"GUI/assets/FragHub_icon.png")).scaled(
                200, 200, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation
            )
            banner.setPixmap(pixmap)
            banner.setAlignment(Qt.AlignmentFlag.AlignCenter)
            main_layout.addWidget(banner)

            self.tabs = QTabWidget()
            self.output_tab = OutputTab()
            self.projects_tab = ProjectsTab()
            self.tabs.addTab(InputTab(), "INPUT")
            self.tabs.addTab(self.output_tab, "OUTPUT")
            self.tabs.addTab(FiltersTab(), "Filters settings")
            self.tabs.addTab(OutputSettingTab(), "Output settings")
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

        def open_progress_window(self):
            if not self.running:
                self.hide()
                self.progress_window = ProgressWindow(parent=None)
                self.progress_window.show()
                self.start_execution()

        def start_execution(self):
            self.running = True
            # Ensure MAIN is accessible here. It should be, due to the deferred import in run_GUI.
            self.thread = Thread(target=self.run_main_function, daemon=True)
            self.thread.start()

        def run_main_function(self):
            try:
                # This call relies on 'MAIN' being imported successfully earlier
                # Ensure the name 'MAIN' matches what was imported (or use MAIN_function if passed)
                from scripts.MAIN import MAIN # Re-import locally or ensure it's accessible
                MAIN(
                    progress_callback=self.progress_window.update_progress_signal.emit,
                    total_items_callback=self.progress_window.update_total_signal.emit,
                    prefix_callback=self.progress_window.update_prefix_signal.emit,
                    item_type_callback=self.progress_window.update_item_type_signal.emit,
                    step_callback=self.progress_window.update_step_signal.emit,
                    completion_callback=self.progress_window.completion_callback.emit,
                    deletion_callback=self.progress_window.deletion_callback.emit,
                )
            except Exception as e:
                traceback.print_exc()
                # Consider showing an error in the progress window or main window
            finally:
                self.show()
                if self.progress_window:
                    self.progress_window.close()
                    self.progress_window = None
                self.running = False

    # Now call run_GUI
    run_GUI()

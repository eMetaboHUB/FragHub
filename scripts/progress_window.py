from PyQt6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QProgressBar, QWidget, QSizePolicy,
    QPushButton, QMainWindow, QSplitter, QTabWidget, QScrollArea, QApplication
)
from PyQt6.QtGui import QFont, QPixmap, QIcon
from PyQt6.QtCore import Qt, pyqtSignal
import sys
import time
import ctypes

ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("FragHub")

def format_time(time_in_seconds):
    """Format the time in seconds to HH:mm:ss"""
    hours, remainder = divmod(time_in_seconds, 3600)  # Convert to hours
    minutes, seconds = divmod(remainder, 60)  # Convert to minutes and seconds
    return f"{int(hours):02}:{int(minutes):02}:{int(seconds):02}"


class ProgressBarWidget(QWidget):
    """
    Widget personnalisé pour une barre de progression mise à jour dynamiquement par des signaux.
    """

    # Signal pour indiquer que la progression a atteint 100 %
    progress_completed_signal = pyqtSignal(str, str)  # (prefix_text, suffix_text)

    def __init__(self, update_progress_signal, update_total_signal, update_prefix_signal, update_item_type_signal,
                 parent=None):
        super().__init__(parent)

        # Layout principal pour organiser les composants de la barre
        layout = QHBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)  # Pas de marges intérieures
        layout.setSpacing(5)  # Espacement entre les éléments

        # Préfixe de la barre de progression
        self.progress_prefix = QLabel("Début...")
        self.progress_prefix.setFont(QFont("Arial", 12))  # Taille de texte ajustée à 12
        self.progress_prefix.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        layout.addWidget(self.progress_prefix)

        # Barre de progression
        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(100)
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(False)  # Texte désactivé pour la barre elle-même
        self.progress_bar.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

        # Stylesheet : rendre la barre épaisse
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

        # Suffixe contenant les statistiques
        self.progress_suffix = QLabel("0.00%")
        self.progress_suffix.setFont(QFont("Arial", 12))  # Taille de texte ajustée à 12
        self.progress_suffix.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.progress_suffix.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        layout.addWidget(self.progress_suffix)

        # Appliquer le layout au widget
        self.setLayout(layout)

        # Variables pour gérer l’état
        self.total_items = 100  # Valeur par défaut
        self.start_time = time.time()
        self.item_type = "items"

        # Connexion des signaux
        update_progress_signal.connect(self.update_progress_bar)
        update_total_signal.connect(self.update_total_items)
        update_prefix_signal.connect(self.update_progress_prefix)
        update_item_type_signal.connect(self.update_item_type)

    def update_item_type(self, item_type):
        """Met à jour dynamiquement l'unité des items pour l'affichage."""
        self.item_type = item_type

    def update_progress_prefix(self, prefix_text):
        """Met à jour dynamiquement le texte du préfixe."""
        self.progress_prefix.setText(prefix_text)

    def update_progress_bar(self, progress):
        """Met à jour la valeur de la barre et les statistiques affichées."""
        self.progress_bar.setValue(progress)

        # Calcul des progrès (pourcentage, temps estimé, vitesse)
        progress_percent = (progress / self.total_items) * 100 if self.total_items > 0 else 0
        elapsed_time = time.time() - self.start_time
        items_per_second = progress / elapsed_time if elapsed_time > 0 else 0
        remaining_items = self.total_items - progress
        estimated_time_left = remaining_items / items_per_second if items_per_second > 0 else 0

        # Mettre à jour le suffixe de progression
        self.progress_suffix.setText(
            f"{progress_percent:.2f}% | {progress}/{self.total_items} {self.item_type} "
            f"[{format_time(elapsed_time)} < {format_time(estimated_time_left)}, {items_per_second:.2f} {self.item_type}/s]"
        )

        # Vérifier si la barre atteint 100 %
        if progress >= self.total_items:
            # Envoyer un signal vers la ProgressWindow pour ajouter un rapport
            self.progress_completed_signal.emit(self.progress_prefix.text(), self.progress_suffix.text())

    def update_total_items(self, total, completed=0):
        """Met à jour le nombre total d'items et réinitialise si nécessaire."""
        self.total_items = total
        self.progress_bar.setMaximum(total)

        # Réinitialiser le temps écoulé
        self.start_time = time.time()  # Temps actuel comme nouvelle référence

        # Met à jour la barre avec les nouveaux totaux et le progrès actuel (completed)
        self.update_progress_bar(completed)


class ProgressWindow(QMainWindow):
    """
    Fenêtre principale avec un onglet Progress (en bas) et un onglet dynamique Report (en haut).
    """

    # Signaux PyQt pour mettre à jour la barre de progression
    update_progress_signal = pyqtSignal(int)  # Progrès actuel
    update_total_signal = pyqtSignal(int, int)  # Total et étape actuelle
    update_prefix_signal = pyqtSignal(str)  # Texte de préfixe
    update_item_type_signal = pyqtSignal(str)  # Type d'items (fichiers, étapes, etc.)
    update_step_signal = pyqtSignal(str)  # Nouveau signal pour afficher une étape dans l'onglet "Report"
    completion_callback = pyqtSignal(str)  # Nouveau signal pour indiquer la fin
    deletion_callback = pyqtSignal(str)  # Nouveau signal pour l'option de suppression

    def __init__(self, parent=None):
        super().__init__(parent)

        # Titre et icône de la fenêtre
        self.setWindowTitle("FragHub 1.3.0")
        self.setWindowIcon(QIcon("./GUI/assets/FragHub_icon.png"))
        self.setGeometry(100, 100, 1280, 720)

        # Configure la fenêtre pour apparaître de manière indépendante dans la barre des tâches
        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowMinimizeButtonHint)
        self.setAttribute(Qt.WidgetAttribute.WA_ShowWithoutActivating,
                          False)  # S'affiche correctement sans voler le focus

        # Ajouter un banner avec l'icône
        banner = QLabel()
        pixmap = QPixmap("./GUI/assets/FragHub_icon.png").scaled(
            200, 200, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation
        )
        banner.setPixmap(pixmap)
        banner.setAlignment(Qt.AlignmentFlag.AlignCenter)

        # Splitter pour diviser l'espace entre les deux widgets
        splitter = QSplitter(Qt.Orientation.Vertical)

        # Premier QTabWidget (haut) avec un onglet "Report"
        self.top_tab_widget = QTabWidget()

        # Création de l'onglet "Report"
        self.report_tab = QWidget()
        self.report_layout = QVBoxLayout()
        self.report_widget = QWidget()
        self.report_content = QVBoxLayout()
        self.report_widget.setLayout(self.report_content)

        # Ajouter un espace extensible pour gérer le contenu dynamique
        self.report_content.addStretch()

        # Régler la taille minimum pour activer le scroll
        self.report_widget.setMinimumWidth(1)

        # Ajout d'une zone défilable pour le rapport
        self.report_scroll = QScrollArea()
        self.report_scroll.setWidgetResizable(True)
        self.report_scroll.setWidget(self.report_widget)
        self.report_layout.addWidget(self.report_scroll)

        self.report_tab.setLayout(self.report_layout)
        self.top_tab_widget.addTab(self.report_tab, "Report")

        splitter.addWidget(self.top_tab_widget)

        # Onglet "Progress" (en bas)
        self.bottom_tab_widget = QTabWidget()
        self.progress_tab = QWidget()
        self.progress_layout = QVBoxLayout()

        # Widget de barre de progression
        self.progress_bar_widget = ProgressBarWidget(
            self.update_progress_signal,
            self.update_total_signal,
            self.update_prefix_signal,
            self.update_item_type_signal
        )
        self.progress_bar_widget.progress_completed_signal.connect(self.add_to_report)  # Connexion au rapport
        self.progress_layout.addWidget(self.progress_bar_widget)

        self.progress_tab.setLayout(self.progress_layout)
        self.bottom_tab_widget.addTab(self.progress_tab, "Progress")

        splitter.addWidget(self.bottom_tab_widget)

        # Configurer les tailles du splitter
        self.bottom_tab_widget.setMinimumHeight(60)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 0)

        # Layout principal
        main_layout = QVBoxLayout()
        main_layout.addWidget(banner)
        main_layout.addWidget(splitter)

        # Bouton STOP
        self.stop_button = QPushButton("STOP")
        self.stop_button.setFixedSize(120, 40)
        self.stop_button.setStyleSheet("background-color: red; color: white; font-weight: bold; font-size: 12px;")
        self.stop_button.clicked.connect(QApplication.quit)  # Quitte entièrement le programme

        # Centrer le bouton
        button_layout = QHBoxLayout()
        button_layout.addStretch()
        button_layout.addWidget(self.stop_button)
        button_layout.addStretch()
        main_layout.addLayout(button_layout)

        # Appliquer le layout principal
        container = QWidget()
        container.setLayout(main_layout)
        self.setCentralWidget(container)

        # Connexion du signal de mise à jour des étapes
        self.update_step_signal.connect(self.add_step_to_report)

        self.completion_callback.connect(self.handle_completion)

        self.deletion_callback.connect(self.add_deletion_to_report)

    def add_to_report(self, prefix_text, suffix_text):
        """
            Ajoute une nouvelle ligne dans le rapport (prefix + barre de progression statique + suffix).
            """
        # Créer un conteneur horizontal pour le préfixe, la barre et le suffixe
        report_layout = QHBoxLayout()
        report_layout.setContentsMargins(10, 5, 10, 5)
        report_layout.setSpacing(10)

        # Préfixe
        prefix_label = QLabel(prefix_text)
        prefix_label.setFont(QFont("Arial", 10, QFont.Weight.Bold))
        prefix_label.setAlignment(Qt.AlignmentFlag.AlignLeft)
        report_layout.addWidget(prefix_label)

        # Fausse barre de progression
        fake_progress_bar = QProgressBar()
        fake_progress_bar.setMinimum(0)
        fake_progress_bar.setMaximum(100)
        fake_progress_bar.setValue(100)  # Définir une valeur fixe (par exemple, 50%)
        fake_progress_bar.setTextVisible(False)  # Désactiver l'affichage du texte

        # Appliquer le même style que votre barre de progression principale
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

        # Suffixe
        suffix_label = QLabel(suffix_text)
        suffix_label.setFont(QFont("Arial", 10))
        suffix_label.setAlignment(Qt.AlignmentFlag.AlignRight)
        report_layout.addWidget(suffix_label)

        # Créer un widget conteneur pour encapsuler le layout
        report_widget = QWidget()
        report_widget.setLayout(report_layout)

        # Ajouter le widget à la mise en page du contenu du rapport
        self.report_content.insertWidget(self.report_content.count() - 1, report_widget)

        # Forcer le scroll vers le bas pour afficher la dernière entrée
        self.report_scroll.verticalScrollBar().setValue(
            self.report_scroll.verticalScrollBar().maximum()
        )

    def add_step_to_report(self, step_message):
        """
            Ajoute une nouvelle étape dans l'onglet "Report".
            """
        new_step = QLabel(step_message)
        new_step.setFont(QFont("Arial", 12, QFont.Weight.Bold))
        new_step.setAlignment(Qt.AlignmentFlag.AlignCenter)  # Centre horizontalement et verticalement

        # Ajouter l'étape avant les éléments extensibles
        self.report_content.insertWidget(self.report_content.count() - 1, new_step)

        # Forcer le scroll vers le bas
        self.report_scroll.verticalScrollBar().setValue(
            self.report_scroll.verticalScrollBar().maximum()
        )

    def handle_completion(self, completion_message):
        """
        Affiche un message final dans l'onglet Progress
        en remplaçant la barre de progression, et change le bouton STOP en FINISH tout en conservant son comportement.
        """
        # Supprimer tous les widgets existants dans progress_layout
        while self.progress_layout.count() > 0:
            item = self.progress_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()  # Supprime proprement chaque widget

        # Créer un QLabel pour afficher le message final
        message_label = QLabel(completion_message)
        message_label.setFont(QFont("Arial", 16, QFont.Weight.Bold))
        message_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.progress_layout.addWidget(message_label)

        # Changer le texte et le style du bouton "STOP"
        self.stop_button.setText("FINISH")  # Change le texte du bouton
        self.stop_button.setStyleSheet(
            "background-color: green; color: white; font-weight: bold; font-size: 14px; padding: 10px; border-radius: 5px;"
        )  # Rend le bouton vert avec du texte blanc

    def add_deletion_to_report(self, deletion_message):
        """
                Ajoute un message de suppression dans l'onglet "Report".
                """
        new_deletion = QLabel(deletion_message)
        new_deletion.setFont(QFont("Arial", 12, QFont.Weight.Normal))  # Texte normal pour les suppressions
        new_deletion.setAlignment(Qt.AlignmentFlag.AlignCenter)  # Centre horizontalement et verticalement
        new_deletion.setStyleSheet("color: red;")  # Par exemple : texte rouge pour indiquer une action de suppression

        # Ajouter le message avant les éléments extensibles
        self.report_content.insertWidget(self.report_content.count() - 1, new_deletion)

        # Forcer le scroll vers le bas
        self.report_scroll.verticalScrollBar().setValue(
            self.report_scroll.verticalScrollBar().maximum()
        )

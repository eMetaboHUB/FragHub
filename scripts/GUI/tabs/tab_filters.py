from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QSpacerItem, QSizePolicy, QPushButton
)
from PyQt6.QtGui import QFont
from PyQt6.QtCore import QSize, pyqtSignal, Qt
from PyQt6.QtGui import QPainter, QColor, QBrush, QPen
from utils.global_vars import parameters_dict


# Classe toggle personnalisée avec "ON/OFF"
class QToggleSwitch(QWidget):
    stateChanged = pyqtSignal(bool)

    def __init__(self, parent=None, initial_state=True):
        super().__init__(parent)
        self.setFixedSize(60, 30)
        self._state = initial_state
        self.setCursor(Qt.CursorShape.PointingHandCursor)

    def mousePressEvent(self, event):
        self._state = not self._state
        self.update()
        self.stateChanged.emit(self._state)

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        rect = self.rect()
        background_color = QColor("#00C853") if self._state else QColor("#f44336")
        painter.setBrush(QBrush(background_color))
        painter.setPen(Qt.PenStyle.NoPen)
        painter.drawRoundedRect(rect, rect.height() // 2, rect.height() // 2)

        font = QFont("Arial", 8, QFont.Weight.Bold)
        painter.setFont(font)
        text_color = QColor("#FFFFFF")

        if self._state:
            text = "ON"
            text_pos = rect.left() + 5
        else:
            text = "OFF"
            text_pos = rect.right() - 25

        painter.setPen(text_color)
        painter.drawText(text_pos, rect.height() // 2 + 5, text)

        button_color = QColor("#FFFFFF")
        button_x = rect.width() - rect.height() + 2 if self._state else 2
        button_y = 2
        painter.setBrush(QBrush(button_color))
        painter.setPen(QPen(QColor("#E0E0E0"), 1))
        painter.drawEllipse(button_x, button_y, rect.height() - 4, rect.height() - 4)

    def sizeHint(self):
        return QSize(60, 30)

    def isChecked(self):
        return self._state

    def setChecked(self, state):
        self._state = state
        self.update()


# Classe principale des filtres
class FiltersTab(QWidget):
    def __init__(self):
        super().__init__()
        self.layout = QVBoxLayout()
        self.create_filters()
        self.setLayout(self.layout)

    def create_filters(self):
        self.add_filters(self.layout)

    def add_filters(self, layout):
        filters = [
            "normalize_intensity",
            "remove_peak_above_precursormz",
            "check_minimum_peak_requiered",
            "reduce_peak_list",
            "remove_spectrum_under_entropy_score",
            "keep_mz_in_range",
            "check_minimum_of_high_peaks_requiered",
        ]

        global parameters_dict

        # Dictionnaire des messages des filtres (à modifier manuellement si besoin)
        filter_messages = {
            "normalize_intensity": "This function normalizes the intensity of all the peaks in a given spectrum to the maximum intensity.",
            "remove_peak_above_precursormz": "This function removes all peaks from the spectrum whose m/z value is greater than the precursor's m/z value + 5 Da.",
            "check_minimum_peak_requiered": "This function checks whether a given mass spectrum contains a minimum number of peaks. If the spectrum contains fewer peaks than the minimum requirement, it delete the spectrum.",
            "reduce_peak_list": "This function reduces the peak list to a specified maximum number of peaks. The peaks to retain are chosen based on their intensity, with peaks of greater intensity being selected.",
            "remove_spectrum_under_entropy_score": "The entropy score of the spectrum is calculated during processing. If a spectrum has an entropy score lower than the minimum required, it is deleted.",
            "keep_mz_in_range": "This function deletes all spectra whose m/z precursor is not between min and max.",
            "check_minimum_of_high_peaks_requiered": "This function is used to check whether a given array containing peak data has a required minimum number of high peaks.\nA high peak is defined as a peak whose intensity is above a certain percentage (intensity_percent) of the maximum intensity.\nIf the array does not contain a sufficient number of high peaks, the function delete the spectrum.",
        }

        # Layout vertical principal
        main_layout = QVBoxLayout()
        main_layout.setSpacing(5)

        for filter_name in filters:
            # Layout horizontal pour une ligne pour chaque filtre
            filter_layout = QHBoxLayout()
            filter_layout.setSpacing(5)  # Réduit l'espace entre les widgets, rendu compact

            # ToggleSwitch (ON/OFF)
            toggle = QToggleSwitch(initial_state=True)
            toggle.stateChanged.connect(
                lambda state, name=filter_name: self.toggle_filter(state, name)
            )
            filter_layout.addWidget(toggle)

            # Texte du filtre (Label avec le nom)
            label = QLabel(filter_name)
            label.setFont(QFont("Arial", 12))
            label.setFixedHeight(30)
            filter_layout.addWidget(label)

            # Bouton Info 🛈 (proche du texte)
            info_button = QPushButton("🛈")
            info_button.setFixedSize(25, 25)

            # Définir le message spécifique pour le bouton
            if filter_name in filter_messages:
                info_button.setToolTip(filter_messages[filter_name])
            else:
                info_button.setToolTip("Pas de message défini pour ce filtre.")

            filter_layout.addWidget(info_button)

            # Ajouter un espace extensible après le bouton (optionnel)
            filter_layout.addStretch()

            # Ajouter des champs particuliers pour certains filtres
            self.add_additional_fields(filter_layout, filter_name)

            # Finaliser la ligne pour ce filtre
            main_layout.addLayout(filter_layout)

        # Appliquer le layout principal
        layout.addLayout(main_layout)

    def add_additional_fields(self, filter_layout, filter_name):
        if filter_name == "check_minimum_peak_requiered":
            n_peaks_layout = QHBoxLayout()
            n_peaks_label = QLabel("N peaks:")
            n_peaks_label.setFont(QFont("Arial", 10))
            n_peaks_label.setFixedWidth(50)
            n_peaks_layout.addWidget(n_peaks_label)

            text_field = QLineEdit("3")
            parameter_key = "check_minimum_peak_requiered_n_peaks"
            parameters_dict[parameter_key] = 3.0
            text_field.setFixedWidth(60)
            n_peaks_layout.addWidget(text_field)

            text_field.textChanged.connect(
                lambda value, key=parameter_key: self.handle_text_change(key, value)
            )
            filter_layout.addLayout(n_peaks_layout)

        elif filter_name == "reduce_peak_list":
            max_peaks_layout = QHBoxLayout()
            max_peaks_label = QLabel("Max peaks:")
            max_peaks_label.setFont(QFont("Arial", 10))
            max_peaks_label.setFixedWidth(70)
            max_peaks_layout.addWidget(max_peaks_label)

            text_field = QLineEdit("500")
            parameter_key = "reduce_peak_list_max_peaks"
            parameters_dict[parameter_key] = 500.0
            text_field.setFixedWidth(60)
            max_peaks_layout.addWidget(text_field)

            text_field.textChanged.connect(
                lambda value, key=parameter_key: self.handle_text_change(key, value)
            )
            filter_layout.addLayout(max_peaks_layout)

        elif filter_name == "keep_mz_in_range":
            mz_range_layout = QHBoxLayout()

            from_label = QLabel("From:")
            from_label.setFont(QFont("Arial", 10))
            from_label.setFixedWidth(40)
            mz_range_layout.addWidget(from_label)

            from_field = QLineEdit("50")
            from_key = "keep_mz_in_range_from_mz"
            parameters_dict[from_key] = 50.0
            from_field.setFixedWidth(60)
            mz_range_layout.addWidget(from_field)
            from_field.textChanged.connect(
                lambda value, key=from_key: self.handle_text_change(key, value)
            )

            # Spacer
            spacer = QSpacerItem(20, 0, QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Minimum)
            mz_range_layout.addSpacerItem(spacer)

            to_label = QLabel("To:")
            to_label.setFont(QFont("Arial", 10))
            to_label.setFixedWidth(40)
            mz_range_layout.addWidget(to_label)

            to_field = QLineEdit("2000")
            to_key = "keep_mz_in_range_to_mz"
            parameters_dict[to_key] = 2000.0
            to_field.setFixedWidth(60)
            mz_range_layout.addWidget(to_field)
            to_field.textChanged.connect(
                lambda value, key=to_key: self.handle_text_change(key, value)
            )

            filter_layout.addLayout(mz_range_layout)

        elif filter_name == "check_minimum_of_high_peaks_requiered":
            high_peaks_layout = QHBoxLayout()

            intensity_label = QLabel("Intensity %:")
            intensity_label.setFont(QFont("Arial", 10))
            intensity_label.setFixedWidth(80)
            high_peaks_layout.addWidget(intensity_label)

            intensity_field = QLineEdit("5")
            intensity_key = "check_minimum_of_high_peaks_requiered_intensity_percent"
            parameters_dict[intensity_key] = 5
            intensity_field.setFixedWidth(60)
            high_peaks_layout.addWidget(intensity_field)
            intensity_field.textChanged.connect(
                lambda value, key=intensity_key: self.handle_text_change(key, value)
            )

            n_peaks_label = QLabel("N peaks:")
            n_peaks_label.setFont(QFont("Arial", 10))
            n_peaks_label.setFixedWidth(70)
            high_peaks_layout.addWidget(n_peaks_label)

            n_peaks_field = QLineEdit("2")
            n_peaks_key = "check_minimum_of_high_peaks_requiered_no_peaks"
            parameters_dict[n_peaks_key] = 2
            n_peaks_field.setFixedWidth(60)
            high_peaks_layout.addWidget(n_peaks_field)
            n_peaks_field.textChanged.connect(
                lambda value, key=n_peaks_key: self.handle_text_change(key, value)
            )

            filter_layout.addLayout(high_peaks_layout)

        elif filter_name == "remove_spectrum_under_entropy_score":
            entropy_layout = QHBoxLayout()

            score_label = QLabel("Score:")
            score_label.setFont(QFont("Arial", 10))
            score_label.setFixedWidth(50)
            entropy_layout.addWidget(score_label)

            score_field = QLineEdit("0.5")
            score_key = "remove_spectrum_under_entropy_score_value"
            parameters_dict[score_key] = 0.5
            score_field.setFixedWidth(60)
            entropy_layout.addWidget(score_field)
            score_field.textChanged.connect(
                lambda value, key=score_key: self.handle_text_change(key, value)
            )

            filter_layout.addLayout(entropy_layout)

    def toggle_filter(self, state, filter_name):
        parameters_dict[filter_name] = 1.0 if state else 0.0

    def handle_text_change(self, key, value):
        try:
            parameters_dict[key] = float(value)
        except ValueError:
            pass

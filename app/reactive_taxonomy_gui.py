"""Small Qt6 wrapper for molecule and reaction featurization."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Literal, Sequence

from PyQt6 import QtCore, QtGui, QtWidgets

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from reactive_taxonomy import featurize_molecule, featurize_reaction  # noqa: E402
from reactive_taxonomy.cli import format_concise_analysis  # noqa: E402


InputKind = Literal["molecule", "reaction"]

REACTION_EXAMPLE = (
    "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
)
MOLECULE_EXAMPLE = "Brc1ccc(N)cc1C#N"


def detect_input_kind(text: str) -> InputKind:
    """Classify input from reaction-SMILES delimiters.

    The ``>`` character is reserved for reaction separators and is not part of
    ordinary molecular SMILES. Detecting any separator also routes malformed
    reaction input to the reaction parser, which can return the useful error.
    """
    return "reaction" if ">" in text else "molecule"


def featurize_text(text: str) -> tuple[InputKind, object]:
    """Featurize stripped text through the appropriate public taxonomy API."""
    value = text.strip()
    if not value:
        raise ValueError("Enter a molecule or reaction SMILES.")
    kind = detect_input_kind(value)
    analysis = (
        featurize_reaction(value)
        if kind == "reaction"
        else featurize_molecule(value)
    )
    return kind, analysis


class ReactiveTaxonomyWindow(QtWidgets.QMainWindow):
    """Auto-detecting molecule and reaction featurization window."""

    def __init__(self) -> None:
        super().__init__()
        self.setObjectName("reactiveTaxonomyWindow")
        self.setFont(QtGui.QFont("Segoe UI", 9))
        self.setWindowTitle("Reactive Taxonomy Featurizer")
        self.resize(900, 620)
        self._build_ui()
        self._connect_signals()
        self._apply_style()
        self._update_detected_kind("")

    def _build_ui(self) -> None:
        central = QtWidgets.QWidget()
        central.setObjectName("centralPanel")
        self.setCentralWidget(central)
        layout = QtWidgets.QVBoxLayout(central)
        layout.setContentsMargins(28, 24, 28, 24)
        layout.setSpacing(14)

        title = QtWidgets.QLabel("Molecule & reaction featurizer")
        title.setObjectName("title")
        layout.addWidget(title)

        description = QtWidgets.QLabel(
            "Paste a molecular SMILES or reaction SMILES. Input containing a "
            "reaction separator (>) is analyzed as a reaction automatically."
        )
        description.setObjectName("description")
        description.setWordWrap(True)
        layout.addWidget(description)

        input_row = QtWidgets.QHBoxLayout()
        input_row.setSpacing(10)
        self.input_edit = QtWidgets.QLineEdit()
        self.input_edit.setObjectName("smilesInput")
        self.input_edit.setClearButtonEnabled(True)
        self.input_edit.setPlaceholderText(
            "SMILES or reactants>>product"
        )
        self.input_edit.setAccessibleName("Molecule or reaction SMILES")
        input_row.addWidget(self.input_edit, 1)

        self.analyze_button = QtWidgets.QPushButton("Analyze")
        self.analyze_button.setObjectName("analyzeInput")
        self.analyze_button.setDefault(True)
        input_row.addWidget(self.analyze_button)
        layout.addLayout(input_row)

        controls = QtWidgets.QHBoxLayout()
        controls.setSpacing(8)
        self.kind_label = QtWidgets.QLabel()
        self.kind_label.setObjectName("detectedKind")
        controls.addWidget(self.kind_label)
        controls.addStretch(1)

        self.reaction_example_button = QtWidgets.QPushButton("Reaction example")
        self.reaction_example_button.setObjectName("reactionExample")
        controls.addWidget(self.reaction_example_button)

        self.molecule_example_button = QtWidgets.QPushButton("Molecule example")
        self.molecule_example_button.setObjectName("moleculeExample")
        controls.addWidget(self.molecule_example_button)

        self.copy_button = QtWidgets.QPushButton("Copy result")
        self.copy_button.setObjectName("copyResult")
        self.copy_button.setEnabled(False)
        controls.addWidget(self.copy_button)
        layout.addLayout(controls)

        self.output = QtWidgets.QPlainTextEdit()
        self.output.setObjectName("analysisOutput")
        self.output.setReadOnly(True)
        self.output.setPlaceholderText("Concise featurization results appear here.")
        self.output.setLineWrapMode(
            QtWidgets.QPlainTextEdit.LineWrapMode.WidgetWidth
        )
        fixed_font = QtGui.QFontDatabase.systemFont(
            QtGui.QFontDatabase.SystemFont.FixedFont
        )
        fixed_font.setPointSize(10)
        self.output.setFont(fixed_font)
        layout.addWidget(self.output, 1)

        self.status_label = QtWidgets.QLabel("Ready")
        self.status_label.setObjectName("statusLabel")
        layout.addWidget(self.status_label)

        QtGui.QShortcut(
            QtGui.QKeySequence("Ctrl+Return"),
            self,
            activated=self.analyze,
        )
        QtGui.QShortcut(
            QtGui.QKeySequence("Ctrl+Enter"),
            self,
            activated=self.analyze,
        )

    def _connect_signals(self) -> None:
        self.input_edit.textChanged.connect(self._update_detected_kind)
        self.input_edit.returnPressed.connect(self.analyze)
        self.analyze_button.clicked.connect(self.analyze)
        self.reaction_example_button.clicked.connect(
            lambda: self._load_example(REACTION_EXAMPLE)
        )
        self.molecule_example_button.clicked.connect(
            lambda: self._load_example(MOLECULE_EXAMPLE)
        )
        self.copy_button.clicked.connect(self.copy_result)

    def _apply_style(self) -> None:
        """Add converter-style accents while preserving the native dark palette."""
        self.setStyleSheet(
            """
            QLabel#title {
                font-size: 22px;
                font-weight: 600;
            }
            QLabel#description {
                color: #d0d7de;
                font-size: 13px;
            }
            QPushButton#analyzeInput {
                background-color: #0078d7;
                color: white;
                font-weight: 700;
                border: none;
                border-radius: 6px;
                padding: 8px 13px;
                padding-left: 20px;
                padding-right: 20px;
            }
            QPushButton#analyzeInput:hover {
                background-color: #1689e5;
            }
            QPushButton#analyzeInput:disabled {
                background-color: #355b78;
                color: #aebdcc;
            }
            QLabel#detectedKind {
                color: #6cb6ff;
                font-weight: 600;
            }
            QLabel#statusLabel {
                color: #c5ced8;
            }
            """
        )

    @QtCore.pyqtSlot(str)
    def _update_detected_kind(self, text: str) -> None:
        value = text.strip()
        if not value:
            self.kind_label.setText("Detected: waiting for input")
            return
        kind = detect_input_kind(value)
        self.kind_label.setText(f"Detected: {kind}")

    def _load_example(self, text: str) -> None:
        self.input_edit.setText(text)
        self.input_edit.setFocus()
        self.input_edit.selectAll()

    @QtCore.pyqtSlot()
    def analyze(self) -> None:
        """Analyze the current input and display the concise CLI-equivalent view."""
        self.analyze_button.setEnabled(False)
        self.status_label.setText("Analyzing…")
        QtWidgets.QApplication.processEvents(
            QtCore.QEventLoop.ProcessEventsFlag.ExcludeUserInputEvents
        )
        try:
            kind, analysis = featurize_text(self.input_edit.text())
            heading = f"{kind.upper()} FEATURIZATION"
            self.output.setPlainText(
                f"{heading}\n\n{format_concise_analysis(analysis)}"
            )
            valid = bool(getattr(analysis, "valid", False))
            state = "valid" if valid else "invalid"
            self.status_label.setText(
                f"Complete · {kind} input · {state}"
            )
            self.copy_button.setEnabled(True)
        except Exception as exc:
            self.output.setPlainText(f"Unable to analyze input.\n\n{exc}")
            self.status_label.setText("Analysis failed")
            self.copy_button.setEnabled(True)
        finally:
            self.analyze_button.setEnabled(True)

    @QtCore.pyqtSlot()
    def copy_result(self) -> None:
        """Copy the visible analysis to the system clipboard."""
        QtWidgets.QApplication.clipboard().setText(self.output.toPlainText())
        self.status_label.setText("Result copied to clipboard")


def main(argv: Sequence[str] | None = None) -> int:
    """Launch the Qt6 featurization application."""
    application = QtWidgets.QApplication(
        list(argv) if argv is not None else sys.argv
    )
    application.setApplicationName("Reactive Taxonomy Featurizer")
    window = ReactiveTaxonomyWindow()
    window.show()
    return application.exec()


if __name__ == "__main__":
    raise SystemExit(main())

"""Small Qt6 wrapper for molecule and reaction featurization."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Literal, Sequence

from PyQt6 import QtCore, QtGui, QtWidgets

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from reactive_taxonomy import (  # noqa: E402
    AtomMappingProvider,
    ExternalMappingAssessment,
    RxnMapperProvider,
    analyze_reaction_with_external_mapping,
    featurize_molecule,
    featurize_reaction,
)
from reactive_taxonomy.cli import format_concise_analysis  # noqa: E402
from visualization import (  # noqa: E402
    build_reaction_core_graphic,
    render_molecule_image_bytes,
    render_reaction_image_bytes,
)
from visualization.qt_widgets import StructureImageLabel  # noqa: E402


InputKind = Literal["molecule", "reaction"]

REACTION_EXAMPLE = (
    "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
)
MOLECULE_EXAMPLE = "Brc1ccc(N)cc1C#N"
REACTION_IMAGE_SIZE = (680, 168)
MOLECULE_IMAGE_SIZE = (480, 300)


def detect_input_kind(text: str) -> InputKind:
    """Classify input from reaction-SMILES delimiters.

    The ``>`` character is reserved for reaction separators and is not part of
    ordinary molecular SMILES. Detecting any separator also routes malformed
    reaction input to the reaction parser, which can return the useful error.
    """
    return "reaction" if ">" in text else "molecule"


def featurize_text(
    text: str,
    *,
    mapping_provider: AtomMappingProvider | None = None,
    force_resolved_mapping: bool = False,
) -> tuple[InputKind, object, ExternalMappingAssessment | None]:
    """Featurize stripped text through the appropriate public taxonomy API."""
    value = text.strip()
    if not value:
        raise ValueError("Enter a molecule or reaction SMILES.")
    kind = detect_input_kind(value)
    if kind == "reaction":
        base_analysis = featurize_reaction(value)
        assessment = (
            analyze_reaction_with_external_mapping(
                value,
                mapping_provider,
                base_analysis=base_analysis,
                force_resolved_shadow=force_resolved_mapping,
            )
            if mapping_provider is not None
            else None
        )
        return (
            kind,
            assessment.analysis if assessment is not None else base_analysis,
            assessment,
        )
    return kind, featurize_molecule(value), None


class ReactiveTaxonomyWindow(QtWidgets.QMainWindow):
    """Auto-detecting molecule and reaction featurization window."""

    def __init__(self) -> None:
        super().__init__()
        self.setObjectName("reactiveTaxonomyWindow")
        self.setFont(QtGui.QFont("Segoe UI", 9))
        self.setWindowTitle("Reactive Taxonomy Featurizer")
        self.resize(900, 620)
        self._mapping_provider: RxnMapperProvider | None = None
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
        self.use_rxnmapper_check = QtWidgets.QCheckBox(
            "Use RXNMapper for unresolved or ambiguous reactions"
        )
        self.use_rxnmapper_check.setObjectName("useRxnMapper")
        self.use_rxnmapper_check.setChecked(True)
        self.use_rxnmapper_check.setToolTip(
            "Checked by default. Supplied maps and resolved internal evidence "
            "still take precedence; generated mapping remains review evidence."
        )
        controls.addWidget(self.use_rxnmapper_check)
        self.force_core_mapping_check = QtWidgets.QCheckBox(
            "Map resolved reactions for minimized graphic"
        )
        self.force_core_mapping_check.setObjectName("forceCoreMapping")
        self.force_core_mapping_check.setChecked(False)
        self.force_core_mapping_check.setToolTip(
            "Optional and slower. Runs RXNMapper even when internal evidence "
            "already resolves the reaction, so a mapped graphical core can be "
            "drawn. The mapped result remains review evidence."
        )
        controls.addWidget(self.force_core_mapping_check)
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

        result_panel = QtWidgets.QWidget()
        result_layout = QtWidgets.QHBoxLayout(result_panel)
        result_layout.setContentsMargins(0, 0, 0, 0)
        result_layout.setSpacing(12)

        analysis_column = QtWidgets.QWidget()
        analysis_layout = QtWidgets.QVBoxLayout(analysis_column)
        analysis_layout.setContentsMargins(0, 0, 0, 0)
        analysis_layout.setSpacing(4)
        analysis_layout.addWidget(QtWidgets.QLabel("Analysis"))

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
        analysis_layout.addWidget(self.output)

        graph_column = QtWidgets.QWidget()
        graph_layout = QtWidgets.QVBoxLayout(graph_column)
        graph_layout.setContentsMargins(0, 0, 0, 0)
        graph_layout.setSpacing(4)
        self.graph_heading = QtWidgets.QLabel("Structure graph")
        self.graph_heading.setObjectName("structureGraphHeading")
        graph_layout.addWidget(self.graph_heading)
        self.graph_tabs = QtWidgets.QTabWidget()
        self.graph_tabs.setObjectName("reactionGraphicTabs")

        self.structure_image_label = StructureImageLabel(
            placeholder="Reaction or compound graph will appear here.",
            object_name="featurizedStructureGraph",
            minimum_height=220,
        )
        self.graph_tabs.addTab(self.structure_image_label, "Full structure")

        minimized_panel = QtWidgets.QWidget()
        minimized_layout = QtWidgets.QVBoxLayout(minimized_panel)
        minimized_layout.setContentsMargins(0, 0, 0, 0)
        minimized_layout.setSpacing(4)
        self.core_image_label = StructureImageLabel(
            placeholder="Mapped minimized reaction will appear here.",
            object_name="minimizedReactionGraphic",
            minimum_height=190,
        )
        minimized_layout.addWidget(self.core_image_label, 1)
        self.core_graphic_note = QtWidgets.QLabel(
            "A mapped reaction core is required."
        )
        self.core_graphic_note.setObjectName("coreGraphicNote")
        self.core_graphic_note.setWordWrap(True)
        minimized_layout.addWidget(self.core_graphic_note)
        self.graph_tabs.addTab(minimized_panel, "Minimized reaction")
        graph_layout.addWidget(self.graph_tabs)

        result_layout.addWidget(analysis_column, stretch=1)
        result_layout.addWidget(graph_column, stretch=1)
        layout.addWidget(result_panel, 1)

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
        self.use_rxnmapper_check.toggled.connect(
            self.force_core_mapping_check.setEnabled
        )
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
            input_text = self.input_edit.text()
            requested_kind = detect_input_kind(input_text.strip())
            mapping_provider = None
            if (
                requested_kind == "reaction"
                and self.use_rxnmapper_check.isChecked()
            ):
                if not RxnMapperProvider.is_available():
                    raise RuntimeError(
                        "RXNMapper is not installed. Run "
                        "'python -m pip install -r requirements-mapping.txt' "
                        "or clear the RXNMapper checkbox."
                    )
                if self._mapping_provider is None:
                    self._mapping_provider = RxnMapperProvider()
                mapping_provider = self._mapping_provider
            kind, analysis, assessment = featurize_text(
                input_text,
                mapping_provider=mapping_provider,
                force_resolved_mapping=(
                    requested_kind == "reaction"
                    and self.force_core_mapping_check.isChecked()
                ),
            )
            heading = f"{kind.upper()} FEATURIZATION"
            mapping_summary = ""
            if kind == "reaction":
                if assessment is None:
                    mapping_summary = "\nRXNMapper: disabled"
                else:
                    result = assessment.mapping_result
                    confidence = (
                        f"{result.mapper_confidence:.3f}"
                        if result is not None
                        and result.mapper_confidence is not None
                        else "not run"
                    )
                    mapping_summary = (
                        f"\nRXNMapper: {assessment.status}"
                        f"\nMapper provider: "
                        f"{assessment.provider_metadata.provider_id}"
                        f"\nMapper confidence: {confidence}"
                    )
            self.output.setPlainText(
                f"{heading}{mapping_summary}\n\n"
                f"{format_concise_analysis(analysis)}"
            )
            self._render_structure(
                kind,
                self.input_edit.text().strip(),
                analysis=analysis,
            )
            valid = bool(getattr(analysis, "valid", False))
            state = "valid" if valid else "invalid"
            self.status_label.setText(
                f"Complete · {kind} input · {state}"
                + (
                    " · RXNMapper on"
                    if kind == "reaction"
                    and self.use_rxnmapper_check.isChecked()
                    else ""
                )
            )
            self.copy_button.setEnabled(True)
        except Exception as exc:
            self.output.setPlainText(f"Unable to analyze input.\n\n{exc}")
            self.graph_heading.setText("Structure graph")
            self.structure_image_label.clear_image(
                "Structure graph unavailable."
            )
            self.core_image_label.clear_image(
                "Minimized reaction graphic unavailable."
            )
            self.core_graphic_note.setText(str(exc))
            self.status_label.setText("Analysis failed")
            self.copy_button.setEnabled(True)
        finally:
            self.analyze_button.setEnabled(True)

    def _render_structure(
        self,
        kind: InputKind,
        text: str,
        *,
        analysis: object,
    ) -> None:
        """Render the full graph and, when available, its minimized core."""
        self.graph_heading.setText(
            "Reaction graph" if kind == "reaction" else "Compound graph"
        )
        try:
            if kind == "reaction":
                drawing = render_reaction_image_bytes(
                    text,
                    size=REACTION_IMAGE_SIZE,
                    image_format="svg",
                )
            else:
                drawing = render_molecule_image_bytes(
                    text,
                    size=MOLECULE_IMAGE_SIZE,
                    image_format="svg",
                )
        except (RuntimeError, ValueError) as exc:
            self.structure_image_label.clear_image(
                f"{self.graph_heading.text()} unavailable."
            )
            self.structure_image_label.setToolTip(str(exc))
            return
        if not self.structure_image_label.set_image_bytes(drawing):
            self.structure_image_label.setToolTip(
                "The renderer returned an unsupported image."
            )
            return
        self.structure_image_label.setToolTip(text)
        if kind != "reaction":
            self.core_image_label.clear_image(
                "Reaction minimization applies only to reactions."
            )
            self.core_graphic_note.setText(
                "Enter a reaction SMILES to generate a minimized graphic."
            )
            self.graph_tabs.setTabEnabled(1, False)
            self.graph_tabs.setCurrentIndex(0)
            return
        self.graph_tabs.setTabEnabled(1, True)
        core = getattr(analysis, "reaction_core", None)
        if core is None:
            self.core_image_label.clear_image(
                "Mapped reaction core unavailable."
            )
            self.core_graphic_note.setText(
                "Supply atom mapping, use RXNMapper for an unresolved "
                "reaction, or enable resolved-reaction mapping."
            )
            self.graph_tabs.setCurrentIndex(0)
            return
        try:
            graphic = build_reaction_core_graphic(
                analysis,
                size=REACTION_IMAGE_SIZE,
                image_format="svg",
            )
        except (RuntimeError, ValueError) as exc:
            self.core_image_label.clear_image(
                "Unable to render minimized reaction."
            )
            self.core_graphic_note.setText(str(exc))
            self.graph_tabs.setCurrentIndex(0)
            return
        if not self.core_image_label.set_image_bytes(graphic.image_bytes):
            self.core_graphic_note.setText(
                "The minimized renderer returned an unsupported image."
            )
            self.graph_tabs.setCurrentIndex(0)
            return
        legend = "; ".join(
            f"{placeholder.label} = {placeholder.fragment_smiles}"
            for placeholder in graphic.placeholders
        ) or "no remote placeholders"
        note = (
            f"{graphic.evidence_status} evidence; "
            f"confidence {graphic.confidence:.3f}; {legend}"
        )
        if graphic.evidence_status == "external":
            note += "; external mapping requires expert review"
        self.core_graphic_note.setText(note)
        self.core_image_label.setToolTip(note)
        self.graph_tabs.setCurrentIndex(1)

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

"""Simple Qt6 interface for type-agnostic condition recommendation."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Tuple

from PyQt6 import QtCore, QtGui, QtWidgets

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DATA_FOLDER = PROJECT_ROOT / "datasets" / "literature"
DEFAULT_INDEX_PATH = DEFAULT_DATA_FOLDER / "generic_index.json.gz"
DEFAULT_MANIFEST_PATH = DEFAULT_DATA_FOLDER / "shard_manifest.json"

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from condition_recommender import (  # noqa: E402
    GenericConditionRecommender,
    GenericRecommendationResult,
)

_RECOMMENDER_CACHE: Dict[
    Tuple[str, int, int],
    GenericConditionRecommender,
] = {}

_RECIPE_ROLE_LABELS = (
    ("catalysts", "Catalyst"),
    ("ligands", "Ligand"),
    ("bases", "Base"),
    ("condensation_agents", "Condensation agent"),
    ("oxidants", "Oxidant"),
    ("reductants", "Reductant"),
    ("acids", "Acid"),
    ("additives", "Additive"),
    ("solvents", "Solvent"),
    ("other_components", "Other"),
)


def default_recommendation_data_path() -> Path:
    """Return the fastest available default recommendation artifact."""
    if DEFAULT_INDEX_PATH.is_file():
        return DEFAULT_INDEX_PATH
    return DEFAULT_MANIFEST_PATH


def _cache_key(path: Path) -> Tuple[str, int, int]:
    resolved = path.resolve()
    stat = resolved.stat()
    return str(resolved), stat.st_size, stat.st_mtime_ns


def _get_cached_recommender(path: str | Path) -> GenericConditionRecommender:
    """Load a validated index once and invalidate it when the file changes."""
    source = Path(path)
    key = _cache_key(source)
    recommender = _RECOMMENDER_CACHE.get(key)
    if recommender is not None:
        return recommender
    resolved = key[0]
    for old_key in tuple(_RECOMMENDER_CACHE):
        if old_key[0] == resolved:
            _RECOMMENDER_CACHE.pop(old_key, None)
    recommender = GenericConditionRecommender.from_path(source)
    _RECOMMENDER_CACHE[key] = recommender
    return recommender


def _component_name(component: Mapping[str, Any]) -> str:
    return str(
        component.get("canonical_name")
        or component.get("raw_identifier")
        or component.get("substance_id")
        or "Unknown"
    )


def _component_names(recipe: Mapping[str, Any], field: str) -> str:
    values = recipe.get(field) or ()
    names = [
        _component_name(component)
        for component in values
        if isinstance(component, Mapping)
    ]
    return ", ".join(names)


def format_recipe_summary(recipe: Mapping[str, Any]) -> str:
    """Render the most useful recipe fields in one compact table cell."""
    parts = []
    for field, label in _RECIPE_ROLE_LABELS:
        names = _component_names(recipe, field)
        if names:
            parts.append(f"{label}: {names}")
    temperature = recipe.get("temperature_c")
    time_h = recipe.get("time_h")
    if temperature is not None:
        parts.append(f"{temperature:g} °C")
    if time_h is not None:
        parts.append(f"{time_h:g} h")
    atmosphere = str(recipe.get("atmosphere") or "").strip()
    if atmosphere:
        parts.append(f"Atmosphere: {atmosphere}")
    return "; ".join(parts) or "No resolved components"


def _display_name(value: Any) -> str:
    text = str(value or "").strip()
    return text.replace("_", " ").title() if text else "Unassigned"


def _friendly_error(error: Any) -> str:
    code = str(error or "RECOMMENDATION_FAILED")
    messages = {
        "EMPTY_GENERIC_INDEX": "The selected recommendation index is empty.",
        "QUERY_HAS_NO_USABLE_REACTION_SIGNATURE": (
            "The molecular graphs did not provide enough evidence for a "
            "recommendable reaction transformation."
        ),
        "NO_CHEMICALLY_COMPATIBLE_PRECEDENT": (
            "No chemically compatible precedent was found in this dataset."
        ),
        "NO_COMPATIBLE_CONDITION_PRECEDENT": (
            "Matching reactions were found, but their condition recipes were "
            "not compatible with the query."
        ),
        "INCOMPATIBLE_REACTION_SIGNATURE_SCHEMA": (
            "The query and saved index use incompatible signature versions. "
            "Rebuild the converted data."
        ),
        "INCOMPATIBLE_REACTION_TAXONOMY_DEFINITIONS": (
            "The saved index uses older chemistry definitions. Rebuild it."
        ),
    }
    return messages.get(code, code.replace("_", " ").title())


class GenericRecommendationWorker(QtCore.QObject):
    """Load the index and recommend outside the Qt event loop."""

    progress = QtCore.pyqtSignal(str)
    finished = QtCore.pyqtSignal(bool, object, str)

    def __init__(
        self,
        data_path: str,
        reaction_smiles: str,
        *,
        top_k: int,
        minimum_pool_size: Optional[int],
    ) -> None:
        super().__init__()
        self.data_path = data_path
        self.reaction_smiles = reaction_smiles
        self.top_k = top_k
        self.minimum_pool_size = minimum_pool_size

    @QtCore.pyqtSlot()
    def run(self) -> None:
        """Load or reuse the index and execute one recommendation."""
        try:
            self.progress.emit("Loading recommendation index…")
            recommender = _get_cached_recommender(self.data_path)
            self.progress.emit("Analyzing reaction and ranking conditions…")
            result = recommender.recommend(
                self.reaction_smiles,
                top_k=self.top_k,
                minimum_pool_size=self.minimum_pool_size,
            )
        except Exception as exc:
            self.finished.emit(
                False,
                None,
                f"{type(exc).__name__}: {exc}",
            )
        else:
            self.finished.emit(True, result, "")


class GenericRecommenderWindow(QtWidgets.QWidget):
    """Simple desktop interface for the clean generic recommender."""

    def __init__(self) -> None:
        super().__init__()
        self.setFont(QtGui.QFont("Segoe UI", 9))
        self.setWindowTitle("Reaction Condition Recommender")
        self.resize(1180, 760)
        self.thread: Optional[QtCore.QThread] = None
        self.worker: Optional[GenericRecommendationWorker] = None
        self.last_result: Optional[GenericRecommendationResult] = None

        self.data_path_edit = QtWidgets.QLineEdit(
            str(default_recommendation_data_path())
        )
        self.data_path_edit.setObjectName("recommendationDataPath")
        self.data_path_edit.setPlaceholderText(
            "generic_index.json.gz or shard_manifest.json"
        )
        self.data_summary = QtWidgets.QLabel()
        self.data_summary.setObjectName("dataSummary")
        self.data_summary.setStyleSheet("color: #52606d;")

        self.reaction_edit = QtWidgets.QPlainTextEdit()
        self.reaction_edit.setObjectName("reactionSmiles")
        self.reaction_edit.setPlaceholderText(
            "Enter reaction SMILES with product, for example:\n"
            "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
        )
        self.reaction_edit.setMaximumHeight(90)

        self.top_k_spin = QtWidgets.QSpinBox()
        self.top_k_spin.setObjectName("topK")
        self.top_k_spin.setRange(1, 50)
        self.top_k_spin.setValue(5)
        self.minimum_pool_spin = QtWidgets.QSpinBox()
        self.minimum_pool_spin.setObjectName("minimumPoolSize")
        self.minimum_pool_spin.setRange(0, 100)
        self.minimum_pool_spin.setSpecialValueText("Default")
        self.minimum_pool_spin.setValue(0)
        self.minimum_pool_spin.setToolTip(
            "Default uses the versioned retrieval definition. A larger value "
            "may force a broader chemistry fallback."
        )

        self.run_button = QtWidgets.QPushButton("Recommend Conditions")
        self.run_button.setObjectName("recommendButton")
        self.run_button.setMinimumHeight(34)
        self.run_button.setStyleSheet(
            "QPushButton { background: #2b6cb0; color: white; "
            "font-weight: 600; border-radius: 4px; padding: 5px 12px; }"
            "QPushButton:hover { background: #2c5282; }"
            "QPushButton:disabled { background: #9bb7d6; }"
        )
        self.example_button = QtWidgets.QPushButton("Load Example")
        self.clear_button = QtWidgets.QPushButton("Clear")
        self.export_button = QtWidgets.QPushButton("Export JSON")
        self.export_button.setObjectName("exportButton")
        self.export_button.setEnabled(False)

        self.summary_box = QtWidgets.QPlainTextEdit()
        self.summary_box.setObjectName("recommendationSummary")
        self.summary_box.setReadOnly(True)
        self.summary_box.setMaximumHeight(130)

        self.results_table = QtWidgets.QTableWidget()
        self.results_table.setObjectName("recommendationTable")
        self.results_table.setColumnCount(9)
        self.results_table.setHorizontalHeaderLabels(
            [
                "Rank",
                "Score",
                "Similarity",
                "Compatibility",
                "Expected yield",
                "Rxn support",
                "Ref support",
                "Conditions",
                "Cautions",
            ]
        )
        self.results_table.setAlternatingRowColors(True)
        self.results_table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
        )
        self.results_table.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows
        )
        self.results_table.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.SingleSelection
        )
        self.results_table.setSortingEnabled(True)
        header = self.results_table.horizontalHeader()
        header.setSectionResizeMode(
            7,
            QtWidgets.QHeaderView.ResizeMode.Stretch,
        )

        self.details_box = QtWidgets.QPlainTextEdit()
        self.details_box.setObjectName("recommendationDetails")
        self.details_box.setReadOnly(True)

        self.progress_bar = QtWidgets.QProgressBar()
        self.progress_bar.setObjectName("recommendationProgress")
        self.progress_bar.setVisible(False)
        self.status_label = QtWidgets.QLabel("Ready")
        self.status_label.setObjectName("recommendationStatus")

        self._build_layout()
        self._bind_signals()
        self._update_data_summary()
        QtGui.QShortcut(
            QtGui.QKeySequence("Ctrl+Return"),
            self,
            activated=self.start_recommendation,
        )

    def _build_layout(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(18, 18, 18, 18)
        layout.setSpacing(10)

        title = QtWidgets.QLabel("Reaction Condition Recommender")
        title.setStyleSheet("font-size: 20px; font-weight: 600;")
        layout.addWidget(title)

        data_row = QtWidgets.QHBoxLayout()
        data_row.addWidget(self.data_path_edit)
        browse = QtWidgets.QPushButton("Browse…")
        browse.clicked.connect(self.choose_data_path)
        data_row.addWidget(browse)
        layout.addWidget(QtWidgets.QLabel("Recommendation data"))
        layout.addLayout(data_row)
        layout.addWidget(self.data_summary)

        form = QtWidgets.QFormLayout()
        form.addRow("Reaction SMILES:", self.reaction_edit)
        options = QtWidgets.QHBoxLayout()
        options.addWidget(QtWidgets.QLabel("Top results"))
        options.addWidget(self.top_k_spin)
        options.addSpacing(18)
        options.addWidget(QtWidgets.QLabel("Minimum precedent pool"))
        options.addWidget(self.minimum_pool_spin)
        options.addStretch()
        form.addRow("Options:", options)
        layout.addLayout(form)

        buttons = QtWidgets.QHBoxLayout()
        buttons.addWidget(self.run_button)
        buttons.addWidget(self.example_button)
        buttons.addWidget(self.clear_button)
        buttons.addWidget(self.export_button)
        buttons.addStretch()
        layout.addLayout(buttons)

        splitter = QtWidgets.QSplitter(
            QtCore.Qt.Orientation.Vertical,
        )
        upper = QtWidgets.QWidget()
        upper_layout = QtWidgets.QVBoxLayout(upper)
        upper_layout.setContentsMargins(0, 0, 0, 0)
        upper_layout.addWidget(QtWidgets.QLabel("Query summary"))
        upper_layout.addWidget(self.summary_box)
        upper_layout.addWidget(QtWidgets.QLabel("Recommended recipes"))
        upper_layout.addWidget(self.results_table)
        splitter.addWidget(upper)

        lower = QtWidgets.QWidget()
        lower_layout = QtWidgets.QVBoxLayout(lower)
        lower_layout.setContentsMargins(0, 0, 0, 0)
        lower_layout.addWidget(QtWidgets.QLabel("Selected recipe details"))
        lower_layout.addWidget(self.details_box)
        splitter.addWidget(lower)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)
        layout.addWidget(splitter, stretch=1)

        status_row = QtWidgets.QHBoxLayout()
        self.progress_bar.setMaximumWidth(220)
        status_row.addWidget(self.progress_bar)
        status_row.addWidget(self.status_label)
        status_row.addStretch()
        layout.addLayout(status_row)

    def _bind_signals(self) -> None:
        self.data_path_edit.editingFinished.connect(self._update_data_summary)
        self.run_button.clicked.connect(self.start_recommendation)
        self.example_button.clicked.connect(self.load_example)
        self.clear_button.clicked.connect(self.clear_results)
        self.export_button.clicked.connect(self.export_json)
        self.results_table.itemSelectionChanged.connect(
            self._show_selected_details
        )

    @QtCore.pyqtSlot()
    def choose_data_path(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Choose recommendation index or manifest",
            self.data_path_edit.text() or str(DEFAULT_DATA_FOLDER),
            (
                "Recommendation data (*.json.gz *.json);;"
                "All files (*)"
            ),
        )
        if path:
            self.data_path_edit.setText(path)
            self._update_data_summary()

    @QtCore.pyqtSlot()
    def _update_data_summary(self) -> None:
        path = Path(self.data_path_edit.text().strip())
        if not path.is_file():
            self.data_summary.setText("Data file not found.")
            return
        report_path = path.parent / "recommendation_artifacts_report.json"
        if report_path.is_file():
            try:
                report = json.loads(report_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                report = {}
            total = report.get("record_count")
            eligible = report.get("eligible_index_record_count")
            if total is not None:
                eligible_text = (
                    f"; {eligible} recommendation-eligible"
                    if eligible is not None
                    else ""
                )
                self.data_summary.setText(
                    f"{total} converted reactions{eligible_text}."
                )
                return
        self.data_summary.setText(
            f"Using {path.name} ({path.stat().st_size / 1024:.1f} KB)."
        )

    @QtCore.pyqtSlot()
    def load_example(self) -> None:
        self.reaction_edit.setPlainText(
            "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
        )

    @QtCore.pyqtSlot()
    def clear_results(self) -> None:
        self.summary_box.clear()
        self.details_box.clear()
        self.results_table.setSortingEnabled(False)
        self.results_table.setRowCount(0)
        self.results_table.setSortingEnabled(True)
        self.last_result = None
        self.export_button.setEnabled(False)
        self.status_label.setText("Ready")

    @QtCore.pyqtSlot()
    def start_recommendation(self) -> None:
        if self.thread is not None:
            return
        data_path = Path(self.data_path_edit.text().strip())
        reaction_smiles = self.reaction_edit.toPlainText().strip()
        if not data_path.is_file():
            QtWidgets.QMessageBox.warning(
                self,
                "Recommendation data required",
                "Choose an existing generic index or shard manifest.",
            )
            return
        if not reaction_smiles:
            QtWidgets.QMessageBox.warning(
                self,
                "Reaction required",
                "Enter a reaction SMILES including its product.",
            )
            return
        if ">>" not in reaction_smiles:
            QtWidgets.QMessageBox.warning(
                self,
                "Product required",
                "Use reaction SMILES in reactants>>product form.",
            )
            return

        self.clear_results()
        self.run_button.setEnabled(False)
        self.status_label.setText("Starting…")
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setVisible(True)

        minimum = self.minimum_pool_spin.value()
        thread = QtCore.QThread(self)
        worker = GenericRecommendationWorker(
            str(data_path),
            reaction_smiles,
            top_k=self.top_k_spin.value(),
            minimum_pool_size=minimum or None,
        )
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.progress.connect(self.status_label.setText)
        worker.finished.connect(self._on_finished)
        worker.finished.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(self._clear_thread)
        self.thread = thread
        self.worker = worker
        thread.start()

    @QtCore.pyqtSlot(bool, object, str)
    def _on_finished(
        self,
        success: bool,
        result: Optional[GenericRecommendationResult],
        error: str,
    ) -> None:
        self.run_button.setEnabled(True)
        self.progress_bar.setVisible(False)
        if not success or result is None:
            self.status_label.setText("Error")
            self.summary_box.setPlainText(error)
            return
        self.last_result = result
        self.export_button.setEnabled(True)
        self._render_result(result)

    def _render_result(self, result: GenericRecommendationResult) -> None:
        if not result.valid:
            self.status_label.setText("No recommendation")
            self.summary_box.setPlainText(_friendly_error(result.error))
            self.results_table.setRowCount(0)
            return
        warnings = ", ".join(result.warnings) if result.warnings else "None"
        self.summary_box.setPlainText(
            "\n".join(
                (
                    f"Detected family: {_display_name(result.named_family)}",
                    (
                        "Transformation: "
                        f"{_display_name(result.transformation_class)}"
                    ),
                    (
                        "Retrieval level: "
                        f"{_display_name(result.retrieval_level)}"
                    ),
                    (
                        f"Candidates: {result.candidate_count}; "
                        f"compatible: {result.compatible_candidate_count}; "
                        f"excluded: {result.excluded_candidate_count}"
                    ),
                    (
                        "Independent compatible precedents: "
                        f"{result.independent_compatible_candidate_count}"
                    ),
                    f"Warnings: {warnings}",
                )
            )
        )
        recommendations = tuple(result.recommendations)
        self.results_table.setSortingEnabled(False)
        self.results_table.setRowCount(len(recommendations))
        for row, recommendation in enumerate(recommendations):
            values = (
                str(recommendation.rank),
                f"{recommendation.score:.3f}",
                f"{recommendation.similarity_score:.3f}",
                f"{recommendation.compatibility_score:.3f}",
                (
                    f"{recommendation.expected_yield_pct:.1f}%"
                    if recommendation.expected_yield_pct is not None
                    else "—"
                ),
                str(recommendation.support),
                str(recommendation.reference_support),
                format_recipe_summary(recommendation.resolved_recipe),
                "; ".join(recommendation.cautions),
            )
            for column, value in enumerate(values):
                item = QtWidgets.QTableWidgetItem(value)
                if column == 0:
                    item.setData(
                        QtCore.Qt.ItemDataRole.UserRole,
                        recommendation,
                    )
                self.results_table.setItem(row, column, item)
        self.results_table.setSortingEnabled(True)
        self.results_table.resizeColumnsToContents()
        self.results_table.horizontalHeader().setSectionResizeMode(
            7,
            QtWidgets.QHeaderView.ResizeMode.Stretch,
        )
        self.status_label.setText(
            f"Done — {len(recommendations)} recipe(s)"
        )
        if recommendations:
            self.results_table.selectRow(0)

    @QtCore.pyqtSlot()
    def _show_selected_details(self) -> None:
        row = self.results_table.currentRow()
        if row < 0:
            self.details_box.clear()
            return
        item = self.results_table.item(row, 0)
        if item is None:
            return
        recommendation = item.data(QtCore.Qt.ItemDataRole.UserRole)
        if recommendation is None:
            return
        recipe = recommendation.resolved_recipe
        lines = [
            f"Rank {recommendation.rank}",
            f"Recipe core: {recommendation.recipe_core_id}",
            (
                "Observed recipe variants: "
                f"{', '.join(recommendation.recipe_variant_ids)}"
            ),
            "",
            "Conditions",
        ]
        for field, label in _RECIPE_ROLE_LABELS:
            names = _component_names(recipe, field)
            if names:
                lines.append(f"{label}: {names}")
        for field, label, unit in (
            ("temperature_c", "Temperature", "°C"),
            ("time_h", "Time", "h"),
            ("concentration_m", "Concentration", "M"),
            ("pressure_bar", "Pressure", "bar"),
        ):
            value = recipe.get(field)
            if value is not None:
                lines.append(f"{label}: {value:g} {unit}")
        atmosphere = str(recipe.get("atmosphere") or "").strip()
        if atmosphere:
            lines.append(f"Atmosphere: {atmosphere}")
        lines.extend(
            (
                "",
                "Evidence",
                (
                    f"Score: {recommendation.score:.3f}; similarity: "
                    f"{recommendation.similarity_score:.3f}; compatibility: "
                    f"{recommendation.compatibility_score:.3f}"
                ),
                (
                    f"Support: {recommendation.support} reaction(s), "
                    f"{recommendation.reference_support} reference(s), "
                    f"{recommendation.dataset_support} dataset(s)"
                ),
            )
        )
        if recommendation.expected_yield_pct is not None:
            lines.append(
                f"Evidence-weighted expected yield: "
                f"{recommendation.expected_yield_pct:.1f}%"
            )
        if recommendation.explanation:
            lines.extend(("", "Why this recipe"))
            lines.extend(f"• {value}" for value in recommendation.explanation)
        if recommendation.compatibility_evidence:
            lines.extend(("", "Compatibility"))
            lines.extend(
                f"• {value}"
                for value in recommendation.compatibility_evidence
            )
        if recommendation.cautions:
            lines.extend(("", "Cautions"))
            lines.extend(f"• {value}" for value in recommendation.cautions)
        if recommendation.precedent_reference_ids:
            lines.extend(
                (
                    "",
                    "Precedent references",
                    *recommendation.precedent_reference_ids,
                )
            )
        if recommendation.precedent_reaction_ids:
            lines.extend(
                (
                    "",
                    "Precedent reactions",
                    *recommendation.precedent_reaction_ids,
                )
            )
        self.details_box.setPlainText("\n".join(lines))

    @QtCore.pyqtSlot()
    def export_json(self) -> None:
        if self.last_result is None:
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Export recommendation result",
            str(PROJECT_ROOT / "results" / "generic_recommendation.json"),
            "JSON files (*.json)",
        )
        if not path:
            return
        destination = Path(path)
        if destination.suffix.casefold() != ".json":
            destination = destination.with_suffix(".json")
        try:
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text(
                json.dumps(
                    self.last_result.to_dict(),
                    ensure_ascii=False,
                    indent=2,
                )
                + "\n",
                encoding="utf-8",
            )
        except OSError as exc:
            QtWidgets.QMessageBox.critical(
                self,
                "Export failed",
                str(exc),
            )
            return
        self.status_label.setText(f"Exported {destination}")

    @QtCore.pyqtSlot()
    def _clear_thread(self) -> None:
        self.thread = None
        self.worker = None

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        if self.thread is not None and self.thread.isRunning():
            QtWidgets.QMessageBox.information(
                self,
                "Recommendation running",
                "Wait for the current recommendation to finish before closing.",
            )
            event.ignore()
            return
        event.accept()


def main() -> None:
    """Launch the generic condition recommender."""
    application = QtWidgets.QApplication(sys.argv)
    window = GenericRecommenderWindow()
    window.show()
    raise SystemExit(application.exec())


if __name__ == "__main__":
    main()


__all__ = [
    "DEFAULT_INDEX_PATH",
    "GenericRecommendationWorker",
    "GenericRecommenderWindow",
    "default_recommendation_data_path",
    "format_recipe_summary",
    "main",
]

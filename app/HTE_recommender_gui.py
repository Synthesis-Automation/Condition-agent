import csv
import json
import html
import math
import os
import sys
import tempfile
import time
from dataclasses import asdict
from pathlib import Path
from typing import Optional, Tuple, Dict, Any, List

from PyQt6 import QtCore, QtGui, QtWidgets


PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def _parse_reaction_smiles(reaction_smiles: str) -> Tuple[str, Optional[str], Optional[str]]:
    text = (reaction_smiles or "").strip()
    if not text:
        return "", None, None
    if ">>" in text:
        reactants_part, product = text.split(">>", 1)
        reactants = [r for r in reactants_part.split(".") if r]
        reactant_a = reactants[0] if reactants else ""
        reactant_b = ".".join(reactants[1:]) if len(reactants) > 1 else None
        product_smiles = product.strip() or None
        return reactant_a, reactant_b, product_smiles
    reactants = [r for r in text.split(".") if r]
    reactant_a = reactants[0] if reactants else ""
    reactant_b = ".".join(reactants[1:]) if len(reactants) > 1 else None
    return reactant_a, reactant_b, None


def _detect_csv_type(path: Path) -> str:
    parts = [part.lower() for part in path.parts]
    for label in ("rules", "literature", "datasets", "experiments", "experiment", "experiements"):
        if label in parts:
            if label in ("literature", "datasets"):
                return "literature"
            if label in ("experiments", "experiment", "experiements"):
                return "experiments"
            return label
    if path.is_dir():
        return "directory"
    try:
        with path.open("r", encoding="utf-8") as handle:
            reader = csv.reader(handle)
            header = next(reader, [])
    except Exception:
        return "unknown"
    header_lower = [h.strip().lower() for h in header]
    if "reaction_type_standardized" in header_lower:
        return "literature"
    if "reaction_id" in header_lower:
        return "literature"
    if "temperature_c" in header_lower:
        return "rules"
    if "reactant_1" in header_lower and "reactant_2" in header_lower:
        return "experiments"
    return "unknown"


def _normalize_source_group_label(value: Any) -> str:
    text = str(value or "").strip().lower()
    if not text or text == "nan":
        return "unknown"
    if text in ("literature", "datasets", "dataset", "lit"):
        return "literature"
    if text in ("experiments", "experiment", "experiements"):
        return "experiments"
    if text == "rules":
        return "rules"
    return text


def _format_float(value: Optional[float]) -> str:
    if value is None:
        return ""
    try:
        return f"{float(value):.2f}".rstrip("0").rstrip(".")
    except (TypeError, ValueError):
        return str(value)


def _safe_text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return str(value)


def _collect_reaction_spectator_groups(reaction_smiles: str) -> List[str]:
    if not reaction_smiles:
        return []
    try:
        from chemtools.featurizers.unified import featurize_reaction
    except Exception:
        return []

    try:
        payload = featurize_reaction(reaction_smiles)
    except Exception:
        return []

    if not isinstance(payload, dict):
        return []
    reaction = payload.get("reaction")
    if not isinstance(reaction, dict):
        return []
    aggregates = reaction.get("aggregates") or {}
    groups = aggregates.get("spectator_groups_ranked") or aggregates.get("spectator_groups_combined") or []
    cleaned = [str(group).strip() for group in groups if str(group).strip()]
    if os.environ.get("HTE_DEBUG_SPECTATORS") == "1":
        print(f"[HTE_DEBUG] reaction_smiles={reaction_smiles}")
        print(f"[HTE_DEBUG] spectator_groups_combined={cleaned}")
    return cleaned


def _reaction_image_path(prefix: str) -> Path:
    output_dir = PROJECT_ROOT / "results" / "visualization"
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir / f"{prefix}_{time.time_ns()}.png"


def _scale_pixmap_half(pixmap: QtGui.QPixmap) -> QtGui.QPixmap:
    if pixmap.isNull():
        return pixmap
    new_width = max(1, int(pixmap.width() * 0.5))
    new_height = max(1, int(pixmap.height() * 0.5))
    return pixmap.scaled(
        new_width,
        new_height,
        QtCore.Qt.AspectRatioMode.KeepAspectRatio,
        QtCore.Qt.TransformationMode.SmoothTransformation,
    )


def _format_nearby_groups(groups: List[str]) -> str:
    return ", ".join(groups) if groups else "None"


def _table_columns_for_type(data_type: str) -> List[Tuple[str, str]]:
    if data_type == "precedent":
        return [
            ("Rank", "rank"),
            ("Similarity", "match_score"),
            ("Yield", "avg_yield"),
            ("Catalyst", "catalyst"),
            ("Ligand", "ligand"),
            ("Base", "base"),
            ("Solvent", "solvent"),
            ("Additive", "additive"),
            ("Reaction Type", "reaction_type"),
            ("Reaction ID", "reaction_id"),
        ]
    base = [
        ("Rank", "rank"),
        ("Avg Z-Score", "avg_z_score"),
        ("Confidence", "confidence_score"),
        ("Success %", "success_rate"),
        ("Avg Yield", "avg_yield"),
        ("Experiments", "num_experiments"),
        ("Catalyst", "catalyst"),
        ("Ligand", "ligand"),
        ("Base", "base"),
        ("Solvent", "solvent"),
        ("Additive", "additive"),
        ("Condensation Agent", "coupling_reagent"),
        ("Spectator Groups", "spectator_groups"),
    ]
    columns = base + [
        ("Reaction Type", "reaction_type"),
        ("Reaction ID", "reaction_id"),
        ("Reactant Types", "reactant_types"),
        ("Match Score", "match_score"),
    ]
    return columns


class RecommendationWorker(QtCore.QObject):
    finished = QtCore.pyqtSignal(bool, object, str, dict)
    progress = QtCore.pyqtSignal(str)

    def __init__(
        self,
        db_path: str,
        reaction_smiles: str,
        top_k: int,
        min_exp: int,
        reaction_filter: str,
        catalyst_filter: str,
        source_group: str,
        use_aryl_weighting: bool,
    ) -> None:
        super().__init__()
        self.db_path = db_path
        self.reaction_smiles = reaction_smiles
        self.top_k = top_k
        self.min_exp = min_exp
        self.reaction_filter = reaction_filter
        self.catalyst_filter = catalyst_filter
        self.source_group = source_group
        self.use_aryl_weighting = use_aryl_weighting

    def run(self) -> None:
        try:
            from chemtools.recommend import HTERecommender

            reactant_a, reactant_b, product = _parse_reaction_smiles(self.reaction_smiles)
            if not reactant_a:
                raise ValueError("Provide a reaction SMILES with at least one reactant.")

            self.progress.emit(f"Loading HTE data from {self.db_path} ...")
            recommender = HTERecommender(self.db_path)
            self.progress.emit("Running recommendation ...")

            result = recommender.recommend(
                reactant_a_smiles=reactant_a,
                reactant_b_smiles=reactant_b,
                product_smiles=product,
                top_k=self.top_k,
                min_experiments=self.min_exp,
                reaction_type_filter=self.reaction_filter or None,
                catalyst_filter=self.catalyst_filter or None,
                source_group=self.source_group or None,
                use_aryl_steric_electronic_weighting=self.use_aryl_weighting,
            )
            stats = recommender.get_statistics()
            self.finished.emit(True, result, "OK", stats)
        except Exception as exc:
            self.finished.emit(False, None, str(exc), {})


class HTERecommenderWindow(QtWidgets.QWidget):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("HTE Recommender (Qt6)")
        self.resize(1300, 850)

        self.db_path_edit = QtWidgets.QLineEdit()
        self.db_path_edit.setPlaceholderText("Select a CSV/JSONL or a folder like data/HTE_db")
        self.data_type_label = QtWidgets.QLabel("Data type: unknown")
        self.data_type_label.setStyleSheet("color: #555;")

        self.reaction_smiles_edit = QtWidgets.QLineEdit()
        self.reaction_smiles_edit.setPlaceholderText("A.B>>P or A.B or A")

        self.top_k_spin = QtWidgets.QSpinBox()
        self.top_k_spin.setRange(1, 200)
        self.top_k_spin.setValue(10)

        self.min_exp_spin = QtWidgets.QSpinBox()
        self.min_exp_spin.setRange(1, 200)
        self.min_exp_spin.setValue(1)

        self.reaction_filter_edit = QtWidgets.QLineEdit()
        self.reaction_filter_edit.setPlaceholderText("Optional reaction type/category filter (e.g., Suzuki)")

        self.catalyst_filter_edit = QtWidgets.QLineEdit()
        self.catalyst_filter_edit.setPlaceholderText("Optional catalyst filter (e.g., Pd, Cu)")

        self.source_group_combo = QtWidgets.QComboBox()
        self.source_group_combo.addItems(["All", "literature", "experiments", "rules"])
        self.source_group_combo.setToolTip(
            "Filter results to a specific HTE source group (rules live in data/HTE_db/rules)."
        )

        self.aryl_weighting_check = QtWidgets.QCheckBox("Aryl steric/electronic weighting")
        self.aryl_weighting_check.setToolTip("Reweight matches by aryl steric/electronic similarity (when available).")

        self.run_button = QtWidgets.QPushButton("Run Recommendation")
        self.clear_button = QtWidgets.QPushButton("Clear Results")
        self.export_json_button = QtWidgets.QPushButton("Export JSON")
        self.run_button.setMinimumSize(160, 34)
        self.clear_button.setMinimumSize(140, 34)
        self.export_json_button.setMinimumSize(140, 34)
        self.export_json_button.setEnabled(False)
        self.run_button.setStyleSheet(
            "QPushButton {"
            " background-color: #2b6cb0;"
            " color: white;"
            " font-weight: 600;"
            " border-radius: 4px;"
            " padding: 4px 10px;"
            "}"
            "QPushButton:hover { background-color: #2c5282; }"
            "QPushButton:pressed { background-color: #2a4365; }"
            "QPushButton:disabled { background-color: #9bb7d6; color: #f2f2f2; }"
        )

        self.summary = QtWidgets.QTextEdit()
        self.summary.setReadOnly(True)

        self.results_tabs = QtWidgets.QTabWidget()
        self._initialize_result_tabs()

        self.status = QtWidgets.QLabel("")
        self.progress_bar = QtWidgets.QProgressBar()
        self.progress_bar.setVisible(False)
        self.progress_bar.setMinimumWidth(160)

        self._setup_layout()
        self._bind_signals()
        self._set_default_path()

        self.thread: Optional[QtCore.QThread] = None
        self.worker: Optional[RecommendationWorker] = None
        self._reaction_dialog: Optional[QtWidgets.QDialog] = None
        self._spectator_groups_summary: str = ""
        self._source_group_label: str = "All"
        self._aryl_weighting_enabled: bool = False
        self._all_json_output: Optional[Dict[str, Any]] = None

    def _setup_layout(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)

        title = QtWidgets.QLabel("HTE Recommender")
        title.setStyleSheet("font-size: 18px; font-weight: bold;")
        title.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title)

        db_row = QtWidgets.QHBoxLayout()
        db_row.addWidget(self.db_path_edit)
        browse_file_btn = QtWidgets.QPushButton("Select CSV")
        browse_dir_btn = QtWidgets.QPushButton("Select Folder")
        db_row.addWidget(browse_file_btn)
        db_row.addWidget(browse_dir_btn)

        browse_file_btn.clicked.connect(self._choose_file)
        browse_dir_btn.clicked.connect(self._choose_dir)

        layout.addLayout(db_row)
        layout.addWidget(self.data_type_label)

        form = QtWidgets.QFormLayout()
        form.addRow("Reaction SMILES:", self.reaction_smiles_edit)

        options_row = QtWidgets.QHBoxLayout()
        options_row.addWidget(QtWidgets.QLabel("Top K:"))
        options_row.addWidget(self.top_k_spin)
        options_row.addSpacing(20)
        options_row.addWidget(QtWidgets.QLabel("Min Experiments:"))
        options_row.addWidget(self.min_exp_spin)
        options_row.addStretch()
        form.addRow("Options:", options_row)

        filters_row = QtWidgets.QHBoxLayout()
        filters_row.addWidget(QtWidgets.QLabel("Reaction:"))
        filters_row.addWidget(self.reaction_filter_edit)
        filters_row.addSpacing(12)
        filters_row.addWidget(QtWidgets.QLabel("Catalyst:"))
        filters_row.addWidget(self.catalyst_filter_edit)
        filters_row.addSpacing(12)
        filters_row.addWidget(QtWidgets.QLabel("Source:"))
        filters_row.addWidget(self.source_group_combo)
        filters_row.addStretch()
        form.addRow("Filters:", filters_row)
        form.addRow("Weighting:", self.aryl_weighting_check)
        layout.addLayout(form)

        button_row = QtWidgets.QHBoxLayout()
        button_row.addWidget(self.run_button)
        button_row.addWidget(self.clear_button)
        button_row.addWidget(self.export_json_button)
        button_row.addStretch()
        layout.addLayout(button_row)

        layout.addWidget(QtWidgets.QLabel("Summary:"))
        layout.addWidget(self.summary, stretch=1)
        layout.addWidget(QtWidgets.QLabel("Recommendations:"))
        layout.addWidget(self.results_tabs, stretch=3)
        status_row = QtWidgets.QHBoxLayout()
        status_row.addWidget(self.progress_bar)
        status_row.addWidget(self.status)
        status_row.addStretch()
        layout.addLayout(status_row)

    def _bind_signals(self) -> None:
        self.run_button.clicked.connect(self._run_recommendation)
        self.clear_button.clicked.connect(self._clear_results)
        self.export_json_button.clicked.connect(self._export_json_output)
        self.db_path_edit.textChanged.connect(self._update_data_type)

    def _set_default_path(self) -> None:
        default_path = PROJECT_ROOT / "data" / "HTE_db"
        if default_path.exists():
            self.db_path_edit.setText(str(default_path))

    def _choose_file(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Select HTE CSV/JSONL",
            str(PROJECT_ROOT),
            "HTE Files (*.csv *.jsonl);;All Files (*)",
        )
        if path:
            self.db_path_edit.setText(path)

    def _choose_dir(self) -> None:
        path = QtWidgets.QFileDialog.getExistingDirectory(
            self,
            "Select HTE Folder",
            str(PROJECT_ROOT),
        )
        if path:
            self.db_path_edit.setText(path)

    def _update_data_type(self) -> None:
        path_text = self.db_path_edit.text().strip()
        if not path_text:
            self.data_type_label.setText("Data type: unknown")
            return
        path = Path(path_text)
        if not path.exists():
            self.data_type_label.setText("Data type: missing path")
            return
        data_type = _detect_csv_type(path)
        self.data_type_label.setText(f"Data type: {data_type}")

    def _create_results_table(self) -> QtWidgets.QTableWidget:
        table = QtWidgets.QTableWidget(0, 0)
        table.setSortingEnabled(True)
        table.horizontalHeader().setStretchLastSection(True)
        table.setAlternatingRowColors(True)
        return table

    def _initialize_result_tabs(self) -> None:
        self.results_tabs.clear()
        self.table = self._create_results_table()
        self.results_tabs.addTab(self.table, "All")
        for label in ("Literature", "Rules", "Experiments", "Precedent"):
            group_table = self._create_results_table()
            self.results_tabs.addTab(group_table, label)

    def _clear_results(self) -> None:
        self.summary.clear()
        self._initialize_result_tabs()
        self.status.setText("")
        self._spectator_groups_summary = ""
        self._all_json_output = None
        self.export_json_button.setEnabled(False)
        self.run_button.setEnabled(True)
        self.run_button.setText("Run Recommendation")
        self.progress_bar.setVisible(False)
        if self._reaction_dialog:
            self._reaction_dialog.close()
            self._reaction_dialog = None

    def _run_recommendation(self) -> None:
        db_path = self.db_path_edit.text().strip()
        reaction_smiles = self.reaction_smiles_edit.text().strip()
        if not db_path:
            QtWidgets.QMessageBox.warning(self, "Missing data", "Select a CSV/JSONL or folder.")
            return
        if not reaction_smiles:
            QtWidgets.QMessageBox.warning(self, "Missing input", "Provide a reaction SMILES.")
            return
        if not Path(db_path).exists():
            QtWidgets.QMessageBox.warning(self, "Invalid path", "The data path does not exist.")
            return

        self._show_reaction_image(reaction_smiles)
        self.summary.clear()
        self.status.setText("Working...")
        self.run_button.setEnabled(False)
        self.run_button.setText("Running...")
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setVisible(True)
        self._source_group_label = self.source_group_combo.currentText().strip()
        self._aryl_weighting_enabled = bool(self.aryl_weighting_check.isChecked())

        self.thread = QtCore.QThread()
        self.worker = RecommendationWorker(
            db_path=db_path,
            reaction_smiles=reaction_smiles,
            top_k=self.top_k_spin.value(),
            min_exp=self.min_exp_spin.value(),
            reaction_filter=self.reaction_filter_edit.text().strip(),
            catalyst_filter=self.catalyst_filter_edit.text().strip(),
            source_group="" if self._source_group_label == "All" else self._source_group_label,
            use_aryl_weighting=self._aryl_weighting_enabled,
        )
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.progress.connect(self._append_status)
        self.worker.finished.connect(self._on_finished)
        self.thread.start()

    def _append_status(self, message: str) -> None:
        self.status.setText(message)

    def _on_finished(self, success: bool, result: object, message: str, stats: Dict[str, Any]) -> None:
        if self.thread:
            self.thread.quit()
            self.thread.wait()
        self.run_button.setEnabled(True)
        self.run_button.setText("Run Recommendation")
        self.progress_bar.setVisible(False)

        if not success:
            self.status.setText("Error")
            QtWidgets.QMessageBox.critical(self, "HTE Recommender", message)
            return

        self.status.setText("Done")
        self._render_result(result, stats)

    def _populate_table(
        self,
        table: QtWidgets.QTableWidget,
        recs: List[Any],
        columns: List[Tuple[str, str]],
    ) -> None:
        table.setColumnCount(len(columns))
        table.setHorizontalHeaderLabels([label for label, _ in columns])
        table.setRowCount(len(recs))

        for row_index, rec in enumerate(recs):
            rec_dict = asdict(rec)
            rec_dict["rank"] = row_index + 1
            reactant_types = []
            for item in list(getattr(rec, "reactant_types", []) or []):
                text = _safe_text(item).strip()
                if text:
                    reactant_types.append(text)
            rec_dict["reactant_types"] = " + ".join(reactant_types)
            for col_index, (_, key) in enumerate(columns):
                value = rec_dict.get(key)
                if isinstance(value, float):
                    cell_text = _format_float(value)
                else:
                    cell_text = str(value) if value is not None else ""
                item = QtWidgets.QTableWidgetItem(cell_text)
                table.setItem(row_index, col_index, item)

        table.resizeColumnsToContents()

    def _render_result(self, result: object, stats: Dict[str, Any]) -> None:
        data_type = _detect_csv_type(Path(self.db_path_edit.text().strip()))
        reaction_smiles = self.reaction_smiles_edit.text().strip()

        reactant_a = getattr(result, "reactant_a_smiles", "")
        reactant_b = getattr(result, "reactant_b_smiles", "") or "None"
        product = getattr(result, "product_smiles", "") or "None"
        type_a = getattr(result, "reactant_a_type", "")
        cat_a = getattr(result, "reactant_a_category", "")
        type_b = getattr(result, "reactant_b_type", "")
        cat_b = getattr(result, "reactant_b_category", "")
        predicted = getattr(result, "predicted_reaction_type", "") or "Unknown"
        confidence = getattr(result, "reaction_type_confidence", 0.0) * 100
        detected = getattr(result, "matched_motifs", None)
        detected_label = ""
        if isinstance(detected, tuple):
            detected_label = ", ".join([item for item in detected if item])
        matches = getattr(result, "total_matching_experiments", 0)
        coverage = getattr(result, "database_coverage", 0.0)
        spectator_summary = _format_nearby_groups(
            _collect_reaction_spectator_groups(reaction_smiles)
        )
        self._spectator_groups_summary = spectator_summary

        source_group = self._source_group_label or "All"
        weighting_label = "on" if self._aryl_weighting_enabled else "off"
        summary_lines = [
            html.escape(
                f"Data: {data_type} | Source: {source_group} | Aryl weighting: {weighting_label} | Matches: {matches} | Coverage: {coverage:.2f}%"
            ),
            (
                f'<span style="color:#ffeb3b; font-weight:700;">'
                f'{html.escape(f"MATCHED: {detected_label or "None"} | PRED: {predicted} ({confidence:.1f}%)")}'
                f"</span>"
            ),
            html.escape(f"A: {reactant_a} | Type: {type_a} ({cat_a})"),
            html.escape(f"B: {reactant_b} | Type: {type_b} ({cat_b})"),
            html.escape(f"Spectator Groups (All): {spectator_summary}"),
            html.escape(f"P: {product}"),
        ]
        if stats:
            stats_line = " | ".join(f"{key}: {stats[key]}" for key in sorted(stats))
            summary_lines.append(html.escape(f"Stats: {stats_line}"))
        self.summary.setHtml("<br>".join(summary_lines))

        recs = list(getattr(result, "recommendations", []) or [])
        source_map = getattr(result, "recommendations_by_source", {}) or {}
        normalized_map: Dict[str, List[Any]] = {}
        for key, items in source_map.items():
            normalized_key = _normalize_source_group_label(key)
            normalized_map.setdefault(normalized_key, []).extend(items)
        source_map = normalized_map

        try:
            from chemtools.formatters import format_hte_output
            reaction_filter = self.reaction_filter_edit.text().strip() or None
            catalyst_filter = self.catalyst_filter_edit.text().strip() or None
            literature_recs = list(source_map.get("literature") or [])
            if not literature_recs:
                literature_recs = list(recs)
            self._all_json_output = format_hte_output(
                result,
                recommendations=literature_recs,
                reaction_smiles=reaction_smiles,
                reaction_type_filter=reaction_filter,
                catalyst_filter=catalyst_filter,
                explanation=None,
            )
        except Exception:
            self._all_json_output = None
        self.export_json_button.setEnabled(bool(self._all_json_output))
        self.results_tabs.clear()

        self.table = self._create_results_table()
        self.results_tabs.addTab(self.table, "All")
        all_columns = _table_columns_for_type(data_type)
        self._populate_table(self.table, recs, all_columns)

        base_groups = ["literature", "rules", "experiments", "precedent"]
        extra_groups = [g for g in sorted(source_map) if g not in base_groups]
        for source_group in base_groups + extra_groups:
            group_recs = list(source_map.get(source_group) or [])
            group_table = self._create_results_table()
            group_columns = _table_columns_for_type(source_group)
            self._populate_table(group_table, group_recs, group_columns)
            label = (source_group or "unknown").replace("_", " ").title()
            self.results_tabs.addTab(group_table, label)

    def _export_json_output(self) -> None:
        if not self._all_json_output:
            QtWidgets.QMessageBox.information(self, "Export JSON", "Run a recommendation first.")
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Export JSON",
            str(PROJECT_ROOT),
            "JSON Files (*.json);;All Files (*)",
        )
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8") as handle:
                json.dump(self._all_json_output, handle, indent=2, ensure_ascii=False)
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "Export JSON", f"Failed to save JSON: {exc}")

    def _show_reaction_image(self, reaction_smiles: str) -> None:
        if not reaction_smiles:
            return
        if self._reaction_dialog:
            try:
                self._reaction_dialog.close()
            except RuntimeError:
                pass
            self._reaction_dialog = None

        dialog = QtWidgets.QDialog(self)
        dialog.setWindowTitle("Reaction Preview")
        dialog.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose)
        dialog.destroyed.connect(lambda: setattr(self, "_reaction_dialog", None))

        layout = QtWidgets.QVBoxLayout(dialog)
        info_label = QtWidgets.QLabel("")
        info_label.setWordWrap(True)
        layout.addWidget(info_label)

        images_widget = QtWidgets.QWidget()
        images_layout = QtWidgets.QHBoxLayout(images_widget)
        images_layout.setContentsMargins(0, 0, 0, 0)
        images_layout.setSpacing(12)
        layout.addWidget(images_widget)

        try:
            from chemtools.visualization import render_molecule_image, render_reaction_image
        except Exception as exc:
            info_label.setText(f"Visualization unavailable: {exc}")
            dialog.show()
            self._reaction_dialog = dialog
            return

        reactant_a, reactant_b, product = _parse_reaction_smiles(reaction_smiles)
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            if product:
                info_label.setText("")
                output_path = temp_path / "reaction.png"
                try:
                    render_reaction_image(reaction_smiles, output_path)
                    pixmap = _scale_pixmap_half(QtGui.QPixmap(str(output_path)))
                except Exception as exc:
                    info_label.setText(f"Unable to render reaction image: {exc}")
                else:
                    image_label = QtWidgets.QLabel()
                    image_label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                    image_label.setPixmap(pixmap)
                    images_layout.addWidget(image_label)
            else:
                info_label.setText("Product SMILES missing; showing reactants only.")
                for label, smiles in (("A", reactant_a), ("B", reactant_b)):
                    if not smiles:
                        continue
                    output_path = temp_path / f"reactant_{label.lower()}.png"
                    try:
                        render_molecule_image(smiles, output_path, legend=f"Reactant {label}")
                        pixmap = _scale_pixmap_half(QtGui.QPixmap(str(output_path)))
                    except Exception as exc:
                        info_label.setText(f"Unable to render reactant images: {exc}")
                        break
                    image_label = QtWidgets.QLabel()
                    image_label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                    image_label.setPixmap(pixmap)
                    images_layout.addWidget(image_label)

        if dialog.layout():
            dialog.layout().setSizeConstraint(QtWidgets.QLayout.SizeConstraint.SetFixedSize)
        dialog.adjustSize()
        dialog.show()
        self._reaction_dialog = dialog


def main() -> None:
    app = QtWidgets.QApplication(sys.argv)
    window = HTERecommenderWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()

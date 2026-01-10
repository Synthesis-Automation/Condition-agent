import csv
import math
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Optional, Tuple, Dict, Any, List

from PyQt6 import QtCore, QtWidgets


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
    for label in ("rules", "datasets", "experiments", "experiment", "experiements"):
        if label in parts:
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
        return "datasets"
    if "temperature_c" in header_lower:
        return "rules"
    if "reactant_1" in header_lower and "reactant_2" in header_lower:
        return "experiments"
    return "unknown"


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


def _table_columns_for_type(data_type: str) -> List[Tuple[str, str]]:
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
    ]
    if data_type == "rules":
        return base
    return base + [
        ("Secondary Solvent", "secondary_solvent"),
        ("Coupling Reagent", "coupling_reagent"),
        ("Reaction Type", "reaction_type"),
        ("Reactant Types", "reactant_types"),
        ("Match Score", "match_score"),
    ]


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
    ) -> None:
        super().__init__()
        self.db_path = db_path
        self.reaction_smiles = reaction_smiles
        self.top_k = top_k
        self.min_exp = min_exp
        self.reaction_filter = reaction_filter
        self.catalyst_filter = catalyst_filter

    def run(self) -> None:
        try:
            from chemtools.HTE import HTERecommender

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
            )
            stats = recommender.get_statistics()
            self.finished.emit(True, result, "OK", stats)
        except Exception as exc:
            self.finished.emit(False, None, str(exc), {})


class HTERecommenderWindow(QtWidgets.QWidget):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("HTE Recommender (Qt6)")
        self.resize(1100, 700)

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
        self.reaction_filter_edit.setPlaceholderText("Optional reaction type filter (e.g., Suzuki)")

        self.catalyst_filter_edit = QtWidgets.QLineEdit()
        self.catalyst_filter_edit.setPlaceholderText("Optional catalyst filter (e.g., Pd, Cu)")

        self.run_button = QtWidgets.QPushButton("Run Recommendation")
        self.clear_button = QtWidgets.QPushButton("Clear Results")

        self.summary = QtWidgets.QPlainTextEdit()
        self.summary.setReadOnly(True)

        self.table = QtWidgets.QTableWidget(0, 0)
        self.table.setSortingEnabled(True)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setAlternatingRowColors(True)

        self.status = QtWidgets.QLabel("")

        self._setup_layout()
        self._bind_signals()
        self._set_default_path()

        self.thread: Optional[QtCore.QThread] = None
        self.worker: Optional[RecommendationWorker] = None

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

        form.addRow("Reaction Filter:", self.reaction_filter_edit)
        form.addRow("Catalyst Filter:", self.catalyst_filter_edit)
        layout.addLayout(form)

        button_row = QtWidgets.QHBoxLayout()
        button_row.addStretch()
        button_row.addWidget(self.run_button)
        button_row.addWidget(self.clear_button)
        layout.addLayout(button_row)

        layout.addWidget(QtWidgets.QLabel("Summary:"))
        layout.addWidget(self.summary, stretch=1)
        layout.addWidget(QtWidgets.QLabel("Recommendations:"))
        layout.addWidget(self.table, stretch=3)
        layout.addWidget(self.status)

    def _bind_signals(self) -> None:
        self.run_button.clicked.connect(self._run_recommendation)
        self.clear_button.clicked.connect(self._clear_results)
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

    def _clear_results(self) -> None:
        self.summary.clear()
        self.table.setRowCount(0)
        self.table.setColumnCount(0)
        self.status.setText("")

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

        self.summary.clear()
        self.status.setText("Working...")
        self.run_button.setEnabled(False)

        self.thread = QtCore.QThread()
        self.worker = RecommendationWorker(
            db_path=db_path,
            reaction_smiles=reaction_smiles,
            top_k=self.top_k_spin.value(),
            min_exp=self.min_exp_spin.value(),
            reaction_filter=self.reaction_filter_edit.text().strip(),
            catalyst_filter=self.catalyst_filter_edit.text().strip(),
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

        if not success:
            self.status.setText("Error")
            QtWidgets.QMessageBox.critical(self, "HTE Recommender", message)
            return

        self.status.setText("Done")
        self._render_result(result, stats)

    def _render_result(self, result: object, stats: Dict[str, Any]) -> None:
        data_type = _detect_csv_type(Path(self.db_path_edit.text().strip()))
        summary_lines = [
            f"Data type: {data_type}",
            f"Reactant A: {getattr(result, 'reactant_a_smiles', '')}",
            f"Reactant B: {getattr(result, 'reactant_b_smiles', '') or 'None'}",
            f"Product: {getattr(result, 'product_smiles', '') or 'None'}",
            f"Type A: {getattr(result, 'reactant_a_type', '')} ({getattr(result, 'reactant_a_category', '')})",
            f"Type B: {getattr(result, 'reactant_b_type', '')} ({getattr(result, 'reactant_b_category', '')})",
            f"Predicted: {getattr(result, 'predicted_reaction_type', '')} ({getattr(result, 'reaction_type_confidence', 0.0) * 100:.1f}%)",
            f"Matches: {getattr(result, 'total_matching_experiments', 0)}",
            f"Coverage: {getattr(result, 'database_coverage', 0.0):.2f}%",
        ]
        if getattr(result, "matched_motifs", None):
            summary_lines.append(f"Matched motifs: {getattr(result, 'matched_motifs')}")
        if stats:
            summary_lines.append("")
            summary_lines.append("Database stats:")
            for key in sorted(stats):
                summary_lines.append(f"  {key}: {stats[key]}")
        self.summary.setPlainText("\n".join(summary_lines))

        recs = list(getattr(result, "recommendations", []) or [])
        columns = _table_columns_for_type(data_type)
        self.table.setColumnCount(len(columns))
        self.table.setHorizontalHeaderLabels([label for label, _ in columns])
        self.table.setRowCount(len(recs))

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
                self.table.setItem(row_index, col_index, item)

        self.table.resizeColumnsToContents()


def main() -> None:
    app = QtWidgets.QApplication(sys.argv)
    window = HTERecommenderWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()

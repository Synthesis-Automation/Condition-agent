"""PyQt6 desktop UI for literature-dataset CAS/SMILES extraction."""

from __future__ import annotations

import os
import sys
import traceback
from pathlib import Path
from typing import Optional

from PyQt6 import QtCore, QtWidgets

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DEFAULT_DATASET_DIR = PROJECT_ROOT / "raw_dataset" / "literature_reaction_dataset"
DEFAULT_OUTPUT_PATH = SCRIPT_DIR / "literature_cas_smiles_pairs.csv"
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cas_tools.cas_smiles_extractor import (  # noqa: E402
    CASSmilesPair,
    discover_csv_files,
    extract_cas_smiles_pairs_from_csv,
    write_cas_smiles_pairs,
)


Signal = QtCore.pyqtSignal
Slot = QtCore.pyqtSlot


class CASSmilesExtractorWorker(QtCore.QObject):
    """Run CSV parsing outside the GUI thread."""

    finished = Signal(bool, str)
    progress = Signal(str)

    def __init__(self, folder_path: str, output_csv_path: str) -> None:
        super().__init__()
        self.folder_path = folder_path
        self.output_csv_path = output_csv_path

    @Slot()
    def run(self) -> None:
        try:
            files = discover_csv_files(
                self.folder_path,
                excluded_paths=[self.output_csv_path],
            )
            if not files:
                self.finished.emit(False, "No CSV files were found in the selected folder.")
                return

            all_pairs: set[CASSmilesPair] = set()
            warnings: list[str] = []
            total_rows = 0
            total_occurrences = 0
            for index, path in enumerate(files, start=1):
                relative = os.path.relpath(path, self.folder_path)
                self.progress.emit(f"[{index}/{len(files)}] Reading {relative}")
                try:
                    result = extract_cas_smiles_pairs_from_csv(path)
                except Exception as exc:
                    warnings.append(f"{relative}: {exc}")
                    self.progress.emit(f"    Warning: failed to read file ({exc})")
                    continue
                all_pairs.update(result.pairs)
                total_rows += result.rows_read
                total_occurrences += result.pair_occurrences
                warnings.extend(f"{relative}: {warning}" for warning in result.warnings)
                self.progress.emit(
                    f"    {result.rows_read:,} rows; {len(result.pairs):,} unique pairs"
                )

            self.progress.emit("Writing deterministic four-column CSV...")
            written = write_cas_smiles_pairs(all_pairs, self.output_csv_path)
            conflict_count = _count_conflicting_cas(all_pairs)
            summary = (
                f"Processed {len(files):,} CSV file(s) and {total_rows:,} row(s).\n"
                f"Observed pair occurrences: {total_occurrences:,}\n"
                f"Unique CAS numbers written (one row each): {written:,}\n"
                f"CAS numbers with multiple observed SMILES: {conflict_count:,}\n"
                f"Warnings: {len(warnings):,}\n\n"
                f"Output CSV: {self.output_csv_path}"
            )
            self.finished.emit(True, summary)
        except Exception as exc:
            self.finished.emit(False, f"Error: {exc}\n\n{traceback.format_exc()}")


class CASSmilesExtractorWindow(QtWidgets.QWidget):
    """Small folder-based extraction window modelled on the CAS extractor."""

    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Literature CAS/SMILES Extractor")
        self.resize(760, 560)

        self.folder_edit = QtWidgets.QLineEdit()
        self.output_edit = QtWidgets.QLineEdit(str(DEFAULT_OUTPUT_PATH))
        self.folder_button = QtWidgets.QPushButton("Browse Folder...")
        self.output_button = QtWidgets.QPushButton("Choose Output...")
        self.file_list = QtWidgets.QListWidget()
        self.file_count_label = QtWidgets.QLabel("No folder selected")
        self.run_button = QtWidgets.QPushButton("Extract CAS/SMILES Pairs")
        self.quit_button = QtWidgets.QPushButton("Quit")
        self.log = QtWidgets.QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumHeight(170)

        self.thread: Optional[QtCore.QThread] = None
        self.worker: Optional[CASSmilesExtractorWorker] = None
        self.candidate_files: list[Path] = []

        self._setup_layout()
        self._wire_events()
        if DEFAULT_DATASET_DIR.is_dir():
            self.folder_edit.setText(str(DEFAULT_DATASET_DIR))
            self._update_file_list()
        else:
            self.run_button.setEnabled(False)

    def _setup_layout(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)
        title = QtWidgets.QLabel("Literature Dataset CAS/SMILES Extractor")
        title.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        title.setStyleSheet("font-size: 16px; font-weight: bold; margin: 10px;")
        layout.addWidget(title)

        form = QtWidgets.QFormLayout()
        folder_row = QtWidgets.QHBoxLayout()
        folder_row.addWidget(self.folder_edit)
        folder_row.addWidget(self.folder_button)
        form.addRow("Input Folder:", folder_row)

        output_row = QtWidgets.QHBoxLayout()
        output_row.addWidget(self.output_edit)
        output_row.addWidget(self.output_button)
        form.addRow("Output CSV:", output_row)

        note = QtWidgets.QLabel(
            "Scans CSV files recursively. CAS/SMILES associations are extracted from "
            "nested JSON objects and matching flat columns. The output contains only "
            "cas_no, compound_smiles, reaction_id, and citation; distinct conflicts "
            "are resolved to one deterministic row per CAS number."
        )
        note.setWordWrap(True)
        note.setStyleSheet("font-style: italic; color: #666; font-size: 10px;")
        form.addRow("", note)
        layout.addLayout(form)

        files_group = QtWidgets.QGroupBox("CSV Files To Process")
        files_layout = QtWidgets.QVBoxLayout(files_group)
        files_layout.addWidget(self.file_count_label)
        files_layout.addWidget(self.file_list)
        layout.addWidget(files_group)

        buttons = QtWidgets.QHBoxLayout()
        buttons.addStretch()
        buttons.addWidget(self.run_button)
        buttons.addWidget(self.quit_button)
        layout.addLayout(buttons)

        log_group = QtWidgets.QGroupBox("Processing Log")
        log_layout = QtWidgets.QVBoxLayout(log_group)
        log_layout.addWidget(self.log)
        layout.addWidget(log_group)

    def _wire_events(self) -> None:
        self.folder_button.clicked.connect(self.pick_folder)
        self.output_button.clicked.connect(self.pick_output)
        self.folder_edit.editingFinished.connect(self._update_file_list)
        self.run_button.clicked.connect(self.run_processing)
        self.quit_button.clicked.connect(self.close)

    def pick_folder(self) -> None:
        start = self.folder_edit.text().strip() or str(PROJECT_ROOT)
        selected = QtWidgets.QFileDialog.getExistingDirectory(
            self,
            "Select folder containing literature CSV files",
            start,
            options=QtWidgets.QFileDialog.Option.ShowDirsOnly,
        )
        if selected:
            self.folder_edit.setText(selected)
            self._update_file_list()

    def pick_output(self) -> None:
        start = self.output_edit.text().strip() or str(DEFAULT_OUTPUT_PATH)
        selected, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save CAS/SMILES CSV",
            start,
            "CSV files (*.csv)",
        )
        if selected:
            if not selected.lower().endswith(".csv"):
                selected += ".csv"
            self.output_edit.setText(selected)
            self._update_file_list()

    def _update_file_list(self) -> None:
        folder = self.folder_edit.text().strip()
        output = self.output_edit.text().strip()
        self.file_list.clear()
        self.candidate_files = []
        if not folder or not Path(folder).is_dir():
            self.file_count_label.setText("No valid folder selected")
            self.run_button.setEnabled(False)
            return

        self.candidate_files = discover_csv_files(folder, excluded_paths=[output] if output else [])
        for path in self.candidate_files:
            self.file_list.addItem(os.path.relpath(path, folder))
        count = len(self.candidate_files)
        self.file_count_label.setText(
            f"Found {count:,} CSV file{'s' if count != 1 else ''} (including subfolders)"
        )
        self.run_button.setEnabled(bool(self.candidate_files and output))

    def validate_inputs(self) -> Optional[str]:
        folder = self.folder_edit.text().strip()
        output = self.output_edit.text().strip()
        if not folder or not Path(folder).is_dir():
            return "Please select a valid input folder."
        if not output:
            return "Please choose an output CSV path."
        if not self.candidate_files:
            return "No CSV files were found in the selected folder."
        return None

    def run_processing(self) -> None:
        error = self.validate_inputs()
        if error:
            QtWidgets.QMessageBox.warning(self, "Invalid Input", error)
            return

        self.setEnabled(False)
        self.log.clear()
        self.log.appendPlainText("Starting CAS/SMILES extraction...")
        self.worker = CASSmilesExtractorWorker(
            self.folder_edit.text().strip(),
            self.output_edit.text().strip(),
        )
        self.thread = QtCore.QThread(self)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.progress.connect(self.log.appendPlainText)
        self.worker.finished.connect(self.on_finished)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        self.thread.finished.connect(self._clear_worker)
        self.thread.start()

    @Slot(bool, str)
    def on_finished(self, success: bool, message: str) -> None:
        self.setEnabled(True)
        self.log.appendPlainText(message)
        if success:
            QtWidgets.QMessageBox.information(self, "Extraction Complete", message)
        else:
            QtWidgets.QMessageBox.critical(self, "Extraction Error", message)

    @Slot()
    def _clear_worker(self) -> None:
        self.worker = None
        self.thread = None


def _count_conflicting_cas(pairs: set[CASSmilesPair]) -> int:
    structures: dict[str, set[str]] = {}
    for pair in pairs:
        structures.setdefault(pair.cas_no, set()).add(pair.compound_smiles)
    return sum(len(smiles) > 1 for smiles in structures.values())


def main() -> None:
    app = QtWidgets.QApplication(sys.argv)
    window = CASSmilesExtractorWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()

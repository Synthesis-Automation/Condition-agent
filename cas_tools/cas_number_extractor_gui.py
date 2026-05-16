from __future__ import annotations

import csv
import os
import sys
import traceback
from pathlib import Path
from typing import List, Optional

from PyQt6 import QtCore, QtWidgets

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cas_tools.cas_number_extractor import (
    discover_candidate_files,
    extract_cas_matches_from_file,
    write_matches_to_csv,
)


if hasattr(QtCore, "Signal") and hasattr(QtCore, "Slot"):
    Signal = QtCore.Signal
    Slot = QtCore.Slot
elif hasattr(QtCore, "pyqtSignal") and hasattr(QtCore, "pyqtSlot"):
    Signal = QtCore.pyqtSignal  # type: ignore[attr-defined]
    Slot = QtCore.pyqtSlot  # type: ignore[attr-defined]
else:  # pragma: no cover
    Signal = None  # type: ignore
    Slot = None  # type: ignore


def default_output_path(folder_path: str) -> str:
    return str(SCRIPT_DIR / "cas_no_all.csv")


class CASExtractorWorker(QtCore.QObject):
    finished = Signal(bool, str) if Signal else None  # type: ignore[misc]
    progress = Signal(str) if Signal else None  # type: ignore[misc]

    def __init__(self, folder_path: str, output_csv_path: str):
        super().__init__()
        self.folder_path = folder_path
        self.output_csv_path = output_csv_path

    def _emit(self, message: str) -> None:
        signal = getattr(self, "progress", None)
        if signal:
            try:
                signal.emit(message)
            except Exception:
                pass

    @Slot() if Slot else (lambda func: func)
    def run(self) -> None:
        try:
            self._emit("Scanning folder for supported files...")
            files = discover_candidate_files(
                self.folder_path,
                excluded_paths=[self.output_csv_path],
            )
            if not files:
                if self.finished:
                    self.finished.emit(False, "No supported files found in the selected folder.")
                return

            self._emit(f"Found {len(files)} supported file(s).")
            all_matches = []
            warnings: List[str] = []

            for index, path in enumerate(files, start=1):
                rel_path = os.path.relpath(path, self.folder_path)
                self._emit(f"[{index}/{len(files)}] Processing {rel_path}")
                matches, file_warnings = extract_cas_matches_from_file(path, base_folder=self.folder_path)
                all_matches.extend(matches)
                warnings.extend(file_warnings)
                self._emit(f"      {len(matches)} CAS match(es)")
                for warning in file_warnings:
                    self._emit(f"      Warning: {warning}")

            self._emit("Writing CSV output...")
            existing_rows = _count_existing_cas_numbers(self.output_csv_path)
            write_matches_to_csv(all_matches, self.output_csv_path)
            updated_rows = _count_existing_cas_numbers(self.output_csv_path)
            added_rows = max(updated_rows - existing_rows, 0)
            self._emit(f"Saved {added_rows} new CAS number(s) to {self.output_csv_path}")

            summary = (
                f"Processed {len(files)} file(s).\n"
                f"CAS matches: {len(all_matches)}\n"
                f"New CAS numbers added: {added_rows}\n"
                f"Warnings: {len(warnings)}\n\n"
                f"Output CSV: {self.output_csv_path}"
            )
            if self.finished:
                self.finished.emit(True, summary)
        except Exception as exc:
            message = f"Error: {exc}\n\n{traceback.format_exc()}"
            if self.finished:
                self.finished.emit(False, message)


class CASExtractorWindow(QtWidgets.QWidget):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("CAS Number Extractor")
        self.resize(700, 500)

        self.folder_edit = QtWidgets.QLineEdit()
        self.btn_folder = QtWidgets.QPushButton("Browse Folder...")

        self.file_list = QtWidgets.QListWidget()
        self.file_count_label = QtWidgets.QLabel("No folder selected")

        self.btn_run = QtWidgets.QPushButton("Extract CAS Numbers")
        self.btn_quit = QtWidgets.QPushButton("Quit")

        self.log = QtWidgets.QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumHeight(150)

        self.thread = None
        self.worker = None
        self.candidate_files: List[Path] = []

        self._setup_layout()
        self._wire_events()
        self.btn_run.setEnabled(False)

    def _setup_layout(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)

        title = QtWidgets.QLabel("Folder-Based CAS Number Extractor")
        title.setStyleSheet("font-size: 16px; font-weight: bold; margin: 10px;")
        title.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title)

        form = QtWidgets.QFormLayout()

        folder_box = QtWidgets.QHBoxLayout()
        folder_box.addWidget(self.folder_edit)
        folder_box.addWidget(self.btn_folder)
        folder_hint = QtWidgets.QLabel("(includes subfolders)")
        folder_hint.setStyleSheet("font-style: italic; color: #666; font-size: 9px;")
        form.addRow("Input Folder:", folder_box)
        form.addRow("", folder_hint)

        note_label = QtWidgets.QLabel(
            "Supported inputs: PDF (text-based), text-like files, CSV/TSV, Word .docx/.docm, Excel .xlsx/.xlsm/.xls\n"
            f"Output: {default_output_path('')} (append-only, unique CAS numbers only)"
        )
        note_label.setStyleSheet("font-style: italic; color: #666; font-size: 10px;")
        form.addRow("", note_label)

        layout.addLayout(form)

        file_group = QtWidgets.QGroupBox("Files To Process")
        file_layout = QtWidgets.QVBoxLayout(file_group)
        file_layout.addWidget(self.file_count_label)
        file_layout.addWidget(self.file_list)
        layout.addWidget(file_group)

        button_layout = QtWidgets.QHBoxLayout()
        button_layout.addStretch()
        button_layout.addWidget(self.btn_run)
        button_layout.addWidget(self.btn_quit)
        layout.addLayout(button_layout)

        log_group = QtWidgets.QGroupBox("Processing Log")
        log_layout = QtWidgets.QVBoxLayout(log_group)
        log_layout.addWidget(self.log)
        layout.addWidget(log_group)

    def _wire_events(self) -> None:
        self.btn_folder.clicked.connect(self.pick_folder)
        self.btn_run.clicked.connect(self.run_processing)
        self.btn_quit.clicked.connect(self.close)

    def log_msg(self, text: str) -> None:
        self.log.appendPlainText(text)

    def pick_folder(self) -> None:
        path = QtWidgets.QFileDialog.getExistingDirectory(
            self,
            "Select folder with files",
            os.getcwd(),
            options=QtWidgets.QFileDialog.Option.ShowDirsOnly,
        )
        if path:
            self.folder_edit.setText(path)
            self._update_file_list()

    def _update_file_list(self) -> None:
        folder_path = self.folder_edit.text().strip()
        self.file_list.clear()
        self.candidate_files = []

        if not folder_path or not os.path.isdir(folder_path):
            self.file_count_label.setText("No valid folder selected")
            self.btn_run.setEnabled(False)
            return

        output_csv = default_output_path(folder_path)
        try:
            self.candidate_files = discover_candidate_files(folder_path, excluded_paths=[output_csv])
            if not self.candidate_files:
                self.file_count_label.setText("No supported files found in this folder or subfolders")
                self.btn_run.setEnabled(False)
                return

            for path in self.candidate_files:
                self.file_list.addItem(os.path.relpath(path, folder_path))

            count = len(self.candidate_files)
            self.file_count_label.setText(f"Found {count} supported file{'s' if count != 1 else ''} (including subfolders)")
            self.btn_run.setEnabled(True)
        except Exception as exc:
            self.file_count_label.setText(f"Error reading folder: {exc}")
            self.btn_run.setEnabled(False)

    def validate_inputs(self) -> Optional[str]:
        folder_path = self.folder_edit.text().strip()
        if not folder_path or not os.path.isdir(folder_path):
            return "Please select a valid folder."
        if not self.candidate_files:
            return "No supported files found in the selected folder."
        return None

    def run_processing(self) -> None:
        error = self.validate_inputs()
        if error:
            QtWidgets.QMessageBox.warning(self, "Invalid Input", error)
            return

        self.setEnabled(False)
        self.log.clear()
        self.log_msg("Starting CAS extraction...")

        folder_path = self.folder_edit.text().strip()
        output_csv = default_output_path(folder_path)

        self.worker = CASExtractorWorker(
            folder_path=folder_path,
            output_csv_path=output_csv,
        )
        self.thread = QtCore.QThread(self)
        self.worker.moveToThread(self.thread)

        self.thread.started.connect(self.worker.run)

        finished_signal = getattr(self.worker, "finished", None)
        if finished_signal:
            finished_signal.connect(self.on_finished)
            finished_signal.connect(self.thread.quit)
            finished_signal.connect(self.worker.deleteLater)

        progress_signal = getattr(self.worker, "progress", None)
        if progress_signal:
            progress_signal.connect(self.log_msg)

        self.thread.finished.connect(self.thread.deleteLater)
        self.thread.finished.connect(lambda: self.setEnabled(True))
        self.thread.finished.connect(lambda: setattr(self, "worker", None))
        self.thread.finished.connect(lambda: setattr(self, "thread", None))
        self.thread.start()

    def on_finished(self, success: bool, message: str) -> None:
        self.setEnabled(True)
        self.log_msg(message)
        if success:
            QtWidgets.QMessageBox.information(self, "Extraction Complete", message)
        else:
            QtWidgets.QMessageBox.critical(self, "Extraction Error", message)


def _count_existing_cas_numbers(output_csv_path: str) -> int:
    path = Path(output_csv_path)
    if not path.exists() or path.stat().st_size == 0:
        return 0

    count = 0
    try:
        with path.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.reader(handle)
            for row in reader:
                if not row:
                    continue
                value = str(row[0]).strip()
                if value and value.lower() != "cas_number":
                    count += 1
    except Exception:
        return 0
    return count


def main() -> None:
    if hasattr(QtWidgets, "QApplication"):
        try:
            QtWidgets.QApplication.setAttribute(QtCore.Qt.ApplicationAttribute.AA_EnableHighDpiScaling, True)
            QtWidgets.QApplication.setAttribute(QtCore.Qt.ApplicationAttribute.AA_UseHighDpiPixmaps, True)
        except Exception:
            pass

    app = QtWidgets.QApplication(sys.argv)
    window = CASExtractorWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()

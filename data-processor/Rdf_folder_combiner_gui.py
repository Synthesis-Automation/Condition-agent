#!/usr/bin/env python3
"""PyQt6 GUI to combine RDF files from a folder into one output file."""
from __future__ import annotations

import os
import shutil
import sys
from typing import List

from PyQt6 import QtCore, QtWidgets


def _sort_key(path: str) -> tuple[str, str]:
    """Sort by folder, then filename for stable ordering."""
    return (os.path.dirname(path).lower(), os.path.basename(path).lower())


def find_rdf_files(folder: str) -> List[str]:
    """Return a sorted list of all *.rdf files under folder (recurses)."""
    matches: List[str] = []
    for root, _, files in os.walk(folder):
        for name in files:
            if name.lower().endswith(".rdf"):
                matches.append(os.path.join(root, name))
    matches.sort(key=_sort_key)
    return matches


def _suggest_output_path(folder: str) -> str:
    """Suggest an output path inside the folder."""
    base = os.path.basename(os.path.abspath(folder)) or "combined"
    return os.path.join(folder, f"{base}_combined.rdf")


class CombineWorker(QtCore.QThread):
    """Background worker that concatenates RDF files."""

    progress = QtCore.pyqtSignal(str)
    progress_step = QtCore.pyqtSignal(int, int, str)  # current, total, relpath
    finished = QtCore.pyqtSignal(str, int)  # output_path, file_count
    failed = QtCore.pyqtSignal(str)

    def __init__(self, folder: str, output_path: str) -> None:
        super().__init__()
        self.folder = folder
        self.output_path = output_path

    def run(self) -> None:  # pragma: no cover - GUI thread
        try:
            folder = self.folder
            files = find_rdf_files(folder)
            if not files:
                self.failed.emit("No .rdf files found in the selected folder.")
                return

            total = len(files)
            self.progress.emit(f"Found {total} RDF file{'s' if total != 1 else ''}.")

            out_dir = os.path.dirname(self.output_path)
            if out_dir:
                os.makedirs(out_dir, exist_ok=True)

            with open(self.output_path, "wb") as out_file:
                for idx, path in enumerate(files, start=1):
                    relpath = os.path.relpath(path, folder)
                    self.progress.emit(f"Appending {relpath}")
                    self.progress_step.emit(idx, total, relpath)
                    with open(path, "rb") as in_file:
                        shutil.copyfileobj(in_file, out_file)
                    if idx != total:
                        out_file.write(b"\n")

            self.finished.emit(self.output_path, total)
        except Exception as exc:  # pragma: no cover - best effort logging
            self.failed.emit(str(exc))


class MainWindow(QtWidgets.QWidget):
    """Simple window to drive the RDF combiner."""

    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("RDF Folder Combiner")
        self.resize(720, 480)

        self.worker: CombineWorker | None = None

        self._build_ui()

    def _build_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)

        row_folder = QtWidgets.QHBoxLayout()
        self.folder_edit = QtWidgets.QLineEdit()
        self.btn_browse = QtWidgets.QPushButton("Select Folder...")
        self.btn_browse.clicked.connect(self._choose_folder)
        row_folder.addWidget(QtWidgets.QLabel("Folder:"))
        row_folder.addWidget(self.folder_edit, 1)
        row_folder.addWidget(self.btn_browse)
        layout.addLayout(row_folder)

        row_output = QtWidgets.QHBoxLayout()
        self.output_edit = QtWidgets.QLineEdit()
        self.btn_output = QtWidgets.QPushButton("Save As...")
        self.btn_output.clicked.connect(self._choose_output)
        row_output.addWidget(QtWidgets.QLabel("Output RDF:"))
        row_output.addWidget(self.output_edit, 1)
        row_output.addWidget(self.btn_output)
        layout.addLayout(row_output)

        self.btn_run = QtWidgets.QPushButton("Combine Files")
        self.btn_run.clicked.connect(self._run)
        layout.addWidget(self.btn_run, alignment=QtCore.Qt.AlignmentFlag.AlignRight)

        self.progress_bar = QtWidgets.QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(1)
        self.progress_bar.setValue(0)
        layout.addWidget(self.progress_bar)

        self.log = QtWidgets.QPlainTextEdit()
        self.log.setReadOnly(True)
        layout.addWidget(self.log, 1)

    def _choose_folder(self) -> None:
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, "Select Folder")
        if folder:
            self.folder_edit.setText(folder)
            if not self.output_edit.text().strip():
                self.output_edit.setText(_suggest_output_path(folder))

    def _choose_output(self) -> None:
        start = self.output_edit.text().strip()
        if not start:
            folder = self.folder_edit.text().strip() or os.getcwd()
            start = _suggest_output_path(folder)
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save Combined RDF",
            start,
            "RDF Files (*.rdf);;All Files (*.*)",
        )
        if path:
            if not path.lower().endswith(".rdf"):
                path += ".rdf"
            self.output_edit.setText(path)

    def _run(self) -> None:
        if self.worker is not None:
            self._log("Already running.")
            return

        folder = self.folder_edit.text().strip()
        if not folder or not os.path.isdir(folder):
            self._log("Please choose a valid folder.")
            return

        output_path = self.output_edit.text().strip()
        if not output_path:
            output_path = _suggest_output_path(folder)
            self.output_edit.setText(output_path)

        self._log(f"Combining RDF files under {folder}")
        self._set_running(True)

        self.worker = CombineWorker(folder, output_path)
        self.worker.progress.connect(self._log)
        self.worker.progress_step.connect(self._on_progress)
        self.worker.finished.connect(self._on_finished)
        self.worker.failed.connect(self._on_failed)
        self.worker.finished.connect(lambda *_: self._cleanup_worker())
        self.worker.failed.connect(lambda *_: self._cleanup_worker())
        self.worker.start()

    def _cleanup_worker(self) -> None:
        self.worker = None
        self._set_running(False)

    def _set_running(self, running: bool) -> None:
        self.btn_run.setEnabled(not running)
        self.btn_browse.setEnabled(not running)
        self.btn_output.setEnabled(not running)
        if running:
            self.progress_bar.setMaximum(0)
            self.progress_bar.setValue(0)
        else:
            self.progress_bar.setMaximum(1)
            self.progress_bar.setValue(0)

    def _log(self, message: str) -> None:
        self.log.appendPlainText(str(message))

    def _on_progress(self, current: int, total: int, relpath: str) -> None:
        try:
            self.progress_bar.setMaximum(total)
        except Exception:
            pass
        self.progress_bar.setValue(current)

    def _on_finished(self, output_path: str, count: int) -> None:
        self._log(f"Success: combined {count} file{'s' if count != 1 else ''} into {output_path}")
        QtWidgets.QMessageBox.information(
            self,
            "Combine Complete",
            f"Combined {count} RDF file{'s' if count != 1 else ''} into:\n{output_path}",
        )

    def _on_failed(self, message: str) -> None:
        self._log(f"Error: {message}")
        QtWidgets.QMessageBox.critical(self, "Combine Failed", message)


def launch() -> None:  # pragma: no cover - GUI runtime
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    launch()

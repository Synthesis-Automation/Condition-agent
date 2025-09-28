#!/usr/bin/env python3
"""PyQt6 GUI to combine RDF files and build/search a CAS index."""
from __future__ import annotations

import importlib.util
import json
import os
import shutil
import sys
import traceback
from pathlib import Path
from typing import Any, List, Optional, Tuple

from PyQt6 import QtCore, QtWidgets

DEFAULT_ENCODING = "latin-1"
PREVIEW_MAX_LINES = 40
PREVIEW_MAX_CHARS = 1500

_INDEXER_MODULE: Optional[Any] = None


def _sort_key(path: str) -> Tuple[str, str]:
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


def _suggest_index_path(output_path: str) -> str:
    base, _ = os.path.splitext(os.path.abspath(output_path))
    return f"{base}_cas_index.jsonl"


def _load_indexer_module() -> Any:
    """Load the indexer helper lazily so imports stay optional."""
    global _INDEXER_MODULE
    if _INDEXER_MODULE is not None:
        return _INDEXER_MODULE

    module_path = Path(__file__).with_name("rdf_reaction_indexer.py")
    if not module_path.exists():
        raise FileNotFoundError(f"Index helper not found next to GUI: {module_path}")

    spec = importlib.util.spec_from_file_location("rdf_reaction_indexer_gui", str(module_path))
    if spec is None or spec.loader is None:
        raise ImportError(f"Unable to load indexer module from {module_path}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module  # ensure visibility during import
    spec.loader.exec_module(module)  # type: ignore[arg-type]
    _INDEXER_MODULE = module
    return module


def _read_reaction_block(
    source: str,
    offset: int,
    length: int,
    *,
    encoding: str = DEFAULT_ENCODING,
    max_lines: Optional[int] = PREVIEW_MAX_LINES,
    max_chars: Optional[int] = PREVIEW_MAX_CHARS,
) -> str:
    """Read a reaction block with optional truncation."""
    try:
        path = Path(source)
        if not path.exists():
            return f"[Source file not found: {source}]"
        with path.open("rb") as handle:
            handle.seek(int(offset))
            data = handle.read(int(length))
        text = data.decode(encoding, errors="replace").strip()
        if not text:
            return ""
        if max_lines is not None or max_chars is not None:
            lines = text.splitlines()
            if max_lines is not None and len(lines) > max_lines:
                text = "\n".join(lines[:max_lines]) + f"\n... (truncated after {max_lines} lines)"
            if max_chars is not None and len(text) > max_chars:
                text = text[:max_chars] + "... (truncated)"
        return text
    except Exception as exc:  # pragma: no cover - best effort preview
        return f"[Error reading block: {exc}]"


class CombineWorker(QtCore.QThread):
    """Background worker that concatenates RDF files (and optionally builds the CAS index)."""

    progress = QtCore.pyqtSignal(str)
    progress_step = QtCore.pyqtSignal(int, int, str)  # current, total, relpath
    finished = QtCore.pyqtSignal(str, int, dict)  # output_path, file_count, extras
    failed = QtCore.pyqtSignal(str)

    def __init__(
        self,
        folder: str,
        output_path: str,
        *,
        build_index: bool,
        index_path: Optional[str],
        metadata_path: Optional[str] = None,
        encoding: str = DEFAULT_ENCODING,
        metadata_only: bool = True,
    ) -> None:
        super().__init__()
        self.folder = folder
        self.output_path = os.path.abspath(output_path)
        self.build_index = build_index
        self.index_path = os.path.abspath(index_path) if index_path else None
        self.metadata_path = os.path.abspath(metadata_path) if metadata_path else None
        self.encoding = encoding
        self.metadata_only = metadata_only

    def run(self) -> None:  # pragma: no cover - GUI thread
        try:
            files = find_rdf_files(self.folder)
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
                    relpath = os.path.relpath(path, self.folder)
                    self.progress.emit(f"Appending {relpath}")
                    self.progress_step.emit(idx, total, relpath)
                    with open(path, "rb") as in_file:
                        shutil.copyfileobj(in_file, out_file)
                    if idx != total:
                        out_file.write(b"\n")

            summary: dict[str, Any] = {
                "index_path": self.index_path if self.build_index else None,
                "metadata_path": self.metadata_path if self.build_index else None,
                "reactions": None,
                "cas_entries": None,
            }

            if self.build_index:
                if not self.index_path and not self.metadata_path:
                    summary["index_error"] = "Index output path not provided."
                    self.progress.emit("Indexing skipped: no output path provided.")
                else:
                    try:
                        self.progress.emit("Building CAS index (this may take a while)...")
                        indexer = _load_indexer_module()
                        if self.index_path:
                            idx_dir = os.path.dirname(self.index_path)
                            if idx_dir:
                                os.makedirs(idx_dir, exist_ok=True)
                        if self.metadata_path:
                            meta_dir = os.path.dirname(self.metadata_path)
                            if meta_dir:
                                os.makedirs(meta_dir, exist_ok=True)

                        blocks = indexer.scan_reactions(
                            self.output_path,
                            encoding=self.encoding,
                        )
                        reactions, cas_entries = indexer.write_outputs(
                            blocks,
                            source_path=self.output_path,
                            jsonl_path=self.metadata_path,
                            index_path=self.index_path,
                            metadata_only=self.metadata_only,
                            quiet=True,
                        )
                        summary["reactions"] = reactions
                        summary["cas_entries"] = cas_entries
                        self.progress.emit(
                            f"CAS index built: {cas_entries} entries across {reactions} reactions."
                        )
                    except Exception as exc:
                        summary["index_error"] = f"{exc}\n{traceback.format_exc()}"
                        self.progress.emit(f"Index build failed: {exc}")

            self.finished.emit(self.output_path, total, summary)
        except Exception as exc:  # pragma: no cover - best effort logging
            self.failed.emit(str(exc))


class RecordDialog(QtWidgets.QDialog):
    """Modal dialog to display full reaction records."""

    def __init__(self, parent: QtWidgets.QWidget, title: str, content: str) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(900, 600)

        layout = QtWidgets.QVBoxLayout(self)
        viewer = QtWidgets.QPlainTextEdit()
        viewer.setReadOnly(True)
        viewer.setPlainText(content)
        layout.addWidget(viewer)

        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Close
        )
        buttons.rejected.connect(self.accept)
        layout.addWidget(buttons)



class MainWindow(QtWidgets.QWidget):
    """Simple window to drive the RDF combiner and CAS index/search."""

    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("RDF Folder Combiner")
        self.resize(820, 560)

        self.worker: Optional[CombineWorker] = None
        self.encoding = DEFAULT_ENCODING
        self.last_index_info: Optional[dict[str, Any]] = None
        self._index_path_dirty = False

        self._build_ui()

    def _build_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)

        row_folder = QtWidgets.QHBoxLayout()
        self.folder_edit = QtWidgets.QLineEdit()
        self.btn_browse = QtWidgets.QPushButton("Select Folder…")
        self.btn_browse.clicked.connect(self._choose_folder)
        row_folder.addWidget(QtWidgets.QLabel("Folder:"))
        row_folder.addWidget(self.folder_edit, 1)
        row_folder.addWidget(self.btn_browse)
        layout.addLayout(row_folder)

        row_output = QtWidgets.QHBoxLayout()
        self.output_edit = QtWidgets.QLineEdit()
        self.output_edit.textChanged.connect(self._handle_output_changed)
        self.btn_output = QtWidgets.QPushButton("Save As…")
        self.btn_output.clicked.connect(self._choose_output)
        row_output.addWidget(QtWidgets.QLabel("Output RDF:"))
        row_output.addWidget(self.output_edit, 1)
        row_output.addWidget(self.btn_output)
        layout.addLayout(row_output)

        index_group = QtWidgets.QGroupBox("CAS Index")
        index_layout = QtWidgets.QVBoxLayout(index_group)
        self.chk_index = QtWidgets.QCheckBox("Build CAS index after combine")
        self.chk_index.setChecked(True)
        index_layout.addWidget(self.chk_index)

        index_row = QtWidgets.QHBoxLayout()
        index_row.addWidget(QtWidgets.QLabel("Index file:"))
        self.index_edit = QtWidgets.QLineEdit()
        self.index_edit.setPlaceholderText("Auto-filled from Output RDF")
        self.index_edit.textEdited.connect(self._mark_index_dirty)
        index_row.addWidget(self.index_edit, 1)
        self.btn_index_browse = QtWidgets.QPushButton("Save As…")
        self.btn_index_browse.clicked.connect(self._choose_index_output)
        self.btn_index_open = QtWidgets.QPushButton("Load Existing…")
        self.btn_index_open.clicked.connect(self._open_index_file)
        self.btn_index_reset = QtWidgets.QPushButton("Reset")
        self.btn_index_reset.clicked.connect(self._reset_index_path)
        index_row.addWidget(self.btn_index_browse)
        index_row.addWidget(self.btn_index_open)
        index_row.addWidget(self.btn_index_reset)
        index_layout.addLayout(index_row)
        layout.addWidget(index_group)

        search_group = QtWidgets.QGroupBox("CAS Search")
        search_layout = QtWidgets.QHBoxLayout(search_group)
        search_layout.addWidget(QtWidgets.QLabel("CAS reaction #:"))
        self.search_edit = QtWidgets.QLineEdit()
        self.search_edit.setPlaceholderText("e.g., 31-352-CAS-5769784")
        self.search_edit.returnPressed.connect(self._search_cas)
        search_layout.addWidget(self.search_edit, 1)
        self.btn_search = QtWidgets.QPushButton("Search")
        self.btn_search.clicked.connect(self._search_cas)
        search_layout.addWidget(self.btn_search)
        layout.addWidget(search_group)

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

    def _choose_index_output(self) -> None:
        start = self.index_edit.text().strip()
        if not start:
            output = self.output_edit.text().strip()
            if output:
                start = _suggest_index_path(output)
            else:
                folder = self.folder_edit.text().strip() or os.getcwd()
                start = os.path.join(folder, "cas_index.jsonl")
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save CAS Index",
            start,
            "JSON Lines (*.jsonl);;All Files (*.*)",
        )
        if path:
            if not path.lower().endswith(".jsonl"):
                path += ".jsonl"
            self.index_edit.setText(path)
            self._index_path_dirty = True

    def _open_index_file(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Open CAS Index",
            self.index_edit.text().strip() or os.getcwd(),
            "JSON Lines (*.jsonl);;All Files (*.*)",
        )
        if path:
            self.index_edit.setText(path)
            self._index_path_dirty = True

    def _reset_index_path(self) -> None:
        output = self.output_edit.text().strip()
        if not output:
            folder = self.folder_edit.text().strip()
            if folder:
                output = _suggest_output_path(folder)
                self.output_edit.setText(output)
        if output:
            self._index_path_dirty = False
            self.index_edit.setText(_suggest_index_path(output))

    def _mark_index_dirty(self) -> None:
        self._index_path_dirty = True

    def _handle_output_changed(self, text: str) -> None:
        if text.strip() and not self._index_path_dirty:
            self.index_edit.setText(_suggest_index_path(text.strip()))

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

        build_index = self.chk_index.isChecked()
        index_path = self.index_edit.text().strip()
        if build_index and not index_path:
            index_path = _suggest_index_path(output_path)
            self.index_edit.setText(index_path)
        if build_index and not index_path:
            self._log("Please provide a path for the CAS index.")
            return

        self._log(f"Combining RDF files under {folder}")
        self._set_running(True)

        self.worker = CombineWorker(
            folder,
            output_path,
            build_index=build_index,
            index_path=index_path if build_index else None,
            encoding=self.encoding,
        )
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
        self.chk_index.setEnabled(not running)
        self.index_edit.setEnabled(not running)
        self.btn_index_browse.setEnabled(not running)
        self.btn_index_open.setEnabled(not running)
        self.btn_index_reset.setEnabled(not running)
        self.btn_search.setEnabled(not running)
        self.search_edit.setEnabled(not running)
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

    def _on_finished(self, output_path: str, count: int, info: dict) -> None:
        self.last_index_info = info or {}
        summary_lines = [
            f"Success: combined {count} file{'s' if count != 1 else ''} into {output_path}",
        ]
        index_path = info.get("index_path")
        index_error = info.get("index_error")
        cas_entries = info.get("cas_entries")
        reactions = info.get("reactions")

        if index_path and not index_error:
            summary_lines.append(
                f"CAS index: {index_path} ({cas_entries or 0} entries across {reactions or 0} reactions)"
            )
            self.index_edit.setText(str(index_path))
            self._index_path_dirty = True
        if index_error:
            summary_lines.append(f"Indexing error: {index_error}")

        for line in summary_lines:
            self._log(line)

        QtWidgets.QMessageBox.information(self, "Combine Complete", "\n".join(summary_lines))

    def _on_failed(self, message: str) -> None:
        self._log(f"Error: {message}")
        QtWidgets.QMessageBox.critical(self, "Combine Failed", message)

    def _search_cas(self) -> None:
        cas = self.search_edit.text().strip()
        if not cas:
            self._log("Enter a CAS reaction number to search.")
            return

        index_path = self.index_edit.text().strip()
        if not index_path:
            self._log("Provide a CAS index file path first.")
            return
        if not os.path.exists(index_path):
            self._log(f"Index file not found: {index_path}")
            QtWidgets.QMessageBox.warning(self, "CAS Search", f"Index file not found:\n{index_path}")
            return

        matches: List[dict[str, Any]] = []
        try:
            with open(index_path, "r", encoding="utf-8") as handle:
                for line_no, line in enumerate(handle, 1):
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        record = json.loads(line)
                    except json.JSONDecodeError as exc:
                        self._log(f"[Index parse error at line {line_no}: {exc}]")
                        continue
                    if record.get("cas") == cas:
                        matches.append(record)
        except Exception as exc:
            self._log(f"Error reading index: {exc}")
            QtWidgets.QMessageBox.critical(self, "CAS Search", f"Cannot read index: {exc}")
            return

        if not matches:
            self._log(f"No entries found for {cas} in {index_path}")
            QtWidgets.QMessageBox.information(self, "CAS Search", f"No entries found for {cas}.")
            return

        full_records: List[str] = []

        for idx, record in enumerate(matches, 1):
            scheme = record.get("scheme") or "-"
            step = record.get("step") or "-"
            variation = record.get("variation")
            sequence = record.get("sequence")
            offset = int(record.get("offset", 0))
            length = int(record.get("length", 0))
            source_file = record.get("source_file") or self.output_edit.text().strip()
            header = (
                f"[{idx}/{len(matches)}] CAS {cas} -> sequence {sequence}, variation {variation}, "
                f"scheme {scheme}, step {step}, offset {offset}, length {length}"
            )
            self._log(header)
            snippet = _read_reaction_block(source_file, offset, length, encoding=self.encoding)
            if snippet:
                self._log(snippet)

            full_text = _read_reaction_block(
                source_file,
                offset,
                length,
                encoding=self.encoding,
                max_lines=None,
                max_chars=None,
            )
            full_records.append(f"{header}\n{full_text or '[No content]'}")

        dialog = RecordDialog(
            self,
            f"CAS Search: {cas}",
            "\n\n".join(full_records),
        )
        dialog.exec()
        self._log(
            f"Displayed full record(s) for {cas} ({len(matches)} match{'es' if len(matches) != 1 else ''})."
        )


def launch() -> None:  # pragma: no cover - GUI runtime
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    launch()

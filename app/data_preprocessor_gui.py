"""PyQt6 app for chemistry-free source-file preprocessing."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict, Iterable, Optional

from PyQt6 import QtCore, QtGui, QtWidgets

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT_FOLDER = PROJECT_ROOT / "datasets" / "intermediate"
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from condition_recommender.ingestion import (  # noqa: E402
    PreprocessingCancelled,
    PreprocessingProgress,
    adapter_ids,
    detect_adapter,
    preprocess_files,
)


def _human_size(size_bytes: int) -> str:
    value = float(size_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024.0 or unit == "TB":
            return f"{value:.0f} {unit}" if unit == "B" else f"{value:.1f} {unit}"
        value /= 1024.0
    return f"{size_bytes} B"


class DataPreprocessingWorker(QtCore.QObject):
    """Preprocess selected files on a worker thread."""

    progress = QtCore.pyqtSignal(object)
    finished = QtCore.pyqtSignal(bool, object, str)

    def __init__(
        self,
        source_files: Iterable[str],
        output_folder: str,
        *,
        adapter_id: str = "auto",
        force: bool = False,
    ) -> None:
        super().__init__()
        self.source_files = tuple(source_files)
        self.output_folder = output_folder
        self.adapter_id = adapter_id
        self.force = force
        self._cancel_requested = False

    def request_cancel(self) -> None:
        """Request a safe stop during or between source files."""
        self._cancel_requested = True

    @QtCore.pyqtSlot()
    def run(self) -> None:
        """Run preprocessing and emit a terminal result."""
        try:
            report = preprocess_files(
                self.source_files,
                self.output_folder,
                adapter_id=self.adapter_id,
                force=self.force,
                progress_callback=self.progress.emit,
                cancel_check=lambda: self._cancel_requested,
            )
        except PreprocessingCancelled as exc:
            self.finished.emit(False, {}, str(exc))
        except Exception as exc:
            self.finished.emit(False, {}, f"{type(exc).__name__}: {exc}")
        else:
            self.finished.emit(True, report, "")


class SourceDataPreprocessorWindow(QtWidgets.QWidget):
    """Desktop controller for file-level source normalization."""

    def __init__(self) -> None:
        super().__init__()
        self.setFont(QtGui.QFont("Segoe UI", 9))
        self.setWindowTitle("Reaction Source Data Preprocessor")
        self.resize(920, 720)
        self.thread: Optional[QtCore.QThread] = None
        self.worker: Optional[DataPreprocessingWorker] = None
        self._completed_output: Optional[Path] = None

        self.source_list = QtWidgets.QListWidget()
        self.source_list.setObjectName("sourceFiles")
        self.source_list.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection
        )
        self.source_list.setAlternatingRowColors(True)

        self.source_summary = QtWidgets.QLabel("No source files selected.")
        self.source_summary.setWordWrap(True)
        self.output_edit = QtWidgets.QLineEdit(str(DEFAULT_OUTPUT_FOLDER))
        self.output_edit.setObjectName("outputFolder")
        self.output_edit.setPlaceholderText(
            "Folder for one intermediate artifact per source file"
        )
        self.adapter_combo = QtWidgets.QComboBox()
        self.adapter_combo.setObjectName("sourceAdapter")
        self.adapter_combo.addItem("Detect from exact columns", "auto")
        for adapter_id in adapter_ids():
            self.adapter_combo.addItem(adapter_id, adapter_id)
        self.force_check = QtWidgets.QCheckBox(
            "Rebuild unchanged files instead of reusing verified artifacts"
        )
        self.force_check.setObjectName("forceRebuild")
        self.force_check.setChecked(False)

        self.start_button = QtWidgets.QPushButton("Preprocess Selected Files")
        self.start_button.setObjectName("preprocessButton")
        self.start_button.setDefault(True)
        self.start_button.setStyleSheet(
            "QPushButton#preprocessButton {"
            "background-color: #0078d7; color: white; font-weight: 700; "
            "padding: 10px 18px; border-radius: 6px;}"
            "QPushButton#preprocessButton:disabled {"
            "background-color: #a6c8f0; color: white;}"
        )
        self.cancel_button = QtWidgets.QPushButton("Cancel")
        self.cancel_button.setObjectName("cancelButton")
        self.cancel_button.setEnabled(False)
        self.open_button = QtWidgets.QPushButton("Open Output Folder")
        self.open_button.setObjectName("openOutputButton")
        self.open_button.setEnabled(False)
        self.progress_bar = QtWidgets.QProgressBar()
        self.progress_bar.setObjectName("preprocessingProgress")
        self.progress_bar.setRange(0, 1)
        self.progress_bar.setValue(0)
        self.progress_bar.setFormat("Waiting")
        self.status_box = QtWidgets.QPlainTextEdit()
        self.status_box.setObjectName("statusBox")
        self.status_box.setReadOnly(True)
        self.status_box.setPlaceholderText(
            "Adapter selection, row progress, cache reuse, and output files "
            "will appear here."
        )

        self._build_layout()
        self.start_button.clicked.connect(self.start_preprocessing)
        self.cancel_button.clicked.connect(self.cancel_preprocessing)
        self.open_button.clicked.connect(self.open_output_folder)

    def _build_layout(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(12)
        title = QtWidgets.QLabel("Reaction Source Data Preprocessor")
        title.setStyleSheet("font-size: 20px; font-weight: 600;")
        layout.addWidget(title)
        description = QtWidgets.QLabel(
            "Normalize each raw source CSV into one reusable, chemistry-free "
            "intermediate file. Molecular analysis and condition-registry "
            "resolution are deliberately deferred to the downstream converter."
        )
        description.setWordWrap(True)
        layout.addWidget(description)

        file_buttons = QtWidgets.QHBoxLayout()
        add_button = QtWidgets.QPushButton("Add CSV Files…")
        add_button.setObjectName("addFilesButton")
        add_button.clicked.connect(self.choose_source_files)
        add_folder_button = QtWidgets.QPushButton("Add Folder…")
        add_folder_button.setObjectName("addFolderButton")
        add_folder_button.clicked.connect(self.choose_source_folder)
        remove_button = QtWidgets.QPushButton("Remove Selected")
        remove_button.setObjectName("removeFilesButton")
        remove_button.clicked.connect(self.remove_selected_files)
        clear_button = QtWidgets.QPushButton("Clear")
        clear_button.setObjectName("clearFilesButton")
        clear_button.clicked.connect(self.clear_source_files)
        file_buttons.addWidget(add_button)
        file_buttons.addWidget(add_folder_button)
        file_buttons.addWidget(remove_button)
        file_buttons.addWidget(clear_button)
        file_buttons.addStretch()
        layout.addLayout(file_buttons)
        layout.addWidget(self.source_list, stretch=1)
        layout.addWidget(self.source_summary)

        form = QtWidgets.QFormLayout()
        form.setFieldGrowthPolicy(
            QtWidgets.QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow
        )
        output_row = QtWidgets.QHBoxLayout()
        output_row.addWidget(self.output_edit)
        output_button = QtWidgets.QPushButton("Browse…")
        output_button.clicked.connect(self.choose_output_folder)
        output_row.addWidget(output_button)
        form.addRow("Output folder:", output_row)
        form.addRow("Source adapter:", self.adapter_combo)
        form.addRow("", self.force_check)
        layout.addLayout(form)

        outputs = QtWidgets.QLabel(
            "Per source: <name>.observations.jsonl.gz. No separate JSON log is "
            "created. Unchanged artifacts are validated from their embedded "
            "provenance before reuse."
        )
        outputs.setWordWrap(True)
        outputs.setStyleSheet(
            "background: #eef4fa; border: 1px solid #ccd9e5; "
            "color: #23313f; padding: 8px; border-radius: 4px;"
        )
        layout.addWidget(outputs)

        button_row = QtWidgets.QHBoxLayout()
        button_row.addWidget(self.start_button)
        button_row.addWidget(self.cancel_button)
        button_row.addWidget(self.open_button)
        button_row.addStretch()
        layout.addLayout(button_row)
        layout.addWidget(self.progress_bar)
        layout.addWidget(QtWidgets.QLabel("Status"))
        layout.addWidget(self.status_box, stretch=1)

    def source_files(self) -> tuple[str, ...]:
        """Return selected source paths in their visible deterministic order."""
        return tuple(
            self.source_list.item(index).data(QtCore.Qt.ItemDataRole.UserRole)
            for index in range(self.source_list.count())
        )

    def add_source_files(self, paths: Iterable[str]) -> None:
        """Add source paths without duplicating an existing selection."""
        existing = {str(Path(value).resolve()) for value in self.source_files()}
        for value in paths:
            path = Path(value)
            resolved = str(path.resolve())
            if resolved in existing:
                continue
            item = QtWidgets.QListWidgetItem(path.name)
            item.setToolTip(resolved)
            item.setData(QtCore.Qt.ItemDataRole.UserRole, resolved)
            self.source_list.addItem(item)
            existing.add(resolved)
        self.refresh_source_summary()

    @QtCore.pyqtSlot()
    def choose_source_files(self) -> None:
        files, _ = QtWidgets.QFileDialog.getOpenFileNames(
            self,
            "Choose source CSV files",
            str(PROJECT_ROOT / "raw_dataset"),
            "CSV files (*.csv *.CSV)",
        )
        if files:
            self.add_source_files(files)

    @QtCore.pyqtSlot()
    def choose_source_folder(self) -> None:
        folder = QtWidgets.QFileDialog.getExistingDirectory(
            self,
            "Choose source CSV folder",
            str(PROJECT_ROOT / "raw_dataset"),
        )
        if folder:
            self.add_source_folder(folder)

    def add_source_folder(self, value: str) -> None:
        """Recursively add a folder's CSV files in deterministic order."""
        root = Path(value)
        files = sorted(
            (
                path
                for path in root.rglob("*")
                if path.is_file() and path.suffix.casefold() == ".csv"
            ),
            key=lambda path: path.relative_to(root).as_posix().casefold(),
        )
        self.add_source_files(str(path) for path in files)

    @QtCore.pyqtSlot()
    def remove_selected_files(self) -> None:
        for item in self.source_list.selectedItems():
            self.source_list.takeItem(self.source_list.row(item))
        self.refresh_source_summary()

    @QtCore.pyqtSlot()
    def clear_source_files(self) -> None:
        self.source_list.clear()
        self.refresh_source_summary()

    @QtCore.pyqtSlot()
    def choose_output_folder(self) -> None:
        folder = QtWidgets.QFileDialog.getExistingDirectory(
            self,
            "Choose intermediate-data output folder",
            self.output_edit.text() or str(PROJECT_ROOT),
        )
        if folder:
            self.output_edit.setText(folder)

    def refresh_source_summary(self) -> None:
        paths = self.source_files()
        if not paths:
            self.source_summary.setText("No source files selected.")
            return
        detected = []
        failures = 0
        for value in paths:
            try:
                detected.append(detect_adapter(value).adapter_id)
            except (OSError, ValueError):
                failures += 1
        counts: Dict[str, int] = {}
        for adapter_id in detected:
            counts[adapter_id] = counts.get(adapter_id, 0) + 1
        detail = ", ".join(
            f"{adapter_id}: {count}" for adapter_id, count in sorted(counts.items())
        )
        if failures:
            detail = (detail + "; " if detail else "") + f"unrecognized: {failures}"
        self.source_summary.setText(
            f"Selected {len(paths)} file(s). " + (detail or "No adapters detected.")
        )

    def _append_status(self, message: str) -> None:
        if message:
            self.status_box.appendPlainText(message)
            scrollbar = self.status_box.verticalScrollBar()
            scrollbar.setValue(scrollbar.maximum())

    @QtCore.pyqtSlot()
    def start_preprocessing(self) -> None:
        if self.thread is not None:
            return
        sources = self.source_files()
        if not sources:
            QtWidgets.QMessageBox.warning(
                self, "Source files required", "Add at least one source CSV file."
            )
            return
        output_text = self.output_edit.text().strip()
        if not output_text:
            QtWidgets.QMessageBox.warning(
                self, "Output required", "Choose an intermediate-data folder."
            )
            return
        output = Path(output_text)
        if any(output.resolve() == Path(value).parent.resolve() for value in sources):
            QtWidgets.QMessageBox.warning(
                self,
                "Separate output folder required",
                "The output folder must be separate from every raw source folder.",
            )
            return
        adapter_id = str(self.adapter_combo.currentData())
        self.status_box.clear()
        self._append_status(f"Selected files: {len(sources)}")
        self._append_status(f"Output folder: {output}")
        self._append_status(f"Adapter: {adapter_id}")
        self.start_button.setEnabled(False)
        self.cancel_button.setEnabled(True)
        self.open_button.setEnabled(False)
        self.progress_bar.setRange(0, len(sources))
        self.progress_bar.setValue(0)
        self.progress_bar.setFormat("Starting…")

        thread = QtCore.QThread(self)
        worker = DataPreprocessingWorker(
            sources,
            str(output),
            adapter_id=adapter_id,
            force=self.force_check.isChecked(),
        )
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.progress.connect(self._on_progress)
        worker.finished.connect(self._on_finished)
        worker.finished.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(self._clear_thread)
        self.thread = thread
        self.worker = worker
        thread.start()

    @QtCore.pyqtSlot()
    def cancel_preprocessing(self) -> None:
        if self.worker is None:
            return
        self.cancel_button.setEnabled(False)
        self._append_status("Cancellation requested; completed files remain reusable.")
        self.worker.request_cancel()

    @QtCore.pyqtSlot(object)
    def _on_progress(self, progress: PreprocessingProgress) -> None:
        if progress.phase in {"completed", "reused"}:
            self.progress_bar.setValue(progress.file_number)
        else:
            self.progress_bar.setValue(max(0, progress.file_number - 1))
        self.progress_bar.setFormat(
            f"File {progress.file_number}/{progress.file_count} • "
            f"{progress.row_count} rows"
        )
        self._append_status(progress.message)

    @QtCore.pyqtSlot(bool, object, str)
    def _on_finished(self, success: bool, report: Dict[str, Any], error: str) -> None:
        self.start_button.setEnabled(True)
        self.cancel_button.setEnabled(False)
        if not success:
            self.progress_bar.setValue(0)
            self.progress_bar.setFormat("Stopped")
            self._append_status(error)
            return
        self.progress_bar.setValue(report["file_count"])
        self.progress_bar.setFormat(
            f"Complete • {report['file_count']} files • {report['row_count']} rows"
        )
        self._completed_output = Path(report["output_dir"])
        self.open_button.setEnabled(True)
        self._append_status(
            f"Ready: {report['file_count']} intermediate file(s); "
            f"{report['reused_file_count']} reused."
        )
        for item in report["files"]:
            self._append_status(
                f"  {Path(item['output_path']).name}: "
                f"{item['output_row_count']} rows, "
                f"{_human_size(item['output_size_bytes'])}"
            )

    @QtCore.pyqtSlot()
    def open_output_folder(self) -> None:
        if self._completed_output is not None:
            QtGui.QDesktopServices.openUrl(
                QtCore.QUrl.fromLocalFile(str(self._completed_output))
            )

    @QtCore.pyqtSlot()
    def _clear_thread(self) -> None:
        self.thread = None
        self.worker = None

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        if self.thread is not None and self.thread.isRunning():
            self.cancel_preprocessing()
            QtWidgets.QMessageBox.information(
                self,
                "Preprocessing still running",
                "Cancellation was requested. Close the app after it stops.",
            )
            event.ignore()
            return
        event.accept()


def main() -> None:
    """Launch the reaction source-data preprocessor."""
    application = QtWidgets.QApplication(sys.argv)
    window = SourceDataPreprocessorWindow()
    window.show()
    raise SystemExit(application.exec())


if __name__ == "__main__":
    main()


__all__ = [
    "DEFAULT_OUTPUT_FOLDER",
    "DataPreprocessingWorker",
    "SourceDataPreprocessorWindow",
    "main",
]

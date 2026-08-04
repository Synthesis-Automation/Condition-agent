"""PyQt6 app for building scalable reaction recommendation artifacts."""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, Optional

from PyQt6 import QtCore, QtGui, QtWidgets

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT_FOLDER = PROJECT_ROOT / "datasets" / "literature"
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from condition_recommender.conversion.artifacts import (  # noqa: E402
    RecommendationArtifactBuildCancelled,
    RecommendationArtifactProgress,
    build_recommendation_artifacts,
)
from condition_recommender.conversion.input_schema import (  # noqa: E402
    discover_conversion_datasets,
)
from reactive_taxonomy import RxnMapperProvider  # noqa: E402


def _human_size(size_bytes: int) -> str:
    value = float(size_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024.0 or unit == "TB":
            return f"{value:.0f} {unit}" if unit == "B" else f"{value:.1f} {unit}"
        value /= 1024.0
    return f"{size_bytes} B"


class ReviewConversionWorker(QtCore.QObject):
    """Run the artifact workflow in a worker thread."""

    progress = QtCore.pyqtSignal(object)
    finished = QtCore.pyqtSignal(bool, object, str)

    def __init__(
        self,
        source_inputs: str | Iterable[str],
        output_folder: str,
        *,
        shard_size: int = 1_000,
        workers: int = 1,
        build_fast_index: bool = True,
        use_rxnmapper: bool = True,
    ) -> None:
        super().__init__()
        self.source_inputs = (
            source_inputs if isinstance(source_inputs, str) else tuple(source_inputs)
        )
        self.output_folder = output_folder
        self.shard_size = shard_size
        self.workers = workers
        self.build_fast_index = build_fast_index
        self.use_rxnmapper = use_rxnmapper
        self._cancel_requested = False

    def request_cancel(self) -> None:
        """Request a safe stop after the active conversion unit."""
        self._cancel_requested = True

    @QtCore.pyqtSlot()
    def run(self) -> None:
        """Build artifacts and emit a terminal result."""
        try:
            report = build_recommendation_artifacts(
                self.source_inputs,
                self.output_folder,
                shard_size=self.shard_size,
                workers=self.workers,
                build_fast_index=self.build_fast_index,
                use_rxnmapper=self.use_rxnmapper,
                progress_callback=self.progress.emit,
                cancel_check=lambda: self._cancel_requested,
            )
        except RecommendationArtifactBuildCancelled as exc:
            self.finished.emit(False, {}, str(exc))
        except Exception as exc:
            self.finished.emit(False, {}, f"{type(exc).__name__}: {exc}")
        else:
            self.finished.emit(True, report, "")


class GenericReactionReviewWindow(QtWidgets.QWidget):
    """Desktop controller for canonical conversion and review export."""

    def __init__(self) -> None:
        super().__init__()
        self.setFont(QtGui.QFont("Segoe UI", 9))
        self.setWindowTitle("Reaction Recommendation Dataset Builder")
        self.resize(880, 680)
        self.thread: Optional[QtCore.QThread] = None
        self.worker: Optional[ReviewConversionWorker] = None
        self._automatic_output = str(DEFAULT_OUTPUT_FOLDER)
        self._completed_output: Optional[Path] = None

        cpu_count = max(1, os.cpu_count() or 1)
        self.source_list = QtWidgets.QListWidget()
        self.source_list.setObjectName("sourceInputs")
        self.source_list.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection
        )
        self.source_list.setAlternatingRowColors(True)
        self.output_edit = QtWidgets.QLineEdit()
        self.output_edit.setObjectName("outputFolder")
        self.output_edit.setPlaceholderText(
            "Folder for canonical data, review CSV, and index"
        )
        self.output_edit.setText(self._automatic_output)
        self.source_summary = QtWidgets.QLabel("No source files or folders selected.")
        self.source_summary.setWordWrap(True)

        self.shard_size_spin = QtWidgets.QSpinBox()
        self.shard_size_spin.setObjectName("shardSize")
        self.shard_size_spin.setRange(100, 100_000)
        self.shard_size_spin.setSingleStep(100)
        self.shard_size_spin.setValue(1_000)
        self.shard_size_spin.setToolTip(
            "Smaller shards checkpoint more often; larger shards reduce "
            "file overhead but use more memory and take longer to cancel."
        )
        self.worker_count_spin = QtWidgets.QSpinBox()
        self.worker_count_spin.setObjectName("workerCount")
        self.worker_count_spin.setRange(1, cpu_count)
        self._parallel_worker_count = min(4, cpu_count)
        self.worker_count_spin.setValue(self._parallel_worker_count)
        self.worker_count_spin.setToolTip(
            "Parallel chemistry workers. More workers can be faster but use "
            "more memory."
        )
        self.build_index_check = QtWidgets.QCheckBox(
            "Build compressed fast-load recommendation index"
        )
        self.build_index_check.setObjectName("buildFastIndex")
        self.build_index_check.setChecked(True)
        self.build_index_check.setToolTip(
            "Recommended for repeated use. It takes extra conversion time and "
            "disk space, but avoids rebuilding lookup maps when recommending."
        )
        self.use_rxnmapper_check = QtWidgets.QCheckBox(
            "Use RXNMapper for unresolved, ambiguous, or missing-core reactions"
        )
        self.use_rxnmapper_check.setObjectName("useRxnMapper")
        self.use_rxnmapper_check.setChecked(True)
        self.use_rxnmapper_check.setToolTip(
            "Checked by default. RXNMapper proposals are reconciled with "
            "internal evidence; resolved reactions retain their internal "
            "signature while consensus mapping supplies only a missing core. "
            "One worker is used to avoid loading multiple model copies."
        )
        self.use_rxnmapper_check.toggled.connect(self._sync_rxnmapper_worker_setting)
        self._sync_rxnmapper_worker_setting(True)

        self.start_button = QtWidgets.QPushButton("Generate Recommendation Data")
        self.start_button.setObjectName("generateButton")
        self.start_button.setDefault(True)
        self.start_button.setStyleSheet(
            "QPushButton#generateButton {"
            "background-color: #0078d7;"
            "color: white;"
            "font-weight: 700;"
            "padding: 10px 18px;"
            "border-radius: 6px;"
            "}"
            "QPushButton#generateButton:disabled {"
            "background-color: #a6c8f0;"
            "color: #ffffff;"
            "}"
        )
        self.cancel_button = QtWidgets.QPushButton("Cancel")
        self.cancel_button.setObjectName("cancelButton")
        self.cancel_button.setEnabled(False)
        self.open_button = QtWidgets.QPushButton("Open Output Folder")
        self.open_button.setObjectName("openOutputButton")
        self.open_button.setEnabled(False)
        self.progress_bar = QtWidgets.QProgressBar()
        self.progress_bar.setObjectName("conversionProgress")
        self.progress_bar.setRange(0, 1)
        self.progress_bar.setValue(0)
        self.progress_bar.setFormat("Waiting")
        self.status_box = QtWidgets.QPlainTextEdit()
        self.status_box.setObjectName("statusBox")
        self.status_box.setReadOnly(True)
        self.status_box.setPlaceholderText(
            "Conversion progress, resume information, and output sizes will "
            "appear here."
        )

        self._build_layout()
        self.start_button.clicked.connect(self.start_conversion)
        self.cancel_button.clicked.connect(self.cancel_conversion)
        self.open_button.clicked.connect(self.open_output_folder)

    def _build_layout(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(12)

        title = QtWidgets.QLabel("Reaction Recommendation Dataset Builder")
        title.setStyleSheet("font-size: 20px; font-weight: 600;")
        layout.addWidget(title)
        description = QtWidgets.QLabel(
            "Convert selected raw CSV or preprocessed observation files, or "
            "every supported file in selected folder trees, then produce "
            "compressed recommendation data and a concise review CSV from the "
            "same canonical records. Interrupted conversions can reuse "
            "completed shards."
        )
        description.setWordWrap(True)
        layout.addWidget(description)

        source_buttons = QtWidgets.QHBoxLayout()
        add_files_button = QtWidgets.QPushButton("Add Files…")
        add_files_button.setObjectName("addSourceFilesButton")
        add_files_button.clicked.connect(self.choose_source_files)
        add_folder_button = QtWidgets.QPushButton("Add Folder…")
        add_folder_button.setObjectName("addSourceFolderButton")
        add_folder_button.clicked.connect(self.choose_source_folder)
        remove_button = QtWidgets.QPushButton("Remove Selected")
        remove_button.setObjectName("removeSourceInputsButton")
        remove_button.clicked.connect(self.remove_selected_inputs)
        clear_button = QtWidgets.QPushButton("Clear")
        clear_button.setObjectName("clearSourceInputsButton")
        clear_button.clicked.connect(self.clear_source_inputs)
        source_buttons.addWidget(add_files_button)
        source_buttons.addWidget(add_folder_button)
        source_buttons.addWidget(remove_button)
        source_buttons.addWidget(clear_button)
        source_buttons.addStretch()
        layout.addLayout(source_buttons)
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

        settings_row = QtWidgets.QHBoxLayout()
        settings_row.addWidget(QtWidgets.QLabel("Rows per shard"))
        settings_row.addWidget(self.shard_size_spin)
        settings_row.addSpacing(16)
        settings_row.addWidget(QtWidgets.QLabel("Parallel workers"))
        settings_row.addWidget(self.worker_count_spin)
        settings_row.addStretch()
        form.addRow("Performance:", settings_row)
        form.addRow("", self.use_rxnmapper_check)
        form.addRow("", self.build_index_check)
        layout.addLayout(form)

        outputs = QtWidgets.QLabel(
            "Outputs: shard_manifest.json + compressed shards (canonical "
            "recommendation data and restart checkpoints) • "
            "reaction_review.csv (concise human review, including spectator "
            "and local steric/electronic context) • "
            "generic_index.sqlite (trusted lazy lookup) + "
            "generic_review_index.sqlite (expert review-core lookup), with "
            "temporary JSON parity indexes. "
            "A duplicate merged "
            "records file is not stored."
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

    @QtCore.pyqtSlot()
    def choose_source_folder(self) -> None:
        folder = QtWidgets.QFileDialog.getExistingDirectory(
            self,
            "Choose raw or preprocessed dataset folder",
            str(PROJECT_ROOT),
        )
        if folder:
            self.add_source_inputs((folder,))

    @QtCore.pyqtSlot()
    def choose_source_files(self) -> None:
        files, _ = QtWidgets.QFileDialog.getOpenFileNames(
            self,
            "Choose raw or preprocessed reaction files",
            str(PROJECT_ROOT),
            (
                "Conversion inputs (*.csv *.CSV *.observations.jsonl.gz);;"
                "CSV files (*.csv *.CSV);;All files (*)"
            ),
        )
        if files:
            self.add_source_inputs(files)

    def source_inputs(self) -> tuple[str, ...]:
        """Return the selected files and folders in visible order."""
        return tuple(
            self.source_list.item(index).data(QtCore.Qt.ItemDataRole.UserRole)
            for index in range(self.source_list.count())
        )

    def add_source_inputs(self, paths: Iterable[str]) -> None:
        """Add files or folders without duplicating an existing selection."""
        existing = {str(Path(value).resolve()).casefold() for value in self.source_inputs()}
        for value in paths:
            path = Path(value)
            resolved = str(path.resolve())
            if resolved.casefold() in existing:
                continue
            label = f"[Folder] {path.name}" if path.is_dir() else path.name
            item = QtWidgets.QListWidgetItem(label)
            item.setToolTip(resolved)
            item.setData(QtCore.Qt.ItemDataRole.UserRole, resolved)
            self.source_list.addItem(item)
            existing.add(resolved.casefold())
        self.refresh_source_summary()

    @QtCore.pyqtSlot()
    def remove_selected_inputs(self) -> None:
        for item in self.source_list.selectedItems():
            self.source_list.takeItem(self.source_list.row(item))
        self.refresh_source_summary()

    @QtCore.pyqtSlot()
    def clear_source_inputs(self) -> None:
        self.source_list.clear()
        self.refresh_source_summary()

    @QtCore.pyqtSlot()
    def choose_output_folder(self) -> None:
        output = QtWidgets.QFileDialog.getExistingDirectory(
            self,
            "Choose recommendation output folder",
            self.output_edit.text() or str(DEFAULT_OUTPUT_FOLDER),
        )
        if output:
            self.output_edit.setText(output)
            self._automatic_output = ""

    @QtCore.pyqtSlot()
    def refresh_source_summary(self) -> None:
        selections = self.source_inputs()
        paths = discover_conversion_datasets(selections)
        if not paths:
            self.source_summary.setText(
                "No raw CSV or preprocessed observation files were found."
            )
            return
        folder_count = sum(Path(value).is_dir() for value in selections)
        file_count = len(selections) - folder_count
        self.source_summary.setText(
            f"Found {len(paths)} conversion input file(s) from "
            f"{file_count} selected file(s) and {folder_count} folder(s)."
        )

    def _append_status(self, message: str) -> None:
        if message:
            self.status_box.appendPlainText(message)
            scrollbar = self.status_box.verticalScrollBar()
            scrollbar.setValue(scrollbar.maximum())

    @QtCore.pyqtSlot(bool)
    def _sync_rxnmapper_worker_setting(self, enabled: bool) -> None:
        """Use one process for RXNMapper and restore parallelism when disabled."""
        if enabled:
            if self.worker_count_spin.isEnabled():
                self._parallel_worker_count = self.worker_count_spin.value()
            self.worker_count_spin.setValue(1)
            self.worker_count_spin.setEnabled(False)
        else:
            self.worker_count_spin.setEnabled(True)
            self.worker_count_spin.setValue(self._parallel_worker_count)

    @QtCore.pyqtSlot()
    def start_conversion(self) -> None:
        if self.thread is not None:
            return
        source_inputs = self.source_inputs()
        output_text = self.output_edit.text().strip()
        if not source_inputs:
            QtWidgets.QMessageBox.warning(
                self,
                "Conversion inputs required",
                "Add at least one supported file or folder.",
            )
            return
        discovered = discover_conversion_datasets(source_inputs)
        if not discovered:
            QtWidgets.QMessageBox.warning(
                self,
                "No conversion inputs",
                "No supported input files were found in the selection.",
            )
            return
        if not output_text:
            QtWidgets.QMessageBox.warning(
                self,
                "Output required",
                "Choose an output folder.",
            )
            return
        if (
            self.use_rxnmapper_check.isChecked()
            and not RxnMapperProvider.is_available()
        ):
            QtWidgets.QMessageBox.warning(
                self,
                "RXNMapper unavailable",
                "RXNMapper is checked but not installed. Run "
                "'python -m pip install -r requirements-mapping.txt' or "
                "clear the RXNMapper checkbox.",
            )
            return
        output = Path(output_text)
        output_resolved = output.resolve()
        invalid_folders = [
            Path(value)
            for value in source_inputs
            if Path(value).is_dir()
            and (
                output_resolved == Path(value).resolve()
                or Path(value).resolve() in output_resolved.parents
            )
        ]
        if invalid_folders:
            QtWidgets.QMessageBox.warning(
                self,
                "Separate output folder required",
                "Choose an output folder outside every selected source folder.",
            )
            return

        self.status_box.clear()
        self._append_status(
            f"Selected inputs: {len(source_inputs)}; discovered files: "
            f"{len(discovered)}"
        )
        self._append_status(f"Output folder: {output}")
        self._append_status(
            f"Settings: {self.shard_size_spin.value()} rows/shard, "
            f"{self.worker_count_spin.value()} worker(s), "
            f"RXNMapper "
            f"{'on' if self.use_rxnmapper_check.isChecked() else 'off'}, "
            f"fast index {'on' if self.build_index_check.isChecked() else 'off'}"
        )
        self.start_button.setEnabled(False)
        self.cancel_button.setEnabled(True)
        self.open_button.setEnabled(False)
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setFormat("Discovering files…")

        thread = QtCore.QThread(self)
        worker = ReviewConversionWorker(
            source_inputs,
            str(output),
            shard_size=self.shard_size_spin.value(),
            workers=self.worker_count_spin.value(),
            build_fast_index=self.build_index_check.isChecked(),
            use_rxnmapper=self.use_rxnmapper_check.isChecked(),
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
    def cancel_conversion(self) -> None:
        if self.worker is None:
            return
        self.cancel_button.setEnabled(False)
        self._append_status(
            "Cancellation requested. The active shard will finish; completed "
            "shards remain reusable when this output folder is selected again."
        )
        self.worker.request_cancel()

    @QtCore.pyqtSlot(object)
    def _on_progress(self, progress: RecommendationArtifactProgress) -> None:
        self.progress_bar.setRange(0, 0)
        if progress.phase == "canonical_shard_completed":
            self.progress_bar.setFormat(
                f"{progress.shard_count} shard(s) • {progress.row_count} reactions"
            )
        elif progress.phase.startswith("review"):
            self.progress_bar.setFormat(
                f"Writing review CSV • {progress.row_count} reactions"
            )
        elif progress.phase.startswith("index"):
            self.progress_bar.setFormat(
                f"Building fast index • {progress.row_count} reactions"
            )
        else:
            self.progress_bar.setFormat(progress.message)
        self._append_status(progress.message)

    @QtCore.pyqtSlot(bool, object, str)
    def _on_finished(
        self,
        success: bool,
        report: Dict[str, Any],
        error: str,
    ) -> None:
        self.start_button.setEnabled(True)
        self.cancel_button.setEnabled(False)
        if success:
            self.progress_bar.setRange(0, 1)
            self.progress_bar.setValue(1)
            self.progress_bar.setFormat(
                f"Complete • {report['record_count']} reactions"
            )
            self._completed_output = Path(report["output_dir"])
            self.open_button.setEnabled(True)
            artifacts = report["artifacts"]
            self._append_status("Ready:")
            self._append_status(
                "  Canonical data manifest: "
                f"{artifacts['canonical_manifest']['path']} "
                f"({_human_size(artifacts['canonical_manifest']['size_bytes'])})"
            )
            self._append_status(
                "  Review CSV: "
                f"{artifacts['review_csv']['path']} "
                f"({_human_size(artifacts['review_csv']['size_bytes'])})"
            )
            if "fast_index" in artifacts:
                self._append_status(
                    "  Fast recommendation index: "
                    f"{artifacts['fast_index']['path']} "
                    f"({_human_size(artifacts['fast_index']['size_bytes'])})"
                )
            self._append_status(
                f"  Shards: {report['shard_count']} "
                f"({_human_size(report['storage']['shard_size_bytes'])}); "
                f"reused this run: {report['reused_shard_count']}"
            )
            if report["eligible_index_record_count"] is not None:
                self._append_status(
                    "  Trusted recommendation precedents: "
                    f"{report['eligible_index_record_count']}"
                )
            if "review_core_index" in artifacts:
                self._append_status(
                    "  Review-core recommendation index: "
                    f"{artifacts['review_core_index']['path']} "
                    f"({_human_size(artifacts['review_core_index']['size_bytes'])})"
                )
            self._append_status(
                "  Qualified review-core precedents: "
                f"{report.get('review_core_precedent_count', 0)}"
            )
            self._append_status(
                "  Query-core eligible reactions: "
                f"{report.get('query_core_eligible_count', 0)}"
            )
            for warning in report.get("warnings") or ():
                self._append_status(f"Warning: {warning}")
        else:
            self.progress_bar.setRange(0, 1)
            self.progress_bar.setValue(0)
            self.progress_bar.setFormat("Stopped")
            self._append_status(error)

    @QtCore.pyqtSlot()
    def open_output_folder(self) -> None:
        if self._completed_output is None:
            return
        QtGui.QDesktopServices.openUrl(
            QtCore.QUrl.fromLocalFile(str(self._completed_output))
        )

    @QtCore.pyqtSlot()
    def _clear_thread(self) -> None:
        self.thread = None
        self.worker = None

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        if self.thread is not None and self.thread.isRunning():
            self.cancel_conversion()
            QtWidgets.QMessageBox.information(
                self,
                "Conversion still running",
                "Cancellation was requested. Close the app after it stops.",
            )
            event.ignore()
            return
        event.accept()


def main() -> None:
    """Launch the recommendation dataset builder."""
    application = QtWidgets.QApplication(sys.argv)
    window = GenericReactionReviewWindow()
    window.show()
    raise SystemExit(application.exec())


if __name__ == "__main__":
    main()


__all__ = [
    "DEFAULT_OUTPUT_FOLDER",
    "GenericReactionReviewWindow",
    "ReviewConversionWorker",
    "main",
]

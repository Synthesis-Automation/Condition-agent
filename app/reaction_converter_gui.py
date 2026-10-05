"""Desktop builder for the complete processed reaction dataset release."""

from __future__ import annotations

import json
import os
from pathlib import Path
import sys
from typing import Any

from PyQt6 import QtCore, QtGui, QtWidgets

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
DEFAULT_INPUT_FOLDER = PROJECT_ROOT / "datasets" / "intermediate_datasets"
DEFAULT_OUTPUT_FOLDER = PROJECT_ROOT / "datasets" / "processed_datasets"

from app.processed_dataset_builder import build_processed_datasets  # noqa: E402
from condition_recommender.conversion.sharded import ShardedConversionCancelled  # noqa: E402
from condition_recommender.processed_release import resolve_processed_release  # noqa: E402


class DatasetBuildWorker(QtCore.QObject):
    """Run the common builder with resumable cancellation in a Qt worker thread."""

    progress = QtCore.pyqtSignal(object)
    finished = QtCore.pyqtSignal(bool, object, str)

    def __init__(self, source: str, output: str, *, workers: int, shard_size: int) -> None:
        super().__init__()
        self.source = source
        self.output = output
        self.workers = workers
        self.shard_size = shard_size
        self.cancelled = False

    def request_cancel(self) -> None:
        """Stop after the current safe unit, retaining completed work."""
        self.cancelled = True

    @QtCore.pyqtSlot()
    def run(self) -> None:
        """Build every required artifact and publish only after validation."""
        try:
            report = build_processed_datasets(
                self.source, self.output, workers=self.workers, shard_size=self.shard_size,
                progress=self.progress.emit, cancel_check=lambda: self.cancelled,
            )
        except (InterruptedError, ShardedConversionCancelled) as exc:
            self.finished.emit(False, {"cancelled": True}, str(exc))
        except Exception as exc:
            self.finished.emit(False, {}, f"{type(exc).__name__}: {exc}")
        else:
            self.finished.emit(True, report, "")


class GenericReactionReviewWindow(QtWidgets.QWidget):
    """Configure one complete corpus build and inspect its progress."""

    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Processed Reaction Dataset Builder")
        self.resize(850, 620)
        self.thread: QtCore.QThread | None = None
        self.worker: DatasetBuildWorker | None = None
        self.source_edit = QtWidgets.QLineEdit(str(DEFAULT_INPUT_FOLDER))
        self.source_edit.setObjectName("sourceFolder")
        self.output_edit = QtWidgets.QLineEdit(str(DEFAULT_OUTPUT_FOLDER))
        self.output_edit.setObjectName("outputFolder")
        self.worker_count_spin = QtWidgets.QSpinBox()
        self.worker_count_spin.setRange(1, max(1, os.cpu_count() or 1))
        self.worker_count_spin.setValue(min(6, os.cpu_count() or 1))
        self.shard_size_spin = QtWidgets.QSpinBox()
        self.shard_size_spin.setRange(100, 5000)
        self.shard_size_spin.setSingleStep(100)
        self.shard_size_spin.setValue(500)
        self.start_button = QtWidgets.QPushButton("Build / Resume All Datasets")
        self.start_button.setObjectName("generateButton")
        self.cancel_button = QtWidgets.QPushButton("Cancel")
        self.cancel_button.setEnabled(False)
        self.validate_button = QtWidgets.QPushButton("Validate Published Release")
        self.open_button = QtWidgets.QPushButton("Open Output Folder")
        self.status_label = QtWidgets.QLabel("Ready. All physical observations will be processed.")
        self.status_label.setWordWrap(True)
        self.log = QtWidgets.QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(3000)
        self.progress_bar = QtWidgets.QProgressBar()
        self.progress_bar.setRange(0, 1)
        self.progress_bar.setValue(0)
        layout = QtWidgets.QVBoxLayout(self)
        explanation = QtWidgets.QLabel(
            "Build one complete release for the app and agent: condition and shared-core "
            "indexes, fragment search on products and starting materials, source evidence, "
            "and forward/retrosynthesis libraries. Completed stages are reused on resume."
        )
        explanation.setWordWrap(True)
        layout.addWidget(explanation)
        form = QtWidgets.QFormLayout()
        for label, edit, default in (("Intermediate input", self.source_edit, DEFAULT_INPUT_FOLDER),
                                     ("Processed output", self.output_edit, DEFAULT_OUTPUT_FOLDER)):
            row = QtWidgets.QWidget()
            row_layout = QtWidgets.QHBoxLayout(row)
            row_layout.setContentsMargins(0, 0, 0, 0)
            browse = QtWidgets.QPushButton("Browse")
            browse.clicked.connect(lambda _, e=edit, d=default: self._browse(e, d))
            row_layout.addWidget(edit)
            row_layout.addWidget(browse)
            form.addRow(label, row)
        form.addRow("Chemistry workers", self.worker_count_spin)
        form.addRow("Observations per shard", self.shard_size_spin)
        layout.addLayout(form)
        actions = QtWidgets.QHBoxLayout()
        for button in (self.start_button, self.cancel_button, self.validate_button, self.open_button):
            actions.addWidget(button)
        layout.addLayout(actions)
        layout.addWidget(self.status_label)
        layout.addWidget(self.progress_bar)
        layout.addWidget(self.log, 1)
        self.start_button.clicked.connect(self.start_conversion)
        self.cancel_button.clicked.connect(self.cancel_conversion)
        self.validate_button.clicked.connect(self.validate_release)
        self.open_button.clicked.connect(self.open_output_folder)

    def _browse(self, edit: QtWidgets.QLineEdit, default: Path) -> None:
        selected = QtWidgets.QFileDialog.getExistingDirectory(self, "Select folder", edit.text() or str(default))
        if selected:
            edit.setText(selected)

    def start_conversion(self) -> None:
        """Start or resume the single production build workflow."""
        if self.thread:
            return
        if not Path(self.source_edit.text()).exists():
            self.status_label.setText("Intermediate input does not exist.")
            return
        self.thread = QtCore.QThread(self)
        self.worker = DatasetBuildWorker(self.source_edit.text(), self.output_edit.text(),
                                          workers=self.worker_count_spin.value(),
                                          shard_size=self.shard_size_spin.value())
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.progress.connect(self._progress)
        self.worker.finished.connect(self._finished)
        self.worker.finished.connect(self.thread.quit)
        self.thread.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self._thread_finished)
        self.start_button.setEnabled(False)
        self.validate_button.setEnabled(False)
        self.cancel_button.setEnabled(True)
        self.source_edit.setEnabled(False)
        self.output_edit.setEnabled(False)
        self.progress_bar.setRange(0, 0)
        self.thread.start()

    def _progress(self, event: dict[str, Any]) -> None:
        message = event.get("message") or str(event.get("phase"))
        if "rows" in event:
            message += f" — {event['rows']:,} observations"
        self.status_label.setText(message)
        self.log.appendPlainText(json.dumps(event, ensure_ascii=False))

    def _finished(self, success: bool, report: dict[str, Any], error: str) -> None:
        self.progress_bar.setRange(0, 1)
        self.progress_bar.setValue(int(success))
        if success:
            count = report["coverage"]["observation_count"]
            self.status_label.setText(f"Published {count:,} observations. Release {report['release_id']}.")
        else:
            self.status_label.setText(error + " Run Build / Resume to continue from completed work.")
        self.log.appendPlainText(self.status_label.text())

    def _thread_finished(self) -> None:
        if self.thread:
            self.thread.deleteLater()
        self.thread = None
        self.worker = None
        self.start_button.setEnabled(True)
        self.validate_button.setEnabled(True)
        self.cancel_button.setEnabled(False)
        self.source_edit.setEnabled(True)
        self.output_edit.setEnabled(True)

    def cancel_conversion(self) -> None:
        """Request cancellation without deleting validated checkpoints."""
        if self.worker:
            self.worker.request_cancel()
            self.cancel_button.setEnabled(False)
            self.status_label.setText("Stopping at the next checkpoint; completed work will be retained.")

    def validate_release(self) -> None:
        """Verify all declared artifacts of the published release."""
        try:
            release = resolve_processed_release(self.output_edit.text(), verify_artifacts=True)
        except Exception as exc:
            self.status_label.setText(f"Validation failed: {exc}")
        else:
            self.status_label.setText(f"Validated release {release.manifest['release_id']}.")

    def open_output_folder(self) -> None:
        """Open the configured processed dataset root."""
        QtGui.QDesktopServices.openUrl(QtCore.QUrl.fromLocalFile(str(Path(self.output_edit.text()).resolve())))

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        if self.thread:
            self.cancel_conversion()
            event.ignore()
        else:
            event.accept()


def main() -> int:
    """Start the standalone desktop builder."""
    application = QtWidgets.QApplication(sys.argv)
    window = GenericReactionReviewWindow()
    window.show()
    return application.exec()


if __name__ == "__main__":
    raise SystemExit(main())

"""PyQt6 app for recursive dataset conversion to concise chemistry review CSV."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict, Optional

from PyQt6 import QtCore, QtGui, QtWidgets

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from condition_recommender.conversion.concise_review import (  # noqa: E402
    ConciseReviewConversionCancelled,
    ConciseReviewProgress,
    convert_dataset_folder_to_concise_review_csv,
)
from condition_recommender.conversion.input_schema import (  # noqa: E402
    discover_csv_datasets,
)


class ReviewConversionWorker(QtCore.QObject):
    """Run chemistry conversion in a worker thread."""

    progress = QtCore.pyqtSignal(object)
    finished = QtCore.pyqtSignal(bool, object, str)

    def __init__(self, source_folder: str, output_path: str) -> None:
        super().__init__()
        self.source_folder = source_folder
        self.output_path = output_path
        self._cancel_requested = False

    def request_cancel(self) -> None:
        """Request cancellation after the current reaction finishes."""
        self._cancel_requested = True

    @QtCore.pyqtSlot()
    def run(self) -> None:
        """Convert the selected folder and emit a terminal result."""
        try:
            report = convert_dataset_folder_to_concise_review_csv(
                self.source_folder,
                self.output_path,
                progress_callback=self.progress.emit,
                cancel_check=lambda: self._cancel_requested,
            )
        except ConciseReviewConversionCancelled as exc:
            self.finished.emit(False, {}, str(exc))
        except Exception as exc:
            self.finished.emit(
                False,
                {},
                f"{type(exc).__name__}: {exc}",
            )
        else:
            self.finished.emit(True, report, "")


class GenericReactionReviewWindow(QtWidgets.QWidget):
    """Folder-to-review-CSV interface for the clean recommendation system."""

    def __init__(self) -> None:
        super().__init__()
        self.setFont(QtGui.QFont("Segoe UI", 9))
        self.setWindowTitle("Reaction Family Review CSV Generator")
        self.resize(820, 560)
        self.thread: Optional[QtCore.QThread] = None
        self.worker: Optional[ReviewConversionWorker] = None
        self._automatic_output = ""

        self.source_edit = QtWidgets.QLineEdit()
        self.source_edit.setObjectName("sourceFolder")
        self.source_edit.setPlaceholderText(
            "Folder containing reaction CSV files; subfolders are included"
        )
        self.output_edit = QtWidgets.QLineEdit()
        self.output_edit.setObjectName("outputCsv")
        self.output_edit.setPlaceholderText("Destination review CSV")
        self.source_summary = QtWidgets.QLabel("No source folder selected.")
        self.source_summary.setWordWrap(True)

        self.start_button = QtWidgets.QPushButton("Generate Review CSV")
        self.start_button.setObjectName("generateButton")
        self.cancel_button = QtWidgets.QPushButton("Cancel")
        self.cancel_button.setObjectName("cancelButton")
        self.cancel_button.setEnabled(False)
        self.progress_bar = QtWidgets.QProgressBar()
        self.progress_bar.setObjectName("conversionProgress")
        self.progress_bar.setRange(0, 1)
        self.progress_bar.setValue(0)
        self.progress_bar.setFormat("Waiting")
        self.status_box = QtWidgets.QPlainTextEdit()
        self.status_box.setObjectName("statusBox")
        self.status_box.setReadOnly(True)
        self.status_box.setPlaceholderText(
            "Conversion progress and errors will appear here."
        )

        self._build_layout()
        self.source_edit.editingFinished.connect(self.refresh_source_summary)
        self.start_button.clicked.connect(self.start_conversion)
        self.cancel_button.clicked.connect(self.cancel_conversion)

    def _build_layout(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(12)

        title = QtWidgets.QLabel("Concise Reaction-Family Review")
        title.setStyleSheet("font-size: 20px; font-weight: 600;")
        layout.addWidget(title)
        description = QtWidgets.QLabel(
            "Select a dataset folder. Every CSV in that folder and its "
            "subfolders will be converted into one five-column review file."
        )
        description.setWordWrap(True)
        layout.addWidget(description)

        form = QtWidgets.QFormLayout()
        form.setFieldGrowthPolicy(
            QtWidgets.QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow
        )
        source_row = QtWidgets.QHBoxLayout()
        source_row.addWidget(self.source_edit)
        source_button = QtWidgets.QPushButton("Browse…")
        source_button.clicked.connect(self.choose_source_folder)
        source_row.addWidget(source_button)
        form.addRow("Dataset folder:", source_row)

        output_row = QtWidgets.QHBoxLayout()
        output_row.addWidget(self.output_edit)
        output_button = QtWidgets.QPushButton("Browse…")
        output_button.clicked.connect(self.choose_output_csv)
        output_row.addWidget(output_button)
        form.addRow("Review CSV:", output_row)
        layout.addLayout(form)
        layout.addWidget(self.source_summary)

        columns = QtWidgets.QLabel(
            "Output columns: canonical reaction SMILES • detailed structural "
            "label • original reaction type • detected reaction family • status"
        )
        columns.setWordWrap(True)
        columns.setStyleSheet(
            "background: #eef4fa; border: 1px solid #ccd9e5; "
            "color: #23313f; padding: 8px; border-radius: 4px;"
        )
        layout.addWidget(columns)

        button_row = QtWidgets.QHBoxLayout()
        button_row.addWidget(self.start_button)
        button_row.addWidget(self.cancel_button)
        button_row.addStretch()
        layout.addLayout(button_row)
        layout.addWidget(self.progress_bar)
        layout.addWidget(QtWidgets.QLabel("Status"))
        layout.addWidget(self.status_box, stretch=1)

    @QtCore.pyqtSlot()
    def choose_source_folder(self) -> None:
        folder = QtWidgets.QFileDialog.getExistingDirectory(
            self,
            "Choose reaction dataset folder",
            self.source_edit.text() or str(PROJECT_ROOT),
        )
        if not folder:
            return
        self.source_edit.setText(folder)
        suggested = str(
            Path(folder).parent / f"{Path(folder).name}_reaction_review.csv"
        )
        if not self.output_edit.text() or (
            self.output_edit.text() == self._automatic_output
        ):
            self.output_edit.setText(suggested)
            self._automatic_output = suggested
        self.refresh_source_summary()

    @QtCore.pyqtSlot()
    def choose_output_csv(self) -> None:
        output, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save concise reaction review",
            self.output_edit.text() or str(PROJECT_ROOT / "reaction_review.csv"),
            "CSV files (*.csv)",
        )
        if output:
            if not output.lower().endswith(".csv"):
                output += ".csv"
            self.output_edit.setText(output)
            self._automatic_output = ""

    @QtCore.pyqtSlot()
    def refresh_source_summary(self) -> None:
        source = Path(self.source_edit.text().strip())
        paths = discover_csv_datasets(source)
        if not paths:
            self.source_summary.setText(
                "No CSV dataset files were found in this folder."
            )
            return
        self.source_summary.setText(
            f"Found {len(paths)} CSV file(s), including subfolders."
        )

    def _append_status(self, message: str) -> None:
        if message:
            self.status_box.appendPlainText(message)
            scrollbar = self.status_box.verticalScrollBar()
            scrollbar.setValue(scrollbar.maximum())

    @QtCore.pyqtSlot()
    def start_conversion(self) -> None:
        if self.thread is not None:
            return
        source = Path(self.source_edit.text().strip())
        output_text = self.output_edit.text().strip()
        if not source.is_dir():
            QtWidgets.QMessageBox.warning(
                self,
                "Dataset folder required",
                "Choose an existing folder containing reaction CSV files.",
            )
            return
        if not discover_csv_datasets(source):
            QtWidgets.QMessageBox.warning(
                self,
                "No CSV files",
                "No CSV files were found in the selected folder or subfolders.",
            )
            return
        if not output_text:
            QtWidgets.QMessageBox.warning(
                self,
                "Output required",
                "Choose where to save the review CSV.",
            )
            return
        output = Path(output_text)
        if output.suffix.lower() != ".csv":
            output = output.with_suffix(".csv")
            self.output_edit.setText(str(output))

        self.status_box.clear()
        self._append_status(f"Source: {source}")
        self._append_status(f"Output: {output}")
        self.start_button.setEnabled(False)
        self.cancel_button.setEnabled(True)
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setFormat("Discovering files…")

        thread = QtCore.QThread(self)
        worker = ReviewConversionWorker(str(source), str(output))
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
            "Cancellation requested; finishing the current reaction…"
        )
        self.worker.request_cancel()

    @QtCore.pyqtSlot(object)
    def _on_progress(self, progress: ConciseReviewProgress) -> None:
        if progress.file_count:
            self.progress_bar.setRange(0, progress.file_count)
            self.progress_bar.setValue(
                min(progress.file_index, progress.file_count)
            )
            self.progress_bar.setFormat(
                f"File %v/%m • {progress.row_count} reactions"
            )
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
            self.progress_bar.setValue(self.progress_bar.maximum())
            self.progress_bar.setFormat(
                f"Complete • {report['row_count']} reactions"
            )
            self._append_status(
                f"Saved {report['row_count']} reaction(s) to "
                f"{report['output_path']}"
            )
        else:
            self.progress_bar.setFormat("Stopped")
            self._append_status(error)

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
    """Launch the concise reaction review application."""
    application = QtWidgets.QApplication(sys.argv)
    window = GenericReactionReviewWindow()
    window.show()
    raise SystemExit(application.exec())


if __name__ == "__main__":
    main()


__all__ = [
    "GenericReactionReviewWindow",
    "ReviewConversionWorker",
    "main",
]

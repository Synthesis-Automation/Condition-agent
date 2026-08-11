"""PyQt6 app for building scalable reaction recommendation artifacts."""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, Optional

from PyQt6 import QtCore, QtGui, QtWidgets

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT_FOLDER = PROJECT_ROOT / "datasets" / "literature"
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from condition_recommender.conversion.artifacts import (  # noqa: E402
    RecommendationArtifactBuildCancelled,
    RecommendationArtifactProgress,
    build_recommendation_artifacts,
    combine_saved_recommendation_batches,
    discover_saved_conversion_batches,
    incomplete_saved_conversion_batches,
    recommendation_library_mode_dir,
    resume_saved_conversion_batch,
    save_recommendation_batch,
)
from condition_recommender.conversion.input_schema import (  # noqa: E402
    discover_conversion_datasets,
)
from core_retrosynthesis_poc import (  # noqa: E402
    FullScaleBuildConfig,
    build_full_scale_operator_library,
)
from core_retrosynthesis_poc.sources import (  # noqa: E402
    COMBINED_BATCH_MANIFEST_FILENAME,
)
from reactive_taxonomy import RxnMapperProvider  # noqa: E402


def _human_size(size_bytes: int) -> str:
    value = float(size_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024.0 or unit == "TB":
            return f"{value:.0f} {unit}" if unit == "B" else f"{value:.1f} {unit}"
        value /= 1024.0
    return f"{size_bytes} B"


DEFAULT_RETROSYNTHESIS_OUTPUT_FOLDER = Path(
    os.environ.get(
        "CORE_RETROSYNTHESIS_LIBRARY_ROOT",
        str(
            PROJECT_ROOT
            / "results"
            / "operator_retrosynthesis_poc"
            / "full_scale_v3"
        ),
    )
)


def build_retrosynthesis_operator_artifacts(
    source_library_folder: str | Path,
    output_root: str | Path,
    *,
    library_mode: str,
    workers: int = 1,
    progress_callback: Optional[
        Callable[[RecommendationArtifactProgress], None]
    ] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
) -> Dict[str, Any]:
    """Build one resumable operator library from every combined record."""

    mode = library_mode.strip().casefold()
    if mode not in {"compact", "full"}:
        raise ValueError(f"Unsupported retrosynthesis library mode: {library_mode}")
    output = Path(output_root) / mode

    def emit(update: Dict[str, Any]) -> None:
        if cancel_check is not None and cancel_check():
            raise RecommendationArtifactBuildCancelled(
                "Retrosynthesis conversion cancelled. Completed operator "
                "shards remain reusable."
            )
        if progress_callback is None:
            return
        phase = str(update.get("phase") or "compile")
        completed = int(update.get("completed_shards") or 0)
        total = int(update.get("total_shards") or 0)
        rows = int(update.get("source_rows") or 0)
        if phase == "compile":
            message = (
                f"Retrosynthesis {mode}: compiled {completed}/{total} "
                f"shard(s), {rows} source reaction(s)."
            )
        elif phase == "merge":
            message = (
                f"Retrosynthesis {mode}: merging {completed}/{total} "
                "operator shard(s)."
            )
        elif phase == "finalize":
            message = f"Retrosynthesis {mode}: finalizing operator library."
        else:
            message = (
                f"Retrosynthesis {mode} operator library ready: "
                f"{update.get('operator_count', 0)} operator(s)."
            )
        progress_callback(
            RecommendationArtifactProgress(
                phase=f"retrosynthesis_{phase}",
                source_file_count=total,
                shard_count=completed,
                row_count=rows,
                message=message,
            )
        )

    _, report = build_full_scale_operator_library(
        source_library_folder,
        output,
        config=FullScaleBuildConfig(),
        workers=workers,
        progress_callback=emit,
    )
    return {
        **report,
        "library_mode": mode,
        "output_dir": str(output.resolve()),
    }


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


class SavedBatchWorker(QtCore.QObject):
    """Save one conversion batch and optionally rebuild the combined index."""

    progress = QtCore.pyqtSignal(object)
    finished = QtCore.pyqtSignal(bool, object, str)

    def __init__(
        self,
        source_inputs: Iterable[str],
        library_folder: str,
        *,
        batch_name: str = "",
        shard_size: int = 1_000,
        workers: int = 1,
        combine_after_save: bool = True,
        use_rxnmapper: bool = True,
        conversion_mode: str = "full",
        build_retrosynthesis: bool = False,
        retrosynthesis_output_folder: str = str(
            DEFAULT_RETROSYNTHESIS_OUTPUT_FOLDER
        ),
    ) -> None:
        super().__init__()
        self.source_inputs = tuple(source_inputs)
        self.library_folder = library_folder
        self.batch_name = batch_name
        self.shard_size = shard_size
        self.workers = workers
        self.combine_after_save = combine_after_save
        self.use_rxnmapper = use_rxnmapper
        self.conversion_mode = conversion_mode
        self.build_retrosynthesis = build_retrosynthesis
        self.retrosynthesis_output_folder = retrosynthesis_output_folder
        self._cancel_requested = False

    def request_cancel(self) -> None:
        """Request a safe stop after the active conversion unit."""
        self._cancel_requested = True

    @QtCore.pyqtSlot()
    def run(self) -> None:
        """Save the selected batch and optionally refresh active artifacts."""
        try:
            batch_report = save_recommendation_batch(
                self.source_inputs,
                self.library_folder,
                batch_name=self.batch_name,
                shard_size=self.shard_size,
                workers=self.workers,
                use_rxnmapper=self.use_rxnmapper,
                conversion_mode=self.conversion_mode,
                progress_callback=self.progress.emit,
                cancel_check=lambda: self._cancel_requested,
            )
            combined_report = None
            if self.combine_after_save or self.build_retrosynthesis:
                combined_report = combine_saved_recommendation_batches(
                    self.library_folder,
                    progress_callback=self.progress.emit,
                    cancel_check=lambda: self._cancel_requested,
                )
            retrosynthesis_report = None
            if self.build_retrosynthesis:
                retrosynthesis_report = build_retrosynthesis_operator_artifacts(
                    self.library_folder,
                    self.retrosynthesis_output_folder,
                    library_mode=self.conversion_mode,
                    workers=self.workers,
                    progress_callback=self.progress.emit,
                    cancel_check=lambda: self._cancel_requested,
                )
        except RecommendationArtifactBuildCancelled as exc:
            self.finished.emit(False, {}, str(exc))
        except Exception as exc:
            self.finished.emit(False, {}, f"{type(exc).__name__}: {exc}")
        else:
            self.finished.emit(
                True,
                {
                    "operation": "save_batch",
                    "output_dir": self.library_folder,
                    "record_count": batch_report["record_count"],
                    "saved_batch": batch_report,
                    "combined": combined_report,
                    "retrosynthesis": retrosynthesis_report,
                },
                "",
            )


class CombineSavedBatchesWorker(QtCore.QObject):
    """Build active recommender artifacts from all saved conversion batches."""

    progress = QtCore.pyqtSignal(object)
    finished = QtCore.pyqtSignal(bool, object, str)

    def __init__(
        self,
        library_folder: str,
        *,
        resume_incomplete: bool = False,
        workers: int = 1,
        library_mode: str = "full",
        build_retrosynthesis: bool = False,
        retrosynthesis_output_folder: str = str(
            DEFAULT_RETROSYNTHESIS_OUTPUT_FOLDER
        ),
    ) -> None:
        super().__init__()
        self.library_folder = library_folder
        self.resume_incomplete = resume_incomplete
        self.workers = workers
        self.library_mode = library_mode
        self.build_retrosynthesis = build_retrosynthesis
        self.retrosynthesis_output_folder = retrosynthesis_output_folder
        self._cancel_requested = False

    def request_cancel(self) -> None:
        """Request cancellation at the next safe checkpoint."""
        self._cancel_requested = True

    @QtCore.pyqtSlot()
    def run(self) -> None:
        """Combine saved batches and emit the resulting active index report."""
        try:
            resumed = []
            if self.resume_incomplete:
                for manifest_path in incomplete_saved_conversion_batches(
                    self.library_folder
                ):
                    self.progress.emit(
                        RecommendationArtifactProgress(
                            phase="resume_started",
                            source_file_count=0,
                            shard_count=0,
                            row_count=0,
                            message=(
                                "Resuming incomplete saved conversion: "
                                f"{manifest_path.parent}"
                            ),
                        )
                    )
                    resumed.append(
                        resume_saved_conversion_batch(
                            manifest_path,
                            workers=self.workers,
                            progress_callback=self.progress.emit,
                            cancel_check=lambda: self._cancel_requested,
                        )
                    )
            report = combine_saved_recommendation_batches(
                self.library_folder,
                progress_callback=self.progress.emit,
                cancel_check=lambda: self._cancel_requested,
            )
            retrosynthesis_report = None
            if self.build_retrosynthesis:
                retrosynthesis_report = build_retrosynthesis_operator_artifacts(
                    self.library_folder,
                    self.retrosynthesis_output_folder,
                    library_mode=self.library_mode,
                    workers=self.workers,
                    progress_callback=self.progress.emit,
                    cancel_check=lambda: self._cancel_requested,
                )
        except RecommendationArtifactBuildCancelled as exc:
            self.finished.emit(False, {}, str(exc))
        except Exception as exc:
            self.finished.emit(False, {}, f"{type(exc).__name__}: {exc}")
        else:
            self.finished.emit(
                True,
                {
                    **report,
                    "operation": "combine_batches",
                    "resumed_batch_count": len(resumed),
                    "retrosynthesis": retrosynthesis_report,
                },
                "",
            )


class RetrosynthesisOnlyWorker(QtCore.QObject):
    """Build operators from an existing combined recommendation corpus."""

    progress = QtCore.pyqtSignal(object)
    finished = QtCore.pyqtSignal(bool, object, str)

    def __init__(
        self,
        library_folder: str,
        *,
        workers: int = 1,
        library_mode: str = "full",
        retrosynthesis_output_folder: str = str(
            DEFAULT_RETROSYNTHESIS_OUTPUT_FOLDER
        ),
    ) -> None:
        super().__init__()
        self.library_folder = library_folder
        self.workers = workers
        self.library_mode = library_mode
        self.retrosynthesis_output_folder = retrosynthesis_output_folder
        self._cancel_requested = False

    def request_cancel(self) -> None:
        """Request cancellation at the next operator-build checkpoint."""
        self._cancel_requested = True

    @QtCore.pyqtSlot()
    def run(self) -> None:
        """Build only retrosynthesis artifacts and emit the result."""
        try:
            report = build_retrosynthesis_operator_artifacts(
                self.library_folder,
                self.retrosynthesis_output_folder,
                library_mode=self.library_mode,
                workers=self.workers,
                progress_callback=self.progress.emit,
                cancel_check=lambda: self._cancel_requested,
            )
        except RecommendationArtifactBuildCancelled as exc:
            self.finished.emit(False, {}, str(exc))
        except Exception as exc:
            self.finished.emit(False, {}, f"{type(exc).__name__}: {exc}")
        else:
            self.finished.emit(
                True,
                {
                    "operation": "retrosynthesis_only",
                    "output_dir": report["output_dir"],
                    "retrosynthesis": report,
                },
                "",
            )


class GenericReactionReviewWindow(QtWidgets.QWidget):
    """Desktop controller for canonical conversion and review export."""

    def __init__(self) -> None:
        super().__init__()
        self.setFont(QtGui.QFont("Segoe UI", 9))
        self.setWindowTitle("Reaction Recommendation Dataset Builder")
        self.resize(880, 680)
        self.thread: Optional[QtCore.QThread] = None
        self.worker: Optional[QtCore.QObject] = None
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
            "Library folder for saved batches and the combined index"
        )
        self.output_edit.setText(self._automatic_output)
        self.output_edit.textChanged.connect(self.refresh_batch_summary)
        self.batch_name_edit = QtWidgets.QLineEdit()
        self.batch_name_edit.setObjectName("batchName")
        self.batch_name_edit.setPlaceholderText(
            "Optional; leave blank for a stable name based on the selected inputs"
        )
        self.source_summary = QtWidgets.QLabel("No source files or folders selected.")
        self.source_summary.setWordWrap(True)
        self.batch_summary = QtWidgets.QLabel()
        self.batch_summary.setObjectName("savedBatchSummary")

        self.conversion_mode_combo = QtWidgets.QComboBox()
        self.conversion_mode_combo.setObjectName("conversionMode")
        self.conversion_mode_combo.addItem("Full", "full")
        self.conversion_mode_combo.addItem("Compact", "compact")
        self.conversion_mode_combo.setCurrentIndex(
            self.conversion_mode_combo.findData("compact")
        )
        self.conversion_mode_combo.setToolTip(
            "Full converts every record. Compact keeps the first 200 records "
            "in each source file plus a deterministic random 15% sample of "
            "the remainder. Each mode has its own saved batches and index."
        )
        self.conversion_mode_combo.currentIndexChanged.connect(
            self.refresh_batch_summary
        )

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
            "Combine all saved batches and rebuild the index after saving"
        )
        self.build_index_check.setObjectName("buildFastIndex")
        self.build_index_check.setChecked(False)
        self.build_index_check.setToolTip(
            "When checked, the active recommender index is refreshed after this "
            "batch is saved. You can also rebuild it later with the Combine / "
            "Build Index button."
        )
        self.build_retrosynthesis_check = QtWidgets.QCheckBox(
            "Build retrosynthesis operators from all saved data"
        )
        self.build_retrosynthesis_check.setObjectName("buildRetrosynthesis")
        self.build_retrosynthesis_check.setChecked(False)
        self.build_retrosynthesis_check.setToolTip(
            "Off by default. When enabled, all saved batches in the selected "
            "Compact or Full library are combined before a resumable "
            "retrosynthesis operator library is built."
        )
        self._combine_before_retrosynthesis = False
        self.build_retrosynthesis_check.toggled.connect(
            self._sync_retrosynthesis_setting
        )
        self.use_rxnmapper_check = QtWidgets.QCheckBox(
            "Use RXNMapper for unresolved, ambiguous, or missing-core reactions"
        )
        self.use_rxnmapper_check.setObjectName("useRxnMapper")
        self.use_rxnmapper_check.setChecked(False)
        self.use_rxnmapper_check.setToolTip(
            "Checked by default. RXNMapper proposals are reconciled with "
            "internal evidence; resolved reactions retain their internal "
            "signature while consensus mapping supplies only a missing core. "
            "One worker is used to avoid loading multiple model copies."
        )
        self.use_rxnmapper_check.toggled.connect(self._sync_rxnmapper_worker_setting)
        self._sync_rxnmapper_worker_setting(False)

        self.options_widget = QtWidgets.QWidget()
        self.options_widget.setObjectName("conversionOptions")
        options_layout = QtWidgets.QHBoxLayout(self.options_widget)
        options_layout.setContentsMargins(0, 0, 0, 0)
        options_layout.addWidget(self.use_rxnmapper_check)
        options_layout.addSpacing(16)
        options_layout.addWidget(self.build_index_check)
        options_layout.addSpacing(16)
        options_layout.addWidget(self.build_retrosynthesis_check)
        options_layout.addStretch()

        self.start_button = QtWidgets.QPushButton("Convert && Save Batch")
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
        self.combine_button = QtWidgets.QPushButton(
            "Combine Saved Batches / Build Index"
        )
        self.combine_button.setObjectName("combineIndexButton")
        self.retrosynthesis_button = QtWidgets.QPushButton(
            "Build Retrosynthesis Only"
        )
        self.retrosynthesis_button.setObjectName("retrosynthesisOnlyButton")
        self.retrosynthesis_button.setToolTip(
            "Build or resume the selected Compact or Full retrosynthesis "
            "operator library from the existing combined corpus. This does "
            "not reconvert data, combine batches, or rebuild the index."
        )
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
        self.combine_button.clicked.connect(self.start_combining)
        self.retrosynthesis_button.clicked.connect(
            self.start_retrosynthesis_only
        )
        self.cancel_button.clicked.connect(self.cancel_conversion)
        self.open_button.clicked.connect(self.open_output_folder)
        self.refresh_batch_summary()

    def _build_layout(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(12)

        title = QtWidgets.QLabel("Reaction Recommendation Dataset Builder")
        title.setStyleSheet("font-size: 20px; font-weight: 600;")
        layout.addWidget(title)
        description = QtWidgets.QLabel(
            "Convert selected raw CSV or preprocessed observation files, or "
            "every supported file in selected folder trees, into a reusable "
            "saved batch. Repeat this for later files, then combine every saved "
            "batch into the active recommender index. Interrupted conversions "
            "can reuse completed shards."
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
        form.addRow("Library mode:", self.conversion_mode_combo)
        form.addRow("Batch name:", self.batch_name_edit)
        form.addRow("Saved data:", self.batch_summary)

        settings_row = QtWidgets.QHBoxLayout()
        settings_row.addWidget(QtWidgets.QLabel("Rows per shard"))
        settings_row.addWidget(self.shard_size_spin)
        settings_row.addSpacing(16)
        settings_row.addWidget(QtWidgets.QLabel("Parallel workers"))
        settings_row.addWidget(self.worker_count_spin)
        settings_row.addStretch()
        form.addRow("Performance:", settings_row)
        form.addRow("Options:", self.options_widget)
        layout.addLayout(form)

        button_row = QtWidgets.QHBoxLayout()
        button_row.addWidget(self.start_button)
        button_row.addWidget(self.combine_button)
        button_row.addWidget(self.retrosynthesis_button)
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

    @QtCore.pyqtSlot()
    @QtCore.pyqtSlot(str)
    @QtCore.pyqtSlot(int)
    def refresh_batch_summary(self, _value: str | int = "") -> None:
        """Show how many independently saved batches can be combined."""
        output_text = self.output_edit.text().strip()
        mode = str(self.conversion_mode_combo.currentData() or "compact")
        library_dir = (
            recommendation_library_mode_dir(output_text, mode)
            if output_text
            else None
        )
        count = (
            len(discover_saved_conversion_batches(library_dir))
            if library_dir is not None
            else 0
        )
        noun = "batch" if count == 1 else "batches"
        self.batch_summary.setText(
            f"{count} saved {noun} available to combine in the "
            f"{mode.title()} library ({library_dir})."
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

    @QtCore.pyqtSlot(bool)
    def _sync_retrosynthesis_setting(self, enabled: bool) -> None:
        """Require a current combined corpus before building operators."""
        if enabled:
            self._combine_before_retrosynthesis = (
                self.build_index_check.isChecked()
            )
            self.build_index_check.setChecked(True)
            self.build_index_check.setEnabled(False)
        else:
            self.build_index_check.setEnabled(True)
            self.build_index_check.setChecked(
                self._combine_before_retrosynthesis
            )

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
        mode = str(self.conversion_mode_combo.currentData() or "compact")
        output_root = Path(output_text)
        output = recommendation_library_mode_dir(output_root, mode)
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
        self._append_status(f"Recommendation library root: {output_root}")
        self._append_status(
            f"Conversion mode: {mode.title()}; artifacts: {output}"
        )
        if self.build_retrosynthesis_check.isChecked():
            self._append_status(
                "Retrosynthesis conversion: on; all available "
                f"{mode} data -> "
                f"{DEFAULT_RETROSYNTHESIS_OUTPUT_FOLDER / mode}"
            )
        if self.batch_name_edit.text().strip():
            self._append_status(
                f"Saved batch name: {self.batch_name_edit.text().strip()}"
            )
        self._append_status(
            f"Settings: {self.shard_size_spin.value()} rows/shard, "
            f"{self.worker_count_spin.value()} worker(s), "
            f"RXNMapper "
            f"{'on' if self.use_rxnmapper_check.isChecked() else 'off'}, "
            "combine after save "
            f"{'on' if self.build_index_check.isChecked() else 'off'}"
        )
        self.start_button.setEnabled(False)
        self.combine_button.setEnabled(False)
        self.retrosynthesis_button.setEnabled(False)
        self.cancel_button.setEnabled(True)
        self.open_button.setEnabled(False)
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setFormat("Discovering files…")

        thread = QtCore.QThread(self)
        worker = SavedBatchWorker(
            source_inputs,
            str(output),
            batch_name=self.batch_name_edit.text().strip(),
            shard_size=self.shard_size_spin.value(),
            workers=self.worker_count_spin.value(),
            combine_after_save=self.build_index_check.isChecked(),
            use_rxnmapper=self.use_rxnmapper_check.isChecked(),
            conversion_mode=mode,
            build_retrosynthesis=(
                self.build_retrosynthesis_check.isChecked()
            ),
            retrosynthesis_output_folder=str(
                DEFAULT_RETROSYNTHESIS_OUTPUT_FOLDER
            ),
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
    def start_retrosynthesis_only(self) -> None:
        """Build operators without rebuilding recommendation artifacts."""
        if self.thread is not None:
            return
        output_text = self.output_edit.text().strip()
        if not output_text:
            QtWidgets.QMessageBox.warning(
                self,
                "Output required",
                "Choose the recommendation library folder.",
            )
            return
        mode = str(self.conversion_mode_combo.currentData() or "compact")
        library_dir = recommendation_library_mode_dir(output_text, mode)
        combined_manifest = library_dir / COMBINED_BATCH_MANIFEST_FILENAME
        if not combined_manifest.is_file():
            QtWidgets.QMessageBox.warning(
                self,
                "Combined library required",
                "No combined recommendation corpus exists for the selected "
                f"{mode.title()} mode. Combine saved batches once before "
                "building retrosynthesis operators.",
            )
            return

        retrosynthesis_output = DEFAULT_RETROSYNTHESIS_OUTPUT_FOLDER / mode
        self.status_box.clear()
        self._append_status(
            f"Building retrosynthesis operators only from {library_dir}"
        )
        self._append_status(
            "Recommendation conversion, batch combination, and index rebuild "
            "will be skipped."
        )
        self._append_status(
            f"Retrosynthesis output: {retrosynthesis_output}; "
            f"workers: {self.worker_count_spin.value()}"
        )
        self.start_button.setEnabled(False)
        self.combine_button.setEnabled(False)
        self.retrosynthesis_button.setEnabled(False)
        self.cancel_button.setEnabled(True)
        self.open_button.setEnabled(False)
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setFormat("Discovering operator source shards…")

        thread = QtCore.QThread(self)
        worker = RetrosynthesisOnlyWorker(
            str(library_dir),
            workers=self.worker_count_spin.value(),
            library_mode=mode,
            retrosynthesis_output_folder=str(
                DEFAULT_RETROSYNTHESIS_OUTPUT_FOLDER
            ),
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
    def start_combining(self) -> None:
        """Combine all saved batches into the active recommender index."""
        if self.thread is not None:
            return
        output_text = self.output_edit.text().strip()
        if not output_text:
            QtWidgets.QMessageBox.warning(
                self,
                "Output required",
                "Choose the recommendation library folder.",
            )
            return
        mode = str(self.conversion_mode_combo.currentData() or "compact")
        library_dir = recommendation_library_mode_dir(output_text, mode)
        manifests = discover_saved_conversion_batches(library_dir)
        if not manifests:
            QtWidgets.QMessageBox.warning(
                self,
                "No saved batches",
                "Convert and save at least one batch before building the index.",
            )
            return
        incomplete = incomplete_saved_conversion_batches(library_dir)
        resume_incomplete = False
        if incomplete:
            listed = "\n".join(f"• {path.parent}" for path in incomplete[:5])
            if len(incomplete) > 5:
                listed += f"\n• …and {len(incomplete) - 5} more"
            choice = QtWidgets.QMessageBox.question(
                self,
                "Resume incomplete conversion?",
                (
                    f"{len(incomplete)} saved conversion(s) are incomplete:\n\n"
                    f"{listed}\n\n"
                    "Resume their checkpointed source files now and combine "
                    "only after coverage is complete? Valid completed shards "
                    "will be reused. This may take some time."
                ),
                (
                    QtWidgets.QMessageBox.StandardButton.Yes
                    | QtWidgets.QMessageBox.StandardButton.No
                ),
                QtWidgets.QMessageBox.StandardButton.Yes,
            )
            if choice != QtWidgets.QMessageBox.StandardButton.Yes:
                return
            resume_incomplete = True
        self.status_box.clear()
        if resume_incomplete:
            self._append_status(
                f"Resuming {len(incomplete)} incomplete conversion(s), then "
                f"combining {len(manifests)} saved {mode} batch(es) in "
                f"{library_dir}"
            )
        else:
            self._append_status(
                f"Combining {len(manifests)} saved {mode} batch(es) in "
                f"{library_dir}"
            )
        self.start_button.setEnabled(False)
        self.combine_button.setEnabled(False)
        self.retrosynthesis_button.setEnabled(False)
        self.cancel_button.setEnabled(True)
        self.open_button.setEnabled(False)
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setFormat("Combining saved batches…")

        thread = QtCore.QThread(self)
        worker = CombineSavedBatchesWorker(
            str(library_dir),
            resume_incomplete=resume_incomplete,
            workers=self.worker_count_spin.value(),
            library_mode=mode,
            build_retrosynthesis=(
                self.build_retrosynthesis_check.isChecked()
            ),
            retrosynthesis_output_folder=str(
                DEFAULT_RETROSYNTHESIS_OUTPUT_FOLDER
            ),
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
        elif progress.phase.startswith("index"):
            self.progress_bar.setFormat(
                f"Building fast index • {progress.row_count} reactions"
            )
        elif progress.phase.startswith("combine"):
            self.progress_bar.setFormat(
                f"Combining batches • {progress.row_count} reactions"
            )
        elif progress.phase.startswith("retrosynthesis"):
            self.progress_bar.setFormat(progress.message)
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
        self.combine_button.setEnabled(True)
        self.retrosynthesis_button.setEnabled(True)
        self.cancel_button.setEnabled(False)
        if success:
            self.progress_bar.setRange(0, 1)
            self.progress_bar.setValue(1)
            operation = str(report.get("operation") or "")
            combined = report.get("combined")
            active_report = (
                combined
                if isinstance(combined, dict)
                else report if operation == "combine_batches" else None
            )
            if operation == "retrosynthesis_only":
                retrosynthesis = report["retrosynthesis"]
                self.progress_bar.setFormat(
                    "Retrosynthesis ready • "
                    f"{retrosynthesis['operator_count']} operators"
                )
            elif operation == "save_batch":
                batch = report["saved_batch"]
                self.progress_bar.setFormat(
                    f"Batch saved • {batch['record_count']} reactions"
                )
                self._append_status(
                    f"Saved batch '{batch['batch_name']}': "
                    f"{batch['record_count']} reaction(s) in {batch['batch_dir']}"
                )
                self._append_status(
                    f"  Canonical manifest: "
                    f"{batch['artifacts']['canonical_manifest']['path']}"
                )
                if active_report is None:
                    self._append_status(
                        "The batch is saved. Use Combine Saved Batches / Build "
                        "Index when you are ready to refresh the recommender."
                    )
            else:
                self.progress_bar.setFormat(
                    f"Combined index ready • {report['record_count']} reactions"
                )
            self._completed_output = Path(report["output_dir"])
            self.open_button.setEnabled(True)
            if active_report is not None:
                artifacts = active_report["artifacts"]
                self._append_status("Active recommender data is ready:")
                self._append_status(
                    "  Combined canonical records: "
                    f"{artifacts['canonical_records']['path']} "
                    f"({_human_size(artifacts['canonical_records']['size_bytes'])})"
                )
                self._append_status(
                    "  Active fast recommendation index: "
                    f"{artifacts['fast_index']['path']} "
                    f"({_human_size(artifacts['fast_index']['size_bytes'])})"
                )
                self._append_status(
                    f"  Saved batches: {active_report['batch_count']}; "
                    f"unique reactions: {active_report['record_count']}; "
                    "duplicates skipped: "
                    f"{active_report['duplicate_record_count']}"
                )
                if int(active_report.get("resumed_batch_count") or 0):
                    self._append_status(
                        "  Resumed incomplete conversions: "
                        f"{active_report['resumed_batch_count']}"
                    )
                self._append_status(
                    "  Trusted recommendation precedents: "
                    f"{active_report['eligible_index_record_count']}"
                )
                for warning in active_report.get("warnings") or ():
                    self._append_status(f"Warning: {warning}")
            elif not operation:
                # Preserve display support for direct, non-batch worker callers.
                artifacts = report["artifacts"]
                self._append_status(
                    f"Ready: {artifacts['canonical_manifest']['path']}"
                )
            retrosynthesis = report.get("retrosynthesis")
            if isinstance(retrosynthesis, dict):
                self._append_status(
                    "Retrosynthesis operator library is ready: "
                    f"{retrosynthesis['library_path']}"
                )
                self._append_status(
                    f"  Mode: {retrosynthesis['library_mode'].title()}; "
                    f"source reactions: {retrosynthesis['source_rows']}; "
                    f"operators: {retrosynthesis['operator_count']}; "
                    f"templates: {retrosynthesis['template_count']}"
                )
            self.refresh_batch_summary()
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
    "CombineSavedBatchesWorker",
    "DEFAULT_OUTPUT_FOLDER",
    "DEFAULT_RETROSYNTHESIS_OUTPUT_FOLDER",
    "GenericReactionReviewWindow",
    "RetrosynthesisOnlyWorker",
    "ReviewConversionWorker",
    "SavedBatchWorker",
    "build_retrosynthesis_operator_artifacts",
    "main",
]

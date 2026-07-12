import sys
import os
import json
import pandas as pd
from pathlib import Path
from typing import Optional, List, Dict, Any

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from PyQt6 import QtWidgets, QtCore
from app.A_convert_to_hte_format import process_reaction_dataset, enrich_reaction_dataset_cas
from chemtools.recommend.reaction_key_utils import (
    build_reaction_events_payload,
    serialize_reaction_events_payload,
)

DEFAULT_LLM_PROVIDER = os.getenv("CHEMTOOLS_LLM_PROVIDER", os.getenv("LLM_PROVIDER", "openai")).strip() or "openai"
DEFAULT_LLM_MODEL = os.getenv("CHEMTOOLS_LLM_MODEL", os.getenv("LLM_MODEL", "")).strip()

try:
    from llmtools.clients import AVAILABLE_MODELS as LLM_AVAILABLE_MODELS
    from llmtools.clients import RECOMMENDED_MODELS as LLM_RECOMMENDED_MODELS
except Exception:
    LLM_AVAILABLE_MODELS = {
        "openai": ["gpt-4o", "gpt-4o-mini", "gpt-5-mini"],
        "aliyun": ["deepseek-v3.2", "deepseek-v3", "deepseek-r1"],
    }
    LLM_RECOMMENDED_MODELS = {
        "openai": {"balanced": "gpt-4o"},
        "aliyun": {"balanced": "deepseek-v3.2"},
    }


class ConversionWorker(QtCore.QObject):
    finished = QtCore.pyqtSignal(bool, str)
    progress = QtCore.pyqtSignal(str)

    def __init__(
        self,
        jobs: List[tuple[str, str]],
        drop_no_catalyst: bool,
        reagent_csv_path: str,
        new_reagents_path: str,
        skip_cas_enrichment: bool,
        num_workers: int,
        llm_assist_options: Optional[Dict[str, Any]] = None,
    ):
        super().__init__()
        self.jobs = jobs
        self.drop_no_catalyst = drop_no_catalyst
        self.reagent_csv_path = reagent_csv_path
        self.new_reagents_path = new_reagents_path
        self.skip_cas_enrichment = skip_cas_enrichment
        self.num_workers = num_workers
        self.llm_assist_options = llm_assist_options or None

    def run(self):
        try:
            # Custom stdout capturing to redirect prints to the GUI log
            class StdoutRedirector:
                def __init__(self, signal):
                    self.signal = signal
                def write(self, text):
                    if text.strip():
                        self.signal.emit(text.strip())
                def flush(self):
                    pass

            old_stdout = sys.stdout
            sys.stdout = StdoutRedirector(self.progress)
            
            try:
                total = len(self.jobs)
                for idx, (input_path, output_path) in enumerate(self.jobs, start=1):
                    print(
                        f"Job {idx}/{total}: {os.path.basename(input_path)} -> "
                        f"{os.path.basename(output_path)}"
                    )
                    # 1. Enrich CAS numbers first
                    if self.skip_cas_enrichment:
                        print("Step 1: Skipping CAS enrichment (fast mode).")
                    else:
                        print(f"Step 1: Enriching CAS numbers in {os.path.basename(input_path)}...")
                        enrich_reaction_dataset_cas(input_path)

                    # 2. Process to HTE format
                    print("Step 2: Converting to HTE format...")
                    process_reaction_dataset(
                        input_path,
                        output_path,
                        drop_no_catalyst=self.drop_no_catalyst,
                        reagent_csv_path=self.reagent_csv_path,
                        new_reagents_path=self.new_reagents_path,
                        num_workers=self.num_workers,
                        llm_assist_options=self.llm_assist_options,
                    )
                    print("")
                self.finished.emit(True, f"Successfully processed {total} file(s).")
            finally:
                sys.stdout = old_stdout
                
        except Exception as e:
            self.finished.emit(False, f"Error: {str(e)}")


def _ensure_reaction_events_column(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure output CSV contains standardized Reaction_Events payloads."""
    if "Reaction_Key" not in df.columns:
        return df
    out = df.copy()
    if "Reaction_Events" not in out.columns:
        out["Reaction_Events"] = ""
    series = out["Reaction_Events"].fillna("").astype(str).str.strip()
    missing_mask = series.eq("") | series.str.lower().eq("nan")
    if not missing_mask.any():
        return out

    key_series = out["Reaction_Key"].fillna("").astype(str)
    out.loc[missing_mask, "Reaction_Events"] = key_series[missing_mask].apply(
        lambda key: serialize_reaction_events_payload(build_reaction_events_payload(key))
    )
    return out

class HTEConverterWindow(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("HTE Dataset Converter (Qt6)")
        self.resize(700, 500)
        
        # UI Elements
        self.dataset_dropdown = QtWidgets.QComboBox()
        self.source_dir_edit = QtWidgets.QLineEdit()
        self.source_dir_edit.setPlaceholderText("Source folder containing *.csv files")
        self.input_edit = QtWidgets.QLineEdit()
        self.input_edit.setPlaceholderText("Select one or more CSV files")
        self.output_edit = QtWidgets.QLineEdit()
        self.drop_catalyst_check = QtWidgets.QCheckBox("Drop reactions without a catalyst")
        self.drop_catalyst_check.setChecked(True)
        self.skip_cas_check = QtWidgets.QCheckBox("Skip CAS enrichment (faster)")
        self.protocol_mode_check = QtWidgets.QCheckBox(
            "Literature+setup CSV mode (skip enrichment, save to HTE_db/literature)"
        )
        self.worker_count_spin = QtWidgets.QSpinBox()
        self.worker_count_spin.setRange(0, 64)
        self.worker_count_spin.setValue(0)
        self.worker_count_spin.setSpecialValueText("Auto")
        self.worker_count_spin.setToolTip("0 = auto (uses a moderate fraction of CPU cores), otherwise use a fixed number of worker processes.")
        self.llm_assist_checkbox = QtWidgets.QCheckBox("LLM assist for reaction typing/key")
        self.llm_provider_combo = QtWidgets.QComboBox()
        self.llm_model_combo = QtWidgets.QComboBox()
        self.llm_model_combo.setEditable(True)
        self.llm_temperature_spin = QtWidgets.QDoubleSpinBox()
        self.llm_temperature_spin.setRange(0.0, 2.0)
        self.llm_temperature_spin.setSingleStep(0.1)
        self.llm_temperature_spin.setDecimals(2)
        self.llm_temperature_spin.setValue(0.0)
        self.llm_max_tokens_spin = QtWidgets.QSpinBox()
        self.llm_max_tokens_spin.setRange(16, 32768)
        self.llm_max_tokens_spin.setValue(700)
        self.llm_timeout_spin = QtWidgets.QSpinBox()
        self.llm_timeout_spin.setRange(5, 600)
        self.llm_timeout_spin.setValue(60)
        self.llm_threshold_spin = QtWidgets.QDoubleSpinBox()
        self.llm_threshold_spin.setRange(0.0, 1.0)
        self.llm_threshold_spin.setSingleStep(0.05)
        self.llm_threshold_spin.setDecimals(2)
        self.llm_threshold_spin.setValue(0.60)
        self.llm_always_checkbox = QtWidgets.QCheckBox("Always run LLM (not only uncertain)")
        self.llm_no_crk_validation_checkbox = QtWidgets.QCheckBox("Disable CRK validation gate")
        
        self.btn_run = QtWidgets.QPushButton("Start Conversion")
        self.btn_quit = QtWidgets.QPushButton("Quit")
        self.log = QtWidgets.QPlainTextEdit()
        self.log.setReadOnly(True)
        
        self.progress_bar = QtWidgets.QProgressBar()
        self.progress_bar.setVisible(False)
        
        self._setup_layout()
        self._populate_llm_providers()
        self._on_llm_assist_toggled(self.llm_assist_checkbox.isChecked())
        self._load_datasets()
        
        # Connections
        self.dataset_dropdown.currentIndexChanged.connect(self._on_dataset_changed)
        self.source_dir_edit.editingFinished.connect(self._load_datasets)
        self.input_edit.textEdited.connect(self._on_input_edited)
        self.protocol_mode_check.stateChanged.connect(self._on_protocol_mode_changed)
        self.llm_assist_checkbox.toggled.connect(self._on_llm_assist_toggled)
        self.llm_provider_combo.currentTextChanged.connect(self._on_llm_provider_changed)
        self.btn_run.clicked.connect(self.run_processing)
        self.btn_quit.clicked.connect(self.close)
        
        # Threading state
        self.thread = None
        self.worker = None
        self.selected_input_paths: List[str] = []

    def _setup_layout(self):
        layout = QtWidgets.QVBoxLayout(self)
        
        title = QtWidgets.QLabel("HTE Dataset Converter")
        title.setStyleSheet("font-size: 18px; font-weight: bold; margin: 10px;")
        title.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title)
        
        form = QtWidgets.QFormLayout()
        
        source_layout = QtWidgets.QHBoxLayout()
        source_layout.addWidget(self.source_dir_edit)
        source_browse_btn = QtWidgets.QPushButton("Browse")
        source_browse_btn.setFixedWidth(80)
        source_browse_btn.clicked.connect(self._choose_source_dir)
        source_layout.addWidget(source_browse_btn)
        form.addRow("Source Folder:", source_layout)

        # Dataset selector
        ds_layout = QtWidgets.QHBoxLayout()
        ds_layout.addWidget(self.dataset_dropdown)
        refresh_btn = QtWidgets.QPushButton("Refresh")
        refresh_btn.setFixedWidth(80)
        refresh_btn.clicked.connect(self._load_datasets)
        ds_layout.addWidget(refresh_btn)
        form.addRow("Available Datasets:", ds_layout)

        input_layout = QtWidgets.QHBoxLayout()
        input_layout.addWidget(self.input_edit)
        input_browse_btn = QtWidgets.QPushButton("Select CSVs")
        input_browse_btn.setFixedWidth(100)
        input_browse_btn.clicked.connect(self._choose_input_files)
        input_layout.addWidget(input_browse_btn)

        form.addRow("Input CSV(s):", input_layout)
        form.addRow("Output CSV/Folder:", self.output_edit)
        form.addRow("", QtWidgets.QLabel("Output includes formed_motifs, spectator_groups, reactant_3, and Reaction_Events columns."))
        form.addRow("", self.drop_catalyst_check)
        form.addRow("", self.skip_cas_check)
        form.addRow("", self.protocol_mode_check)
        form.addRow("Workers:", self.worker_count_spin)

        llm_toggle_row = QtWidgets.QHBoxLayout()
        llm_toggle_row.addWidget(self.llm_assist_checkbox)
        llm_toggle_row.addWidget(self.llm_always_checkbox)
        llm_toggle_row.addWidget(self.llm_no_crk_validation_checkbox)
        llm_toggle_row.addStretch()
        form.addRow("LLM:", llm_toggle_row)

        llm_model_row = QtWidgets.QHBoxLayout()
        llm_model_row.addWidget(QtWidgets.QLabel("Provider"))
        llm_model_row.addWidget(self.llm_provider_combo)
        llm_model_row.addSpacing(12)
        llm_model_row.addWidget(QtWidgets.QLabel("Model"))
        llm_model_row.addWidget(self.llm_model_combo)
        llm_model_row.addStretch()
        form.addRow("LLM Model:", llm_model_row)

        llm_params_row = QtWidgets.QHBoxLayout()
        llm_params_row.addWidget(QtWidgets.QLabel("Temperature"))
        llm_params_row.addWidget(self.llm_temperature_spin)
        llm_params_row.addSpacing(12)
        llm_params_row.addWidget(QtWidgets.QLabel("Max tokens"))
        llm_params_row.addWidget(self.llm_max_tokens_spin)
        llm_params_row.addSpacing(12)
        llm_params_row.addWidget(QtWidgets.QLabel("Timeout (s)"))
        llm_params_row.addWidget(self.llm_timeout_spin)
        llm_params_row.addSpacing(12)
        llm_params_row.addWidget(QtWidgets.QLabel("Confidence threshold"))
        llm_params_row.addWidget(self.llm_threshold_spin)
        llm_params_row.addStretch()
        form.addRow("LLM Params:", llm_params_row)
        
        layout.addLayout(form)
        
        layout.addWidget(QtWidgets.QLabel("Processing Log:"))
        layout.addWidget(self.log)
        layout.addWidget(self.progress_bar)
        
        btns = QtWidgets.QHBoxLayout()
        btns.addStretch()
        btns.addWidget(self.btn_run)
        btns.addWidget(self.btn_quit)
        layout.addLayout(btns)

    def _load_datasets(self):
        self.dataset_dropdown.clear()
        self.dataset_dropdown.addItem("-- Select or enter custom path --", None)

        source_text = self.source_dir_edit.text().strip()
        dataset_dir = Path(source_text) if source_text else PROJECT_ROOT / "data" / "reaction_dataset"
        self.source_dir_edit.setText(str(dataset_dir))

        if dataset_dir.exists():
            files = sorted(p for p in dataset_dir.rglob("*.csv") if p.is_file())
            for f in files:
                self.dataset_dropdown.addItem(f.name, str(f))

    def _on_dataset_changed(self, index):
        if index <= 0:
            return
            
        data_path = self.dataset_dropdown.currentData()
        if not data_path:
            return

        input_path = Path(data_path)
        output_path = PROJECT_ROOT / "data" / "HTE_db" / "literature" / f"{input_path.stem}_canonical.csv"

        self.input_edit.setText(str(input_path))
        self.output_edit.setText(str(output_path))
        self.selected_input_paths = [str(input_path)]

    def _choose_source_dir(self):
        path = QtWidgets.QFileDialog.getExistingDirectory(
            self,
            "Select Source Folder",
            str(PROJECT_ROOT),
        )
        if path:
            self.source_dir_edit.setText(path)
            self._load_datasets()

    def _choose_input_files(self):
        paths, _ = QtWidgets.QFileDialog.getOpenFileNames(
            self,
            "Select CSV Files",
            str(PROJECT_ROOT),
            "CSV Files (*.csv)",
        )
        if not paths:
            return

        self.dataset_dropdown.setCurrentIndex(0)
        self.selected_input_paths = paths
        self.input_edit.setText("; ".join(paths))

        if len(paths) > 1:
            default_out_dir = PROJECT_ROOT / "data" / "HTE_db" / "literature"
            self.output_edit.setText(str(default_out_dir))
        else:
            input_path = Path(paths[0])
            output_path = PROJECT_ROOT / "data" / "HTE_db" / "literature" / f"{input_path.stem}_canonical.csv"
            self.output_edit.setText(str(output_path))

    def _on_input_edited(self):
        if self.dataset_dropdown.currentIndex() > 0:
            self.dataset_dropdown.setCurrentIndex(0)
        self.selected_input_paths = []

    def _on_protocol_mode_changed(self):
        """Update output path when literature+setup mode is toggled."""
        if not self.protocol_mode_check.isChecked():
            return
        
        # Auto-set output to literature directory (protocol-like datasets are folded into literature)
        input_paths = self._get_input_paths()
        if input_paths:
            if len(input_paths) > 1:
                output_path = PROJECT_ROOT / "data" / "HTE_db" / "literature"
                self.output_edit.setText(str(output_path))
            else:
                input_path = Path(input_paths[0])
                output_path = PROJECT_ROOT / "data" / "HTE_db" / "literature" / f"{input_path.stem}_hte.csv"
                self.output_edit.setText(str(output_path))
        
        # Enable skip CAS enrichment automatically
        self.skip_cas_check.setChecked(True)

    def _populate_llm_providers(self) -> None:
        providers = sorted(LLM_AVAILABLE_MODELS.keys()) or ["openai", "aliyun"]
        self.llm_provider_combo.blockSignals(True)
        self.llm_provider_combo.clear()
        self.llm_provider_combo.addItems(providers)
        fallback_provider = (
            DEFAULT_LLM_PROVIDER if DEFAULT_LLM_PROVIDER in providers else providers[0]
        )
        self.llm_provider_combo.setCurrentText(fallback_provider)
        self.llm_provider_combo.blockSignals(False)
        self._populate_llm_models(fallback_provider, prefer=DEFAULT_LLM_MODEL)

    def _populate_llm_models(self, provider: str, *, prefer: Optional[str] = None) -> None:
        provider_key = (provider or "openai").strip().lower()
        presets = list(LLM_AVAILABLE_MODELS.get(provider_key) or [])
        current_text = (prefer or self.llm_model_combo.currentText() or "").strip()
        self.llm_model_combo.blockSignals(True)
        try:
            self.llm_model_combo.clear()
            self.llm_model_combo.addItems(presets)
            if current_text and current_text in presets:
                self.llm_model_combo.setCurrentText(current_text)
            else:
                recommended = (LLM_RECOMMENDED_MODELS.get(provider_key) or {}).get("balanced")
                fallback = (recommended or (presets[0] if presets else "")).strip()
                if provider_key == DEFAULT_LLM_PROVIDER.lower() and DEFAULT_LLM_MODEL:
                    fallback = DEFAULT_LLM_MODEL
                self.llm_model_combo.setCurrentText(fallback)
        finally:
            self.llm_model_combo.blockSignals(False)

    def _on_llm_provider_changed(self, provider: str) -> None:
        self._populate_llm_models(provider)

    def _on_llm_assist_toggled(self, enabled: bool) -> None:
        self.llm_provider_combo.setEnabled(enabled)
        self.llm_model_combo.setEnabled(enabled)
        self.llm_temperature_spin.setEnabled(enabled)
        self.llm_max_tokens_spin.setEnabled(enabled)
        self.llm_timeout_spin.setEnabled(enabled)
        self.llm_threshold_spin.setEnabled(enabled)
        self.llm_always_checkbox.setEnabled(enabled)
        self.llm_no_crk_validation_checkbox.setEnabled(enabled)

    def _gather_llm_assist_options(self) -> Optional[Dict[str, Any]]:
        if not self.llm_assist_checkbox.isChecked():
            return None
        llm_model = self.llm_model_combo.currentText().strip()
        if not llm_model:
            raise ValueError("LLM model is required when LLM assist is enabled.")
        return {
            "enabled": True,
            "provider": self.llm_provider_combo.currentText().strip().lower() or "openai",
            "model": llm_model,
            "temperature": float(self.llm_temperature_spin.value()),
            "max_tokens": int(self.llm_max_tokens_spin.value()),
            "timeout": int(self.llm_timeout_spin.value()),
            "only_on_uncertain": not self.llm_always_checkbox.isChecked(),
            "confidence_threshold": float(self.llm_threshold_spin.value()),
            "require_crk_validation": not self.llm_no_crk_validation_checkbox.isChecked(),
        }

    def _get_input_paths(self) -> List[str]:
        if self.selected_input_paths:
            return list(self.selected_input_paths)
        text = self.input_edit.text().strip()
        if not text:
            return []
        return [p.strip() for p in text.split(";") if p.strip()]

    def _build_jobs(self, input_paths: List[str], output_text: str) -> List[tuple[str, str]]:
        output_path = Path(output_text)
        if len(input_paths) > 1:
            if output_path.suffix.lower() == ".csv":
                QtWidgets.QMessageBox.critical(
                    self,
                    "Error",
                    "Output must be a folder when converting multiple CSV files.",
                )
                return []
            output_path.mkdir(parents=True, exist_ok=True)
            return [
                (path, str(output_path / f"{Path(path).stem}_canonical.csv"))
                for path in input_paths
            ]

        if output_path.suffix.lower() == ".csv" or (output_path.exists() and output_path.is_file()):
            output_path.parent.mkdir(parents=True, exist_ok=True)
            return [(input_paths[0], str(output_path))]

        output_path.mkdir(parents=True, exist_ok=True)
        return [
            (input_paths[0], str(output_path / f"{Path(input_paths[0]).stem}_canonical.csv"))
        ]

    def log_msg(self, text: str):
        self.log.appendPlainText(text)

    def run_processing(self):
        input_paths = self._get_input_paths()
        output_text = self.output_edit.text().strip()

        if not input_paths:
            QtWidgets.QMessageBox.critical(self, "Error", "Input path is required.")
            return

        missing_paths = [path for path in input_paths if not os.path.exists(path)]
        if missing_paths:
            QtWidgets.QMessageBox.critical(
                self,
                "Error",
                f"Invalid input path: {missing_paths[0]}",
            )
            return

        if not output_text:
            QtWidgets.QMessageBox.critical(self, "Error", "Output path is required.")
            return

        jobs = self._build_jobs(input_paths, output_text)
        if not jobs:
            return
        try:
            llm_assist_options = self._gather_llm_assist_options()
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "Error", str(exc))
            return

        # Check if literature+setup mode
        if self.protocol_mode_check.isChecked():
            self._run_protocol_conversion(jobs, llm_assist_options=llm_assist_options)
            return

        self.setEnabled(False)
        self.log.clear()
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0) # Indeterminate
        
        # Start worker thread
        self.thread = QtCore.QThread()
        self.worker = ConversionWorker(
            jobs,
            self.drop_catalyst_check.isChecked(),
            str(PROJECT_ROOT / "condition_registry" / "definitions" / "substances.v1.csv"),
            str(PROJECT_ROOT / "condition_registry" / "definitions" / "pending_substances.csv"),
            self.skip_cas_check.isChecked(),
            self.worker_count_spin.value(),
            llm_assist_options=llm_assist_options,
        )
        self.worker.moveToThread(self.thread)
        
        self.thread.started.connect(self.worker.run)
        self.worker.progress.connect(self.log_msg)
        self.worker.finished.connect(self.on_finished)
        
        self.thread.start()

    def _run_protocol_conversion(
        self,
        jobs: List[tuple[str, str]],
        *,
        llm_assist_options: Optional[Dict[str, Any]],
    ):
        """Convert literature+setup CSV to HTE format with full processing."""
        try:
            self.log.clear()
            self.log_msg("Literature+Setup CSV Conversion Mode")
            self.log_msg("="*50)
            self.log_msg("Running full HTE processing pipeline:")
            self.log_msg("  - Reaction type detection")
            self.log_msg("  - Motif extraction")
            self.log_msg("  - Spectator group ranking")
            worker_setting = self.worker_count_spin.value()
            self.log_msg(
                f"  - Worker processes: {'auto' if worker_setting == 0 else worker_setting}"
            )
            if llm_assist_options:
                self.log_msg(
                    "  - LLM assist: "
                    f"{llm_assist_options.get('provider')}/{llm_assist_options.get('model')}"
                )
            self.log_msg("")
            
            # Redirect stdout to capture processing logs
            class StdoutRedirector:
                def __init__(self, log_func):
                    self.log_func = log_func
                def write(self, text):
                    if text.strip():
                        self.log_func(text.strip())
                def flush(self):
                    pass
            
            old_stdout = sys.stdout
            sys.stdout = StdoutRedirector(self.log_msg)
            
            try:
                total = len(jobs)
                for idx, (input_path, output_path) in enumerate(jobs, start=1):
                    print(f"Job {idx}/{total}: {os.path.basename(input_path)} -> {os.path.basename(output_path)}")
                    
                    # Read literature+setup CSV and prepare for processing
                    df_original = pd.read_csv(input_path)
                    print(f"Loaded {len(df_original)} literature/setup rows")
                    
                    # Backup original columns to restore after processing
                    has_reaction_id = 'reaction_id' in df_original.columns
                    has_setup_json = 'reaction_setup_json' in df_original.columns
                    has_reference = 'reference' in df_original.columns
                    
                    original_ids = None
                    setup_json_backup = None
                    reference_backup = None
                    
                    if has_reaction_id:
                        original_ids = df_original['reaction_id'].copy()
                        print(f"Preserving {len(original_ids)} original reaction IDs")
                    
                    if has_setup_json:
                        setup_json_backup = df_original['reaction_setup_json'].copy()
                        print(f"Preserving reaction_setup_json column ({len(setup_json_backup)} rows)")
                    
                    if has_reference:
                        reference_backup = df_original['reference'].copy()
                        print(f"Preserving reference column ({len(reference_backup)} rows)")
                    
                    # Process through full HTE pipeline
                    print("Processing with full HTE pipeline...")
                    process_reaction_dataset(
                        input_path,
                        output_path,
                        drop_no_catalyst=self.drop_catalyst_check.isChecked(),
                        reagent_csv_path=str(PROJECT_ROOT / "condition_registry" / "definitions" / "substances.v1.csv"),
                        new_reagents_path=str(PROJECT_ROOT / "condition_registry" / "definitions" / "pending_substances.csv"),
                        num_workers=self.worker_count_spin.value(),
                        llm_assist_options=llm_assist_options,
                    )
                    
                    # Restore original columns
                    print("Restoring original setup metadata...")
                    output_df = pd.read_csv(output_path)
                    print(f"Output loaded: {len(output_df)} rows, {len(output_df.columns)} columns")
                    
                    if len(df_original) != len(output_df):
                        print(f"⚠ Row count mismatch: input={len(df_original)}, output={len(output_df)}")
                        print("  Some rows may have been filtered by processing pipeline (e.g., drop_no_catalyst)")
                        print("  Mapping metadata by reaction_smiles...")
                        
                        # Create mapping dictionaries using reaction_smiles as key
                        if original_ids is not None:
                            id_map = dict(zip(df_original['reaction_smiles'], original_ids))
                            output_df['reaction_id'] = output_df['reaction_smiles'].map(id_map)
                            restored = output_df['reaction_id'].notna().sum()
                            print(f"✓ Restored {restored}/{len(output_df)} reaction IDs by mapping")
                        
                        if reference_backup is not None:
                            ref_map = dict(zip(df_original['reaction_smiles'], reference_backup))
                            output_df['reference'] = output_df['reaction_smiles'].map(ref_map)
                            restored = output_df['reference'].notna().sum()
                            print(f"✓ Restored {restored}/{len(output_df)} references by mapping")
                        
                        if setup_json_backup is not None:
                            setup_map = dict(zip(df_original['reaction_smiles'], setup_json_backup))
                            output_df['reaction_setup_json'] = output_df['reaction_smiles'].map(setup_map)
                            restored = output_df['reaction_setup_json'].notna().sum()
                            print(f"✓ Added reaction_setup_json to {restored}/{len(output_df)} rows by mapping")
                    else:
                        # Same length - can use direct assignment
                        if original_ids is not None:
                            output_df['reaction_id'] = original_ids
                            print(f"✓ Restored {len(original_ids)} original reaction IDs")
                        
                        if reference_backup is not None:
                            output_df['reference'] = reference_backup
                            print(f"✓ Restored reference information")
                        
                        if setup_json_backup is not None:
                            output_df['reaction_setup_json'] = setup_json_backup
                            print(f"✓ Added reaction_setup_json to {len(output_df)} rows")
                    
                    output_df = _ensure_reaction_events_column(output_df)
                    output_df.to_csv(output_path, index=False)
                    print(f"✓ Saved enhanced output with {len(output_df)} rows and {len(output_df.columns)} columns")
                    
                    print("")
                
                print(f"✓ Successfully processed {total} file(s)")
                
            finally:
                sys.stdout = old_stdout
            
            QtWidgets.QMessageBox.information(
                self, 
                "Success", 
                f"Successfully processed {len(jobs)} literature/setup file(s) with full HTE pipeline"
            )
            
        except Exception as e:
            sys.stdout = old_stdout
            error_msg = f"Error during literature/setup conversion: {str(e)}"
            self.log_msg(f"\n✗ {error_msg}")
            import traceback
            self.log_msg(traceback.format_exc())
            QtWidgets.QMessageBox.critical(self, "Error", error_msg)

    def on_finished(self, success: bool, message: str):
        self.thread.quit()
        self.thread.wait()
        
        self.setEnabled(True)
        self.progress_bar.setVisible(False)
        
        if success:
            QtWidgets.QMessageBox.information(self, "Success", message)
            self.log_msg("\n--- Processing Complete ---")
        else:
            QtWidgets.QMessageBox.warning(self, "Error", message)

def main():
    app = QtWidgets.QApplication(sys.argv)
    window = HTEConverterWindow()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()

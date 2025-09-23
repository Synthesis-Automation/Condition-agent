#!/usr/bin/env python3

"""PyQt6 UI for the reagent taxonomy generator workflow."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from PyQt6.QtCore import QObject, QRunnable, QThreadPool, pyqtSignal
from PyQt6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

MODULE_DIR = Path(__file__).resolve().parent
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

from reagent_taxonomy_generator import (  # type: ignore
    DEFAULT_RESOLVER_TIMEOUT,
    ROLE_FILES,
    RoleHeuristics,
    TaxonomyStore,
    build_entry,
    dedupe_synonyms,
    normalize_cas,
    resolve_identity_from_cas,
)

DEFAULT_TAXONOMY_DIR = (MODULE_DIR.parent / "data" / "compound_taxonomy").resolve()


class TaxonomyGenerationError(RuntimeError):
    """Raised when taxonomy generation fails for the requested input."""


def generate_taxonomy_entry(
    *,
    cas: str,
    role: Optional[str],
    taxonomy_dir: Path,
    allow_default_family: bool,
    dry_run: bool,
    resolver_timeout: float = DEFAULT_RESOLVER_TIMEOUT,
    name_override: Optional[str] = None,
    smiles_override: Optional[str] = None,
) -> Dict[str, Any]:
    """Execute the taxonomy generation workflow and return a JSON-friendly result."""
    if not cas:
        raise TaxonomyGenerationError("CAS number is required.")
    normalized_cas = normalize_cas(cas)

    if not taxonomy_dir.exists():
        raise TaxonomyGenerationError(f"Taxonomy directory '{taxonomy_dir}' does not exist.")

    store = TaxonomyStore(taxonomy_dir)
    heuristics = RoleHeuristics(store)

    resolved_identity = resolve_identity_from_cas(normalized_cas, timeout=resolver_timeout)
    name = name_override or (resolved_identity.get("name") if resolved_identity else None)
    if not name:
        raise TaxonomyGenerationError(
            "Unable to determine reagent name. Provide a name override or verify the CAS number."
        )

    abbr = name
    resolved_synonyms: Sequence[str] = resolved_identity.get("synonyms", []) if resolved_identity else []
    synonyms = dedupe_synonyms([name, abbr, *resolved_synonyms])

    family_id: Optional[str] = None
    family_reason: Optional[List[str]] = None
    role_reason: Optional[str] = None
    used_default = False
    debug_log: List[str] = []

    auto_resolve_source = resolved_identity.get("source") if resolved_identity else None
    resolved_smiles = resolved_identity.get("smiles") if resolved_identity else None
    smiles = smiles_override or resolved_smiles

    inference = heuristics.infer_family(name, synonyms)
    if inference:
        inferred_role, inferred_family, reason_tokens = inference
        if role and inferred_role != role:
            debug_log.append(
                f"Inferred family '{inferred_family}' mapped to role '{inferred_role}' "
                f"but UI selection forced role '{role}'."
            )
        else:
            role = role or inferred_role
            family_id = inferred_family
            family_reason = reason_tokens

    if not role:
        role_inference = heuristics.infer_role(name, synonyms)
        if role_inference:
            role, pattern = role_inference
            role_reason = pattern
        else:
            raise TaxonomyGenerationError("Unable to infer role; please select one manually.")

    if not family_id:
        if allow_default_family:
            default_family = heuristics.default_family_for_role(role)
            if default_family:
                family_id = default_family
                used_default = True
        if not family_id:
            raise TaxonomyGenerationError(
                "Unable to determine a family for this reagent. Enable default family fallback "
                "or provide a more specific CAS number."
            )

    family_role = store.role_for_family(family_id)
    if family_role and family_role != role:
        raise TaxonomyGenerationError(
            f"Family '{family_id}' belongs to role '{family_role}', which conflicts with role '{role}'."
        )

    existing = store.find_by_cas(normalized_cas)
    if existing:
        existing_role, existing_family, data = existing
        result = {
            "cas": normalized_cas,
            "name": data.get("name"),
            "role": existing_role,
            "family_id": existing_family,
            "status": "exists",
            "taxonomy_file": str(store.file_for_role(existing_role)),
        }
        if debug_log:
            result["debug_log"] = debug_log
        return result

    family_data = store.family_data(role, family_id)
    numeric_features = store.numeric_baseline(role, family_id)
    entry = build_entry(role, family_data, normalized_cas, name, abbr, synonyms, smiles, numeric_features)

    result: Dict[str, Any] = {
        "cas": normalized_cas,
        "name": name,
        "role": role,
        "family_id": family_id,
        "taxonomy_file": str(store.file_for_role(role)),
        "dry_run": dry_run,
        "used_default_family": used_default,
    }
    if auto_resolve_source:
        result["auto_resolve_source"] = auto_resolve_source
    if resolved_smiles and not smiles_override:
        result["smiles_source"] = auto_resolve_source or "resolver"
    if family_reason:
        result["family_tokens"] = family_reason
    if role_reason:
        result["role_pattern"] = role_reason
    if debug_log:
        result["debug_log"] = debug_log

    if dry_run:
        result["status"] = "dry_run"
        result["entry_preview"] = entry
        return result

    store.add_entry(role, family_id, entry)
    path = store.save_role(role)
    result["status"] = "written"
    result["written_to"] = str(path)
    return result


class GenerationWorkerSignals(QObject):
    """Qt signals emitted by the taxonomy generation worker."""

    finished = pyqtSignal(dict)
    error = pyqtSignal(str)


class GenerationWorker(QRunnable):
    """Run taxonomy generation logic on a background thread."""

    def __init__(self, params: Dict[str, Any]) -> None:
        super().__init__()
        self.params = params
        self.signals = GenerationWorkerSignals()

    def run(self) -> None:  # pragma: no cover - UI worker thread
        try:
            result = generate_taxonomy_entry(**self.params)
        except Exception as exc:  # noqa: BLE001 - we want to surface all failures
            self.signals.error.emit(str(exc))
        else:
            self.signals.finished.emit(result)


class TaxonomyGeneratorWindow(QMainWindow):
    """Main window that wraps the reagent taxonomy workflow."""

    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Reagent Taxonomy Generator")
        self.thread_pool = QThreadPool.globalInstance()
        self._current_worker: Optional[GenerationWorker] = None

        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)

        form_layout = QFormLayout()
        layout.addLayout(form_layout)

        self.cas_input = QLineEdit()
        self.cas_input.setPlaceholderText("e.g. 100-00-5")
        form_layout.addRow("CAS number", self.cas_input)

        self.name_input = QLineEdit()
        self.name_input.setPlaceholderText("Optional; overrides resolver result")
        form_layout.addRow("Name override", self.name_input)

        self.role_combo = QComboBox()
        self.role_combo.addItem("Select role", userData=None)
        for role_key in sorted(ROLE_FILES.keys()):
            label = role_key.replace("_", " ").title()
            self.role_combo.addItem(label, userData=role_key)
        form_layout.addRow("Reagent role", self.role_combo)

        path_layout = QHBoxLayout()
        self.taxonomy_path_input = QLineEdit(str(DEFAULT_TAXONOMY_DIR))
        browse_button = QPushButton("Browse")
        browse_button.clicked.connect(self.on_browse_taxonomy_dir)
        path_layout.addWidget(self.taxonomy_path_input)
        path_layout.addWidget(browse_button)
        form_layout.addRow("Taxonomy directory", path_layout)

        options_layout = QHBoxLayout()
        self.allow_default_checkbox = QCheckBox("Allow default family fallback")
        self.allow_default_checkbox.setChecked(True)
        self.dry_run_checkbox = QCheckBox("Dry run (preview only)")
        self.dry_run_checkbox.setChecked(True)
        options_layout.addWidget(self.allow_default_checkbox)
        options_layout.addWidget(self.dry_run_checkbox)
        options_layout.addStretch()
        layout.addLayout(options_layout)

        button_layout = QHBoxLayout()
        self.generate_button = QPushButton("Generate")
        self.generate_button.clicked.connect(self.on_generate_clicked)
        self.clear_button = QPushButton("Clear")
        self.clear_button.clicked.connect(self.clear_output)
        button_layout.addWidget(self.generate_button)
        button_layout.addWidget(self.clear_button)
        button_layout.addStretch()
        layout.addLayout(button_layout)

        self.status_label = QLabel()
        layout.addWidget(self.status_label)

        self.output_view = QPlainTextEdit()
        self.output_view.setReadOnly(True)
        self.output_view.setPlaceholderText("Generation results will appear here.")
        layout.addWidget(self.output_view, stretch=1)

        self.resize(820, 600)

    def on_browse_taxonomy_dir(self) -> None:
        """Let the user choose a taxonomy directory."""
        directory = QFileDialog.getExistingDirectory(self, "Select taxonomy directory", str(DEFAULT_TAXONOMY_DIR))
        if directory:
            self.taxonomy_path_input.setText(directory)

    def on_generate_clicked(self) -> None:
        """Kick off taxonomy generation with the current form values."""
        cas = self.cas_input.text().strip()
        if not cas:
            self.show_error("CAS number is required.")
            return

        role = self.role_combo.currentData()
        if not role:
            self.show_error("Select a reagent role before generating.")
            return

        taxonomy_dir_text = self.taxonomy_path_input.text().strip()
        if not taxonomy_dir_text:
            self.show_error("Taxonomy directory is required.")
            return

        taxonomy_dir = Path(taxonomy_dir_text).expanduser()
        allow_default = self.allow_default_checkbox.isChecked()
        dry_run = self.dry_run_checkbox.isChecked()
        name_override = self.name_input.text().strip() or None

        self.status_label.setText("Running taxonomy generation...")
        self.output_view.clear()
        self.set_inputs_enabled(False)

        params = {
            "cas": cas,
            "role": role,
            "taxonomy_dir": taxonomy_dir,
            "allow_default_family": allow_default,
            "dry_run": dry_run,
            "resolver_timeout": DEFAULT_RESOLVER_TIMEOUT,
            "name_override": name_override,
        }
        worker = GenerationWorker(params)
        worker.signals.finished.connect(self.on_generation_success)
        worker.signals.error.connect(self.on_generation_failure)
        self._current_worker = worker
        self.thread_pool.start(worker)

    def on_generation_success(self, result: Dict[str, Any]) -> None:
        """Handle successful completion of the taxonomy worker."""
        self.status_label.setText(f"Status: {result.get('status', 'ok')}")
        self.output_view.setPlainText(json.dumps(result, indent=2, ensure_ascii=False))
        self.set_inputs_enabled(True)
        self._current_worker = None

    def on_generation_failure(self, message: str) -> None:
        """Display worker errors to the user."""
        self.set_inputs_enabled(True)
        self.status_label.setText("Generation failed.")
        QMessageBox.critical(self, "Reagent taxonomy generation failed", message)
        self._current_worker = None

    def clear_output(self) -> None:
        """Clear the output preview."""
        self.output_view.clear()
        self.status_label.clear()

    def set_inputs_enabled(self, enabled: bool) -> None:
        """Toggle input widgets while the worker is running."""
        for widget in (
            self.cas_input,
            self.name_input,
            self.role_combo,
            self.taxonomy_path_input,
            self.allow_default_checkbox,
            self.dry_run_checkbox,
            self.generate_button,
            self.clear_button,
        ):
            widget.setEnabled(enabled)

    def show_error(self, message: str) -> None:
        """Show a modal error dialog."""
        QMessageBox.warning(self, "Input error", message)


def main() -> None:
    """Launch the PyQt6 application."""
    app = QApplication(sys.argv)
    window = TaxonomyGeneratorWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()

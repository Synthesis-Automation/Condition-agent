from __future__ import annotations

import io
import json
import os
import sys
from contextlib import redirect_stdout
from concurrent.futures import Future, ThreadPoolExecutor
from pathlib import Path
from typing import Any, Dict, List, Optional

from PyQt6 import QtCore, QtGui, QtWidgets

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from chemtools.featurizers.unified import featurize_molecule, featurize_reaction  # noqa: E402
from app.Cpd_rxn_featurization_cli import _print_readable  # noqa: E402

DEFAULT_LLM_PROVIDER = (
    os.getenv("CHEMTOOLS_LLM_PROVIDER", os.getenv("LLM_PROVIDER", "openai")).strip()
    or "openai"
)
DEFAULT_LLM_MODEL = os.getenv(
    "CHEMTOOLS_LLM_MODEL",
    os.getenv("LLM_MODEL", ""),
).strip()

try:
    from llmtools.clients import AVAILABLE_MODELS as LLM_AVAILABLE_MODELS  # noqa: E402
    from llmtools.clients import RECOMMENDED_MODELS as LLM_RECOMMENDED_MODELS  # noqa: E402
except Exception:  # pragma: no cover - optional dependency
    LLM_AVAILABLE_MODELS = {
        "openai": ["gpt-4o", "gpt-4o-mini", "gpt-5-mini"],
        "aliyun": ["deepseek-v3.2", "deepseek-v3", "deepseek-r1"],
    }
    LLM_RECOMMENDED_MODELS = {
        "openai": {"balanced": "gpt-4o"},
        "aliyun": {"balanced": "deepseek-v3.2"},
    }


class FeaturizationWindow(QtWidgets.QWidget):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Compound / Reaction Featurization")
        self.resize(1120, 840)

        self.mode_combo = QtWidgets.QComboBox()
        self.smiles_input = QtWidgets.QPlainTextEdit()
        self.smiles_input.setPlaceholderText(
            "Enter molecule SMILES or reaction SMILES (A.B>>P)."
        )
        self.smiles_input.setFixedHeight(64)
        self.target_groups_input = QtWidgets.QLineEdit()
        self.target_groups_input.setPlaceholderText("Optional, comma-separated (e.g. Br,CN)")

        self.include_ar_h_checkbox = QtWidgets.QCheckBox("Include Ar-H motifs")
        self.discovery_checkbox = QtWidgets.QCheckBox("Discovery mode")
        self.reactant_coverage_guard_checkbox = QtWidgets.QCheckBox(
            "Reactant coverage guard"
        )
        self.reactant_coverage_guard_checkbox.setChecked(True)
        self.show_roles_checkbox = QtWidgets.QCheckBox("Show role details in summary")
        self.show_rdkit_checkbox = QtWidgets.QCheckBox("Show RDKit properties in summary")

        self.llm_assist_checkbox = QtWidgets.QCheckBox("LLM assist (reaction only)")
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
        self.llm_no_crk_validation_checkbox = QtWidgets.QCheckBox(
            "Disable CRK validation gate"
        )

        self.run_button = QtWidgets.QPushButton("Run Featurization")
        self.clear_button = QtWidgets.QPushButton("Clear")
        self.copy_json_button = QtWidgets.QPushButton("Copy JSON")

        self.output_tabs = QtWidgets.QTabWidget()
        self.summary_output = QtWidgets.QPlainTextEdit()
        self.summary_output.setReadOnly(True)
        self.json_output = QtWidgets.QPlainTextEdit()
        self.json_output.setReadOnly(True)
        self.output_tabs.addTab(self.summary_output, "Summary")
        self.output_tabs.addTab(self.json_output, "JSON")

        self.status_label = QtWidgets.QLabel("")
        self.status_label.setStyleSheet("color: #5DADE2;")
        self._status_anim_timer = QtCore.QTimer(self)
        self._status_anim_timer.setInterval(180)
        self._status_anim_timer.timeout.connect(self._tick_status_animation)
        self._status_anim_base = "Running featurization"
        self._status_anim_step = 0

        self._executor = ThreadPoolExecutor(max_workers=1)
        self._future: Optional[Future[Dict[str, Any]]] = None
        self._poll_timer = QtCore.QTimer(self)
        self._poll_timer.setInterval(100)
        self._poll_timer.timeout.connect(self._on_poll)

        self._build_layout()
        self._set_defaults()
        self._wire_events()

    def _build_layout(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)

        title = QtWidgets.QLabel("Compound / Reaction Featurization Tool")
        title.setStyleSheet("font-size: 18px; font-weight: bold; margin: 6px;")
        title.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title)

        form = QtWidgets.QFormLayout()
        self.mode_combo.addItem("Auto-detect", "auto")
        self.mode_combo.addItem("Reaction", "reaction")
        self.mode_combo.addItem("Compound", "compound")
        form.addRow("Mode:", self.mode_combo)
        form.addRow("SMILES / Reaction:", self.smiles_input)
        form.addRow("Target groups:", self.target_groups_input)

        options_row = QtWidgets.QHBoxLayout()
        options_row.addWidget(self.include_ar_h_checkbox)
        options_row.addWidget(self.discovery_checkbox)
        options_row.addWidget(self.reactant_coverage_guard_checkbox)
        options_row.addWidget(self.show_roles_checkbox)
        options_row.addWidget(self.show_rdkit_checkbox)
        options_row.addStretch()
        form.addRow("Options:", options_row)

        llm_toggle_row = QtWidgets.QHBoxLayout()
        llm_toggle_row.addWidget(self.llm_assist_checkbox)
        llm_toggle_row.addWidget(self.llm_always_checkbox)
        llm_toggle_row.addWidget(self.llm_no_crk_validation_checkbox)
        llm_toggle_row.addStretch()
        form.addRow("LLM:", llm_toggle_row)

        llm_row = QtWidgets.QHBoxLayout()
        llm_row.addWidget(QtWidgets.QLabel("Provider"))
        llm_row.addWidget(self.llm_provider_combo)
        llm_row.addSpacing(12)
        llm_row.addWidget(QtWidgets.QLabel("Model"))
        llm_row.addWidget(self.llm_model_combo)
        llm_row.addStretch()
        form.addRow("LLM Model:", llm_row)

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

        actions = QtWidgets.QHBoxLayout()
        actions.addStretch()
        actions.addWidget(self.run_button)
        actions.addWidget(self.copy_json_button)
        actions.addWidget(self.clear_button)
        layout.addLayout(actions)

        layout.addWidget(self.output_tabs)
        layout.addWidget(self.status_label)

    def _set_defaults(self) -> None:
        self._populate_llm_providers()
        self._on_llm_assist_toggled(self.llm_assist_checkbox.isChecked())

    def _wire_events(self) -> None:
        self.run_button.clicked.connect(self._on_run)
        self.clear_button.clicked.connect(self._on_clear)
        self.copy_json_button.clicked.connect(self._on_copy_json)
        self.llm_assist_checkbox.toggled.connect(self._on_llm_assist_toggled)
        self.llm_provider_combo.currentTextChanged.connect(self._on_llm_provider_changed)

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

    def _gather_options(self) -> Dict[str, Any]:
        target_groups = [
            tg.strip()
            for tg in self.target_groups_input.text().split(",")
            if tg.strip()
        ]
        options: Dict[str, Any] = {
            "include_ar_h": self.include_ar_h_checkbox.isChecked(),
            "target_groups": target_groups or None,
            "discovery_mode": self.discovery_checkbox.isChecked(),
            "confirm_coupling_products": True,
            "reactant_coverage_guard": self.reactant_coverage_guard_checkbox.isChecked(),
        }
        if self.llm_assist_checkbox.isChecked():
            llm_model = self.llm_model_combo.currentText().strip()
            if not llm_model:
                raise ValueError("LLM model is required when LLM assist is enabled.")
            options["llm_assist"] = {
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
        return options

    def _on_run(self) -> None:
        if self._future and not self._future.done():
            return
        smiles = self.smiles_input.toPlainText().strip()
        if not smiles:
            self._show_error("SMILES / reaction input is required.")
            return
        try:
            options = self._gather_options()
        except Exception as exc:
            self._show_error(str(exc))
            return

        mode = str(self.mode_combo.currentData() or "auto").strip().lower()
        show_roles = self.show_roles_checkbox.isChecked()
        show_rdkit = self.show_rdkit_checkbox.isChecked()

        self.run_button.setEnabled(False)
        self._start_status_animation("Running featurization")
        self._future = self._executor.submit(
            self._run_featurization_job,
            smiles,
            mode,
            options,
            show_roles,
            show_rdkit,
        )
        self._poll_timer.start()

    def _run_featurization_job(
        self,
        smiles: str,
        mode: str,
        options: Dict[str, Any],
        show_roles: bool,
        show_rdkit: bool,
    ) -> Dict[str, Any]:
        is_reaction = ">" in smiles
        treat_as_reaction = (
            mode == "reaction" or (mode == "auto" and is_reaction)
        )
        if treat_as_reaction:
            reaction_options = dict(options)
            reaction_options.setdefault("motif_site_filter", "substituent")
            reaction_options["detailed"] = True
            payload = featurize_reaction(smiles, options=reaction_options)
        else:
            molecule_options = dict(options)
            molecule_options.pop("llm_assist", None)
            molecule_options["detailed"] = True
            payload = featurize_molecule(smiles, options=molecule_options)

        summary_text = self._render_summary(payload, show_roles=show_roles, show_rdkit=show_rdkit)
        json_text = json.dumps(payload, indent=2, sort_keys=True)
        return {"payload": payload, "summary": summary_text, "json": json_text}

    def _render_summary(
        self,
        payload: Dict[str, Any],
        *,
        show_roles: bool,
        show_rdkit: bool,
    ) -> str:
        buffer = io.StringIO()
        try:
            with redirect_stdout(buffer):
                _print_readable(
                    payload,
                    show_roles=show_roles,
                    show_rdkit=show_rdkit,
                    show_extended=True,
                )
            text = buffer.getvalue().strip()
            evidence_text = self._format_detection_evidence_summary(payload)
            if evidence_text:
                if text:
                    text = f"{text}\n\n{evidence_text}"
                else:
                    text = evidence_text
            return text or "(No summary output)"
        except Exception:
            return json.dumps(payload, indent=2, sort_keys=True)

    def _format_detection_evidence_summary(self, payload: Dict[str, Any]) -> str:
        detection = payload.get("detection") if isinstance(payload, dict) else None
        if not isinstance(detection, dict):
            return ""

        lines: List[str] = []
        evidence = detection.get("evidence")
        if isinstance(evidence, dict):
            lines.append("Detection Evidence")
            lines.append("-" * 72)
            matcher = str(evidence.get("matcher") or "").strip()
            if matcher:
                lines.append(f"matcher: {matcher}")

            selected = evidence.get("selected")
            if isinstance(selected, dict):
                sel_rid = str(selected.get("reaction_id") or "").strip() or "Unknown"
                lines.append(f"selected: {sel_rid}")
                sel_slots = selected.get("matched_reactant_slots") or []
                if sel_slots:
                    lines.append(f"selected_reactant_slots: {', '.join(str(s) for s in sel_slots)}")
                sel_prod_slots = selected.get("matched_product_slots")
                if sel_prod_slots is not None:
                    lines.append(f"selected_product_slots: {sel_prod_slots}")
                sel_score = selected.get("score")
                if isinstance(sel_score, list):
                    lines.append(f"selected_score: {self._format_score_vector(sel_score)}")

            top_candidates = evidence.get("top_candidates") or []
            if isinstance(top_candidates, list) and top_candidates:
                lines.append("top_candidates:")
                for idx, candidate in enumerate(top_candidates[:5], start=1):
                    if not isinstance(candidate, dict):
                        continue
                    rid = str(candidate.get("reaction_id") or "Unknown")
                    reactant_support = candidate.get("reactant_support")
                    product_support = candidate.get("product_support")
                    matched_slots = candidate.get("matched_reactant_slots") or []
                    score = candidate.get("score")
                    detail_parts: List[str] = []
                    if matched_slots:
                        detail_parts.append(f"slots={len(matched_slots)}")
                    if reactant_support is not None:
                        detail_parts.append(f"reactant_support={reactant_support}")
                    if product_support is not None:
                        detail_parts.append(f"product_support={product_support}")
                    if isinstance(score, list):
                        detail_parts.append(f"score={self._format_score_vector(score)}")
                    if detail_parts:
                        lines.append(f"  {idx}. {rid} ({'; '.join(detail_parts)})")
                    else:
                        lines.append(f"  {idx}. {rid}")

        llm_assist = detection.get("llm_assist")
        if isinstance(llm_assist, dict):
            if lines:
                lines.append("")
            lines.append("LLM Assist")
            lines.append("-" * 72)
            for key in ("used", "status", "decision"):
                if key in llm_assist:
                    lines.append(f"{key}: {llm_assist.get(key)}")

        return "\n".join(lines).strip()

    @staticmethod
    def _format_score_vector(score: List[Any]) -> str:
        return "[" + ", ".join(str(v) for v in score) + "]"

    def _on_poll(self) -> None:
        if not self._future:
            self._poll_timer.stop()
            self._stop_status_animation()
            self.run_button.setEnabled(True)
            return
        if not self._future.done():
            return

        self._poll_timer.stop()
        self._stop_status_animation()
        self.run_button.setEnabled(True)
        future = self._future
        self._future = None

        try:
            result = future.result()
        except Exception as exc:
            self._show_error(f"Featurization failed: {exc}")
            return

        self.summary_output.setPlainText(str(result.get("summary") or ""))
        self.json_output.setPlainText(str(result.get("json") or ""))
        self.status_label.setText("Done.")

    def _on_copy_json(self) -> None:
        text = self.json_output.toPlainText().strip()
        if not text:
            return
        QtWidgets.QApplication.clipboard().setText(text)
        self.status_label.setText("JSON copied to clipboard.")

    def _on_clear(self) -> None:
        self.smiles_input.clear()
        self.target_groups_input.clear()
        self.summary_output.clear()
        self.json_output.clear()
        self.status_label.clear()

    def _show_error(self, message: str) -> None:
        self._stop_status_animation()
        QtWidgets.QMessageBox.warning(self, "Featurization", message)

    def _start_status_animation(self, base_text: str) -> None:
        self._status_anim_base = base_text
        self._status_anim_step = 0
        self.status_label.setText(base_text)
        self._status_anim_timer.start()

    def _tick_status_animation(self) -> None:
        dots = "." * ((self._status_anim_step % 3) + 1)
        self.status_label.setText(f"{self._status_anim_base}{dots}")
        self._status_anim_step += 1

    def _stop_status_animation(self) -> None:
        if self._status_anim_timer.isActive():
            self._status_anim_timer.stop()

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        self._poll_timer.stop()
        self._stop_status_animation()
        self._executor.shutdown(wait=False, cancel_futures=True)
        super().closeEvent(event)


def main() -> None:
    app = QtWidgets.QApplication(sys.argv)
    window = FeaturizationWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()

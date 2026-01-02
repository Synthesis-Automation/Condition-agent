"""Reusable Qt dialogs for the ChemTools assistant UI."""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Optional

from PyQt6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QVBoxLayout,
)

DEFAULT_LLM_PROVIDER = os.getenv("LLM_PROVIDER", "openai").strip() or "openai"
DEFAULT_LLM_MODEL = os.getenv("LLM_MODEL", "gpt-4o").strip() or "gpt-4o"

try:
    from llmtools.clients import AVAILABLE_MODELS as _AVAILABLE_LLM_MODELS
    from llmtools.clients import RECOMMENDED_MODELS as _RECOMMENDED_LLM_MODELS
except Exception:  # pragma: no cover - optional dependency
    _AVAILABLE_LLM_MODELS = {
        "openai": ["gpt-4o", "gpt-4o-mini"],
        "aliyun": ["deepseek-v3", "deepseek-r1"],
    }
    _RECOMMENDED_LLM_MODELS = {
        "openai": {"balanced": "gpt-4o"},
        "aliyun": {"balanced": "deepseek-v3"},
    }


@dataclass(frozen=True)
class LLMConfig:
    provider: str
    model: str
    temperature: float


class LLMConfigDialog(QDialog):
    """Dialog for choosing the assistant LLM provider/model settings."""

    def __init__(
        self,
        *,
        provider: Optional[str] = None,
        model: Optional[str] = None,
        temperature: float = 0.0,
        parent: Optional[object] = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("LLM Mode")
        self.setModal(True)

        self._model_presets = dict(_AVAILABLE_LLM_MODELS)
        self._recommended_models = dict(_RECOMMENDED_LLM_MODELS)

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.provider_combo = QComboBox()
        self.provider_combo.addItems(sorted(self._model_presets.keys()) or ["openai", "aliyun"])
        self.provider_combo.setCurrentText((provider or DEFAULT_LLM_PROVIDER) or "openai")
        form.addRow("Provider:", self.provider_combo)

        self.model_combo = QComboBox()
        self.model_combo.setEditable(True)
        self._populate_models_for_provider(self.provider_combo.currentText(), prefer=model)
        form.addRow("Model:", self.model_combo)

        self.temperature_spin = QDoubleSpinBox()
        self.temperature_spin.setMinimum(0.0)
        self.temperature_spin.setMaximum(2.0)
        self.temperature_spin.setSingleStep(0.1)
        self.temperature_spin.setValue(float(temperature or 0.0))
        form.addRow("Temperature:", self.temperature_spin)

        layout.addLayout(form)

        self.provider_combo.currentTextChanged.connect(self._on_provider_changed)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _populate_models_for_provider(self, provider: str, *, prefer: Optional[str] = None) -> None:
        provider_key = (provider or "openai").strip().lower()
        presets = list(self._model_presets.get(provider_key) or [])

        current_text = (prefer or self.model_combo.currentText() or "").strip()
        self.model_combo.blockSignals(True)
        try:
            self.model_combo.clear()
            self.model_combo.addItems(presets)
            if current_text and current_text in presets:
                self.model_combo.setCurrentText(current_text)
            else:
                recommended = (self._recommended_models.get(provider_key) or {}).get("balanced")
                fallback = (recommended or (presets[0] if presets else "")).strip()
                if not prefer and provider_key == DEFAULT_LLM_PROVIDER.lower():
                    fallback = (DEFAULT_LLM_MODEL or fallback).strip()
                self.model_combo.setCurrentText(fallback)
        finally:
            self.model_combo.blockSignals(False)

    def _on_provider_changed(self, provider: str) -> None:
        self._populate_models_for_provider(provider)

    def get_config(self) -> LLMConfig:
        return LLMConfig(
            provider=self.provider_combo.currentText().strip() or "openai",
            model=self.model_combo.currentText().strip() or "gpt-4o",
            temperature=float(self.temperature_spin.value()),
        )

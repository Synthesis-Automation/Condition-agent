import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from PyQt6 import QtCore, QtWidgets

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from chemtools.reagent.addition_service import (  # noqa: E402
    ReagentAdditionError,
    ReagentAdditionService,
)

DEFAULT_LLM_PROVIDER = os.getenv("LLM_PROVIDER", "openai").strip() or "openai"
DEFAULT_LLM_MODEL = os.getenv("LLM_MODEL", "gpt-4o").strip() or "gpt-4o"

try:
    from llmtools.clients import AVAILABLE_MODELS as LLM_AVAILABLE_MODELS
    from llmtools.clients import RECOMMENDED_MODELS as LLM_RECOMMENDED_MODELS
except Exception:  # pragma: no cover - optional dependency
    LLM_AVAILABLE_MODELS = {
        "openai": ["gpt-4o", "gpt-4o-mini"],
        "aliyun": ["deepseek-v3", "deepseek-r1"],
    }
    LLM_RECOMMENDED_MODELS = {
        "openai": {"balanced": "gpt-4o"},
        "aliyun": {"balanced": "deepseek-v3"},
    }


class ReagentAdditionWindow(QtWidgets.QWidget):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Reagent Addition (reagents.csv)")
        self.resize(900, 860)

        self.service: Optional[ReagentAdditionService] = None

        self.registry_dir_input = QtWidgets.QLineEdit()

        self.cas_input = QtWidgets.QLineEdit()
        self.name_input = QtWidgets.QLineEdit()
        self.abbr_input = QtWidgets.QLineEdit()
        self.smiles_input = QtWidgets.QLineEdit()
        self.tag_input = QtWidgets.QLineEdit()
        self.role_combo = QtWidgets.QComboBox()
        self.family_combo = QtWidgets.QComboBox()
        self.family_combo.setEditable(True)
        self.family_combo.setInsertPolicy(QtWidgets.QComboBox.InsertPolicy.NoInsert)

        self.allow_default_checkbox = QtWidgets.QCheckBox("Allow default family fallback")
        self.auto_resolve_checkbox = QtWidgets.QCheckBox("Auto-resolve name/smiles from CAS")
        self.auto_resolve_checkbox.setChecked(True)
        self.ai_assist_checkbox = QtWidgets.QCheckBox("AI-assist")
        self.llm_provider_combo = QtWidgets.QComboBox()
        self.llm_model_combo = QtWidgets.QComboBox()
        self.llm_model_combo.setEditable(True)

        self.preview_button = QtWidgets.QPushButton("Preview (dry run)")
        self.save_button = QtWidgets.QPushButton("Save to registry")
        self.clear_button = QtWidgets.QPushButton("Clear")

        self.preview_table = QtWidgets.QTableWidget()
        self.preview_table.setColumnCount(2)
        self.preview_table.setHorizontalHeaderLabels(["Field", "Value"])
        header = self.preview_table.horizontalHeader()
        header.setStretchLastSection(True)
        header.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Stretch)
        self.preview_table.verticalHeader().setVisible(False)
        self.preview_table.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.preview_table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.DoubleClicked
            | QtWidgets.QAbstractItemView.EditTrigger.SelectedClicked
            | QtWidgets.QAbstractItemView.EditTrigger.EditKeyPressed
        )
        self.status_label = QtWidgets.QLabel("")
        self._preview_header: List[str] = []

        self._build_layout()
        self._set_defaults()
        self._wire_events()
        self._reload_taxonomy()

    def _build_layout(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)

        title = QtWidgets.QLabel("Reagent Addition Tool")
        title.setStyleSheet("font-size: 18px; font-weight: bold; margin: 6px;")
        title.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title)

        form = QtWidgets.QFormLayout()

        registry_row = QtWidgets.QHBoxLayout()
        registry_row.addWidget(self.registry_dir_input)
        registry_btn = QtWidgets.QPushButton("Browse")
        registry_btn.setFixedWidth(90)
        registry_btn.clicked.connect(self._choose_registry_dir)
        registry_row.addWidget(registry_btn)
        form.addRow("Registry folder:", registry_row)

        self.cas_input.setPlaceholderText("e.g. 14221-01-3")
        self.name_input.setPlaceholderText("Optional (auto-resolve if enabled)")
        self.smiles_input.setPlaceholderText("Optional SMILES")
        self.tag_input.setPlaceholderText("Optional tag")
        input_grid = QtWidgets.QGridLayout()
        input_grid.setHorizontalSpacing(12)
        input_grid.setVerticalSpacing(8)
        input_grid.setColumnStretch(1, 1)
        input_grid.setColumnStretch(3, 1)

        input_grid.addWidget(QtWidgets.QLabel("CAS (required):"), 0, 0)
        input_grid.addWidget(self.cas_input, 0, 1)
        input_grid.addWidget(QtWidgets.QLabel("Tag:"), 0, 2)
        input_grid.addWidget(self.tag_input, 0, 3)

        input_grid.addWidget(QtWidgets.QLabel("Name:"), 1, 0)
        input_grid.addWidget(self.name_input, 1, 1)
        input_grid.addWidget(QtWidgets.QLabel("Role:"), 1, 2)
        input_grid.addWidget(self.role_combo, 1, 3)

        input_grid.addWidget(QtWidgets.QLabel("SMILES:"), 2, 0)
        input_grid.addWidget(self.smiles_input, 2, 1)
        input_grid.addWidget(QtWidgets.QLabel("Family:"), 2, 2)
        input_grid.addWidget(self.family_combo, 2, 3)
        form.addRow(input_grid)

        self.abbr_input.setPlaceholderText("Optional abbreviation")
        form.addRow("Abbreviation:", self.abbr_input)
        checkbox_row = QtWidgets.QHBoxLayout()
        checkbox_row.addWidget(self.allow_default_checkbox)
        checkbox_row.addWidget(self.auto_resolve_checkbox)
        checkbox_row.addWidget(self.ai_assist_checkbox)
        checkbox_row.addStretch()
        form.addRow("Options:", checkbox_row)

        llm_row = QtWidgets.QHBoxLayout()
        llm_row.addWidget(QtWidgets.QLabel("Provider"))
        llm_row.addWidget(self.llm_provider_combo)
        llm_row.addSpacing(12)
        llm_row.addWidget(QtWidgets.QLabel("Model"))
        llm_row.addWidget(self.llm_model_combo)
        llm_row.addStretch()
        form.addRow("LLM:", llm_row)

        layout.addLayout(form)

        actions = QtWidgets.QHBoxLayout()
        actions.addStretch()
        actions.addWidget(self.preview_button)
        actions.addWidget(self.save_button)
        actions.addWidget(self.clear_button)
        layout.addLayout(actions)

        layout.addWidget(QtWidgets.QLabel("CSV preview (editable table):"))
        layout.addWidget(self.preview_table)
        layout.addWidget(self.status_label)

    def _set_defaults(self) -> None:
        default_registry = PROJECT_ROOT / "data" / "reagent_db"
        self.registry_dir_input.setText(str(default_registry))
        self.save_button.setEnabled(False)
        self._populate_llm_providers()
        self._on_ai_assist_toggled(self.ai_assist_checkbox.isChecked())

    def _wire_events(self) -> None:
        self.role_combo.currentIndexChanged.connect(self._on_role_changed)
        self.preview_button.clicked.connect(self._on_preview)
        self.save_button.clicked.connect(self._on_save)
        self.clear_button.clicked.connect(self._on_clear)
        self.ai_assist_checkbox.toggled.connect(self._on_ai_assist_toggled)
        self.llm_provider_combo.currentTextChanged.connect(self._on_llm_provider_changed)

    def _choose_registry_dir(self) -> None:
        directory = QtWidgets.QFileDialog.getExistingDirectory(
            self,
            "Select registry directory",
            str(PROJECT_ROOT),
        )
        if directory:
            self.registry_dir_input.setText(directory)

    def _reload_taxonomy(self) -> None:
        try:
            self.service = ReagentAdditionService(
                registry_dir=self._registry_dir(),
            )
        except Exception as exc:
            self.service = None
            self.role_combo.clear()
            self.family_combo.clear()
            QtWidgets.QMessageBox.critical(
                self,
                "Taxonomy load failed",
                f"Unable to load taxonomy data.\n\n{exc}",
            )
            return

        self._populate_roles()
        self._populate_families(None)

    def _populate_roles(self) -> None:
        if not self.service:
            return
        self.role_combo.blockSignals(True)
        self.role_combo.clear()
        self.role_combo.addItem("Auto-detect", None)
        for role in self.service.list_roles():
            label = role.get("label") or role["role"]
            display = f"{label} ({role['role']})"
            self.role_combo.addItem(display, role["role"])
        self.role_combo.blockSignals(False)

    def _populate_families(self, role: Optional[str]) -> None:
        if not self.service:
            return
        self.family_combo.blockSignals(True)
        self.family_combo.clear()
        self.family_combo.addItem("Auto-detect", None)
        families = self.service.list_families(role)
        for entry in families:
            label = entry.get("label") or entry["family_id"]
            if role:
                display = f"{label} ({entry['family_id']})"
            else:
                display = f"{label} ({entry['family_id']}) [{entry['role']}]"
            self.family_combo.addItem(display, entry["family_id"])
        self.family_combo.blockSignals(False)

    def _populate_llm_providers(self) -> None:
        providers = sorted(LLM_AVAILABLE_MODELS.keys()) or ["openai", "aliyun"]
        self.llm_provider_combo.blockSignals(True)
        self.llm_provider_combo.clear()
        self.llm_provider_combo.addItems(providers)
        self.llm_provider_combo.setCurrentText(DEFAULT_LLM_PROVIDER)
        self.llm_provider_combo.blockSignals(False)
        self._populate_llm_models(self.llm_provider_combo.currentText(), prefer=DEFAULT_LLM_MODEL)

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
                if provider_key == DEFAULT_LLM_PROVIDER.lower():
                    fallback = (DEFAULT_LLM_MODEL or fallback).strip()
                self.llm_model_combo.setCurrentText(fallback)
        finally:
            self.llm_model_combo.blockSignals(False)

    def _on_llm_provider_changed(self, provider: str) -> None:
        self._populate_llm_models(provider)

    def _on_ai_assist_toggled(self, enabled: bool) -> None:
        self.llm_provider_combo.setEnabled(enabled)
        self.llm_model_combo.setEnabled(enabled)

    def _on_role_changed(self) -> None:
        role = self.role_combo.currentData()
        self._populate_families(role)

    def _registry_dir(self) -> str:
        return self.registry_dir_input.text().strip()

    def _taxonomy_dir(self) -> Optional[str]:
        return None

    def _family_value(self) -> Optional[str]:
        data = self.family_combo.currentData()
        if data:
            return str(data).strip() or None
        text = self.family_combo.currentText().strip()
        if not text:
            return None
        if text.lower() == "auto-detect":
            return None
        if "(" in text and ")" in text:
            candidate = text.split("(")[-1].split(")")[0].strip()
            if candidate:
                return candidate
        return text

    def _gather_payload(self) -> Dict[str, Any]:
        return {
            "cas": self.cas_input.text().strip(),
            "name": self.name_input.text().strip() or None,
            "abbreviation": self.abbr_input.text().strip() or None,
            "smiles": self.smiles_input.text().strip() or None,
            "tag": self.tag_input.text().strip() or None,
            "role": self.role_combo.currentData(),
            "family_id": self._family_value(),
            "allow_default_family": self.allow_default_checkbox.isChecked(),
            "auto_resolve": self.auto_resolve_checkbox.isChecked(),
            "ai_assist": self.ai_assist_checkbox.isChecked(),
            "llm_provider": self.llm_provider_combo.currentText().strip() or None,
            "llm_model": self.llm_model_combo.currentText().strip() or None,
        }

    def _ensure_service(self) -> bool:
        if self.service is None:
            self._reload_taxonomy()
        if self.service is None:
            return False
        try:
            self.service.update_paths(registry_dir=self._registry_dir())
        except Exception:
            return False
        return True

    def _on_preview(self) -> None:
        if not self._ensure_service():
            return
        payload = self._gather_payload()
        try:
            result = self.service.preview_entry(**payload)
        except (ReagentAdditionError, ValueError) as exc:
            self._show_error(str(exc))
            return
        except Exception as exc:
            self._show_error(f"Unexpected error: {exc}")
            return

        self._populate_preview_table(result)
        self.status_label.setText(f"Status: {result.get('status', 'ok')}")
        self.save_button.setEnabled("entry_preview" in result or result.get("status") == "exists")

    def _on_save(self) -> None:
        if not self._ensure_service():
            return
        if self._table_has_data():
            try:
                row = self._row_from_table()
                result = self.service.save_csv_row(row)
            except (ReagentAdditionError, ValueError) as exc:
                self._show_error(str(exc))
                return
            except Exception as exc:
                self._show_error(f"Unexpected error: {exc}")
                return
            self._populate_preview_table(result)
            self.status_label.setText(self._format_status(result, prefix="Saved"))
            self.save_button.setEnabled(False)
            return

        payload = self._gather_payload()
        try:
            result = self.service.save_entry(**payload)
        except (ReagentAdditionError, ValueError) as exc:
            self._show_error(str(exc))
            return
        except Exception as exc:
            self._show_error(f"Unexpected error: {exc}")
            return

        self._populate_preview_table(result)
        self.status_label.setText(self._format_status(result, prefix="Saved"))
        self.save_button.setEnabled(False)

    def _on_clear(self) -> None:
        self.cas_input.clear()
        self.name_input.clear()
        self.abbr_input.clear()
        self.smiles_input.clear()
        self.tag_input.clear()
        self.role_combo.setCurrentIndex(0)
        self.family_combo.setCurrentIndex(0)
        self._clear_preview_table()
        self.status_label.clear()
        self.save_button.setEnabled(False)

    def _populate_preview_table(self, result: Dict[str, Any]) -> None:
        if not self.service:
            self._clear_preview_table()
            return

        entry = None
        role = result.get("role")
        status = result.get("status", "")
        if status == "exists":
            cas = result.get("cas")
            if cas:
                try:
                    existing = self.service.find_existing(cas)
                    if existing and existing.get("entry"):
                        entry = existing["entry"]
                        role = existing.get("role")
                except Exception:
                    entry = None

        if entry is None:
            candidate = result.get("entry_preview")
            if isinstance(candidate, dict):
                entry = candidate

        if not entry:
            self._clear_preview_table()
            return

        header, row = self.service.build_csv_row(entry, role=role)
        self._preview_header = header
        self.preview_table.setRowCount(len(header))
        for idx, field in enumerate(header):
            field_item = QtWidgets.QTableWidgetItem(field)
            field_item.setFlags(field_item.flags() & ~QtCore.Qt.ItemFlag.ItemIsEditable)
            self.preview_table.setItem(idx, 0, field_item)
            value_item = QtWidgets.QTableWidgetItem(row.get(field, ""))
            self.preview_table.setItem(idx, 1, value_item)

        self.preview_table.resizeColumnsToContents()
        self._resize_preview_table()

    def _clear_preview_table(self) -> None:
        self.preview_table.clearContents()
        self.preview_table.setRowCount(0)
        self._preview_header = []
        self.preview_table.setMinimumHeight(0)
        self.preview_table.setMaximumHeight(16777215)

    def _table_has_data(self) -> bool:
        return self.preview_table.rowCount() > 0

    def _row_from_table(self) -> Dict[str, str]:
        header = self._preview_header
        if not header:
            header = []
            for row in range(self.preview_table.rowCount()):
                item = self.preview_table.item(row, 0)
                if item:
                    header.append(item.text().strip())
        row_data: Dict[str, str] = {}
        for idx, field in enumerate(header):
            value_item = self.preview_table.item(idx, 1)
            row_data[field] = value_item.text() if value_item else ""
        return row_data

    def _resize_preview_table(self) -> None:
        header = self.preview_table.horizontalHeader()
        total_height = header.height()
        for row in range(self.preview_table.rowCount()):
            total_height += self.preview_table.rowHeight(row)
        total_height += self.preview_table.frameWidth() * 2
        self.preview_table.setMinimumHeight(total_height)
        self.preview_table.setMaximumHeight(total_height)
        base = max(self.height() - self.preview_table.height(), 240)
        desired = total_height + base
        if desired > self.height():
            self.resize(self.width(), desired)

    def _format_status(self, result: Dict[str, Any], *, prefix: str = "Status") -> str:
        status = result.get("status", "ok")
        name = None
        entry = result.get("entry_preview")
        if isinstance(entry, dict):
            name = entry.get("name")
        if not name:
            name = result.get("name")
        if name:
            return f"{prefix}: {name} ({status})"
        return f"{prefix}: {status}"

    def _show_error(self, message: str) -> None:
        QtWidgets.QMessageBox.warning(self, "Reagent addition", message)


def main() -> None:
    app = QtWidgets.QApplication(sys.argv)
    window = ReagentAdditionWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()

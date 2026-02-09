import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from concurrent.futures import Future, ThreadPoolExecutor

from PyQt6 import QtCore, QtGui, QtWidgets

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
        self.resize(920, 760)

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
        self.ai_assist_checkbox.setChecked(True)
        self.llm_provider_combo = QtWidgets.QComboBox()
        self.llm_model_combo = QtWidgets.QComboBox()
        self.llm_model_combo.setEditable(True)

        self.preview_button = QtWidgets.QPushButton("Preview (dry run)")
        self.save_button = QtWidgets.QPushButton("Save to registry")
        self.clear_button = QtWidgets.QPushButton("Clear")

        self.preview_table = QtWidgets.QTableWidget()
        self.preview_table.setColumnCount(4)
        self.preview_table.setHorizontalHeaderLabels(["Field", "Value", "Field", "Value"])
        header = self.preview_table.horizontalHeader()
        header.setStretchLastSection(False)
        for col in range(4):
            header.setSectionResizeMode(col, QtWidgets.QHeaderView.ResizeMode.Interactive)
        self.preview_table.verticalHeader().setVisible(False)
        self.preview_table.setHorizontalScrollBarPolicy(
            QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        self.preview_table.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.preview_table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.DoubleClicked
            | QtWidgets.QAbstractItemView.EditTrigger.SelectedClicked
            | QtWidgets.QAbstractItemView.EditTrigger.EditKeyPressed
        )
        self.status_label = QtWidgets.QLabel("")
        self.status_label.setStyleSheet("color: #5DADE2;")
        self._preview_header: List[str] = []
        self._preview_executor = ThreadPoolExecutor(max_workers=1)
        self._preview_future: Optional[Future[Dict[str, Any]]] = None
        self._preview_poll_timer = QtCore.QTimer(self)
        self._preview_poll_timer.setInterval(100)
        self._preview_poll_timer.timeout.connect(self._on_preview_poll)
        self._status_anim_timer = QtCore.QTimer(self)
        self._status_anim_timer.setInterval(180)
        self._status_anim_timer.timeout.connect(self._tick_status_animation)
        self._status_anim_base = "Dry run in progress"
        self._status_anim_step = 0

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
        if self._preview_future and not self._preview_future.done():
            return
        if not self._ensure_service():
            return
        payload = self._gather_payload()
        self.preview_button.setEnabled(False)
        self.save_button.setEnabled(False)
        self._start_status_animation("Dry run in progress")
        try:
            self._preview_future = self._preview_executor.submit(self.service.preview_entry, **payload)
            self._preview_poll_timer.start()
        except Exception as exc:
            self._stop_status_animation()
            self.preview_button.setEnabled(True)
            self._show_error(f"Unexpected error: {exc}")
            return

    def _on_preview_poll(self) -> None:
        future = self._preview_future
        if future is None:
            self._preview_poll_timer.stop()
            self._stop_status_animation()
            self.preview_button.setEnabled(True)
            return
        if not future.done():
            return

        self._preview_poll_timer.stop()
        self._stop_status_animation()
        self.preview_button.setEnabled(True)
        self._preview_future = None

        try:
            result = future.result()
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
        field_positions: Dict[str, Tuple[int, int]] = {}
        row_count = (len(header) + 1) // 2
        self.preview_table.setRowCount(row_count)
        for idx in range(row_count):
            left_field_idx = idx * 2
            right_field_idx = left_field_idx + 1

            left_field = header[left_field_idx]
            field_positions[left_field] = (idx, 1)
            field_item = QtWidgets.QTableWidgetItem(left_field)
            field_item.setFlags(field_item.flags() & ~QtCore.Qt.ItemFlag.ItemIsEditable)
            self.preview_table.setItem(idx, 0, field_item)
            value_item = QtWidgets.QTableWidgetItem(row.get(left_field, ""))
            self.preview_table.setItem(idx, 1, value_item)

            if right_field_idx < len(header):
                right_field = header[right_field_idx]
                field_positions[right_field] = (idx, 3)
                field_item_2 = QtWidgets.QTableWidgetItem(right_field)
                field_item_2.setFlags(field_item_2.flags() & ~QtCore.Qt.ItemFlag.ItemIsEditable)
                self.preview_table.setItem(idx, 2, field_item_2)
                value_item_2 = QtWidgets.QTableWidgetItem(row.get(right_field, ""))
                self.preview_table.setItem(idx, 3, value_item_2)
            else:
                empty_field_item = QtWidgets.QTableWidgetItem("")
                empty_field_item.setFlags(
                    empty_field_item.flags() & ~QtCore.Qt.ItemFlag.ItemIsEditable
                )
                self.preview_table.setItem(idx, 2, empty_field_item)
                self.preview_table.setItem(idx, 3, QtWidgets.QTableWidgetItem(""))

        self._attach_preview_role_family_selectors(row, field_positions)
        self._adjust_preview_column_widths()
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
                for field_col in (0, 2):
                    item = self.preview_table.item(row, field_col)
                    if item:
                        field = item.text().strip()
                        if field:
                            header.append(field)
        row_data: Dict[str, str] = {}
        for row in range(self.preview_table.rowCount()):
            left_field_item = self.preview_table.item(row, 0)
            if left_field_item:
                left_field = left_field_item.text().strip()
                if left_field:
                    row_data[left_field] = self._preview_value_from_cell(row, 1)

            right_field_item = self.preview_table.item(row, 2)
            if right_field_item:
                right_field = right_field_item.text().strip()
                if right_field:
                    row_data[right_field] = self._preview_value_from_cell(row, 3)
        return row_data

    def _preview_value_from_cell(self, row: int, value_col: int) -> str:
        widget = self.preview_table.cellWidget(row, value_col)
        if isinstance(widget, QtWidgets.QComboBox):
            data = widget.currentData()
            if data is None:
                return ""
            return str(data).strip()
        value_item = self.preview_table.item(row, value_col)
        return value_item.text() if value_item else ""

    def _attach_preview_role_family_selectors(
        self,
        row_data: Dict[str, str],
        field_positions: Dict[str, Tuple[int, int]],
    ) -> None:
        if not self.service:
            return
        role_pos = field_positions.get("role_1")
        family_pos = field_positions.get("family_1")
        if not role_pos and not family_pos:
            return

        role_combo: Optional[QtWidgets.QComboBox] = None
        family_combo: Optional[QtWidgets.QComboBox] = None
        role_value = str(row_data.get("role_1") or "").strip()
        family_value = str(row_data.get("family_1") or "").strip()

        if role_pos:
            role_combo = QtWidgets.QComboBox()
            role_combo.addItem("Auto-detect", "")
            role_ids: set[str] = set()
            for entry in self.service.list_roles():
                role_id = str(entry.get("role") or "").strip()
                if not role_id or role_id in role_ids:
                    continue
                role_ids.add(role_id)
                label = str(entry.get("label") or role_id).strip()
                role_combo.addItem(f"{label} ({role_id})", role_id)
            if role_value and role_combo.findData(role_value) < 0:
                role_combo.addItem(role_value, role_value)
            role_combo.setCurrentIndex(max(role_combo.findData(role_value), 0))
            self.preview_table.setCellWidget(role_pos[0], role_pos[1], role_combo)

        if family_pos:
            family_combo = QtWidgets.QComboBox()
            self._populate_preview_family_combo(
                family_combo,
                role_value or "",
                preferred_family=family_value,
            )
            self.preview_table.setCellWidget(family_pos[0], family_pos[1], family_combo)

        if role_combo and family_combo:
            role_combo.currentIndexChanged.connect(
                lambda _idx: self._populate_preview_family_combo(
                    family_combo,
                    str(role_combo.currentData() or "").strip(),
                    preferred_family="",
                )
            )

    def _populate_preview_family_combo(
        self,
        combo: QtWidgets.QComboBox,
        role_value: str,
        *,
        preferred_family: str,
    ) -> None:
        if not self.service:
            return
        combo.blockSignals(True)
        try:
            combo.clear()
            combo.addItem("Auto-detect", "")
            family_ids: set[str] = set()
            role_arg = role_value or None
            for entry in self.service.list_families(role_arg):
                family_id = str(entry.get("family_id") or "").strip()
                if not family_id or family_id in family_ids:
                    continue
                family_ids.add(family_id)
                label = str(entry.get("label") or family_id).strip()
                combo.addItem(f"{label} ({family_id})", family_id)
            if preferred_family and combo.findData(preferred_family) < 0:
                combo.addItem(preferred_family, preferred_family)
            combo.setCurrentIndex(max(combo.findData(preferred_family), 0))
        finally:
            combo.blockSignals(False)

    def _resize_preview_table(self) -> None:
        self._adjust_preview_column_widths()
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

    def _adjust_preview_column_widths(self) -> None:
        if self.preview_table.columnCount() != 4:
            return
        viewport_width = self.preview_table.viewport().width()
        if viewport_width <= 0:
            return

        field_width = int(viewport_width * 0.18)
        value_width = int(viewport_width * 0.32)
        widths = [field_width, value_width, field_width, value_width]
        remainder = viewport_width - sum(widths)
        widths[-1] += remainder

        for col, width in enumerate(widths):
            self.preview_table.setColumnWidth(col, max(60, width))

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:  # type: ignore[name-defined]
        super().resizeEvent(event)
        self._adjust_preview_column_widths()

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
        self._stop_status_animation()
        QtWidgets.QMessageBox.warning(self, "Reagent addition", message)

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
        self._preview_poll_timer.stop()
        self._stop_status_animation()
        self._preview_executor.shutdown(wait=False, cancel_futures=True)
        super().closeEvent(event)


def main() -> None:
    app = QtWidgets.QApplication(sys.argv)
    window = ReagentAdditionWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()

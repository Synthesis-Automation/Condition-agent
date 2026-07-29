"""PyQt6 wrapper for the single-event reaction-template registry."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Optional

from PyQt6 import QtCore, QtGui, QtWidgets

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from reactive_taxonomy.reaction_templates import (  # noqa: E402
    DEFAULT_REACTION_TEMPLATE_REGISTRY_PATH,
    ReactionTemplate,
    ReactionTemplateError,
    derive_reaction_template,
    load_reaction_template_registry,
    match_reaction_templates,
    upsert_reaction_template,
    validate_reaction_template_registry,
)


class ReactionTemplateRegistryWindow(QtWidgets.QMainWindow):
    """Author, inspect, validate, and query a reaction-template registry."""

    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Reaction Template Registry")
        self.resize(1120, 820)
        self.setFont(QtGui.QFont("Segoe UI", 9))
        self._templates: tuple[ReactionTemplate, ...] = ()

        self.registry_edit = QtWidgets.QLineEdit(
            str(DEFAULT_REACTION_TEMPLATE_REGISTRY_PATH)
        )
        self.registry_edit.setObjectName("registryPath")
        self.template_id_edit = QtWidgets.QLineEdit()
        self.template_id_edit.setObjectName("templateId")
        self.template_id_edit.setPlaceholderText("carbonyl_to_dialkoxy")
        self.name_edit = QtWidgets.QLineEdit()
        self.name_edit.setObjectName("templateName")
        self.name_edit.setPlaceholderText("Carbonyl to dialkoxy")
        self.family_edit = QtWidgets.QLineEdit()
        self.family_edit.setObjectName("familyId")
        self.family_edit.setPlaceholderText("acetalization (optional)")
        self.aliases_edit = QtWidgets.QLineEdit()
        self.aliases_edit.setObjectName("aliases")
        self.aliases_edit.setPlaceholderText(
            "comma-separated aliases (optional)"
        )
        self.transformation_edit = QtWidgets.QLineEdit()
        self.transformation_edit.setObjectName("transformationClass")
        self.transformation_edit.setPlaceholderText(
            "generic transformation class (optional)"
        )
        self.status_combo = QtWidgets.QComboBox()
        self.status_combo.setObjectName("templateStatus")
        self.status_combo.addItems(("draft", "active", "retired"))
        self.mapped_reaction_edit = QtWidgets.QPlainTextEdit()
        self.mapped_reaction_edit.setObjectName("mappedReferenceReaction")
        self.mapped_reaction_edit.setPlaceholderText(
            "Paste one fully atom-mapped, single-event reaction SMILES"
        )
        self.mapped_reaction_edit.setMinimumHeight(110)
        self.notes_edit = QtWidgets.QPlainTextEdit()
        self.notes_edit.setObjectName("templateNotes")
        self.notes_edit.setPlaceholderText(
            "Provenance or generalization decisions (optional)"
        )
        self.notes_edit.setMaximumHeight(80)
        self.replace_check = QtWidgets.QCheckBox(
            "Replace an existing template with the same ID"
        )
        self.replace_check.setObjectName("replaceExisting")

        self.import_button = QtWidgets.QPushButton("Import mapped reference")
        self.import_button.setObjectName("importTemplate")
        self.import_button.setDefault(True)
        self.validate_button = QtWidgets.QPushButton("Validate registry")
        self.validate_button.setObjectName("validateRegistry")
        self.refresh_button = QtWidgets.QPushButton("Refresh")
        self.refresh_button.setObjectName("refreshRegistry")

        self.table = QtWidgets.QTableWidget(0, 6)
        self.table.setObjectName("templateTable")
        self.table.setHorizontalHeaderLabels(
            ("ID", "Status", "Family", "Name", "Edits", "Fingerprint")
        )
        self.table.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows
        )
        self.table.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.SingleSelection
        )
        self.table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
        )
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.horizontalHeader().setSectionResizeMode(
            3, QtWidgets.QHeaderView.ResizeMode.Stretch
        )

        self.query_edit = QtWidgets.QPlainTextEdit()
        self.query_edit.setObjectName("queryReaction")
        self.query_edit.setPlaceholderText(
            "Paste a query reaction SMILES to derive its signature and match "
            "registered edit templates"
        )
        self.query_edit.setMaximumHeight(90)
        self.include_drafts_check = QtWidgets.QCheckBox(
            "Include draft templates"
        )
        self.include_drafts_check.setChecked(True)
        self.match_button = QtWidgets.QPushButton("Match query")
        self.match_button.setObjectName("matchQuery")

        self.details = QtWidgets.QPlainTextEdit()
        self.details.setObjectName("templateDetails")
        self.details.setReadOnly(True)
        self.status_label = QtWidgets.QLabel()
        self.status_label.setObjectName("registryStatus")
        self.status_label.setWordWrap(True)

        self._build_layout()
        self._connect_signals()
        self.refresh_registry()

    def _build_layout(self) -> None:
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QVBoxLayout(central)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(10)

        title = QtWidgets.QLabel("Expandable Reaction-Template Registry")
        title.setStyleSheet("font-size: 20px; font-weight: 600;")
        layout.addWidget(title)
        description = QtWidgets.QLabel(
            "Templates are compiled from fully mapped, single connected edit "
            "events. They provide structural interpretation candidates; query "
            "signatures are always generated from query evidence."
        )
        description.setWordWrap(True)
        layout.addWidget(description)

        registry_row = QtWidgets.QHBoxLayout()
        registry_row.addWidget(QtWidgets.QLabel("Registry"))
        registry_row.addWidget(self.registry_edit, stretch=1)
        browse = QtWidgets.QPushButton("Browse…")
        browse.clicked.connect(self.choose_registry)
        registry_row.addWidget(browse)
        registry_row.addWidget(self.refresh_button)
        registry_row.addWidget(self.validate_button)
        layout.addLayout(registry_row)

        splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        layout.addWidget(splitter, stretch=1)

        authoring = QtWidgets.QWidget()
        authoring_layout = QtWidgets.QVBoxLayout(authoring)
        form = QtWidgets.QFormLayout()
        identity_row = QtWidgets.QHBoxLayout()
        identity_row.addWidget(self.template_id_edit, stretch=2)
        identity_row.addWidget(QtWidgets.QLabel("Name"))
        identity_row.addWidget(self.name_edit, stretch=3)
        form.addRow("Template ID", identity_row)
        family_row = QtWidgets.QHBoxLayout()
        family_row.addWidget(self.family_edit)
        family_row.addWidget(QtWidgets.QLabel("Transformation"))
        family_row.addWidget(self.transformation_edit)
        family_row.addWidget(QtWidgets.QLabel("Status"))
        family_row.addWidget(self.status_combo)
        form.addRow("Family", family_row)
        form.addRow("Aliases", self.aliases_edit)
        form.addRow("Mapped reference", self.mapped_reaction_edit)
        form.addRow("Notes", self.notes_edit)
        authoring_layout.addLayout(form)
        action_row = QtWidgets.QHBoxLayout()
        action_row.addWidget(self.import_button)
        action_row.addWidget(self.replace_check)
        action_row.addStretch()
        authoring_layout.addLayout(action_row)
        splitter.addWidget(authoring)

        browser = QtWidgets.QWidget()
        browser_layout = QtWidgets.QVBoxLayout(browser)
        browser_layout.addWidget(QtWidgets.QLabel("Registered templates"))
        browser_layout.addWidget(self.table, stretch=2)
        query_row = QtWidgets.QHBoxLayout()
        query_row.addWidget(self.query_edit, stretch=1)
        query_actions = QtWidgets.QVBoxLayout()
        query_actions.addWidget(self.include_drafts_check)
        query_actions.addWidget(self.match_button)
        query_actions.addStretch()
        query_row.addLayout(query_actions)
        browser_layout.addLayout(query_row)
        browser_layout.addWidget(QtWidgets.QLabel("Details / query result"))
        browser_layout.addWidget(self.details, stretch=1)
        splitter.addWidget(browser)
        splitter.setSizes((330, 440))

        self.status_label.setStyleSheet(
            "background: #eef4fa; border: 1px solid #ccd9e5; "
            "color: #23313f; padding: 7px; border-radius: 4px;"
        )
        layout.addWidget(self.status_label)

    def _connect_signals(self) -> None:
        self.import_button.clicked.connect(self.import_template)
        self.validate_button.clicked.connect(self.validate_registry)
        self.refresh_button.clicked.connect(self.refresh_registry)
        self.match_button.clicked.connect(self.match_query)
        self.table.itemSelectionChanged.connect(self.show_selected_template)
        self.registry_edit.editingFinished.connect(self.refresh_registry)

    def _registry_path(self) -> Path:
        text = self.registry_edit.text().strip()
        if not text:
            raise ReactionTemplateError("Choose a registry JSON path")
        return Path(text)

    @QtCore.pyqtSlot()
    def choose_registry(self) -> None:
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Choose reaction-template registry",
            self.registry_edit.text()
            or str(DEFAULT_REACTION_TEMPLATE_REGISTRY_PATH),
            "JSON files (*.json);;All files (*)",
        )
        if filename:
            self.registry_edit.setText(filename)
            self.refresh_registry()

    def _show_error(self, title: str, exc: Exception) -> None:
        self.status_label.setText(f"{title}: {exc}")
        QtWidgets.QMessageBox.warning(self, title, str(exc))

    @QtCore.pyqtSlot()
    def import_template(self) -> None:
        aliases = tuple(
            value.strip()
            for value in self.aliases_edit.text().split(",")
            if value.strip()
        )
        try:
            template = derive_reaction_template(
                self.mapped_reaction_edit.toPlainText().strip(),
                template_id=self.template_id_edit.text().strip(),
                display_name=self.name_edit.text().strip(),
                family_id=self.family_edit.text().strip() or None,
                aliases=aliases,
                transformation_class=(
                    self.transformation_edit.text().strip() or None
                ),
                status=self.status_combo.currentText(),  # type: ignore[arg-type]
                provenance="qt6_manual_mapped_reference",
                notes=self.notes_edit.toPlainText().strip(),
            )
            path = upsert_reaction_template(
                template,
                self._registry_path(),
                replace_existing=self.replace_check.isChecked(),
            )
        except (OSError, ReactionTemplateError) as exc:
            self._show_error("Template import failed", exc)
            return
        self.status_label.setText(
            f"Saved {template.template_id} to {path.resolve()} — "
            f"{len(template.edits)} edits, {template.edit_archetype}."
        )
        self.refresh_registry(select_id=template.template_id)

    @QtCore.pyqtSlot()
    def validate_registry(self) -> None:
        try:
            errors = validate_reaction_template_registry(self._registry_path())
        except OSError as exc:
            self._show_error("Registry validation failed", exc)
            return
        if errors:
            self.details.setPlainText("\n".join(errors))
            self.status_label.setText(
                f"Registry invalid: {len(errors)} error(s)."
            )
        else:
            self.status_label.setText(
                f"Registry valid: {len(self._templates)} template(s)."
            )

    def refresh_registry(self, *, select_id: Optional[str] = None) -> None:
        try:
            self._templates = load_reaction_template_registry(
                self._registry_path()
            )
        except (OSError, ReactionTemplateError) as exc:
            self._templates = ()
            self.table.setRowCount(0)
            self.status_label.setText(f"Cannot load registry: {exc}")
            return
        self.table.setRowCount(len(self._templates))
        selected_row = None
        for row, template in enumerate(self._templates):
            values = (
                template.template_id,
                template.status,
                template.family_id or "",
                template.display_name,
                str(len(template.edits)),
                template.edit_fingerprint,
            )
            for column, value in enumerate(values):
                item = QtWidgets.QTableWidgetItem(value)
                item.setData(
                    QtCore.Qt.ItemDataRole.UserRole, template.template_id
                )
                self.table.setItem(row, column, item)
            if template.template_id == select_id:
                selected_row = row
        self.table.resizeColumnsToContents()
        self.table.horizontalHeader().setSectionResizeMode(
            3, QtWidgets.QHeaderView.ResizeMode.Stretch
        )
        self.table.horizontalHeader().setSectionResizeMode(
            5, QtWidgets.QHeaderView.ResizeMode.Stretch
        )
        self.status_label.setText(
            f"Loaded {len(self._templates)} template(s) from "
            f"{self._registry_path().resolve()}."
        )
        if selected_row is not None:
            self.table.selectRow(selected_row)

    @QtCore.pyqtSlot()
    def show_selected_template(self) -> None:
        selected = self.table.selectedItems()
        if not selected:
            return
        template_id = selected[0].data(QtCore.Qt.ItemDataRole.UserRole)
        template = next(
            (
                item
                for item in self._templates
                if item.template_id == template_id
            ),
            None,
        )
        if template is not None:
            self.details.setPlainText(
                json.dumps(
                    template.to_dict(),
                    ensure_ascii=False,
                    sort_keys=True,
                    indent=2,
                )
            )

    @QtCore.pyqtSlot()
    def match_query(self) -> None:
        reaction = self.query_edit.toPlainText().strip()
        if not reaction:
            self._show_error(
                "Query required",
                ReactionTemplateError("Paste a query reaction SMILES"),
            )
            return
        try:
            result = match_reaction_templates(
                reaction,
                path=self._registry_path(),
                include_drafts=self.include_drafts_check.isChecked(),
            )
        except (OSError, ReactionTemplateError) as exc:
            self._show_error("Query match failed", exc)
            return
        self.details.setPlainText(
            json.dumps(
                result.to_dict(),
                ensure_ascii=False,
                sort_keys=True,
                indent=2,
            )
        )
        self.status_label.setText(
            f"Query: signature {result.signature_id or 'unavailable'}; "
            f"{len(result.matches)} template match(es); evidence "
            f"{result.evidence}."
        )


def main() -> None:
    """Launch the reaction-template registry Qt6 application."""
    application = QtWidgets.QApplication(sys.argv)
    window = ReactionTemplateRegistryWindow()
    window.show()
    raise SystemExit(application.exec())


if __name__ == "__main__":
    main()


__all__ = ["ReactionTemplateRegistryWindow", "main"]

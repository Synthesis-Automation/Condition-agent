"""PyQt6 app for validated condition-registry compound additions."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

from PyQt6 import QtCore, QtGui, QtWidgets

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from condition_registry import (  # noqa: E402
    CompoundAdditionError,
    CompoundAdditionRequest,
    CompoundAliasInput,
    RoleAssignment,
    Substance,
    add_compound,
    resolve_substance,
    update_compound,
)
from condition_registry.loader import load_taxonomy  # noqa: E402
from condition_registry.models import (  # noqa: E402
    CONDITION_NAME_IDENTIFIER_TYPES,
)
from cas_tools import (  # noqa: E402
    CompoundLookupResult,
    is_valid_cas_number,
    lookup_compound_by_cas,
)

ALIAS_TYPES = tuple(
    value
    for value in CONDITION_NAME_IDENTIFIER_TYPES
    if value not in {"canonical_name", "abbreviation"}
) + ("inchi_key", "database_id")


class CompoundLookupWorker(QtCore.QObject):
    """Resolve web metadata without blocking the Qt event loop."""

    finished = QtCore.pyqtSignal(object, str)

    def __init__(self, cas: str) -> None:
        super().__init__()
        self.cas = cas

    @QtCore.pyqtSlot()
    def run(self) -> None:
        try:
            result = lookup_compound_by_cas(self.cas)
        except Exception as error:
            self.finished.emit(None, f"{type(error).__name__}: {error}")
        else:
            self.finished.emit(result, "")


class CompoundRegistryWindow(QtWidgets.QWidget):
    """Curator-facing form for adding one validated registry compound."""

    def __init__(self) -> None:
        super().__init__()
        self.setFont(QtGui.QFont("Segoe UI", 9))
        self.setWindowTitle("Condition Registry — Add Compound")
        self.resize(940, 820)
        self.setWindowState(
            self.windowState() | QtCore.Qt.WindowState.WindowMaximized
        )
        self._taxonomy = load_taxonomy()
        self._families_by_role = self._build_family_map()
        self.lookup_thread: Optional[QtCore.QThread] = None
        self.lookup_worker: Optional[CompoundLookupWorker] = None
        self.last_lookup_result: Optional[CompoundLookupResult] = None
        self.editing_substance_id: Optional[str] = None

        self.cas_edit = QtWidgets.QLineEdit()
        self.cas_edit.setObjectName("casNumber")
        self.cas_edit.setPlaceholderText("e.g. 64-17-5")
        self.cas_check_button = QtWidgets.QPushButton("Check CAS")
        self.cas_check_button.setObjectName("checkCasButton")
        self.lookup_button = QtWidgets.QPushButton("Auto-fill from Web")
        self.lookup_button.setObjectName("autoFillCompoundButton")
        self.lookup_button.setToolTip(
            "Look up identity and physical fields using PubChem with an NCI fallback"
        )
        self.cas_status = QtWidgets.QLabel("Not checked")
        self.cas_status.setObjectName("casStatus")

        self.name_edit = QtWidgets.QLineEdit()
        self.name_edit.setObjectName("canonicalName")
        self.name_edit.setPlaceholderText("Preferred canonical compound name")
        self.abbreviation_edit = QtWidgets.QLineEdit()
        self.abbreviation_edit.setObjectName("abbreviation")
        self.abbreviation_edit.setPlaceholderText("Optional, e.g. EtOH")
        self.smiles_edit = QtWidgets.QLineEdit()
        self.smiles_edit.setObjectName("smiles")
        self.smiles_edit.setPlaceholderText("Optional; validated and canonicalized with RDKit")
        self.formula_edit = QtWidgets.QLineEdit()
        self.formula_edit.setObjectName("formula")
        self.formula_edit.setPlaceholderText("Derived from SMILES when omitted")
        self.molecular_weight_edit = self._numeric_edit(
            "molecularWeight", "Derived from SMILES when omitted", positive=True
        )
        self.kind_combo = QtWidgets.QComboBox()
        self.kind_combo.setObjectName("substanceKind")
        self.kind_combo.setEditable(True)
        self.kind_combo.addItems(
            (
                "",
                "solid",
                "powder",
                "liquid",
                "gas",
                "solution",
                "mixture",
                "polymer",
                "supported catalyst",
            )
        )
        self.density_edit = self._numeric_edit("density", "Optional g/mL", positive=True)
        self.boiling_point_edit = self._numeric_edit(
            "boilingPoint", "Optional °C", positive=False
        )
        self.melting_point_edit = self._numeric_edit(
            "meltingPoint", "Optional °C", positive=False
        )

        self.primary_role_combo = self._role_combo("primaryRole", optional=False)
        self.primary_family_combo = QtWidgets.QComboBox()
        self.primary_family_combo.setObjectName("primaryFamily")
        self.primary_tag_edit = QtWidgets.QLineEdit()
        self.primary_tag_edit.setObjectName("primaryRoleTag")
        self.primary_tag_edit.setPlaceholderText("Optional chemistry note")
        self.secondary_role_combo = self._role_combo("secondaryRole", optional=True)
        self.secondary_family_combo = QtWidgets.QComboBox()
        self.secondary_family_combo.setObjectName("secondaryFamily")
        self.secondary_tag_edit = QtWidgets.QLineEdit()
        self.secondary_tag_edit.setObjectName("secondaryRoleTag")
        self.secondary_tag_edit.setPlaceholderText("Optional chemistry note")

        self.alias_table = QtWidgets.QTableWidget(0, 5)
        self.alias_table.setObjectName("aliasTable")
        self.alias_table.setHorizontalHeaderLabels(
            (
                "Identifier type",
                "Alias / identifier",
                "Language",
                "Preferred",
                "Shared",
            )
        )
        self.alias_table.horizontalHeader().setSectionResizeMode(
            0, QtWidgets.QHeaderView.ResizeMode.ResizeToContents
        )
        self.alias_table.horizontalHeader().setSectionResizeMode(
            1, QtWidgets.QHeaderView.ResizeMode.Stretch
        )
        self.alias_table.horizontalHeader().setSectionResizeMode(
            2, QtWidgets.QHeaderView.ResizeMode.ResizeToContents
        )
        self.alias_table.horizontalHeader().setSectionResizeMode(
            3, QtWidgets.QHeaderView.ResizeMode.ResizeToContents
        )
        self.alias_table.horizontalHeader().setSectionResizeMode(
            4, QtWidgets.QHeaderView.ResizeMode.ResizeToContents
        )
        self.alias_table.verticalHeader().setVisible(False)
        self.alias_table.setAlternatingRowColors(True)
        self.alias_table.setMinimumHeight(150)

        self.source_edit = QtWidgets.QLineEdit()
        self.source_edit.setObjectName("provenanceSource")
        self.source_edit.setPlaceholderText(
            "Required: DOI, database record, vendor catalogue, or curator reference"
        )
        self.notes_edit = QtWidgets.QPlainTextEdit()
        self.notes_edit.setObjectName("curatorNotes")
        self.notes_edit.setPlaceholderText("Optional curation rationale or cautions")
        self.notes_edit.setMaximumHeight(90)

        self.save_button = QtWidgets.QPushButton("Validate and Add Compound")
        self.save_button.setObjectName("addCompoundButton")
        self.save_button.setDefault(True)
        self.save_button.setStyleSheet(
            "QPushButton#addCompoundButton {"
            "background-color: #0078d7; color: white; font-weight: 700; "
            "padding: 10px 18px; border-radius: 6px;}"
            "QPushButton#addCompoundButton:disabled {"
            "background-color: #a6c8f0; color: white;}"
        )
        self.clear_button = QtWidgets.QPushButton("Clear Form")
        self.clear_button.setObjectName("clearCompoundFormButton")
        self.status_box = QtWidgets.QPlainTextEdit()
        self.status_box.setObjectName("compoundAdditionStatus")
        self.status_box.setReadOnly(True)
        self.status_box.setMaximumHeight(145)
        self.status_box.setPlaceholderText(
            "Validation errors or the new stable substance ID will appear here."
        )

        self._build_layout()
        self._connect_signals()
        self._set_role(self.primary_role_combo, "other_reagent")
        self._refresh_family_combo(
            self.primary_role_combo, self.primary_family_combo
        )
        self._refresh_family_combo(
            self.secondary_role_combo, self.secondary_family_combo
        )
        self.add_alias_row()

    @staticmethod
    def _numeric_edit(
        object_name: str,
        placeholder: str,
        *,
        positive: bool,
    ) -> QtWidgets.QLineEdit:
        edit = QtWidgets.QLineEdit()
        edit.setObjectName(object_name)
        edit.setPlaceholderText(placeholder)
        minimum = 0.0 if positive else -1_000_000.0
        edit.setValidator(QtGui.QDoubleValidator(minimum, 1_000_000.0, 8))
        return edit

    def _build_family_map(self) -> dict[str, tuple[tuple[str, str], ...]]:
        values: dict[str, list[tuple[str, str]]] = {}
        for family in self._taxonomy.get("families", ()):
            role_id = str(family["role_id"])
            family_id = str(family["id"])
            description = str(family.get("description") or "").strip()
            label = family_id if not description else f"{family_id} — {description}"
            values.setdefault(role_id, []).append((label, family_id))
        return {
            role_id: tuple(sorted(families, key=lambda item: item[1]))
            for role_id, families in values.items()
        }

    def _role_combo(self, object_name: str, *, optional: bool) -> QtWidgets.QComboBox:
        combo = QtWidgets.QComboBox()
        combo.setObjectName(object_name)
        if optional:
            combo.addItem("No secondary role", "")
        roles = sorted(
            self._taxonomy.get("roles", ()),
            key=lambda item: (int(item.get("priority", 100)), str(item["id"])),
        )
        for role in roles:
            role_id = str(role["id"])
            combo.addItem(role_id.replace("_", " ").title(), role_id)
        return combo

    @staticmethod
    def _set_role(combo: QtWidgets.QComboBox, role_id: str) -> None:
        index = combo.findData(role_id)
        if index >= 0:
            combo.setCurrentIndex(index)

    def _build_layout(self) -> None:
        root = QtWidgets.QVBoxLayout(self)
        root.setContentsMargins(20, 20, 20, 20)
        root.setSpacing(12)
        title = QtWidgets.QLabel("Add a Condition-Registry Compound")
        title.setStyleSheet("font-size: 20px; font-weight: 600;")
        root.addWidget(title)
        description = QtWidgets.QLabel(
            "Create one stable substance identity from curator-supplied evidence. "
            "The app validates identifiers, molecular structure, role taxonomy, "
            "duplicates, and aliases before changing the versioned definitions."
        )
        description.setWordWrap(True)
        root.addWidget(description)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        content = QtWidgets.QWidget()
        content_layout = QtWidgets.QVBoxLayout(content)
        content_layout.setContentsMargins(0, 0, 8, 0)
        content_layout.setSpacing(12)

        identity_group = QtWidgets.QGroupBox("Identity")
        identity_group.setObjectName("identityGroup")
        identity_form = QtWidgets.QFormLayout(identity_group)
        cas_row = QtWidgets.QHBoxLayout()
        cas_row.addWidget(self.cas_edit, stretch=1)
        cas_row.addWidget(self.cas_check_button)
        cas_row.addWidget(self.lookup_button)
        cas_row.addWidget(self.cas_status)
        identity_form.addRow("CAS number *", cas_row)
        identity_form.addRow("Canonical name *", self.name_edit)
        identity_form.addRow("Abbreviation", self.abbreviation_edit)
        identity_form.addRow("SMILES", self.smiles_edit)
        identity_form.addRow("Formula", self.formula_edit)
        identity_form.addRow("Molecular weight", self.molecular_weight_edit)
        identity_form.addRow("Substance kind", self.kind_combo)
        identity_form.addRow("Density", self.density_edit)
        identity_form.addRow("Boiling point", self.boiling_point_edit)
        identity_form.addRow("Melting point", self.melting_point_edit)
        roles_group = QtWidgets.QGroupBox("Condition-role capabilities")
        roles_group.setObjectName("conditionRoleCapabilitiesGroup")
        roles_form = QtWidgets.QFormLayout(roles_group)
        primary = QtWidgets.QHBoxLayout()
        primary.addWidget(self.primary_role_combo)
        primary.addWidget(self.primary_family_combo, stretch=1)
        primary.addWidget(self.primary_tag_edit, stretch=1)
        roles_form.addRow("Primary role", primary)
        secondary = QtWidgets.QHBoxLayout()
        secondary.addWidget(self.secondary_role_combo)
        secondary.addWidget(self.secondary_family_combo, stretch=1)
        secondary.addWidget(self.secondary_tag_edit, stretch=1)
        roles_form.addRow("Secondary role", secondary)
        role_note = QtWidgets.QLabel(
            "These are curated capabilities. Recipe construction still resolves "
            "the actual role from reaction context."
        )
        role_note.setWordWrap(True)
        roles_form.addRow("", role_note)

        top_columns = QtWidgets.QHBoxLayout()
        top_columns.setObjectName("identityAndRolesColumns")
        top_columns.setSpacing(12)
        top_columns.addWidget(identity_group, stretch=1)
        top_columns.addWidget(roles_group, stretch=1)
        content_layout.addLayout(top_columns)

        alias_group = QtWidgets.QGroupBox("Additional aliases and identifiers")
        alias_group.setObjectName("additionalAliasesGroup")
        alias_layout = QtWidgets.QVBoxLayout(alias_group)
        alias_layout.addWidget(self.alias_table)
        alias_buttons = QtWidgets.QHBoxLayout()
        add_alias_button = QtWidgets.QPushButton("Add Alias Row")
        add_alias_button.setObjectName("addAliasRowButton")
        add_alias_button.clicked.connect(self.add_alias_row)
        remove_alias_button = QtWidgets.QPushButton("Remove Selected Rows")
        remove_alias_button.setObjectName("removeAliasRowsButton")
        remove_alias_button.clicked.connect(self.remove_selected_alias_rows)
        alias_buttons.addWidget(add_alias_button)
        alias_buttons.addWidget(remove_alias_button)
        alias_buttons.addStretch()
        alias_layout.addLayout(alias_buttons)
        alias_help = QtWidgets.QLabel(
            "Use one row per alias. Check Shared only when the same normalized "
            "identifier intentionally belongs to more than one substance."
        )
        alias_help.setWordWrap(True)
        alias_layout.addWidget(alias_help)

        provenance_group = QtWidgets.QGroupBox("Provenance")
        provenance_group.setObjectName("provenanceGroup")
        provenance_form = QtWidgets.QFormLayout(provenance_group)
        provenance_form.addRow("Source *", self.source_edit)
        provenance_form.addRow("Curator notes", self.notes_edit)

        lower_columns = QtWidgets.QHBoxLayout()
        lower_columns.setObjectName("aliasesAndProvenanceColumns")
        lower_columns.setSpacing(12)
        lower_columns.addWidget(alias_group, stretch=1)
        lower_columns.addWidget(provenance_group, stretch=1)
        content_layout.addLayout(lower_columns)
        scroll.setWidget(content)
        root.addWidget(scroll, stretch=1)

        buttons = QtWidgets.QHBoxLayout()
        buttons.addWidget(self.save_button)
        buttons.addWidget(self.clear_button)
        buttons.addStretch()
        root.addLayout(buttons)
        root.addWidget(QtWidgets.QLabel("Curation status"))
        root.addWidget(self.status_box)

    def _connect_signals(self) -> None:
        self.cas_check_button.clicked.connect(self.check_cas)
        self.lookup_button.clicked.connect(self.start_web_lookup)
        self.cas_edit.editingFinished.connect(self.check_cas)
        self.primary_role_combo.currentIndexChanged.connect(
            lambda: self._refresh_family_combo(
                self.primary_role_combo, self.primary_family_combo
            )
        )
        self.secondary_role_combo.currentIndexChanged.connect(
            lambda: self._refresh_family_combo(
                self.secondary_role_combo, self.secondary_family_combo
            )
        )
        self.save_button.clicked.connect(self.save_compound)
        self.clear_button.clicked.connect(self.clear_form)

    def _refresh_family_combo(
        self,
        role_combo: QtWidgets.QComboBox,
        family_combo: QtWidgets.QComboBox,
    ) -> None:
        current = family_combo.currentData()
        role_id = str(role_combo.currentData() or "")
        family_combo.clear()
        family_combo.addItem("No family annotation", "")
        for label, family_id in self._families_by_role.get(role_id, ()):
            family_combo.addItem(label, family_id)
        index = family_combo.findData(current)
        if index >= 0:
            family_combo.setCurrentIndex(index)
        family_combo.setEnabled(bool(role_id))

    @QtCore.pyqtSlot()
    def add_alias_row(
        self,
        identifier_type: str = "common_name",
        value: str = "",
        language: str = "en",
        is_preferred: bool = False,
        allow_ambiguous: bool = False,
    ) -> None:
        row = self.alias_table.rowCount()
        self.alias_table.insertRow(row)
        type_combo = QtWidgets.QComboBox()
        for alias_type in ALIAS_TYPES:
            type_combo.addItem(alias_type.replace("_", " ").title(), alias_type)
        index = type_combo.findData(identifier_type)
        if index >= 0:
            type_combo.setCurrentIndex(index)
        value_edit = QtWidgets.QLineEdit(value)
        value_edit.setPlaceholderText("Alias or external identifier")
        language_edit = QtWidgets.QLineEdit(language)
        language_edit.setMaximumWidth(80)
        preferred_check = QtWidgets.QCheckBox()
        preferred_check.setChecked(is_preferred)
        preferred_check.setToolTip("Mark this as the preferred identifier of its type")
        shared_check = QtWidgets.QCheckBox()
        shared_check.setChecked(allow_ambiguous)
        shared_check.setToolTip("Allow this identifier to resolve ambiguously")
        self.alias_table.setCellWidget(row, 0, type_combo)
        self.alias_table.setCellWidget(row, 1, value_edit)
        self.alias_table.setCellWidget(row, 2, language_edit)
        self.alias_table.setCellWidget(row, 3, preferred_check)
        self.alias_table.setCellWidget(row, 4, shared_check)

    @QtCore.pyqtSlot()
    def remove_selected_alias_rows(self) -> None:
        rows = sorted(
            {index.row() for index in self.alias_table.selectedIndexes()},
            reverse=True,
        )
        for row in rows:
            self.alias_table.removeRow(row)

    def alias_inputs(self) -> tuple[CompoundAliasInput, ...]:
        aliases = []
        for row in range(self.alias_table.rowCount()):
            type_combo = self.alias_table.cellWidget(row, 0)
            value_edit = self.alias_table.cellWidget(row, 1)
            language_edit = self.alias_table.cellWidget(row, 2)
            preferred_check = self.alias_table.cellWidget(row, 3)
            shared_check = self.alias_table.cellWidget(row, 4)
            if not isinstance(type_combo, QtWidgets.QComboBox):
                continue
            if not isinstance(value_edit, QtWidgets.QLineEdit):
                continue
            value = value_edit.text().strip()
            if not value:
                continue
            language = (
                language_edit.text().strip()
                if isinstance(language_edit, QtWidgets.QLineEdit)
                else ""
            )
            aliases.append(
                CompoundAliasInput(
                    identifier_type=str(type_combo.currentData()),
                    value=value,
                    language=language or None,
                    is_preferred=(
                        preferred_check.isChecked()
                        if isinstance(preferred_check, QtWidgets.QCheckBox)
                        else False
                    ),
                    allow_ambiguous=(
                        shared_check.isChecked()
                        if isinstance(shared_check, QtWidgets.QCheckBox)
                        else False
                    ),
                )
            )
        return tuple(aliases)

    @staticmethod
    def _optional_number(edit: QtWidgets.QLineEdit, label: str) -> Optional[float]:
        text = edit.text().strip()
        if not text:
            return None
        try:
            return float(text)
        except ValueError as error:
            raise CompoundAdditionError((f"INVALID_NUMBER:{label}",)) from error

    @staticmethod
    def _role_assignment(
        role_combo: QtWidgets.QComboBox,
        family_combo: QtWidgets.QComboBox,
        tag_edit: QtWidgets.QLineEdit,
    ) -> Optional[RoleAssignment]:
        role_id = str(role_combo.currentData() or "")
        if not role_id:
            return None
        return RoleAssignment(
            role_id=role_id,
            family_id=str(family_combo.currentData() or "") or None,
            tag=tag_edit.text().strip() or None,
            evidence="curator_gui",
        )

    def compound_request(self) -> CompoundAdditionRequest:
        """Return typed form data without changing registry definitions."""
        roles = tuple(
            assignment
            for assignment in (
                self._role_assignment(
                    self.primary_role_combo,
                    self.primary_family_combo,
                    self.primary_tag_edit,
                ),
                self._role_assignment(
                    self.secondary_role_combo,
                    self.secondary_family_combo,
                    self.secondary_tag_edit,
                ),
            )
            if assignment is not None
        )
        return CompoundAdditionRequest(
            canonical_name=self.name_edit.text().strip(),
            cas=self.cas_edit.text().strip(),
            source=self.source_edit.text().strip(),
            smiles=self.smiles_edit.text().strip() or None,
            abbreviation=self.abbreviation_edit.text().strip() or None,
            formula=self.formula_edit.text().strip() or None,
            molecular_weight=self._optional_number(
                self.molecular_weight_edit, "molecular_weight"
            ),
            substance_kind=str(self.kind_combo.currentText()).strip() or None,
            density=self._optional_number(self.density_edit, "density"),
            boiling_point_c=self._optional_number(
                self.boiling_point_edit, "boiling_point_c"
            ),
            melting_point_c=self._optional_number(
                self.melting_point_edit, "melting_point_c"
            ),
            roles=roles,
            aliases=self.alias_inputs(),
            curator_notes=self.notes_edit.toPlainText().strip() or None,
        )

    @QtCore.pyqtSlot()
    def check_cas(self) -> None:
        cas = self.cas_edit.text().strip()
        if not cas:
            self.cas_status.setText("CAS required")
            self.cas_status.setStyleSheet("color: #a4262c;")
            return
        result = resolve_substance(cas=cas)
        if result.status == "unresolved":
            self.editing_substance_id = None
            self.cas_edit.setReadOnly(False)
            self.save_button.setText("Validate and Add Compound")
            self.cas_status.setText("Available")
            self.cas_status.setStyleSheet("color: #107c10; font-weight: 600;")
        elif result.status == "resolved" and result.substance is not None:
            self.load_existing_substance(result.substance)
        else:
            self.cas_status.setText("Invalid or ambiguous")
            self.cas_status.setStyleSheet("color: #a4262c; font-weight: 600;")

    def _load_role_assignment(
        self,
        role_combo: QtWidgets.QComboBox,
        family_combo: QtWidgets.QComboBox,
        tag_edit: QtWidgets.QLineEdit,
        assignment: Optional[RoleAssignment],
    ) -> None:
        if assignment is None:
            role_combo.setCurrentIndex(-1)
            self._refresh_family_combo(role_combo, family_combo)
            tag_edit.clear()
            return
        self._set_role(role_combo, assignment.role_id)
        self._refresh_family_combo(role_combo, family_combo)
        family_index = family_combo.findData(assignment.family_id or "")
        if family_index >= 0:
            family_combo.setCurrentIndex(family_index)
        tag_edit.setText(assignment.tag or "")

    def load_existing_substance(self, substance: Substance) -> None:
        """Load one resolved registry substance into the editable form."""
        self.editing_substance_id = substance.substance_id
        self.last_lookup_result = None
        self.cas_edit.setText(substance.cas or "")
        self.cas_edit.setReadOnly(True)
        self.name_edit.setText(substance.canonical_name)
        abbreviations = sorted(
            (
                identifier
                for identifier in substance.identifiers
                if identifier.status == "active"
                and identifier.identifier_type == "abbreviation"
            ),
            key=lambda identifier: (not identifier.is_preferred, identifier.value),
        )
        self.abbreviation_edit.setText(
            abbreviations[0].value if abbreviations else ""
        )
        self.smiles_edit.setText(substance.smiles or "")
        properties = substance.properties
        self.formula_edit.setText(str(properties.get("formula") or ""))
        self.molecular_weight_edit.setText(str(properties.get("mw") or ""))
        self.kind_combo.setCurrentText(str(properties.get("type") or ""))
        self.density_edit.setText(str(properties.get("density") or ""))
        self.boiling_point_edit.setText(str(properties.get("bp") or ""))
        self.melting_point_edit.setText(str(properties.get("mp") or ""))

        roles = (*substance.roles, None, None)[:2]
        self._load_role_assignment(
            self.primary_role_combo,
            self.primary_family_combo,
            self.primary_tag_edit,
            roles[0],
        )
        self._load_role_assignment(
            self.secondary_role_combo,
            self.secondary_family_combo,
            self.secondary_tag_edit,
            roles[1],
        )

        self.alias_table.setRowCount(0)
        for identifier in substance.identifiers:
            if (
                identifier.status != "active"
                or identifier.identifier_type not in ALIAS_TYPES
            ):
                continue
            self.add_alias_row(
                identifier.identifier_type,
                identifier.value,
                identifier.language or "",
                identifier.is_preferred,
                identifier.allow_ambiguous,
            )
        if self.alias_table.rowCount() == 0:
            self.add_alias_row()

        self.source_edit.setText(str(properties.get("source") or ""))
        self.notes_edit.setPlainText(str(properties.get("curator_notes") or ""))
        self.save_button.setText("Validate and Save Changes")
        self.cas_status.setText("Loaded for editing")
        self.cas_status.setStyleSheet("color: #005a9e; font-weight: 600;")
        source_note = (
            "Review the loaded fields and save changes."
            if self.source_edit.text().strip()
            else "Enter a provenance source before saving this legacy record."
        )
        self.status_box.setPlainText(
            f"Loaded {substance.canonical_name} ({substance.substance_id}).\n"
            f"{source_note}\n"
            "The CAS number is locked to preserve the stable substance identity."
        )

    @QtCore.pyqtSlot()
    def start_web_lookup(self) -> None:
        """Start a background web lookup for the current CAS number."""
        cas = self.cas_edit.text().strip()
        if not is_valid_cas_number(cas):
            self.cas_status.setText("Invalid CAS")
            self.cas_status.setStyleSheet("color: #a4262c; font-weight: 600;")
            self.status_box.setPlainText(
                "Auto-fill was not started because the CAS checksum is invalid."
            )
            return
        if self.lookup_thread is not None and self.lookup_thread.isRunning():
            return
        self.lookup_button.setEnabled(False)
        self.save_button.setEnabled(False)
        self.cas_status.setText("Looking up…")
        self.cas_status.setStyleSheet("color: #005a9e; font-weight: 600;")
        self.status_box.setPlainText(
            f"Looking up {cas} using PubChem; NCI/CADD will be used as a fallback."
        )
        thread = QtCore.QThread(self)
        worker = CompoundLookupWorker(cas)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.finished.connect(self._web_lookup_finished)
        worker.finished.connect(thread.quit)
        worker.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(self._clear_lookup_thread)
        self.lookup_thread = thread
        self.lookup_worker = worker
        thread.start()

    @QtCore.pyqtSlot(object, str)
    def _web_lookup_finished(
        self,
        result: Optional[CompoundLookupResult],
        error: str,
    ) -> None:
        self.lookup_button.setEnabled(True)
        self.save_button.setEnabled(True)
        if error:
            self.cas_status.setText("Lookup failed")
            self.cas_status.setStyleSheet("color: #a4262c; font-weight: 600;")
            self.status_box.setPlainText(f"Web lookup failed: {error}")
            return
        if result is None:
            self.cas_status.setText("Lookup failed")
            self.status_box.setPlainText("Web lookup returned no result object.")
            return
        if self.cas_edit.text().strip() != result.cas:
            self.status_box.setPlainText(
                f"Ignored lookup result for {result.cas} because the CAS field changed."
            )
            self.cas_status.setText("CAS changed")
            return
        self.last_lookup_result = result
        if not result.found:
            self.cas_status.setText("No web record")
            self.cas_status.setStyleSheet("color: #a4262c; font-weight: 600;")
            details = [f"No compound metadata was found for {result.cas}."]
            details.extend(f"Warning: {warning}" for warning in result.warnings)
            self.status_box.setPlainText("\n".join(details))
            return
        self.apply_lookup_result(result)

    @QtCore.pyqtSlot()
    def _clear_lookup_thread(self) -> None:
        self.lookup_thread = None
        self.lookup_worker = None

    def _set_optional_text(
        self,
        edit: QtWidgets.QLineEdit,
        value: Optional[str | float],
    ) -> None:
        if value is None:
            return
        edit.setText(str(value))

    def _existing_alias_keys(self) -> set[tuple[str, str]]:
        return {
            (alias.identifier_type, alias.value.casefold())
            for alias in self.alias_inputs()
        }

    def apply_lookup_result(self, result: CompoundLookupResult) -> None:
        """Populate supported fields and retain lookup provenance for review."""
        self._set_optional_text(self.name_edit, result.canonical_name)
        self._set_optional_text(self.abbreviation_edit, result.abbreviation)
        self._set_optional_text(self.smiles_edit, result.smiles)
        self._set_optional_text(self.formula_edit, result.formula)
        self._set_optional_text(
            self.molecular_weight_edit, result.molecular_weight
        )
        self._set_optional_text(self.density_edit, result.density)
        self._set_optional_text(
            self.boiling_point_edit, result.boiling_point_c
        )
        self._set_optional_text(
            self.melting_point_edit, result.melting_point_c
        )
        if result.substance_kind:
            self.kind_combo.setCurrentText(result.substance_kind)
        if result.source_ids:
            self.source_edit.setText(";".join(result.source_ids))

        existing = self._existing_alias_keys()
        canonical_key = (result.canonical_name or "").casefold()
        abbreviation_key = (result.abbreviation or "").casefold()

        def add_if_new(identifier_type: str, value: Optional[str]) -> None:
            if not value:
                return
            key = (identifier_type, value.casefold())
            if key in existing:
                return
            if value.casefold() in {canonical_key, abbreviation_key}:
                return
            self.add_alias_row(identifier_type, value, "en")
            existing.add(key)

        add_if_new("systematic_name", result.iupac_name)
        for synonym in result.synonyms:
            add_if_new("common_name", synonym)
        add_if_new("inchi_key", result.inchi_key)
        if result.pubchem_cid is not None:
            add_if_new("database_id", f"pubchem:cid:{result.pubchem_cid}")
        empty_rows = []
        for row in range(self.alias_table.rowCount()):
            value_edit = self.alias_table.cellWidget(row, 1)
            if isinstance(value_edit, QtWidgets.QLineEdit) and not value_edit.text().strip():
                empty_rows.append(row)
        for row in reversed(empty_rows):
            self.alias_table.removeRow(row)

        notes = self.notes_edit.toPlainText().strip()
        lookup_notes = ["Web auto-fill sources:", *result.source_urls]
        lookup_notes.extend(f"Lookup warning: {item}" for item in result.warnings)
        lookup_notes.append(
            "Condition roles and families were not inferred; review them manually."
        )
        self.notes_edit.setPlainText(
            "\n".join(part for part in (notes, "\n".join(lookup_notes)) if part)
        )
        self.cas_status.setText("Web data loaded")
        self.cas_status.setStyleSheet("color: #107c10; font-weight: 600;")
        filled = [
            label
            for label, value in (
                ("name", result.canonical_name),
                ("SMILES", result.smiles),
                ("formula", result.formula),
                ("molecular weight", result.molecular_weight),
                ("density", result.density),
                ("boiling point", result.boiling_point_c),
                ("melting point", result.melting_point_c),
                ("physical state", result.substance_kind),
            )
            if value is not None
        ]
        status_lines = [
            f"Auto-filled {result.cas} ({result.status}).",
            f"Fields: {', '.join(filled) or 'aliases only'}.",
            f"Aliases/identifiers offered: {len(self.alias_inputs())}.",
            "Roles and families require curator review.",
        ]
        status_lines.extend(f"Warning: {item}" for item in result.warnings)
        self.status_box.setPlainText("\n".join(status_lines))

    @QtCore.pyqtSlot()
    def save_compound(self) -> None:
        self.save_button.setEnabled(False)
        self.status_box.clear()
        editing_substance_id = self.editing_substance_id
        try:
            request = self.compound_request()
            result = (
                update_compound(editing_substance_id, request)
                if editing_substance_id is not None
                else add_compound(request)
            )
        except CompoundAdditionError as error:
            action = "saved" if editing_substance_id is not None else "added"
            self.status_box.appendPlainText(f"Compound was not {action}:")
            for item in error.errors:
                self.status_box.appendPlainText(f"  • {item}")
            self.cas_status.setText("Validation failed")
            self.cas_status.setStyleSheet("color: #a4262c; font-weight: 600;")
        except Exception as error:
            action = "saved" if editing_substance_id is not None else "added"
            self.status_box.appendPlainText(
                f"Compound was not {action}: {type(error).__name__}: {error}"
            )
            self.cas_status.setText("Write failed")
            self.cas_status.setStyleSheet("color: #a4262c; font-weight: 600;")
        else:
            was_updated = editing_substance_id is not None
            action = "Updated" if was_updated else "Added"
            self.editing_substance_id = (
                result.substance.substance_id if was_updated else None
            )
            self.cas_status.setText(action)
            self.cas_status.setStyleSheet("color: #107c10; font-weight: 600;")
            self.status_box.appendPlainText(
                f"{action} {result.substance.canonical_name} as "
                f"{result.substance.substance_id}."
            )
            self.status_box.appendPlainText(
                f"Canonical SMILES: {result.canonical_smiles or '-'}"
            )
            self.status_box.appendPlainText(f"Formula: {result.formula or '-'}")
            molecular_weight = (
                "-"
                if result.molecular_weight is None
                else f"{result.molecular_weight:.6g}"
            )
            self.status_box.appendPlainText(
                f"Molecular weight: {molecular_weight}"
            )
            self.status_box.appendPlainText(
                f"Supplemental identifiers: {result.alias_count}"
            )
        finally:
            self.save_button.setEnabled(True)

    @QtCore.pyqtSlot()
    def clear_form(self) -> None:
        for edit in (
            self.cas_edit,
            self.name_edit,
            self.abbreviation_edit,
            self.smiles_edit,
            self.formula_edit,
            self.molecular_weight_edit,
            self.density_edit,
            self.boiling_point_edit,
            self.melting_point_edit,
            self.primary_tag_edit,
            self.secondary_tag_edit,
            self.source_edit,
        ):
            edit.clear()
        self.kind_combo.setCurrentIndex(0)
        self._set_role(self.primary_role_combo, "other_reagent")
        self.secondary_role_combo.setCurrentIndex(0)
        self.alias_table.setRowCount(0)
        self.add_alias_row()
        self.notes_edit.clear()
        self.cas_status.setText("Not checked")
        self.cas_status.setStyleSheet("")
        self.cas_edit.setReadOnly(False)
        self.save_button.setText("Validate and Add Compound")
        self.status_box.clear()
        self.last_lookup_result = None
        self.editing_substance_id = None

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        if self.lookup_thread is not None and self.lookup_thread.isRunning():
            self.status_box.setPlainText(
                "A web lookup is still running. Close the window after it finishes."
            )
            event.ignore()
            return
        event.accept()


def main() -> None:
    """Launch the condition-registry compound curator."""
    application = QtWidgets.QApplication(sys.argv)
    window = CompoundRegistryWindow()
    window.showMaximized()
    raise SystemExit(application.exec())


if __name__ == "__main__":
    main()


__all__ = [
    "ALIAS_TYPES",
    "CompoundLookupWorker",
    "CompoundRegistryWindow",
    "main",
]

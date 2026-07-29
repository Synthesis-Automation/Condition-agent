"""PyQt6 wrapper for the single-event reaction-template registry."""

from __future__ import annotations

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
from reactive_taxonomy.reaction_api import featurize_reaction  # noqa: E402


def _parse_key_values(text: str, *, field: str) -> dict[str, str]:
    """Parse comma-separated KEY=VALUE authoring options."""
    result = {}
    for raw_value in (value.strip() for value in text.split(",")):
        if not raw_value:
            continue
        if "=" not in raw_value:
            raise ReactionTemplateError(
                f"{field} entries must use KEY=VALUE"
            )
        key, value = (
            part.strip() for part in raw_value.split("=", 1)
        )
        if not key or not value:
            raise ReactionTemplateError(
                f"{field} entries must use non-empty KEY=VALUE"
            )
        if key in result:
            raise ReactionTemplateError(f"{field} repeats key {key}")
        result[key] = value
    return result


def _format_template(template: ReactionTemplate) -> str:
    """Render one stored template as a compact authoring summary."""
    lines = [
        "REACTION TEMPLATE",
        "",
        f"ID: {template.template_id}",
        f"Name: {template.display_name}",
        f"Status: {template.status}",
        f"Family: {template.family_id or 'unassigned'}",
        f"Reaction label: {template.reaction_label}",
        f"Product label: {template.product_label}",
        f"Transformation: {template.transformation_class or 'unassigned'}",
        f"Edit archetype: {template.edit_archetype}",
        f"Edit fingerprint: {template.edit_fingerprint}",
        f"Definition hash: {template.definition_hash}",
    ]
    if template.aliases:
        lines.append(f"Aliases: {', '.join(template.aliases)}")
    lines.extend(("", "SEMANTIC ROLES"))
    for role in template.roles:
        maps = ", ".join(str(value) for value in role.atom_map_numbers)
        label = (
            f"; display {role.display_label}"
            if role.display_label is not None
            else ""
        )
        tokens = (
            "; requires " + ", ".join(role.required_context_tokens)
            if role.required_context_tokens
            else ""
        )
        lines.append(
            f"  {role.role_id}: {role.site_type}; reference maps {maps}"
            f"{label}{tokens}"
        )
    if template.atom_element_alternatives:
        lines.extend(("", "ELEMENT ALTERNATIVES"))
        for item in template.atom_element_alternatives:
            lines.append(
                f"  map {item.atom_map_number}: "
                + ", ".join(item.elements)
            )
    lines.extend(("", "PARTICIPANTS"))
    for participant in template.participants:
        lines.append(
            f"  {participant.side} × {participant.explicit_count}: "
            f"{participant.canonical_smiles}"
        )
    lines.extend(("", "NORMALIZED EDITS"))
    for edit in template.edits:
        atom_1 = (
            f"{edit.atom_1.element}(map {edit.atom_1.atom_map_number})"
        )
        atom_2 = (
            f"{edit.atom_2.element}(map {edit.atom_2.atom_map_number})"
            if edit.atom_2 is not None
            else "H"
        )
        lines.append(
            f"  {edit.edit_type}: {atom_1}–{atom_2}, "
            f"{edit.old_order or 'none'} → {edit.new_order or 'none'}"
        )
    lines.extend(
        (
            "",
            "MAPPED REFERENCE",
            f"  {template.mapped_reference_reaction}",
        )
    )
    if template.notes:
        lines.extend(("", "NOTES", f"  {template.notes}"))
    return "\n".join(lines)


def _format_reaction_analysis(analysis: object) -> str:
    """Render a full reaction analysis without exposing raw nested JSON."""
    completeness = getattr(analysis, "reaction_completeness", None)
    signature = getattr(analysis, "reaction_signature", None)
    display_label = getattr(analysis, "display_label", None)
    lines = [
        "REACTION FEATURIZATION",
        "",
        f"Input: {getattr(analysis, 'input_reaction_smiles', '')}",
        f"Status: {'valid' if getattr(analysis, 'valid', False) else 'invalid'}",
        f"Evidence: {getattr(analysis, 'evidence_quality', 'unresolved')}",
        f"Reaction: {getattr(analysis, 'reaction_label', None) or 'unavailable'}",
        (
            "Structural label: "
            f"{getattr(display_label, 'structural_label', None) or 'unavailable'}"
        ),
        (
            "Transformation: "
            f"{getattr(analysis, 'transformation_class', None) or 'unresolved'}"
        ),
        f"Named family: {getattr(analysis, 'named_family', None) or 'unassigned'}",
    ]
    compatible = tuple(
        getattr(analysis, "compatible_named_families", ()) or ()
    )
    if compatible:
        lines.append(f"Compatible families: {', '.join(compatible)}")
    if completeness is not None:
        coverage = getattr(completeness, "product_heavy_atom_coverage", None)
        lines.extend(
            (
                "",
                "COMPLETENESS",
                f"  Status: {completeness.status}",
                f"  Evidence: {completeness.evidence}",
                (
                    "  Product heavy-atom coverage: "
                    + (
                        f"{100.0 * float(coverage):.1f}%"
                        if coverage is not None
                        else "unavailable"
                    )
                ),
            )
        )
        product_excess = dict(
            getattr(completeness, "product_element_excess", {}) or {}
        )
        if product_excess:
            lines.append(
                "  Unaccounted product atoms: "
                + ", ".join(
                    f"{element} × {count}"
                    for element, count in sorted(product_excess.items())
                )
            )
    lines.extend(("", "SIGNATURE"))
    if signature is None:
        lines.append("  RS3 signature: unavailable")
    else:
        lines.extend(
            (
                f"  ID: {signature.signature_id}",
                f"  Event scope: {signature.event_scope}",
                f"  Event count: {signature.event_count}",
                f"  Edit archetype: {signature.edit_archetype}",
                f"  Bond-edit key: {signature.bond_edit_signature_key}",
            )
        )
        if signature.edits:
            lines.append("  Edits:")
            for edit in signature.edits:
                atom_2 = edit.atom_2.element if edit.atom_2 is not None else "H"
                lines.append(
                    f"    {edit.edit_type}: {edit.atom_1.element}–{atom_2}, "
                    f"{edit.old_order or 'none'} → "
                    f"{edit.new_order or 'none'}"
                )
    candidates = tuple(getattr(analysis, "candidates", ()) or ())
    lines.extend(("", f"Candidate interpretations: {len(candidates)}"))
    warnings = tuple(getattr(analysis, "warnings", ()) or ())
    if warnings:
        lines.extend(("", "WARNINGS"))
        lines.extend(f"  • {warning}" for warning in warnings)
    error = getattr(analysis, "error", None)
    if error:
        lines.extend(("", "ERROR", f"  {error}"))
    return "\n".join(lines)


def _format_match_result(result: object) -> str:
    """Render template matches and provisional evidence clearly."""
    lines = [
        "REACTION TEMPLATE MATCH",
        "",
        f"Input: {getattr(result, 'reaction_smiles', '')}",
        f"Status: {'valid' if getattr(result, 'valid', False) else 'invalid'}",
        f"Evidence: {getattr(result, 'evidence', 'unresolved')}",
        f"RS3 signature: {getattr(result, 'signature_id', None) or 'unavailable'}",
        (
            "Generic edit key: "
            f"{getattr(result, 'edit_fingerprint', None) or 'unavailable'}"
        ),
    ]
    lines.extend(("", _format_match_section(result)))
    warnings = tuple(getattr(result, "warnings", ()) or ())
    if warnings:
        lines.extend(("", "WARNINGS"))
        lines.extend(f"  • {warning}" for warning in warnings)
    error = getattr(result, "error", None)
    if error:
        lines.extend(("", "ERROR", f"  {error}"))
    return "\n".join(lines)


def _format_match_section(result: object) -> str:
    """Render only the registry interpretation section for combined analysis."""
    matches = tuple(getattr(result, "matches", ()) or ())
    lines = [f"TEMPLATE REGISTRY MATCHES ({len(matches)})"]
    if not matches:
        lines.append("  None")
    for match in matches:
        if (
            match.evidence
            == "exact_template_reconstruction_with_inferred_multiplicity"
        ):
            mode = "multiplicity-assisted exact reconstruction"
        elif match.evidence == "exact_template_reconstruction":
            mode = "exact template reconstruction"
        elif match.provisional:
            mode = "provisional centre-transition match"
        else:
            mode = "verified edit match"
        lines.extend(
            (
                f"  {match.display_name}",
                f"    ID: {match.template_id}",
                f"    Family: {match.family_id or 'unassigned'}",
                f"    Status: {match.status}",
                f"    Match: {mode}",
                f"    Evidence: {match.evidence}",
                f"    Confidence: {match.confidence:.2f}",
            )
        )
        if match.predicted_product_smiles:
            lines.append(
                f"    Reconstructed product: "
                f"{match.predicted_product_smiles}"
            )
        if match.inferred_multiplicity:
            lines.append("    Inferred repeated participant: yes")
        interpretation = match.interpretation
        if interpretation is not None:
            lines.extend(
                (
                    f"    Reaction label: {interpretation.reaction_label}",
                    f"    Structural label: "
                    f"{interpretation.structural_label}",
                    "    Bound roles:",
                )
            )
            for binding in interpretation.roles:
                count = (
                    f" × {binding.multiplicity}"
                    if binding.multiplicity > 1
                    else ""
                )
                steric = (
                    f"{binding.steric_class} "
                    f"({binding.steric_score:.2f})"
                    if binding.steric_class is not None
                    and binding.steric_score is not None
                    else "unavailable"
                )
                lines.extend(
                    (
                        f"      {binding.role_id}{count}: "
                        f"{binding.chemist_label}",
                        f"        Site: {binding.site_type}; component "
                        f"{binding.component_index}; atoms "
                        + ", ".join(
                            str(value) for value in binding.atom_indices
                        ),
                        f"        Steric: {steric}; electronic: "
                        f"{binding.electronic_class or 'unavailable'}",
                    )
                )
                if binding.nearby_groups:
                    lines.append(
                        "        Nearby groups: "
                        + ", ".join(binding.nearby_groups)
                    )
    return "\n".join(lines)


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
        self.reaction_label_edit = QtWidgets.QLineEdit()
        self.reaction_label_edit.setObjectName("reactionLabel")
        self.reaction_label_edit.setPlaceholderText(
            "Acetal formation (defaults to name)"
        )
        self.product_label_edit = QtWidgets.QLineEdit()
        self.product_label_edit.setObjectName("productLabel")
        self.product_label_edit.setPlaceholderText(
            "acetal (defaults to product)"
        )
        self.role_labels_edit = QtWidgets.QLineEdit()
        self.role_labels_edit.setObjectName("roleLabels")
        self.role_labels_edit.setPlaceholderText(
            "activated_sp3_carbon=α-halo ester (optional)"
        )
        self.role_tokens_edit = QtWidgets.QLineEdit()
        self.role_tokens_edit.setObjectName("roleTokens")
        self.role_tokens_edit.setPlaceholderText(
            "activated_sp3_carbon=alpha_to:ester (optional)"
        )
        self.atom_alternatives_edit = QtWidgets.QLineEdit()
        self.atom_alternatives_edit.setObjectName("atomAlternatives")
        self.atom_alternatives_edit.setPlaceholderText(
            "1=Cl|Br|I (optional)"
        )
        self.transformation_edit = QtWidgets.QLineEdit()
        self.transformation_edit.setObjectName("transformationClass")
        self.transformation_edit.setPlaceholderText(
            "generic transformation class (optional)"
        )
        self.status_combo = QtWidgets.QComboBox()
        self.status_combo.setObjectName("templateStatus")
        self.status_combo.addItems(("draft", "active", "retired"))
        self.mapped_reaction_edit = QtWidgets.QLineEdit()
        self.mapped_reaction_edit.setObjectName("mappedReferenceReaction")
        self.mapped_reaction_edit.setPlaceholderText(
            "Paste one fully atom-mapped, single-event reaction SMILES"
        )
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

        self.query_edit = QtWidgets.QLineEdit()
        self.query_edit.setObjectName("queryReaction")
        self.query_edit.setPlaceholderText(
            "Paste a query reaction SMILES to derive its signature and match "
            "registered edit templates"
        )
        self.include_drafts_check = QtWidgets.QCheckBox(
            "Include draft templates"
        )
        self.include_drafts_check.setChecked(True)
        self.featurize_button = QtWidgets.QPushButton("Featurize reaction")
        self.featurize_button.setObjectName("featurizeQuery")
        self.match_button = QtWidgets.QPushButton("Match query")
        self.match_button.setObjectName("matchQuery")

        self.details = QtWidgets.QPlainTextEdit()
        self.details.setObjectName("templateDetails")
        self.details.setReadOnly(True)
        self.details.setMinimumHeight(180)
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
        label_row = QtWidgets.QHBoxLayout()
        label_row.addWidget(self.reaction_label_edit)
        label_row.addWidget(QtWidgets.QLabel("Product"))
        label_row.addWidget(self.product_label_edit)
        form.addRow("Reaction label", label_row)
        generalization_row = QtWidgets.QHBoxLayout()
        generalization_row.addWidget(self.role_labels_edit)
        generalization_row.addWidget(QtWidgets.QLabel("Atom alternatives"))
        generalization_row.addWidget(self.atom_alternatives_edit)
        form.addRow("Role labels", generalization_row)
        form.addRow("Required role tokens", self.role_tokens_edit)
        form.addRow("Mapped reference", self.mapped_reaction_edit)
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
        query_row.addWidget(self.include_drafts_check)
        query_row.addWidget(self.featurize_button)
        query_row.addWidget(self.match_button)
        browser_layout.addLayout(query_row)
        browser_layout.addWidget(QtWidgets.QLabel("Details / query result"))
        browser_layout.addWidget(self.details, stretch=1)
        splitter.addWidget(browser)
        splitter.setSizes((220, 590))

        self.status_label.setStyleSheet(
            "background: #eef4fa; border: 1px solid #ccd9e5; "
            "color: #23313f; padding: 7px; border-radius: 4px;"
        )
        layout.addWidget(self.status_label)

    def _connect_signals(self) -> None:
        self.import_button.clicked.connect(self.import_template)
        self.validate_button.clicked.connect(self.validate_registry)
        self.refresh_button.clicked.connect(self.refresh_registry)
        self.featurize_button.clicked.connect(self.featurize_query)
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
            role_labels = _parse_key_values(
                self.role_labels_edit.text(),
                field="Role labels",
            )
            raw_alternatives = _parse_key_values(
                self.atom_alternatives_edit.text(),
                field="Atom alternatives",
            )
            raw_role_tokens = _parse_key_values(
                self.role_tokens_edit.text(),
                field="Required role tokens",
            )
            try:
                atom_element_alternatives = {
                    int(map_number): tuple(
                        element.strip()
                        for element in elements.split("|")
                        if element.strip()
                    )
                    for map_number, elements in raw_alternatives.items()
                }
            except ValueError as exc:
                raise ReactionTemplateError(
                    "Atom-alternative keys must be integer map numbers"
                ) from exc
            template = derive_reaction_template(
                self.mapped_reaction_edit.text().strip(),
                template_id=self.template_id_edit.text().strip(),
                display_name=self.name_edit.text().strip(),
                family_id=self.family_edit.text().strip() or None,
                aliases=aliases,
                reaction_label=(
                    self.reaction_label_edit.text().strip() or None
                ),
                product_label=(
                    self.product_label_edit.text().strip() or None
                ),
                role_labels=role_labels,
                role_required_tokens={
                    role: tuple(
                        token.strip()
                        for token in tokens.split("|")
                        if token.strip()
                    )
                    for role, tokens in raw_role_tokens.items()
                },
                atom_element_alternatives=atom_element_alternatives,
                transformation_class=(
                    self.transformation_edit.text().strip() or None
                ),
                status=self.status_combo.currentText(),  # type: ignore[arg-type]
                provenance="qt6_manual_mapped_reference",
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
            self.details.setPlainText(_format_template(template))

    @QtCore.pyqtSlot()
    def featurize_query(self) -> None:
        """Run full standalone reaction featurization for the test input."""
        reaction = self.query_edit.text().strip()
        if not reaction:
            self._show_error(
                "Reaction required",
                ReactionTemplateError("Paste a reaction SMILES to featurize"),
            )
            return
        analysis = featurize_reaction(reaction)
        try:
            template_result = match_reaction_templates(
                reaction,
                path=self._registry_path(),
                include_drafts=self.include_drafts_check.isChecked(),
            )
        except (OSError, ReactionTemplateError) as exc:
            template_result = None
            template_summary = (
                "TEMPLATE REGISTRY\n"
                f"  Match unavailable: {exc}"
            )
        else:
            template_summary = _format_match_section(template_result)
        self.details.setPlainText(
            _format_reaction_analysis(analysis)
            + "\n\n"
            + template_summary
        )
        completeness = analysis.reaction_completeness
        completeness_status = (
            completeness.status if completeness is not None else "unavailable"
        )
        signature_id = (
            analysis.reaction_signature.signature_id
            if analysis.reaction_signature is not None
            else "unavailable"
        )
        match_status = "template matching unavailable"
        if template_result is not None:
            if template_result.matches:
                first_match = template_result.matches[0]
                if (
                    first_match.evidence
                    == "exact_template_reconstruction_with_inferred_multiplicity"
                ):
                    match_kind = "multiplicity-assisted exact"
                elif first_match.evidence == "exact_template_reconstruction":
                    match_kind = "exact reconstruction"
                elif first_match.provisional:
                    match_kind = "provisional centre match"
                else:
                    match_kind = "verified-edit"
                match_status = (
                    f"{len(template_result.matches)} template match(es), "
                    f"{first_match.family_id or first_match.template_id} "
                    f"({match_kind})"
                )
            else:
                match_status = "no template match"
        self.status_label.setText(
            f"Featurization: {'valid' if analysis.valid else 'invalid'}; "
            f"evidence {analysis.evidence_quality}; completeness "
            f"{completeness_status}; signature {signature_id}; "
            f"{match_status}."
        )

    @QtCore.pyqtSlot()
    def match_query(self) -> None:
        reaction = self.query_edit.text().strip()
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
        self.details.setPlainText(_format_match_result(result))
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

"""Reusable Qt dialogs for the ChemTools assistant UI."""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Dict, List, Optional

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QPlainTextEdit,
    QVBoxLayout,
    QWidget,
)

from chem_assistant.constraint_parser import ConstraintSpec, build_constraint_spec
from chemtools.rule import RuleBuilder
from chem_assistant.chemtools_wrapper import (
    RuleBuilderAutoInput,
    run_rule_builder_autofill,
)


def _split_csv(value: str) -> List[str]:
    return [item.strip() for item in (value or "").split(",") if item.strip()]


class ConstraintDialog(QDialog):
    """Dialog for editing constraint settings."""

    def __init__(
        self,
        spec: ConstraintSpec,
        constraint_text: str = "",
        parent: Optional[QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Manage Constraints")
        self.setModal(True)

        layout = QVBoxLayout(self)

        form = QFormLayout()
        self.text_edit = QPlainTextEdit()
        self.text_edit.setPlaceholderText("Natural language preferences...")
        self.text_edit.setPlainText(constraint_text.strip())
        form.addRow("Constraint text:", self.text_edit)

        self.allow_edit = QLineEdit(", ".join(sorted(spec.allow_metals)))
        self.exclude_edit = QLineEdit(", ".join(sorted(spec.exclude_metals)))
        self.prefer_edit = QLineEdit(", ".join(sorted(spec.prefer_metals)))
        form.addRow("Allow metals (comma separated):", self.allow_edit)
        form.addRow("Exclude metals:", self.exclude_edit)
        form.addRow("Prefer metals:", self.prefer_edit)

        self.cross_checkbox = QCheckBox("Enable cross-family search")
        self.cross_checkbox.setChecked(spec.search_all_families)
        form.addRow(self.cross_checkbox)

        self.rules_edit = QLineEdit(", ".join(sorted(spec.constraint_rules.keys())))
        form.addRow("Optional rules (comma separated):", self.rules_edit)

        layout.addLayout(form)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def build_spec(self) -> ConstraintSpec:
        """Return the updated constraint spec."""
        rules = {key: True for key in _split_csv(self.rules_edit.text())}
        spec = build_constraint_spec(
            text=self.text_edit.toPlainText().strip() or None,
            allow_metals=_split_csv(self.allow_edit.text()),
            exclude_metals=_split_csv(self.exclude_edit.text()),
            prefer_metals=_split_csv(self.prefer_edit.text()),
            search_all_families=self.cross_checkbox.isChecked(),
            base_constraint_rules=rules or None,
        )
        return spec

    def get_constraint_text(self) -> str:
        return self.text_edit.toPlainText().strip()


class RuleBuilderDialog(QDialog):
    """Form-based editor for rule database JSON."""

    def __init__(
        self,
        parent: Optional[QWidget] = None,
        initial_data: Optional[Dict[str, object]] = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Rule Builder Editor")
        self.resize(900, 700)
        self._loaded_path: Optional[str] = None

        self.family_edit = QLineEdit()
        self.metadata_id_edit = QLineEdit()
        self.metadata_name_edit = QLineEdit()
        self.metadata_version_edit = QLineEdit()
        self.metadata_created_edit = QLineEdit()
        self.metadata_status_edit = QLineEdit("draft")
        self.metadata_tags_edit = QLineEdit()

        self.references_edit = QPlainTextEdit()
        self.references_edit.setPlaceholderText("One reaction SMILES per line.")
        self.notes_edit = QPlainTextEdit()

        self.scope_type_edit = QLineEdit()
        self.scope_compatible_edit = QLineEdit()
        self.scope_incompatible_edit = QLineEdit()

        self.applies_all_edit = QLineEdit()
        self.applies_any_edit = QLineEdit()

        self.default_rule_edit = QPlainTextEdit()
        self.base_rules_edit = QPlainTextEdit()
        self.modifiers_edit = QPlainTextEdit()

        layout = QVBoxLayout(self)

        layout.addWidget(self._build_metadata_group())
        layout.addWidget(self._build_reaction_group())
        layout.addWidget(self._build_scope_group())
        layout.addWidget(self._build_applicability_group())
        layout.addWidget(self._build_rule_group())
        layout.addWidget(self._build_modifier_group())

        button_row = QHBoxLayout()
        load_btn = QPushButton("Load JSON...")
        load_btn.clicked.connect(self.load_from_file)
        button_row.addWidget(load_btn)

        validate_btn = QPushButton("Validate")
        validate_btn.clicked.connect(self.validate_data)
        button_row.addWidget(validate_btn)

        save_btn = QPushButton("Save As...")
        save_btn.clicked.connect(self.save_to_file)
        button_row.addWidget(save_btn)

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        button_row.addWidget(close_btn)
        button_row.addStretch()

        layout.addLayout(button_row)

        if initial_data:
            self.apply_data(initial_data)

    def _build_metadata_group(self) -> QGroupBox:
        group = QGroupBox("Metadata")
        form = QFormLayout(group)
        form.addRow("Family:", self.family_edit)
        form.addRow("Metadata ID:", self.metadata_id_edit)
        form.addRow("Name:", self.metadata_name_edit)
        form.addRow("Version:", self.metadata_version_edit)
        form.addRow("Created date:", self.metadata_created_edit)
        form.addRow("Status:", self.metadata_status_edit)
        form.addRow("Tags (comma separated):", self.metadata_tags_edit)
        return group

    def _build_reaction_group(self) -> QGroupBox:
        group = QGroupBox("Reaction Notes")
        layout = QVBoxLayout(group)
        layout.addWidget(QLabel("Reference reactions:"))
        layout.addWidget(self.references_edit)
        layout.addWidget(QLabel("Notes / protocol summary:"))
        layout.addWidget(self.notes_edit)
        return group

    def _build_scope_group(self) -> QGroupBox:
        group = QGroupBox("Scope")
        form = QFormLayout(group)
        form.addRow("Scope type:", self.scope_type_edit)
        form.addRow("Compatible groups:", self.scope_compatible_edit)
        form.addRow("Incompatible groups:", self.scope_incompatible_edit)
        return group

    def _build_applicability_group(self) -> QGroupBox:
        group = QGroupBox("Applicability")
        form = QFormLayout(group)
        form.addRow("Requires ALL features:", self.applies_all_edit)
        form.addRow("ANY-of features:", self.applies_any_edit)
        return group

    def _build_rule_group(self) -> QGroupBox:
        group = QGroupBox("Rules")
        layout = QVBoxLayout(group)
        layout.addWidget(QLabel("Default rule (JSON dict):"))
        self.default_rule_edit.setPlaceholderText('{"id": "...", "description": "...", "conditions": {...}}')
        layout.addWidget(self.default_rule_edit)
        layout.addWidget(QLabel("Base rules (JSON list):"))
        self.base_rules_edit.setPlaceholderText("[{...}, {...}]")
        layout.addWidget(self.base_rules_edit)
        return group

    def _build_modifier_group(self) -> QGroupBox:
        group = QGroupBox("Modifiers")
        layout = QVBoxLayout(group)
        layout.addWidget(QLabel("Modifiers (JSON list):"))
        self.modifiers_edit.setPlaceholderText("[{...}]")
        layout.addWidget(self.modifiers_edit)
        return group

    def _parse_lines(self, edit: QPlainTextEdit) -> List[str]:
        return [line.strip() for line in edit.toPlainText().splitlines() if line.strip()]

    def _json_from_edit(self, edit: QPlainTextEdit, default):
        text = edit.toPlainText().strip()
        if not text:
            return default
        return json.loads(text)

    def _build_rule_builder(self) -> RuleBuilder:
        family = self.family_edit.text().strip()
        if not family:
            raise ValueError("Reaction family is required.")

        builder = RuleBuilder.new(family)
        builder.set_metadata(
            id=self.metadata_id_edit.text().strip() or "draft_rule",
            name=self.metadata_name_edit.text().strip() or "Draft Rule",
            version=self.metadata_version_edit.text().strip() or "v0",
            created_date=self.metadata_created_edit.text().strip() or "",
            status=self.metadata_status_edit.text().strip() or "draft",
            tags=_split_csv(self.metadata_tags_edit.text()),
        )
        references = self._parse_lines(self.references_edit)
        if not references:
            raise ValueError("At least one reference reaction is required.")
        builder.add_reference_reactions(references)
        if self.notes_edit.toPlainText().strip():
            builder.set_notes(self.notes_edit.toPlainText().strip())

        builder.set_scope(
            scope_type=self.scope_type_edit.text().strip() or None,
            compatible=_split_csv(self.scope_compatible_edit.text()),
            incompatible=_split_csv(self.scope_incompatible_edit.text()),
        )

        builder.set_applies_if(
            all_features=_split_csv(self.applies_all_edit.text()),
            any_features=_split_csv(self.applies_any_edit.text()),
        )

        default_rule = self._json_from_edit(self.default_rule_edit, {})
        if not default_rule:
            raise ValueError("Default rule JSON is required.")
        builder.set_default_rule(
            rule_id=default_rule.get("id"),
            description=default_rule.get("description"),
            conditions=default_rule.get("conditions"),
        )

        base_rules = self._json_from_edit(self.base_rules_edit, [])
        for rule in base_rules:
            builder.upsert_base_rule(
                rule.get("id", "rule"),
                name=rule.get("name", "Unnamed rule"),
                description=rule.get("description", ""),
                reactant_features=rule.get("reactant_features") or {},
                conditions=rule.get("conditions") or {},
                priority=rule.get("priority"),
            )

        modifiers = self._json_from_edit(self.modifiers_edit, [])
        for modifier in modifiers:
            builder.upsert_modifier(
                modifier.get("id", "modifier"),
                when=modifier.get("when") or [],
                suggestion=modifier.get("suggest") or modifier.get("suggestion", ""),
                rationale=modifier.get("rationale"),
            )

        return builder

    def validate_data(self) -> None:
        try:
            builder = self._build_rule_builder()
            issues = builder.validate(strict=False)
            if issues:
                msg = "\n".join(f"- [{issue.severity}] {issue.field}: {issue.message}" for issue in issues)
                QMessageBox.warning(self, "Validation issues", msg)
            else:
                QMessageBox.information(self, "Validation", "No issues detected.")
        except Exception as exc:
            QMessageBox.critical(self, "Validation error", str(exc))

    def save_to_file(self) -> None:
        try:
            builder = self._build_rule_builder()
            issues = builder.validate(strict=False)
            if issues:
                msg = "\n".join(f"- [{issue.severity}] {issue.field}: {issue.message}" for issue in issues)
                answer = QMessageBox.question(
                    self,
                    "Validation issues",
                    f"Issues detected:\n{msg}\n\nSave anyway?",
                )
                if answer != QMessageBox.StandardButton.Yes:
                    return
            path, _ = QFileDialog.getSaveFileName(
                self,
                "Save rule database",
                self._loaded_path or "rule_db.json",
                "JSON files (*.json)",
            )
            if not path:
                return
            builder.save(path)
            self._loaded_path = path
            QMessageBox.information(self, "Saved", f"Rule database saved to {path}.")
        except Exception as exc:
            QMessageBox.critical(self, "Save failed", str(exc))

    def load_from_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Open rule database",
            self._loaded_path or "",
            "JSON files (*.json)",
        )
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8") as handle:
                data = json.load(handle)
            self.apply_data(data)
            self._loaded_path = path
            QMessageBox.information(self, "Loaded", f"Loaded {path}")
        except Exception as exc:
            QMessageBox.critical(self, "Load failed", str(exc))

    def apply_data(self, data: Dict[str, object]) -> None:
        metadata = data.get("metadata", {})
        reaction = data.get("reaction", {})
        scope = reaction.get("scope", {})
        applies_if = data.get("applies_if", {})

        self.family_edit.setText(str(reaction.get("family", "")))
        self.metadata_id_edit.setText(str(metadata.get("id", "")))
        self.metadata_name_edit.setText(str(metadata.get("name", "")))
        self.metadata_version_edit.setText(str(metadata.get("version", "")))
        self.metadata_created_edit.setText(str(metadata.get("created_date", "")))
        self.metadata_status_edit.setText(str(metadata.get("status", "")))
        self.metadata_tags_edit.setText(", ".join(metadata.get("tags", [])))

        references = reaction.get("reference_reactions", [])
        self.references_edit.setPlainText("\n".join(references or []))
        self.notes_edit.setPlainText(str(reaction.get("notes", "")))

        self.scope_type_edit.setText(str(scope.get("scope_type", "")))
        self.scope_compatible_edit.setText(", ".join(scope.get("compatible_functional_groups", [])))
        self.scope_incompatible_edit.setText(", ".join(scope.get("incompatible_functional_groups", [])))

        self.applies_all_edit.setText(", ".join(applies_if.get("all", [])))
        self.applies_any_edit.setText(", ".join(applies_if.get("any", [])))

        self.default_rule_edit.setPlainText(json.dumps(data.get("default_rule", {}), indent=2))
        self.base_rules_edit.setPlainText(json.dumps(data.get("base_rules", []), indent=2))
        self.modifiers_edit.setPlainText(json.dumps(data.get("modifiers", []), indent=2))


class RuleBuilderAutofillDialog(QDialog):
    """Dialog for running the LLM-assisted rule builder autofill."""

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Rule Builder Autofill")
        self.resize(800, 700)
        self.latest_result: Optional[Dict[str, object]] = None
        self.accepted_data: Optional[Dict[str, object]] = None

        layout = QVBoxLayout(self)

        self.family_edit = QLineEdit()
        self.metadata_id_edit = QLineEdit()
        self.metadata_name_edit = QLineEdit()
        self.metadata_version_edit = QLineEdit()
        self.created_date_edit = QLineEdit()
        self.status_edit = QLineEdit("draft")
        self.tags_edit = QLineEdit()

        metadata_group = QGroupBox("Metadata")
        meta_form = QFormLayout(metadata_group)
        meta_form.addRow("Family:", self.family_edit)
        meta_form.addRow("Metadata ID:", self.metadata_id_edit)
        meta_form.addRow("Name:", self.metadata_name_edit)
        meta_form.addRow("Version:", self.metadata_version_edit)
        meta_form.addRow("Created date:", self.created_date_edit)
        meta_form.addRow("Status:", self.status_edit)
        meta_form.addRow("Tags:", self.tags_edit)
        layout.addWidget(metadata_group)

        ref_group = QGroupBox("Reference reactions")
        ref_layout = QVBoxLayout(ref_group)
        self.references_edit = QPlainTextEdit()
        self.references_edit.setPlaceholderText("Enter reaction SMILES, one per line.")
        ref_layout.addWidget(self.references_edit)
        load_refs_btn = QPushButton("Load references from JSON...")
        load_refs_btn.clicked.connect(self.load_references_from_file)
        ref_layout.addWidget(load_refs_btn)
        layout.addWidget(ref_group)

        proto_group = QGroupBox("Protocol / notes")
        proto_layout = QVBoxLayout(proto_group)
        self.protocol_edit = QPlainTextEdit()
        self.protocol_edit.setPlaceholderText("Paste HTE notes or protocol text here.")
        proto_layout.addWidget(self.protocol_edit)
        proto_btn = QPushButton("Load protocol text from file...")
        proto_btn.clicked.connect(self.load_protocol_from_file)
        proto_layout.addWidget(proto_btn)
        layout.addWidget(proto_group)

        options_group = QGroupBox("Options")
        options_form = QFormLayout(options_group)
        self.notes_edit = QPlainTextEdit()
        self.notes_edit.setPlaceholderText("(Optional) Override notes for reaction block.")
        options_form.addRow("Override notes:", self.notes_edit)

        self.focus_edit = QLineEdit()
        self.applies_hint_edit = QLineEdit()
        self.modifier_hint_edit = QLineEdit()
        self.max_rules_edit = QLineEdit("4")

        options_form.addRow("Focus hint:", self.focus_edit)
        options_form.addRow("Applies-if hints:", self.applies_hint_edit)
        options_form.addRow("Modifier hints:", self.modifier_hint_edit)
        options_form.addRow("Max base rules:", self.max_rules_edit)
        layout.addWidget(options_group)

        self.result_display = QPlainTextEdit()
        self.result_display.setReadOnly(True)
        layout.addWidget(QLabel("Result preview:"))
        layout.addWidget(self.result_display)

        button_row = QHBoxLayout()
        run_btn = QPushButton("Run Autofill")
        run_btn.clicked.connect(self.run_autofill)
        button_row.addWidget(run_btn)

        save_btn = QPushButton("Save Result...")
        save_btn.clicked.connect(self.save_result)
        button_row.addWidget(save_btn)

        accept_btn = QPushButton("Open in Builder")
        accept_btn.clicked.connect(self.accept_draft)
        button_row.addWidget(accept_btn)

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.reject)
        button_row.addWidget(close_btn)
        button_row.addStretch()

        layout.addLayout(button_row)

    def _gather_references(self) -> List[str]:
        refs = [line.strip() for line in self.references_edit.toPlainText().splitlines() if line.strip()]
        if not refs:
            raise ValueError("Provide at least one reference reaction.")
        return refs

    def _load_text_file(self) -> Optional[str]:
        path, _ = QFileDialog.getOpenFileName(self, "Select text file", "", "Text/Markdown (*.txt *.md);;All files (*.*)")
        if not path:
            return None
        with open(path, "r", encoding="utf-8") as handle:
            return handle.read()

    def load_references_from_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Load rule database", "", "JSON files (*.json);;All files (*.*)"
        )
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8") as handle:
                data = json.load(handle)
            refs = data.get("reaction", {}).get("reference_reactions", [])
            if not refs:
                raise ValueError("No reference_reactions found in file.")
            self.references_edit.setPlainText("\n".join(refs))
        except Exception as exc:
            QMessageBox.critical(self, "Load failed", str(exc))

    def load_protocol_from_file(self) -> None:
        text = self._load_text_file()
        if text is not None:
            self.protocol_edit.setPlainText(text)

    def run_autofill(self) -> None:
        try:
            params = RuleBuilderAutoInput(
                family=self.family_edit.text().strip(),
                metadata_id=self.metadata_id_edit.text().strip(),
                metadata_name=self.metadata_name_edit.text().strip(),
                metadata_version=self.metadata_version_edit.text().strip(),
                created_date=self.created_date_edit.text().strip() or None,
                status=self.status_edit.text().strip() or "draft",
                tags=_split_csv(self.tags_edit.text()),
                reference_reactions=self._gather_references(),
                protocol_text=self.protocol_edit.toPlainText().strip(),
                notes=self.notes_edit.toPlainText().strip() or None,
                desired_focus=self.focus_edit.text().strip() or None,
                applies_if_hints=_split_csv(self.applies_hint_edit.text()),
                modifier_hints=_split_csv(self.modifier_hint_edit.text()),
                max_base_rules=int(self.max_rules_edit.text().strip() or "4"),
            )
        except Exception as exc:
            QMessageBox.critical(self, "Invalid input", str(exc))
            return

        try:
            result = run_rule_builder_autofill(params)
        except Exception as exc:
            QMessageBox.critical(self, "Autofill failed", str(exc))
            return

        if not result.get("success"):
            QMessageBox.critical(self, "Autofill error", result.get("error", "Unknown error"))
            return

        self.latest_result = result
        rule_db = result.get("rule_database", {})
        issues = result.get("issues", [])
        issue_lines = "\n".join(f"- [{issue.get('severity')}] {issue.get('field')}: {issue.get('message')}" for issue in issues)
        display_text = json.dumps(rule_db, indent=2)
        if issue_lines:
            display_text = f"Validation issues:\n{issue_lines}\n\n{display_text}"
        self.result_display.setPlainText(display_text)
        QMessageBox.information(self, "Autofill complete", "Draft generated. Review the preview pane.")

    def save_result(self) -> None:
        if not self.latest_result:
            QMessageBox.information(self, "No result", "Run autofill first.")
            return
        rule_db = self.latest_result.get("rule_database", {})
        path, _ = QFileDialog.getSaveFileName(self, "Save rule database", "rule_db_autofill.json", "JSON files (*.json)")
        if not path:
            return
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(rule_db, handle, indent=2)
        QMessageBox.information(self, "Saved", f"Draft saved to {path}")

    def accept_draft(self) -> None:
        if not self.latest_result:
            QMessageBox.information(self, "No draft", "Run autofill first.")
            return
        self.accepted_data = self.latest_result.get("rule_database")
        if not self.accepted_data:
            QMessageBox.critical(self, "Error", "Draft missing rule_database payload.")
            return
        self.accept()


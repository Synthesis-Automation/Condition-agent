"""
Interactive session helpers for the rule-builder workflow.

This module sits between the deterministic :class:`RuleBuilder` utilities and
higher-level UX (CLI, LLM agent). It guides a user through collecting the
required fields for a rule database, while keeping the prompt logic testable
and re-usable.
"""

from __future__ import annotations

from typing import Any, Callable, Dict, Iterable, List, Optional

from chemtools.rule import RuleBuilder, ValidationIssue

PromptFn = Callable[[str], str]
OutputFn = Callable[[str], None]


class RuleBuilderSession:
    """Guide users through assembling a rule database via prompts."""

    def __init__(
        self,
        builder: RuleBuilder,
        *,
        prompt_fn: Optional[PromptFn] = None,
        output_fn: Optional[OutputFn] = None,
    ):
        self.builder = builder
        self._prompt_fn = prompt_fn or input  # type: ignore[assignment]
        self._output_fn = output_fn or print  # type: ignore[assignment]

    # ------------------------------------------------------------------ #
    # Public API
    # ------------------------------------------------------------------ #

    def run_wizard(self) -> None:
        """Execute the full interactive wizard."""
        self._output("\n=== Rule Builder Wizard ===")
        self.collect_metadata()
        self.collect_reference_reactions()
        self.collect_scope_and_notes()
        self.collect_applies_if()
        self.collect_default_rule()
        self.collect_base_rules()
        self.collect_modifiers()
        self._output("Wizard complete. Run validation before saving.\n")

    def validate(self, strict: bool = False) -> List[ValidationIssue]:
        """Proxy to the underlying builder validation."""
        return self.builder.validate(strict=strict)

    # ------------------------------------------------------------------ #
    # Collection steps
    # ------------------------------------------------------------------ #

    def collect_metadata(self) -> None:
        self._output("\n--- Metadata ---")
        meta = self.builder.metadata
        fields = [
            ("id", "Metadata id"),
            ("name", "Metadata name"),
            ("version", "Metadata version"),
            ("created_date", "Metadata created date (YYYY-MM-DD)"),
            ("status", "Metadata status"),
        ]
        updates: Dict[str, Any] = {}
        for field, label in fields:
            default = self._default_str(meta.get(field))
            value = self._ask(label, default=default)
            updates[field] = value.strip()
        tags_default = ", ".join(meta.get("tags", []))
        tags_text = self._ask(
            "Metadata tags (comma separated)", default=tags_default or None
        )
        updates["tags"] = self._split_list(tags_text)
        self.builder.set_metadata(**updates)

    def collect_reference_reactions(self) -> None:
        self._output("\n--- Reference Reactions ---")
        refs = list(self.builder.reaction.get("reference_reactions", []))
        if refs:
            self._output(
                f"Existing references ({len(refs)}):\n  "
                + "\n  ".join(refs)
            )
        self._output(
            "Enter reaction SMILES (format: reactants>>products). "
            "Leave blank to finish."
        )
        while True:
            rxn = self._ask("Reference reaction SMILES (blank to finish)")
            if not rxn:
                if self.builder.reaction.get("reference_reactions"):
                    break
                self._output("At least one reference reaction is required.")
                continue
            self.builder.add_reference_reactions([rxn.strip()])

    def collect_scope_and_notes(self) -> None:
        self._output("\n--- Scope & Notes ---")
        scope = self.builder.reaction.get("scope", {})
        scope_type = self._ask(
            "Scope type (optional)",
            default=self._default_str(scope.get("scope_type")),
        )
        compat = self._ask(
            "Compatible groups (comma separated, optional)",
            default=", ".join(scope.get("compatible_functional_groups", [])) or None,
        )
        incompatible = self._ask(
            "Incompatible groups (comma separated, optional)",
            default=", ".join(scope.get("incompatible_functional_groups", [])) or None,
        )
        self.builder.set_scope(
            scope_type=scope_type.strip() or None,
            compatible=self._split_list(compat),
            incompatible=self._split_list(incompatible),
        )
        notes = self._ask(
            "Notes / protocol summary (optional)",
            default=self._default_str(self.builder.reaction.get("notes")),
        )
        if notes.strip():
            self.builder.set_notes(notes.strip())

    def collect_applies_if(self) -> None:
        self._output("\n--- Applicability ---")
        applies = self.builder.data.get("applies_if", {})
        all_feats = self._ask(
            "Applicability – required ALL features (comma separated)",
            default=", ".join(applies.get("all", [])) or None,
        )
        any_feats = self._ask(
            "Applicability – ANY-of features (comma separated, optional)",
            default=", ".join(applies.get("any", [])) or None,
        )
        self.builder.set_applies_if(
            all_features=self._split_list(all_feats),
            any_features=self._split_list(any_feats),
        )

    def collect_default_rule(self) -> None:
        self._output("\n--- Default Rule ---")
        default_rule = self.builder.data.get("default_rule", {})
        rule_id = self._ask(
            "Default rule id", default=self._default_str(default_rule.get("id"))
        ).strip()
        description = self._ask(
            "Default rule description",
            default=self._default_str(default_rule.get("description")),
        ).strip()
        while True:
            conditions = self._collect_conditions(
                default_rule.get("conditions", {}),
                prompt="Default conditions key=value (blank to finish)",
            )
            if conditions:
                break
            self._output("Default rule requires at least one condition.")
        self.builder.set_default_rule(
            rule_id=rule_id or None,
            description=description or None,
            conditions=conditions,
        )

    def collect_base_rules(self) -> None:
        self._output("\n--- Base Rules ---")
        while True:
            has_rules = bool(self.builder.data.get("base_rules"))
            default_choice = "n" if has_rules else "y"
            choice = self._ask(
                "Add or update base rule? (y/n)", default=default_choice
            ).strip().lower()
            if choice not in {"y", "yes"}:
                if not has_rules:
                    self._output("At least one base rule is required.")
                    continue
                break
            rule_id = self._ask("Base rule id").strip()
            if not rule_id:
                self._output("Base rule id cannot be empty.")
                continue
            existing = self._find_base_rule(rule_id)
            name = self._ask(
                "Base rule name",
                default=self._default_str(existing.get("name") if existing else None),
            ).strip()
            description = self._ask(
                "Base rule description",
                default=self._default_str(
                    existing.get("description") if existing else None
                ),
            ).strip()
            all_feats = self._ask(
                "Base rule ALL features (comma separated, optional)",
                default=", ".join(
                    (existing or {}).get("reactant_features", {}).get("all", [])
                )
                or None,
            )
            any_feats = self._ask(
                "Base rule ANY features (comma separated, optional)",
                default=", ".join(
                    (existing or {}).get("reactant_features", {}).get("any", [])
                )
                or None,
            )
            priority_text = self._ask(
                "Base rule priority (integer, optional)",
                default=self._default_str(
                    (existing or {}).get("priority") if existing else None
                ),
            ).strip()
            priority = None
            if priority_text:
                try:
                    priority = int(priority_text)
                except ValueError:
                    self._output("Invalid priority; ignoring.")
            while True:
                conditions = self._collect_conditions(
                    (existing or {}).get("conditions", {}),
                    prompt="Base rule conditions key=value (blank to finish)",
                )
                if conditions:
                    break
                self._output("Base rule conditions cannot be empty.")
            reactant_features: Dict[str, Any] = {}
            all_list = self._split_list(all_feats)
            any_list = self._split_list(any_feats)
            if all_list:
                reactant_features["all"] = all_list
            if any_list:
                reactant_features["any"] = any_list
            self.builder.upsert_base_rule(
                rule_id,
                name=name or rule_id,
                description=description or "",
                reactant_features=reactant_features or {},
                conditions=conditions,
                priority=priority,
            )

    def collect_modifiers(self) -> None:
        self._output("\n--- Modifiers (optional) ---")
        while True:
            choice = self._ask("Add or update modifier? (y/n)", default="n").strip().lower()
            if choice not in {"y", "yes"}:
                break
            mod_id = self._ask("Modifier id").strip()
            if not mod_id:
                self._output("Modifier id cannot be empty.")
                continue
            existing = self._find_modifier(mod_id)
            default_triggers = (
                ", ".join(existing.get("when", [])) if existing else None
            )
            when_text = self._ask(
                "Modifier triggers (comma separated)",
                default=default_triggers or None,
            )
            when = self._split_list(when_text)
            if not when:
                self._output("At least one trigger is required.")
                continue
            suggestion = self._ask(
                "Modifier suggestion text",
                default=self._default_str(existing.get("suggest") if existing else None),
            ).strip()
            if not suggestion:
                self._output("Suggestion cannot be empty.")
                continue
            rationale = self._ask(
                "Modifier rationale (optional)",
                default=self._default_str(existing.get("rationale") if existing else None),
            ).strip()
            self.builder.upsert_modifier(
                mod_id, when=when, suggestion=suggestion, rationale=rationale or None
            )

    # ------------------------------------------------------------------ #
    # Helpers
    # ------------------------------------------------------------------ #

    def _ask(self, prompt: str, default: Optional[str] = None) -> str:
        suffix = f" [{default}]" if default else ""
        response = self._prompt_fn(f"{prompt}{suffix}: ").strip()
        if not response and default is not None:
            return default
        return response

    def _split_list(self, value: Optional[str]) -> List[str]:
        if not value:
            return []
        return [item.strip() for item in value.split(",") if item.strip()]

    def _collect_conditions(
        self,
        existing: Dict[str, Any],
        *,
        prompt: str,
    ) -> Dict[str, Any]:
        self._output(
            "Enter conditions as key=value pairs. Leave blank when finished."
        )
        if existing:
            self._output(
                "Current values:\n  " + "\n  ".join(f"{k} = {v}" for k, v in existing.items())
            )
        entries = dict(existing)
        while True:
            pair = self._ask(prompt)
            if not pair:
                break
            if "=" not in pair:
                self._output("  Please use key=value format.")
                continue
            key, value = pair.split("=", 1)
            entries[key.strip()] = value.strip()
        return {k: v for k, v in entries.items() if v}

    def _find_base_rule(self, rule_id: str) -> Dict[str, Any]:
        for rule in self.builder.data.get("base_rules", []):
            if rule.get("id") == rule_id:
                return rule
        return {}

    def _find_modifier(self, modifier_id: str) -> Dict[str, Any]:
        for modifier in self.builder.data.get("modifiers", []):
            if modifier.get("id") == modifier_id:
                return modifier
        return {}

    def _default_str(self, value: Any) -> Optional[str]:
        if value is None:
            return None
        if isinstance(value, list):
            return ", ".join(value)
        return str(value)

    def _output(self, message: str) -> None:
        self._output_fn(message)


__all__ = ["RuleBuilderSession"]

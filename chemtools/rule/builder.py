"""
Rule Database Builder Utilities.

Provides a deterministic API for assembling and editing rule database JSON
files (schema v2.0). The builder keeps data in-memory as Python dictionaries,
performs structural validation, and can emit changes back to disk for use by
the rule engine and unified recommender.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional
import copy
import difflib
import json

from .database import RuleDatabase


@dataclass(frozen=True)
class ValidationIssue:
    """Represents a validation issue detected by the builder."""

    field: str
    message: str
    severity: str = "error"  # either "error" or "warning"


class RuleBuilder:
    """
    Imperative helper for constructing rule databases programmatically.

    Typical usage:

        builder = RuleBuilder.new(family="Suzuki_Miyaura")
        builder.set_metadata(id="suzuki_v3", name="Suzuki", version="v3")
        builder.add_reference_reactions([...])
        builder.set_default_rule(...)
        builder.upsert_base_rule(...)
        builder.validate()
        builder.save("data/rule_db_v2/Suzuki_db.json")
    """

    DEFAULT_SCHEMA_VERSION = "2.0"
    DEFAULT_SOURCE_TYPE = "rule"

    def __init__(self, data: Optional[Dict[str, Any]] = None):
        self._data: Dict[str, Any] = copy.deepcopy(data) if data else {}
        if not self._data:
            self._data = {
                "schema_version": self.DEFAULT_SCHEMA_VERSION,
                "source_type": self.DEFAULT_SOURCE_TYPE,
                "metadata": {},
                "reaction": {
                    "family": "",
                    "reference_reactions": [],
                    "scope": {
                        "scope_type": "",
                        "compatible_functional_groups": [],
                        "incompatible_functional_groups": [],
                    },
                    "notes": "",
                },
                "applies_if": {},
                "default_rule": {"conditions": {}},
                "base_rules": [],
                "modifiers": [],
            }
        self._original_snapshot = copy.deepcopy(self._data)

    # ------------------------------------------------------------------ #
    # Constructors
    # ------------------------------------------------------------------ #

    @classmethod
    def from_file(cls, path: str | Path) -> "RuleBuilder":
        """Load an existing rule database from disk."""
        path = Path(path)
        data = json.loads(path.read_text(encoding="utf-8"))
        builder = cls(data=data)
        builder._source_path = path
        return builder

    @classmethod
    def new(cls, family: str, metadata: Optional[Dict[str, Any]] = None) -> "RuleBuilder":
        """Create a builder pre-populated with minimal required fields."""
        metadata = metadata or {}
        builder = cls()
        builder._data["reaction"]["family"] = family
        builder._data["metadata"].update(metadata)
        return builder

    # ------------------------------------------------------------------ #
    # Convenience getters
    # ------------------------------------------------------------------ #

    @property
    def data(self) -> Dict[str, Any]:
        """Return the raw mutable dictionary."""
        return self._data

    @property
    def metadata(self) -> Dict[str, Any]:
        return self._data.setdefault("metadata", {})

    @property
    def reaction(self) -> Dict[str, Any]:
        return self._data.setdefault("reaction", {})

    # ------------------------------------------------------------------ #
    # Mutation helpers
    # ------------------------------------------------------------------ #

    def set_metadata(self, **fields: Any) -> None:
        """Update metadata fields."""
        self.metadata.update({k: v for k, v in fields.items() if v is not None})

    def set_reaction_family(self, family: str) -> None:
        self.reaction["family"] = family

    def add_reference_reactions(self, reactions: Iterable[str]) -> None:
        refs = self.reaction.setdefault("reference_reactions", [])
        for rxn in reactions:
            if rxn and rxn not in refs:
                refs.append(rxn)

    def set_scope(
        self,
        scope_type: Optional[str] = None,
        compatible: Optional[Iterable[str]] = None,
        incompatible: Optional[Iterable[str]] = None,
    ) -> None:
        scope = self.reaction.setdefault("scope", {})
        if scope_type is not None:
            scope["scope_type"] = scope_type
        if compatible is not None:
            scope["compatible_functional_groups"] = list(dict.fromkeys(compatible))
        if incompatible is not None:
            scope["incompatible_functional_groups"] = list(
                dict.fromkeys(incompatible)
            )

    def set_notes(self, notes: str) -> None:
        self.reaction["notes"] = notes

    def set_applies_if(self, *, all_features=None, any_features=None, raw=None) -> None:
        if raw is not None:
            self._data["applies_if"] = raw
            return
        applies: Dict[str, Any] = {}
        if all_features:
            applies["all"] = list(dict.fromkeys(all_features))
        if any_features:
            applies["any"] = list(dict.fromkeys(any_features))
        self._data["applies_if"] = applies

    def set_default_rule(
        self,
        *,
        rule_id: Optional[str] = None,
        description: Optional[str] = None,
        conditions: Optional[Dict[str, Any]] = None,
    ) -> None:
        default_rule = self._data.setdefault("default_rule", {})
        if rule_id is not None:
            default_rule["id"] = rule_id
        if description is not None:
            default_rule["description"] = description
        if conditions is not None:
            default_rule["conditions"] = conditions

    def upsert_base_rule(
        self,
        rule_id: str,
        *,
        name: str,
        description: str,
        reactant_features: Dict[str, Any],
        conditions: Dict[str, Any],
        priority: Optional[int] = None,
    ) -> None:
        base_rules = self._data.setdefault("base_rules", [])
        payload = {
            "id": rule_id,
            "name": name,
            "description": description,
            "reactant_features": reactant_features,
            "conditions": conditions,
        }
        if priority is not None:
            payload["priority"] = priority
        for idx, rule in enumerate(base_rules):
            if rule.get("id") == rule_id:
                base_rules[idx] = payload
                break
        else:
            base_rules.append(payload)

    def remove_base_rule(self, rule_id: str) -> bool:
        base_rules = self._data.get("base_rules", [])
        for idx, rule in enumerate(base_rules):
            if rule.get("id") == rule_id:
                base_rules.pop(idx)
                return True
        return False

    def upsert_modifier(
        self,
        modifier_id: str,
        *,
        when: Iterable[str],
        suggestion: str,
        rationale: Optional[str] = None,
    ) -> None:
        modifiers = self._data.setdefault("modifiers", [])
        payload = {
            "id": modifier_id,
            "when": list(dict.fromkeys(when)),
            "suggest": suggestion,
        }
        if rationale:
            payload["rationale"] = rationale
        for idx, mod in enumerate(modifiers):
            if mod.get("id") == modifier_id:
                modifiers[idx] = payload
                break
        else:
            modifiers.append(payload)

    def remove_modifier(self, modifier_id: str) -> bool:
        modifiers = self._data.get("modifiers", [])
        for idx, mod in enumerate(modifiers):
            if mod.get("id") == modifier_id:
                modifiers.pop(idx)
                return True
        return False

    # ------------------------------------------------------------------ #
    # Validation & persistence
    # ------------------------------------------------------------------ #

    def to_dict(self) -> Dict[str, Any]:
        return copy.deepcopy(self._data)

    def save(self, path: Optional[str | Path] = None, *, indent: int = 2) -> Path:
        """Write the current builder state to disk."""
        target = Path(path) if path else getattr(self, "_source_path", None)
        if target is None:
            raise ValueError("No path specified for save().")
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(self._data, indent=indent), encoding="utf-8")
        self._original_snapshot = copy.deepcopy(self._data)
        self._source_path = target
        return target

    def diff(self) -> str:
        """Return a unified diff between the original snapshot and current state."""
        old = json.dumps(self._original_snapshot or {}, indent=2, sort_keys=True).splitlines()
        new = json.dumps(self._data, indent=2, sort_keys=True).splitlines()
        return "\n".join(
            difflib.unified_diff(old, new, fromfile="original", tofile="current")
        )

    def validate(self, *, strict: bool = True) -> List[ValidationIssue]:
        """Validate the builder contents, optionally raising if errors exist."""
        issues = self._collect_structural_issues()
        try:
            rule_db = RuleDatabase.from_dict(self.to_dict())
            for msg in rule_db.validate():
                issues.append(ValidationIssue(field="rule_database", message=msg))
        except Exception as exc:  # pragma: no cover - defensive
            issues.append(ValidationIssue(field="rule_database", message=str(exc)))

        if strict:
            errors = [issue for issue in issues if issue.severity == "error"]
            if errors:
                details = "\n".join(f"- [{i.field}] {i.message}" for i in errors)
                raise ValueError(f"RuleBuilder validation failed:\n{details}")
        return issues

    # ------------------------------------------------------------------ #
    # Internal helpers
    # ------------------------------------------------------------------ #

    def _collect_structural_issues(self) -> List[ValidationIssue]:
        issues: List[ValidationIssue] = []
        meta = self.metadata
        required_meta = ["id", "name", "version", "created_date", "status", "tags"]
        for field in required_meta:
            value = meta.get(field)
            if not value:
                issues.append(ValidationIssue(field=f"metadata.{field}", message="Missing value"))
        tags = meta.get("tags")
        if tags is not None and not isinstance(tags, list):
            issues.append(ValidationIssue(field="metadata.tags", message="Must be a list"))

        reaction = self.reaction
        if not reaction.get("family"):
            issues.append(ValidationIssue(field="reaction.family", message="Required"))

        refs = reaction.get("reference_reactions", [])
        if not refs:
            issues.append(ValidationIssue(field="reaction.reference_reactions", message="Provide at least one reaction"))
        else:
            for idx, rxn in enumerate(refs):
                if ">>" not in rxn:
                    issues.append(
                        ValidationIssue(
                            field=f"reaction.reference_reactions[{idx}]",
                            message="Should look like reaction SMILES (contains '>>')",
                            severity="warning",
                        )
                    )

        applies = self._data.get("applies_if", {})
        if applies and not any(k in applies for k in ("all", "any")):
            issues.append(
                ValidationIssue(
                    field="applies_if",
                    message="Use 'all' or 'any' keys to describe applicability",
                )
            )

        default_rule = self._data.get("default_rule", {})
        if not default_rule.get("conditions"):
            issues.append(
                ValidationIssue(field="default_rule.conditions", message="Provide baseline conditions")
            )

        base_rules = self._data.get("base_rules", [])
        if not base_rules:
            issues.append(
                ValidationIssue(field="base_rules", message="At least one base rule is required")
            )
        else:
            seen_ids: set[str] = set()
            for idx, rule in enumerate(base_rules):
                rid = rule.get("id")
                if not rid:
                    issues.append(
                        ValidationIssue(
                            field=f"base_rules[{idx}].id", message="Each rule needs an id"
                        )
                    )
                elif rid in seen_ids:
                    issues.append(
                        ValidationIssue(
                            field=f"base_rules[{idx}].id", message="Duplicate rule id"
                        )
                    )
                seen_ids.add(rid)
                if not rule.get("conditions"):
                    issues.append(
                        ValidationIssue(
                            field=f"base_rules[{idx}].conditions",
                            message="Conditions dict cannot be empty",
                        )
                    )

        modifiers = self._data.get("modifiers", [])
        for idx, mod in enumerate(modifiers):
            if not mod.get("when"):
                issues.append(
                    ValidationIssue(
                        field=f"modifiers[{idx}].when",
                        message="Provide at least one trigger",
                    )
                )
            if not mod.get("suggest"):
                issues.append(
                    ValidationIssue(
                        field=f"modifiers[{idx}].suggest",
                        message="Suggestion text required",
                    )
                )

        return issues


__all__ = ["RuleBuilder", "ValidationIssue"]

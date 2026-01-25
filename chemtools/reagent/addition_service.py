"""
UI-friendly helpers for adding reagents to the flattened registry CSV.

This module keeps the UI thin by centralizing path handling, taxonomy
lookups, and input normalization for reagent additions.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence
import re

from .registry_addition import add_reagent_entry, ReagentAdditionError
from .registry_store import RegistryStore
from .taxonomy_store import TaxonomyStore
from .taxonomy_utils import dedupe_synonyms, normalize_cas, resolve_identity_from_cas

try:  # Optional LLM support
    from llmtools.clients import LLMClient
    from llmtools.reagent_classifier import classify_role, assign_fields
    _LLM_AVAILABLE = True
except Exception:
    LLMClient = None  # type: ignore[assignment]
    classify_role = None  # type: ignore[assignment]
    assign_fields = None  # type: ignore[assignment]
    _LLM_AVAILABLE = False

_SPLIT_PATTERN = re.compile(r"[,\n;]+")


def _coerce_path(value: Optional[str | Path]) -> Optional[Path]:
    if value is None:
        return None
    if isinstance(value, Path):
        return value.expanduser()
    text = str(value).strip()
    return Path(text).expanduser() if text else None


def parse_synonym_text(text: str) -> List[str]:
    """Split a free-form synonym field into a list."""
    if not text:
        return []
    parts = [item.strip() for item in _SPLIT_PATTERN.split(text) if item and item.strip()]
    if not parts:
        return []
    return dedupe_synonyms(parts)


class ReagentAdditionService:
    """Convenience wrapper for UI clients."""

    def __init__(
        self,
        *,
        registry_dir: Optional[str | Path] = None,
        taxonomy_dir: Optional[str | Path] = None,
    ) -> None:
        self.registry_dir = _coerce_path(registry_dir) or Path("data/reagent_db")
        self.taxonomy_dir = _coerce_path(taxonomy_dir)
        self.taxonomy = TaxonomyStore(self.taxonomy_dir)

    def update_paths(
        self,
        *,
        registry_dir: Optional[str | Path] = None,
        taxonomy_dir: Optional[str | Path] = None,
    ) -> None:
        """Update registry/taxonomy directories and reload taxonomy."""
        if registry_dir is not None:
            resolved = _coerce_path(registry_dir)
            if resolved is not None:
                self.registry_dir = resolved
        if taxonomy_dir is not None:
            self.taxonomy_dir = _coerce_path(taxonomy_dir)
        self.taxonomy = TaxonomyStore(self.taxonomy_dir)

    def list_roles(self) -> List[Dict[str, Any]]:
        """Return roles with labels/descriptions from the taxonomy."""
        roles: List[Dict[str, Any]] = []
        for role_id, data in self.taxonomy.role_data.items():
            priority = data.get("priority")
            try:
                priority_val = int(priority) if priority is not None else 99
            except (TypeError, ValueError):
                priority_val = 99
            roles.append(
                {
                    "role": role_id,
                    "label": data.get("name") or role_id,
                    "description": data.get("description") or "",
                    "priority": priority_val,
                }
            )
        roles.sort(key=lambda item: (item.get("priority", 99), item.get("role", "")))
        return roles

    def list_families(self, role: Optional[str] = None) -> List[Dict[str, Any]]:
        """Return family metadata for a role or for all roles."""
        families = self.taxonomy.list_families(role)
        return [
            {"role": role_id, "family_id": family_id, "label": label}
            for role_id, family_id, label in families
        ]

    def find_existing(self, cas: str) -> Optional[Dict[str, Any]]:
        """Lookup an existing reagent entry by CAS."""
        normalized = normalize_cas(cas)
        registry = RegistryStore(self.registry_dir)
        found = registry.find_by_cas(normalized)
        if not found:
            return None
        role, family_id, entry = found
        return {"role": role, "family_id": family_id, "entry": entry}

    def build_csv_row(self, entry: Dict[str, Any], role: Optional[str] = None) -> tuple[List[str], Dict[str, str]]:
        """Return the csv header and row mapping for a registry entry."""
        registry = RegistryStore(self.registry_dir)
        row = registry.entry_to_csv_row(entry, role=role)
        return registry.csv_header(), row

    def format_csv_preview(self, entry: Dict[str, Any], role: Optional[str] = None) -> str:
        """Format a registry entry as edit-friendly text."""
        registry = RegistryStore(self.registry_dir)
        row = registry.entry_to_csv_row(entry, role=role)
        header = registry.csv_header()
        width = max((len(field) for field in header), default=8)
        lines = ["Editable fields (reagents.csv order):"]
        for field in header:
            value = row.get(field, "")
            lines.append(f"{field.ljust(width)} : {value}")
        return "\n".join(lines)

    def save_csv_row(self, row: Dict[str, Any]) -> Dict[str, Any]:
        """Persist a CSV row mapping to the registry."""
        registry = RegistryStore(self.registry_dir)
        entry = registry.row_to_entry(row)
        role = entry.get("role") or ""
        if not role:
            role = str(row.get("role_1") or "").strip()
        if not role:
            raise ReagentAdditionError("role_1 is required to save a reagent entry.")
        registry.add_entry(role, entry)
        path = registry.save_role(role)
        return {
            "status": "written",
            "written_to": str(path),
            "entry_preview": entry,
            "role": role,
            "family_id": entry.get("family_id") or "",
        }

    def preview_entry(
        self,
        *,
        cas: str,
        name: Optional[str] = None,
        synonyms: Optional[Sequence[str]] = None,
        abbreviation: Optional[str] = None,
        tag: Optional[str] = None,
        role: Optional[str] = None,
        family_id: Optional[str] = None,
        smiles: Optional[str] = None,
        allow_default_family: bool = False,
        auto_resolve: bool = True,
        ai_assist: bool = False,
        llm_provider: Optional[str] = None,
        llm_model: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Run the addition workflow without writing to disk."""
        if ai_assist:
            return self._ai_assist_entry(
                cas=cas,
                name=name,
                synonyms=synonyms,
                abbreviation=abbreviation,
                tag=tag,
                role=role,
                family_id=family_id,
                smiles=smiles,
                allow_default_family=allow_default_family,
                auto_resolve=auto_resolve,
                llm_provider=llm_provider,
                llm_model=llm_model,
                dry_run=True,
            )
        return add_reagent_entry(
            cas=cas,
            taxonomy_dir=self.taxonomy_dir,
            registry_dir=self.registry_dir,
            name=name,
            synonyms=synonyms,
            abbreviation=abbreviation,
            tag=tag,
            role=role,
            family_id=family_id,
            smiles=smiles,
            allow_default_family=allow_default_family,
            dry_run=True,
            auto_resolve=auto_resolve,
        )

    def save_entry(
        self,
        *,
        cas: str,
        name: Optional[str] = None,
        synonyms: Optional[Sequence[str]] = None,
        abbreviation: Optional[str] = None,
        tag: Optional[str] = None,
        role: Optional[str] = None,
        family_id: Optional[str] = None,
        smiles: Optional[str] = None,
        allow_default_family: bool = False,
        auto_resolve: bool = True,
        ai_assist: bool = False,
        llm_provider: Optional[str] = None,
        llm_model: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Run the addition workflow and write to the registry CSV."""
        if ai_assist:
            return self._ai_assist_entry(
                cas=cas,
                name=name,
                synonyms=synonyms,
                abbreviation=abbreviation,
                tag=tag,
                role=role,
                family_id=family_id,
                smiles=smiles,
                allow_default_family=allow_default_family,
                auto_resolve=auto_resolve,
                llm_provider=llm_provider,
                llm_model=llm_model,
                dry_run=False,
            )
        return add_reagent_entry(
            cas=cas,
            taxonomy_dir=self.taxonomy_dir,
            registry_dir=self.registry_dir,
            name=name,
            synonyms=synonyms,
            abbreviation=abbreviation,
            tag=tag,
            role=role,
            family_id=family_id,
            smiles=smiles,
            allow_default_family=allow_default_family,
            dry_run=False,
            auto_resolve=auto_resolve,
        )

    def _ai_assist_entry(
        self,
        *,
        cas: str,
        name: Optional[str],
        synonyms: Optional[Sequence[str]],
        abbreviation: Optional[str],
        tag: Optional[str],
        role: Optional[str],
        family_id: Optional[str],
        smiles: Optional[str],
        allow_default_family: bool,
        auto_resolve: bool,
        llm_provider: Optional[str],
        llm_model: Optional[str],
        dry_run: bool,
    ) -> Dict[str, Any]:
        if not _LLM_AVAILABLE:
            raise ReagentAdditionError("AI-assist requires llmtools and configured API keys.")

        normalized_cas = normalize_cas(cas)
        identity: Dict[str, Any] = {"cas": normalized_cas}
        resolved: Dict[str, Any] = {}
        if auto_resolve:
            resolved = resolve_identity_from_cas(normalized_cas) or {}

        identity["name"] = (name or resolved.get("name") or "").strip()
        if not identity["name"]:
            if auto_resolve:
                raise ReagentAdditionError(
                    f"AI-assist could not auto-resolve a name for CAS {normalized_cas}. "
                    "Provide a name or check the CAS/connection."
                )
            raise ReagentAdditionError("AI-assist needs a name. Provide it or enable auto-resolve.")
        identity["smiles"] = (smiles or resolved.get("smiles") or "").strip() or None
        resolved_synonyms = resolved.get("synonyms") or []
        resolved_formula = resolved.get("molecular_formula")
        resolved_mw = resolved.get("molecular_weight")
        resolved_density = resolved.get("density")
        resolved_bp = resolved.get("boiling_point")
        resolved_mp = resolved.get("melting_point")
        provided_synonyms = dedupe_synonyms(list(synonyms or []))
        identity["synonyms"] = dedupe_synonyms([*provided_synonyms, *resolved_synonyms])
        if resolved_formula:
            identity["molecular_formula"] = resolved_formula
        if resolved_mw is not None:
            identity["molecular_weight"] = resolved_mw
        if resolved_density:
            identity["density"] = resolved_density
        if resolved_bp:
            identity["boiling_point"] = resolved_bp
        if resolved_mp:
            identity["melting_point"] = resolved_mp

        provider = (llm_provider or os.getenv("LLM_PROVIDER") or "openai").strip().lower()
        model = (llm_model or os.getenv("LLM_MODEL") or "").strip() or None
        model_type = (os.getenv("LLM_MODEL_TYPE") or "balanced").strip()

        try:
            if model:
                llm_client = LLMClient(provider=provider, model=model)
            else:
                llm_client = LLMClient.from_env(provider=provider, model_type=model_type)
        except Exception as exc:
            raise ReagentAdditionError(f"AI-assist LLM client error: {exc}") from exc

        llm_role = role
        if not llm_role:
            role_result = classify_role(identity, llm_client)
            if role_result.get("status") != "success":
                llm_role = "other_reagent"
            else:
                llm_role = role_result.get("role")
        if not llm_role:
            llm_role = "other_reagent"

        llm_family = family_id
        if not llm_family:
            fields_result = assign_fields(identity, llm_role, self.registry_dir, llm_client)
            if fields_result.get("status") != "success":
                raise ReagentAdditionError(fields_result.get("error") or "AI family assignment failed.")
            llm_family = fields_result.get("family")

        return add_reagent_entry(
            cas=normalized_cas,
            taxonomy_dir=self.taxonomy_dir,
            registry_dir=self.registry_dir,
            name=identity.get("name"),
            formula=identity.get("molecular_formula"),
            molecular_weight=identity.get("molecular_weight"),
            density=identity.get("density"),
            boiling_point=identity.get("boiling_point"),
            melting_point=identity.get("melting_point"),
            synonyms=identity.get("synonyms"),
            abbreviation=abbreviation,
            tag=tag,
            role=llm_role,
            family_id=llm_family,
            smiles=identity.get("smiles"),
            allow_default_family=allow_default_family,
            dry_run=dry_run,
            auto_resolve=False,
        )


__all__ = [
    "ReagentAdditionService",
    "ReagentAdditionError",
    "parse_synonym_text",
]

"""Taxonomy-owned normalization of external reactive-handle labels."""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Literal, Optional

from .api import featurize_molecule
from .validation import validate_taxonomy


_PATH = Path(__file__).with_name("definitions") / "source_label_crosswalk.v1.json"
MappingStatus = Literal["canonical", "qualified", "unresolved"]
SignatureScope = Literal["exact", "family"]


@dataclass(frozen=True)
class SourceLabelMapping:
    """One source label resolved to a handle plus optional constraints."""

    source_label: str
    base_label: str
    display_label: str
    representative_smiles: str
    canonical_signature: str
    mapping_status: MappingStatus
    signature_scope: SignatureScope = "exact"
    center_substitution_class: Optional[str] = None
    attachment_carbon_class: Optional[str] = None
    alpha_branched: Optional[bool] = None
    qualifier_scope: Optional[str] = None

    def to_columns(self, prefix: str) -> Dict[str, str]:
        """Return consistently named flat-CSV fields for one reactive site."""
        alpha_branched = (
            "" if self.alpha_branched is None else str(self.alpha_branched).lower()
        )
        return {
            f"{prefix}_source_label": self.source_label,
            f"{prefix}_normalized_label": self.base_label,
            f"{prefix}_display_label": self.display_label,
            f"{prefix}_signature": self.canonical_signature,
            f"{prefix}_center_class": self.center_substitution_class or "",
            f"{prefix}_attachment_class": self.attachment_carbon_class or "",
            f"{prefix}_alpha_branched": alpha_branched,
            f"{prefix}_qualifier_scope": self.qualifier_scope or "",
            f"{prefix}_mapping_status": self.mapping_status,
        }


@lru_cache(maxsize=1)
def load_source_label_mappings() -> Dict[str, SourceLabelMapping]:
    """Load the versioned source-label crosswalk keyed by source label."""
    with _PATH.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    mappings: Dict[str, SourceLabelMapping] = {}
    for raw in payload.get("mappings") or []:
        mapping = SourceLabelMapping(**raw)
        if mapping.source_label in mappings:
            raise ValueError(f"Duplicate source label: {mapping.source_label}")
        mappings[mapping.source_label] = mapping
    return mappings


def resolve_source_label(source_label: str) -> SourceLabelMapping:
    """Resolve a source label, preserving unsupported labels without inference."""
    source = str(source_label or "").strip()
    mapping = load_source_label_mappings().get(source)
    if mapping is not None:
        return mapping
    return SourceLabelMapping(
        source_label=source,
        base_label=source,
        display_label=source,
        representative_smiles="",
        canonical_signature="",
        mapping_status="unresolved",
    )


def validate_source_label_mappings() -> list[str]:
    """Validate mapping signatures, labels, and declared qualifiers."""
    errors = [f"taxonomy:{error}" for error in validate_taxonomy()]
    for source, mapping in load_source_label_mappings().items():
        analysis = featurize_molecule(
            mapping.representative_smiles,
            label_style="hte_legacy",
        )
        environments = {
            environment.site_id: environment
            for environment in analysis.site_environments
        }
        matches = [
            (site, environments[site.site_id])
            for site in analysis.sites
            if (
                site.canonical_signature == mapping.canonical_signature
                if mapping.signature_scope == "exact"
                else site.canonical_signature.startswith(
                    f"{mapping.canonical_signature}|"
                )
            )
            and site.site_id in environments
        ]
        if not analysis.valid or len(matches) != 1:
            errors.append(f"{source}:signature")
            continue
        site, environment = matches[0]
        if site.chemist_label != mapping.base_label:
            errors.append(f"{source}:base_label")
        steric: Dict[str, Any] = dict(environment.steric)
        if (
            mapping.center_substitution_class
            and steric.get("center_substitution_class")
            != mapping.center_substitution_class
        ):
            errors.append(f"{source}:center_substitution_class")
        attached = list(steric.get("attached_groups") or [])
        if mapping.attachment_carbon_class and not any(
            item.get("attachment_carbon_class") == mapping.attachment_carbon_class
            for item in attached
        ):
            errors.append(f"{source}:attachment_carbon_class")
        if mapping.alpha_branched is not None and not any(
            item.get("alpha_branched") is mapping.alpha_branched for item in attached
        ):
            errors.append(f"{source}:alpha_branched")
        if mapping.mapping_status == "qualified" and not mapping.qualifier_scope:
            errors.append(f"{source}:qualifier_scope")
        if (
            mapping.signature_scope == "family"
            and len(mapping.canonical_signature.split("|")) < 2
        ):
            errors.append(f"{source}:signature_scope")
    return errors


__all__ = [
    "MappingStatus",
    "SourceLabelMapping",
    "SignatureScope",
    "load_source_label_mappings",
    "resolve_source_label",
    "validate_source_label_mappings",
]

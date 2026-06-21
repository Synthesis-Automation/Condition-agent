"""Molecule domain data models."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional


@dataclass(frozen=True)
class MoleculeFeatures:
    """Structured molecule feature bundle."""

    smiles: str
    motifs: List[Mapping[str, Any]] = field(default_factory=list)
    context_motifs: List[Mapping[str, Any]] = field(default_factory=list)
    ranked_motifs: List[str] = field(default_factory=list)
    steric: Mapping[str, Any] = field(default_factory=dict)
    electronics: Mapping[str, Any] = field(default_factory=dict)
    nearby: List[Mapping[str, Any]] = field(default_factory=list)
    aryl_analysis: Mapping[str, Any] = field(default_factory=dict)
    analyses: List[Mapping[str, Any]] = field(default_factory=list)
    meta: Mapping[str, Any] = field(default_factory=dict)
    schema_version: str = "v2"

    @classmethod
    def from_payload(cls, payload: Mapping[str, Any]) -> "MoleculeFeatures":
        """Build a model from the existing molecule feature payload shape."""
        return cls(
            schema_version=str(payload.get("schema_version") or "v2"),
            smiles=str(payload.get("smiles") or ""),
            motifs=list(payload.get("motifs") or []),
            context_motifs=list(payload.get("context_motifs") or []),
            ranked_motifs=[str(item) for item in payload.get("ranked_motifs") or []],
            steric=payload.get("steric") or {},
            electronics=payload.get("electronics") or {},
            nearby=list(payload.get("nearby") or []),
            aryl_analysis=payload.get("aryl_analysis") or {},
            analyses=list(payload.get("analyses") or []),
            meta=payload.get("meta") or {},
        )

    def to_payload(self) -> Dict[str, Any]:
        """Return the stable dictionary payload used by current callers."""
        return {
            "schema_version": self.schema_version,
            "smiles": self.smiles,
            "motifs": list(self.motifs),
            "context_motifs": list(self.context_motifs),
            "ranked_motifs": list(self.ranked_motifs),
            "steric": dict(self.steric),
            "electronics": dict(self.electronics),
            "nearby": list(self.nearby),
            "aryl_analysis": dict(self.aryl_analysis),
            "analyses": list(self.analyses),
            "meta": dict(self.meta),
        }


@dataclass(frozen=True)
class MoleculeParseResult:
    """Parsed molecule result independent of feature extraction."""

    input_smiles: str
    canonical_smiles: Optional[str]
    rdkit_available: bool
    valid: bool
    error: Optional[str] = None


__all__ = ["MoleculeFeatures", "MoleculeParseResult"]

from __future__ import annotations

"""Typed data structures for rule-driven reaction condition matchers."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, Literal, Mapping, Sequence, Tuple

from rdkit.Chem.rdchem import Mol


@dataclass(slots=True)
class CompiledSmarts:
    """Container for a compiled SMARTS pattern used during matching."""

    pattern: str
    mol: Mol
    atom_count: int


@dataclass(slots=True)
class SchemeEntry:
    """Representation of a scheme or default-condition row from the database."""

    id: str
    type: Literal["scheme", "default_condition"]
    name: str | None
    reactant_smarts: Tuple[str, ...]
    applies_to: Mapping[str, Any] | None
    feature_requirements: Mapping[str, Any] | None
    priority: int
    conditions: Mapping[str, Any] | None
    default_condition: Mapping[str, Any] | None
    notes: Tuple[str, ...]
    compiled_smarts: Tuple[CompiledSmarts, ...]
    raw: Mapping[str, Any] = field(repr=False)

    def __post_init__(self) -> None:
        if self.type not in {"scheme", "default_condition"}:
            raise ValueError(f"Unsupported entry type '{self.type}' for {self.id!r}")

    @property
    def specificity_score(self) -> int:
        """Heuristic specificity score based on SMARTS size."""

        return sum(item.atom_count for item in self.compiled_smarts)

    def as_dict(self) -> Dict[str, Any]:
        """Return a JSON-friendly representation without the compiled SMARTS."""

        payload = dict(self.raw)
        payload.pop("_compiled", None)
        return payload


@dataclass(slots=True)
class RuleDB:
    """Base metadata shared by all condition rule databases."""

    path: Path
    schema_version: str
    reaction_type: str
    metadata: Dict[str, Any]


@dataclass(slots=True)
class SchemeConditionDB(RuleDB):
    """Top-level container for SMARTS-driven SchemeConditionDB payloads."""

    updated_at: str | None
    entries: Tuple[SchemeEntry, ...]

    def __iter__(self) -> Iterable[SchemeEntry]:
        return iter(self.entries)

    def schemes(self) -> Tuple[SchemeEntry, ...]:
        return tuple(entry for entry in self.entries if entry.type == "scheme")

    def defaults(self) -> Tuple[SchemeEntry, ...]:
        return tuple(entry for entry in self.entries if entry.type == "default_condition")


@dataclass(slots=True)
class SelectorRule:
    """Representation of a feature-driven selector rule."""

    id: str
    conditions: Mapping[str, Any]
    payload: Mapping[str, Any]
    rank_hint: float | None
    raw: Mapping[str, Any] = field(repr=False)


@dataclass(slots=True)
class SelectorRuleDB(RuleDB):
    """Top-level container for feature-based selector databases."""

    selectors: Tuple[SelectorRule, ...]
    defaults: Mapping[str, Any]
    feature_schema: Mapping[str, Any] | None
    guardrails: Tuple[Mapping[str, Any], ...]
    priors: Tuple[Mapping[str, Any], ...]
    repairs: Tuple[Mapping[str, Any], ...]


@dataclass(slots=True)
class AnchorHit:
    """Record of an anchor SMARTS match during normalization."""

    anchor_smarts: str
    reactant_index: int
    atom_indices: Tuple[int, ...]


@dataclass(slots=True)
class EssentialCoreNormalizationResult:
    """Artifacts produced by the essential-core normalization pipeline."""

    reactant_smiles: Tuple[str, ...]
    kept_reactants: Tuple[str, ...]
    dropped_reactants: Tuple[str, ...]
    sanitized_mols: Tuple[Mol, ...]
    masked_mols: Tuple[Mol, ...]
    masked_source_indices: Tuple[int, ...]
    anchor_hits: Tuple[AnchorHit, ...]
    smarts_bag: Tuple[str, ...]


@dataclass(slots=True)
class MatchResult:
    """Final outcome of running the matcher for a reaction."""

    reaction_type: str
    match_type: Literal["scheme", "default", "global_default", "selector", "rule_default"]
    entry_id: str
    entry_name: str | None
    priority: int
    conditions: Mapping[str, Any]
    trace: Dict[str, Any]
    reagents: Any = None  # New field for standardized reagents array

    def to_json_dict(self) -> Dict[str, Any]:
        """Return a JSON-serialisable dictionary."""

        result = {
            "reaction_type": self.reaction_type,
            "match_type": self.match_type,
            "entry_id": self.entry_id,
            "entry_name": self.entry_name,
            "priority": self.priority,
            "conditions": self._to_plain(self.conditions),
            "trace": self.trace,
        }
        
        # Include reagents if present (new standardized format)
        if self.reagents is not None:
            result["reagents"] = self.reagents
        
        return result

    @staticmethod
    def _to_plain(value: Mapping[str, Any]) -> Dict[str, Any]:
        return dict(value)

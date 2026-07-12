"""Typed public contracts for condition substance resolution."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, Optional, Tuple


@dataclass(frozen=True)
class RoleAssignment:
    role_id: str
    family_id: Optional[str]
    tag: Optional[str] = None
    confidence: float = 1.0
    evidence: str = "curated_registry"


@dataclass(frozen=True)
class Substance:
    substance_id: str
    canonical_name: str
    cas: Optional[str]
    smiles: Optional[str]
    aliases: Tuple[str, ...] = ()
    roles: Tuple[RoleAssignment, ...] = ()
    properties: Dict[str, Any] = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        if self.properties is None:
            object.__setattr__(self, "properties", {})

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ResolutionResult:
    query: str
    status: str
    substance: Optional[Substance] = None
    match_kind: Optional[str] = None
    candidates: Tuple[str, ...] = ()
    warnings: Tuple[str, ...] = ()

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

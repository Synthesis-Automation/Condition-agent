"""Reusable, evidence-based heuristic for precursor realism.

The scorer deliberately consumes normalized evidence flags rather than opening
particular databases.  Callers may therefore obtain the three evidence signals
from SQLite indexes, in-memory registries, APIs, or future evidence providers
without changing the scoring contract.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable, Mapping

from rdkit import Chem

from .molecule_index import MoleculeIdentity, molecule_identity


PRECURSOR_REALISM_SCHEMA_VERSION = "1.0"
PRECURSOR_REALISM_DEFINITION_PATH = (
    Path(__file__).with_name("definitions") / "precursor_realism.v1.json"
)

_EVIDENCE_TIERS = (
    "buyable_registry_literature",
    "buyable_corroborated",
    "buyable_only",
    "registry_literature",
    "registry_only",
    "literature_only",
    "unsupported",
)


@dataclass(frozen=True)
class PrecursorEvidence:
    """Exact normalized structure matches from the three trusted sources."""

    buyable: bool
    in_compound_registry: bool
    in_literature: bool

    def __post_init__(self) -> None:
        for field, value in asdict(self).items():
            if not isinstance(value, bool):
                raise TypeError(f"{field} must be boolean")


@dataclass(frozen=True)
class PrecursorRealismPolicy:
    """Validated immutable scoring policy loaded from a versioned definition."""

    definition_id: str
    schema_version: str
    maximum_smallness_da: float
    no_smallness_da: float
    tier_scores: tuple[tuple[str, float, float], ...]

    def tier(self, tier_id: str) -> tuple[float, float]:
        """Return the base score and maximum MW penalty for an evidence tier."""

        tiers = {
            name: (base_score, maximum_penalty)
            for name, base_score, maximum_penalty in self.tier_scores
        }
        try:
            return tiers[tier_id]
        except KeyError as error:
            raise ValueError(f"unsupported precursor evidence tier: {tier_id}") from error


@dataclass(frozen=True)
class PrecursorRealismAssessment:
    """Auditable precursor-realism result with every scoring contribution."""

    canonical_smiles: str
    inchi_key: str | None
    molecular_weight: float
    evidence: PrecursorEvidence
    evidence_tier: str
    base_score: float
    molecular_weight_smallness: float
    molecular_weight_penalty: float
    score: float
    definition_id: str
    schema_version: str = PRECURSOR_REALISM_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible score trace."""

        return asdict(self)


def _finite_score(value: object, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{field} must be a finite number between zero and one")
    result = float(value)
    if not math.isfinite(result) or not 0.0 <= result <= 1.0:
        raise ValueError(f"{field} must be a finite number between zero and one")
    return result


def _positive_float(value: object, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{field} must be a positive finite number")
    result = float(value)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"{field} must be a positive finite number")
    return result


def validate_precursor_realism_definition(value: Mapping[str, Any]) -> None:
    """Validate a precursor-realism definition without mutating global state."""

    if value.get("definition_id") != "precursor_realism.v1":
        raise ValueError("unexpected precursor-realism definition ID")
    if value.get("schema_version") != PRECURSOR_REALISM_SCHEMA_VERSION:
        raise ValueError("unsupported precursor-realism schema")
    molecular_weight = value.get("molecular_weight")
    tiers = value.get("evidence_tiers")
    if not isinstance(molecular_weight, Mapping) or not isinstance(tiers, Mapping):
        raise ValueError("precursor-realism definition requires MW and tier rules")
    maximum_smallness = _positive_float(
        molecular_weight.get("maximum_smallness_da"),
        "maximum-smallness molecular weight",
    )
    no_smallness = _positive_float(
        molecular_weight.get("no_smallness_da"),
        "no-smallness molecular weight",
    )
    if no_smallness <= maximum_smallness:
        raise ValueError("no-smallness MW must exceed maximum-smallness MW")
    if set(tiers) != set(_EVIDENCE_TIERS):
        raise ValueError("precursor-realism evidence tiers are incomplete")
    for tier_id in _EVIDENCE_TIERS:
        tier = tiers[tier_id]
        if not isinstance(tier, Mapping):
            raise ValueError(f"precursor-realism tier {tier_id} must be an object")
        _finite_score(tier.get("base_score"), f"{tier_id} base score")
        _finite_score(
            tier.get("maximum_molecular_weight_penalty"),
            f"{tier_id} maximum molecular-weight penalty",
        )


@lru_cache(maxsize=1)
def load_precursor_realism_policy() -> PrecursorRealismPolicy:
    """Load and validate the canonical precursor-realism policy."""

    value = json.loads(
        PRECURSOR_REALISM_DEFINITION_PATH.read_text(encoding="utf-8")
    )
    validate_precursor_realism_definition(value)
    molecular_weight = value["molecular_weight"]
    tiers = value["evidence_tiers"]
    return PrecursorRealismPolicy(
        definition_id=str(value["definition_id"]),
        schema_version=str(value["schema_version"]),
        maximum_smallness_da=float(molecular_weight["maximum_smallness_da"]),
        no_smallness_da=float(molecular_weight["no_smallness_da"]),
        tier_scores=tuple(
            (
                tier_id,
                float(tiers[tier_id]["base_score"]),
                float(tiers[tier_id]["maximum_molecular_weight_penalty"]),
            )
            for tier_id in _EVIDENCE_TIERS
        ),
    )


def _evidence_tier(evidence: PrecursorEvidence) -> str:
    if evidence.buyable:
        corroboration = int(evidence.in_compound_registry) + int(
            evidence.in_literature
        )
        if corroboration == 2:
            return "buyable_registry_literature"
        if corroboration == 1:
            return "buyable_corroborated"
        return "buyable_only"
    if evidence.in_compound_registry and evidence.in_literature:
        return "registry_literature"
    if evidence.in_compound_registry:
        return "registry_only"
    if evidence.in_literature:
        return "literature_only"
    return "unsupported"


def _smallness(molecular_weight: float, policy: PrecursorRealismPolicy) -> float:
    span = policy.no_smallness_da - policy.maximum_smallness_da
    value = (policy.no_smallness_da - molecular_weight) / span
    return min(1.0, max(0.0, value))


def assess_precursor_realism(
    precursor: MoleculeIdentity | str,
    evidence: PrecursorEvidence,
    *,
    policy: PrecursorRealismPolicy | None = None,
) -> PrecursorRealismAssessment:
    """Score one valid precursor from exact source matches and molecular weight.

    The result is a ranking heuristic, not a probability of chemical existence
    or experimental success. Invalid molecular structures are rejected because
    graph validity belongs upstream of evidence-based ranking.
    """

    identity = (
        molecule_identity(precursor) if isinstance(precursor, str) else precursor
    )
    if identity is None:
        raise ValueError("precursor must be a valid molecular structure")
    resolved = policy or load_precursor_realism_policy()
    tier_id = _evidence_tier(evidence)
    base_score, maximum_penalty = resolved.tier(tier_id)
    smallness = _smallness(identity.molecular_weight, resolved)
    penalty = maximum_penalty * smallness
    score = min(1.0, max(0.0, base_score - penalty))
    return PrecursorRealismAssessment(
        canonical_smiles=identity.canonical_smiles,
        inchi_key=identity.inchi_key,
        molecular_weight=round(identity.molecular_weight, 6),
        evidence=evidence,
        evidence_tier=tier_id,
        base_score=round(base_score, 6),
        molecular_weight_smallness=round(smallness, 6),
        molecular_weight_penalty=round(penalty, 6),
        score=round(score, 6),
        definition_id=resolved.definition_id,
        schema_version=resolved.schema_version,
    )


def aggregate_precursor_realism(
    assessments: tuple[PrecursorRealismAssessment, ...],
) -> float:
    """Return the weakest required precursor score for one reaction candidate."""

    if not assessments:
        raise ValueError("at least one precursor assessment is required")
    return min(assessment.score for assessment in assessments)


def assess_precursor_components(
    precursor_smiles: str,
    evidence_provider: Callable[[MoleculeIdentity], PrecursorEvidence],
    *,
    policy: PrecursorRealismPolicy | None = None,
) -> tuple[PrecursorRealismAssessment, ...]:
    """Assess every connected precursor component in canonical order.

    Retrosynthesis candidates encode required precursor molecules as a
    dot-separated SMILES string. Evidence collection stays caller-owned: the
    provider receives each strict molecular identity and returns exact match
    flags for the configured databases.
    """

    molecule = Chem.MolFromSmiles(str(precursor_smiles or "").strip())
    if molecule is None:
        raise ValueError("precursors must be valid molecular structures")
    component_smiles = sorted(
        Chem.MolToSmiles(component, canonical=True, isomericSmiles=True)
        for component in Chem.GetMolFrags(
            molecule,
            asMols=True,
            sanitizeFrags=True,
        )
    )
    if not component_smiles:
        raise ValueError("at least one precursor component is required")
    assessments = []
    for smiles in component_smiles:
        identity = molecule_identity(smiles)
        if identity is None:
            raise ValueError("precursors must be valid molecular structures")
        assessments.append(
            assess_precursor_realism(
                identity,
                evidence_provider(identity),
                policy=policy,
            )
        )
    return tuple(assessments)


__all__ = [
    "PRECURSOR_REALISM_DEFINITION_PATH",
    "PRECURSOR_REALISM_SCHEMA_VERSION",
    "PrecursorEvidence",
    "PrecursorRealismAssessment",
    "PrecursorRealismPolicy",
    "aggregate_precursor_realism",
    "assess_precursor_components",
    "assess_precursor_realism",
    "load_precursor_realism_policy",
    "validate_precursor_realism_definition",
]

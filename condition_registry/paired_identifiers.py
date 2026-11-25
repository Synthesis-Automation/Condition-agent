"""Independent registry evidence for source-paired names and molecular graphs."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Mapping

from .api import resolve_identifier
from .models import ResolutionResult
from .normalization import normalize_smiles


@dataclass(frozen=True)
class PairedIdentifierEvidence:
    """Independent candidate sets and the identity corroborated by both channels."""

    status: str
    name_candidates: tuple[str, ...]
    structure_candidates: tuple[str, ...]
    corroborated_substance_id: str | None = None
    definition_version: str = "paired_identifiers.v1"


def _candidates(result: ResolutionResult) -> tuple[str, ...]:
    return (result.substance.substance_id,) if result.substance else result.candidates


@lru_cache(maxsize=8192)
def _compare_pair(name: str, smiles: str) -> PairedIdentifierEvidence:
    by_name = resolve_identifier(name, identifier_type="name")
    by_structure = resolve_identifier(smiles, identifier_type="smiles")
    names, structures = _candidates(by_name), _candidates(by_structure)
    shared = set(names).intersection(structures)
    corroborated = None
    if names and structures:
        if len(shared) == 1 and (by_name.substance or by_structure.substance):
            corroborated = next(iter(shared))
            status = "corroborated"
        else:
            status = "ambiguous" if shared else "conflict"
    elif by_structure.status == "invalid_identifier":
        status = "invalid_structure"
    elif by_name.substance and by_name.substance.smiles and (
        normalize_smiles(by_name.substance.smiles) != normalize_smiles(smiles)
    ):
        status = "conflict"
    else:
        status = "uncorroborated"
    return PairedIdentifierEvidence(status, names, structures, corroborated)


def paired_identifier_evidence(
    provenance: Mapping[str, Any] | None,
) -> PairedIdentifierEvidence | None:
    """Compare co-reported identifiers without pairing separate claims by position."""
    identifiers = (provenance or {}).get("source_identifiers") or ()
    names = {i["value"] for i in identifiers if i.get("identifier_type") == "name"}
    structures = {i["value"] for i in identifiers if i.get("identifier_type") == "smiles"}
    if not names or not structures:
        return None
    if len(names) != 1 or len(structures) != 1:
        return PairedIdentifierEvidence("ambiguous", (), ())
    return _compare_pair(next(iter(names)), next(iter(structures)))

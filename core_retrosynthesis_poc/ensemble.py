"""Baseline-first composition of exact and core-derived retrosynthesis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Literal, Sequence

from retrosynthesis_poc import RetrosynthesisLibrary
from retrosynthesis_poc import disconnect_target as baseline_disconnect

from .models import CoreTemplateLibrary
from .search import disconnect_target as core_disconnect


@dataclass(frozen=True)
class EnsembleCandidate:
    """One baseline or core-fallback proposal with method provenance."""

    source_method: Literal["rdchiral_baseline", "core_l1_context"]
    candidate: Any

    @property
    def proposed_reaction_smiles(self) -> str:
        """Return the proposed precursor-to-target reaction SMILES."""

        return str(self.candidate.proposed_reaction_smiles)

    def to_dict(self) -> dict[str, Any]:
        """Return candidate fields plus explicit source-method provenance."""

        return {
            "source_method": self.source_method,
            "candidate": self.candidate.to_dict(),
        }


def merge_baseline_first(
    baseline: Sequence[Any],
    fallback: Sequence[Any],
    top_k: int,
) -> tuple[Any, ...]:
    """Keep baseline order and fill unused ranks with unique fallback proposals."""

    if top_k < 1:
        raise ValueError("top-k must be positive")
    selected = []
    seen = set()
    for candidate in (*baseline, *fallback):
        if candidate.precursor_smiles in seen:
            continue
        selected.append(candidate)
        seen.add(candidate.precursor_smiles)
        if len(selected) >= top_k:
            break
    return tuple(selected)


def disconnect_ensemble(
    target_smiles: str,
    baseline_library: RetrosynthesisLibrary,
    core_library: CoreTemplateLibrary,
    *,
    allowed_bonds: Iterable[str] = ("C-N", "C-O", "C-S"),
    top_k: int = 20,
    max_candidates_to_validate: int = 50,
    validate_forward: bool = True,
) -> tuple[EnsembleCandidate, ...]:
    """Return exact baseline candidates followed by contextual L1 fallbacks."""

    bonds = tuple(allowed_bonds)
    baseline = baseline_disconnect(
        target_smiles,
        baseline_library,
        allowed_bonds=bonds,
        top_k=top_k,
        validate_forward=validate_forward,
    )
    fallback = core_disconnect(
        target_smiles,
        core_library,
        allowed_bonds=bonds,
        levels=("L1",),
        top_k=top_k,
        max_candidates_to_validate=max_candidates_to_validate,
        validate_forward=validate_forward,
        use_context=True,
    )
    baseline_ids = {id(candidate) for candidate in baseline}
    return tuple(
        EnsembleCandidate(
            source_method=(
                "rdchiral_baseline"
                if id(candidate) in baseline_ids
                else "core_l1_context"
            ),
            candidate=candidate,
        )
        for candidate in merge_baseline_first(baseline, fallback, top_k)
    )


__all__ = ["EnsembleCandidate", "disconnect_ensemble", "merge_baseline_first"]

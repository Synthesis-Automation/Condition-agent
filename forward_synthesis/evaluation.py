"""Leakage-guarded product-hidden replay evaluation for forward prediction."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Iterable, Mapping, Optional, Tuple

from reactive_taxonomy import canonical_molecule_collection

from .models import ForwardOperatorLibrary
from .prediction import predict_products


FORWARD_REPLAY_EVALUATION_SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class ForwardReplayCaseResult:
    """Target-hidden outcome for one held-out observed reaction."""

    reaction_id: str
    reference_id: str
    starting_materials: str
    expected_product: str
    expected_product_rank: Optional[int]
    exact_recovered: bool
    predicted_product_count: int
    status: str


@dataclass(frozen=True)
class ForwardReplayEvaluationReport:
    """Aggregate exact-product retrieval metrics without yield claims."""

    case_count: int
    exact_top_1_count: int
    exact_top_5_count: int
    exact_top_10_count: int
    exact_any_count: int
    no_product_count: int
    reciprocal_rank_mean: float
    cases: Tuple[ForwardReplayCaseResult, ...]
    split_requirement: str = "reference_disjoint"
    schema_version: str = FORWARD_REPLAY_EVALUATION_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _reaction_sides(reaction_smiles: str) -> tuple[str, str]:
    if reaction_smiles.count(">>") == 1:
        left, right = reaction_smiles.split(">>")
    else:
        parts = reaction_smiles.split(">")
        if len(parts) != 3:
            raise ValueError("evaluation reaction must have two or three parts")
        left, right = parts[0], parts[2]
    if not left or not right:
        raise ValueError("evaluation reaction requires reactants and product")
    return left, right


def _library_reference_ids(library: ForwardOperatorLibrary) -> set[str]:
    return {
        precedent.reference_id
        for operator in library.operators
        for precedent in operator.precedents
        if precedent.reference_id
    }


def evaluate_product_hidden_replay(
    rows: Iterable[Mapping[str, Any]],
    library: ForwardOperatorLibrary,
    *,
    top_k: int = 20,
    require_reference_disjoint: bool = True,
) -> ForwardReplayEvaluationReport:
    """Hide each observed product and measure exact product recovery.

    Enumerated alternatives are not labelled as false products.  The metric is
    therefore recall/rank of the observed product, not selectivity accuracy.
    """

    library_references = _library_reference_ids(library)
    cases = []
    for index, row in enumerate(rows):
        reaction_id = str(row.get("reaction_id") or f"row-{index}")
        reference_id = str(
            row.get("reference_id")
            or (row.get("reference_identity") or {}).get("reference_id")
            or ""
        )
        if (
            require_reference_disjoint
            and reference_id
            and reference_id in library_references
        ):
            raise ValueError(
                f"evaluation reference overlaps operator library: {reference_id}"
            )
        left, right = _reaction_sides(str(row.get("reaction_smiles") or ""))
        expected = canonical_molecule_collection(right)
        if expected is None or "." in expected:
            cases.append(
                ForwardReplayCaseResult(
                    reaction_id=reaction_id,
                    reference_id=reference_id,
                    starting_materials=left,
                    expected_product=right,
                    expected_product_rank=None,
                    exact_recovered=False,
                    predicted_product_count=0,
                    status="invalid_or_multi_product_expected",
                )
            )
            continue
        prediction = predict_products(left, library, top_k=top_k)
        ranks = tuple(
            candidate.rank
            for candidate in prediction.candidates
            if candidate.product_smiles == expected
        )
        rank = ranks[0] if ranks else None
        cases.append(
            ForwardReplayCaseResult(
                reaction_id=reaction_id,
                reference_id=reference_id,
                starting_materials=prediction.canonical_starting_materials or left,
                expected_product=expected,
                expected_product_rank=rank,
                exact_recovered=rank is not None,
                predicted_product_count=len(prediction.candidates),
                status=("recovered" if rank is not None else prediction.status),
            )
        )
    values = tuple(cases)
    ranks = tuple(
        item.expected_product_rank
        for item in values
        if item.expected_product_rank is not None
    )
    return ForwardReplayEvaluationReport(
        case_count=len(values),
        exact_top_1_count=sum(rank <= 1 for rank in ranks),
        exact_top_5_count=sum(rank <= 5 for rank in ranks),
        exact_top_10_count=sum(rank <= 10 for rank in ranks),
        exact_any_count=len(ranks),
        no_product_count=sum(item.predicted_product_count == 0 for item in values),
        reciprocal_rank_mean=round(
            sum(1.0 / rank for rank in ranks) / len(values),
            8,
        )
        if values
        else 0.0,
        cases=values,
    )


__all__ = [
    "FORWARD_REPLAY_EVALUATION_SCHEMA_VERSION",
    "ForwardReplayCaseResult",
    "ForwardReplayEvaluationReport",
    "evaluate_product_hidden_replay",
]

"""Nested chemical-space expansion from observed two-step routes.

This proof of concept treats a precedent route as executable structural evidence:
each observed step is compiled into source-round-tripped forward operators, and
only products that pass the normal forward/reverse graph checks are propagated.
The expansion levels are cumulative and make increasing extrapolation explicit.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Literal, Mapping, Optional

from cas_tools import open_stock_lookup
from reactive_taxonomy import canonical_molecule_collection

from forward_synthesis import build_forward_library, predict_products

from .chemistry import digest, split_reaction_smiles
from .generic_library import build_generic_library


PRECEDENT_ROUTE_EXPANSION_SCHEMA_VERSION = "1.1"
SUPPORTED_PRECEDENT_ROUTE_EXPANSION_DEFINITION_SCHEMAS = frozenset(
    {"1.0", "1.1"}
)
PRECEDENT_ROUTE_EXPANSION_ALGORITHM_VERSION = (
    "two_step_precedent_route_expansion.v1.1"
)
DEFAULT_PRECEDENT_ROUTE_EXPANSION_DEFINITION = (
    Path(__file__).resolve().parent
    / "definitions"
    / "two_step_precedent_route_expansion_poc.v1.json"
)

ExpansionLevel = Literal["R0", "R1", "R2"]
ExpansionInputSource = Literal["observed", "curated_stock", "declared"]
EXPANSION_LEVELS: tuple[ExpansionLevel, ...] = ("R0", "R1", "R2")
_LEVEL_RANK = {level: index for index, level in enumerate(EXPANSION_LEVELS)}
_OPERATOR_LEVELS: dict[ExpansionLevel, tuple[str, ...]] = {
    "R0": ("L2",),
    "R1": ("L2",),
    "R2": ("L2", "L1"),
}
_LEVEL_CLAIMS: dict[ExpansionLevel, str] = {
    "R0": "exact observed-route replay",
    "R1": "route expansion under exact-context L2 operators",
    "R2": "cumulative expansion allowing local-context L1 operators",
}


@dataclass(frozen=True)
class PrecedentRouteExpansionInput:
    """One starting-material or partner set admitted from a declared level."""

    smiles: str
    minimum_level: ExpansionLevel
    source_kind: ExpansionInputSource = "declared"
    stock_components: tuple[str, ...] = ()
    label: str = ""

    def __post_init__(self) -> None:
        if self.minimum_level not in EXPANSION_LEVELS:
            raise ValueError(f"unsupported expansion level: {self.minimum_level}")
        if self.source_kind not in {"observed", "curated_stock", "declared"}:
            raise ValueError(f"unsupported input source kind: {self.source_kind}")
        if self.source_kind == "curated_stock" and not self.stock_components:
            raise ValueError("curated-stock inputs require stock components")
        if any(not component.strip() for component in self.stock_components):
            raise ValueError("stock components must be nonempty molecules")


@dataclass(frozen=True)
class TwoStepPrecedentRouteDefinition:
    """One linear two-step route plus bounded building-block choices."""

    route_id: str
    reference_id: str
    first_step_id: str
    first_reaction_smiles: str
    second_step_id: str
    second_reaction_smiles: str
    intermediate_smiles: str
    target_smiles: str
    first_step_inputs: tuple[PrecedentRouteExpansionInput, ...]
    second_step_partners: tuple[PrecedentRouteExpansionInput, ...]
    source_dataset_id: str = ""
    source_route_id: str = ""
    patent_id: str = ""
    source_tree_id: str = ""
    source_route_core_id: str = ""
    source_lineage_id: str = ""
    first_source_reaction_id: str = ""
    second_source_reaction_id: str = ""

    def __post_init__(self) -> None:
        if not self.route_id or not self.reference_id:
            raise ValueError("route and reference IDs are required")
        if self.first_step_id == self.second_step_id:
            raise ValueError("two-step route step IDs must be distinct")
        exact_first = tuple(
            item for item in self.first_step_inputs if item.minimum_level == "R0"
        )
        if not exact_first:
            raise ValueError("the first step requires an R0 input")
        exact_partners = tuple(
            item
            for item in self.second_step_partners
            if item.minimum_level == "R0"
        )
        if not exact_partners:
            raise ValueError("the second step requires an R0 partner")
        if any(item.source_kind != "observed" for item in exact_first):
            raise ValueError("R0 first-step inputs must be observed")
        if any(item.source_kind != "observed" for item in exact_partners):
            raise ValueError("R0 second-step partners must be observed")
        first = split_reaction_smiles(self.first_reaction_smiles)
        second = split_reaction_smiles(self.second_reaction_smiles)
        intermediate = canonical_molecule_collection(self.intermediate_smiles)
        target = canonical_molecule_collection(self.target_smiles)
        if first is None or second is None or intermediate is None or target is None:
            raise ValueError("route reactions and declared products must parse")
        if canonical_molecule_collection(first[1]) != intermediate:
            raise ValueError("first-step product contradicts declared intermediate")
        if canonical_molecule_collection(second[1]) != target:
            raise ValueError("second-step product contradicts declared target")
        second_components = set(
            (canonical_molecule_collection(second[0]) or "").split(".")
        )
        if not set(intermediate.split(".")).issubset(second_components):
            raise ValueError("declared intermediate is absent from the second step")


@dataclass(frozen=True)
class PrecedentRouteExpansionDefinition:
    """Versioned panel of two-step precedent routes."""

    definition_id: str
    routes: tuple[TwoStepPrecedentRouteDefinition, ...]
    schema_version: str

    def __post_init__(self) -> None:
        if self.schema_version not in (
            SUPPORTED_PRECEDENT_ROUTE_EXPANSION_DEFINITION_SCHEMAS
        ):
            raise ValueError("unsupported precedent-route expansion schema")
        route_ids = tuple(route.route_id for route in self.routes)
        if not self.routes or len(route_ids) != len(set(route_ids)):
            raise ValueError("precedent-route IDs must be nonempty and unique")


@dataclass(frozen=True)
class PrecedentRouteExpansionPathway:
    """One fully graph-validated two-step path to a generated product."""

    pathway_id: str
    first_discovered_level: ExpansionLevel
    product_smiles: str
    intermediate_smiles: str
    first_step_starting_materials: str
    second_step_partner_smiles: str
    second_step_starting_materials: str
    first_operator_id: str
    first_forward_operator_id: str
    first_abstraction_level: str
    first_reverse_round_trip: bool
    first_operator_edit_agreement: bool
    second_operator_id: str
    second_forward_operator_id: str
    second_abstraction_level: str
    second_reverse_round_trip: bool
    second_operator_edit_agreement: bool
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible pathway record."""

        value = asdict(self)
        value["warnings"] = list(self.warnings)
        return value


@dataclass(frozen=True)
class PrecedentRouteStockComponentEvidence:
    """Exact supplier-stock evidence for one declared building block."""

    smiles: str
    available: bool
    canonical_smiles: str = ""
    offer_count: int = 0
    source_records: tuple[dict[str, str], ...] = ()

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible exact-stock assessment."""

        return {
            "smiles": self.smiles,
            "available": self.available,
            "canonical_smiles": self.canonical_smiles,
            "offer_count": self.offer_count,
            "source_records": [dict(record) for record in self.source_records],
        }


@dataclass(frozen=True)
class PrecedentRouteInputEvidence:
    """Provenance and optional stock verification for one expansion input."""

    role: Literal["first_step", "second_step_partner"]
    input_smiles: str
    minimum_level: ExpansionLevel
    source_kind: ExpansionInputSource
    label: str
    stock_evidence_complete: Optional[bool]
    stock_components: tuple[PrecedentRouteStockComponentEvidence, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible input evidence record."""

        return {
            "role": self.role,
            "input_smiles": self.input_smiles,
            "minimum_level": self.minimum_level,
            "source_kind": self.source_kind,
            "label": self.label,
            "stock_evidence_complete": self.stock_evidence_complete,
            "stock_components": [
                component.to_dict() for component in self.stock_components
            ],
        }


@dataclass(frozen=True)
class PrecedentRouteExpansionLevelResult:
    """Cumulative products and pathways available at one expansion level."""

    level: ExpansionLevel
    claim: str
    allowed_operator_levels: tuple[str, ...]
    first_step_input_count: int
    second_step_partner_count: int
    intermediate_count: int
    product_count: int
    new_product_count: int
    product_smiles: tuple[str, ...]
    pathways: tuple[PrecedentRouteExpansionPathway, ...]
    input_evidence: tuple[PrecedentRouteInputEvidence, ...]
    rejection_counts: dict[str, int]

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible cumulative level result."""

        return {
            "level": self.level,
            "claim": self.claim,
            "allowed_operator_levels": list(self.allowed_operator_levels),
            "first_step_input_count": self.first_step_input_count,
            "second_step_partner_count": self.second_step_partner_count,
            "intermediate_count": self.intermediate_count,
            "product_count": self.product_count,
            "new_product_count": self.new_product_count,
            "product_smiles": list(self.product_smiles),
            "pathways": [pathway.to_dict() for pathway in self.pathways],
            "input_evidence": [
                evidence.to_dict() for evidence in self.input_evidence
            ],
            "rejection_counts": dict(sorted(self.rejection_counts.items())),
        }


@dataclass(frozen=True)
class PrecedentRouteExpansionResult:
    """Auditable nested chemical space generated from one observed route."""

    expansion_id: str
    route_id: str
    reference_id: str
    source_dataset_id: str
    source_route_id: str
    patent_id: str
    source_tree_id: str
    source_route_core_id: str
    source_lineage_id: str
    first_source_reaction_id: str
    second_source_reaction_id: str
    intermediate_smiles: str
    target_smiles: str
    exact_replay_valid: bool
    compiled_template_count: int
    admitted_forward_operator_count: int
    first_step_operator_ids: tuple[str, ...]
    second_step_operator_ids: tuple[str, ...]
    levels: tuple[PrecedentRouteExpansionLevelResult, ...]
    warnings: tuple[str, ...]
    definition_id: str
    algorithm_version: str = PRECEDENT_ROUTE_EXPANSION_ALGORITHM_VERSION
    schema_version: str = PRECEDENT_ROUTE_EXPANSION_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible route-expansion report."""

        return {
            "expansion_id": self.expansion_id,
            "route_id": self.route_id,
            "reference_id": self.reference_id,
            "source_dataset_id": self.source_dataset_id,
            "source_route_id": self.source_route_id,
            "patent_id": self.patent_id,
            "source_tree_id": self.source_tree_id,
            "source_route_core_id": self.source_route_core_id,
            "source_lineage_id": self.source_lineage_id,
            "first_source_reaction_id": self.first_source_reaction_id,
            "second_source_reaction_id": self.second_source_reaction_id,
            "intermediate_smiles": self.intermediate_smiles,
            "target_smiles": self.target_smiles,
            "exact_replay_valid": self.exact_replay_valid,
            "compiled_template_count": self.compiled_template_count,
            "admitted_forward_operator_count": (
                self.admitted_forward_operator_count
            ),
            "first_step_operator_ids": list(self.first_step_operator_ids),
            "second_step_operator_ids": list(self.second_step_operator_ids),
            "levels": [level.to_dict() for level in self.levels],
            "warnings": list(self.warnings),
            "definition_id": self.definition_id,
            "algorithm_version": self.algorithm_version,
            "schema_version": self.schema_version,
        }


def _input(value: Mapping[str, Any]) -> PrecedentRouteExpansionInput:
    minimum_level = str(value["minimum_level"])
    return PrecedentRouteExpansionInput(
        smiles=str(value["smiles"]),
        minimum_level=minimum_level,  # type: ignore[arg-type]
        source_kind=str(
            value.get("source_kind")
            or ("observed" if minimum_level == "R0" else "declared")
        ),  # type: ignore[arg-type]
        stock_components=tuple(
            str(item) for item in value.get("stock_components") or ()
        ),
        label=str(value.get("label") or ""),
    )


def load_precedent_route_expansion_definition(
    path: str | Path = DEFAULT_PRECEDENT_ROUTE_EXPANSION_DEFINITION,
) -> PrecedentRouteExpansionDefinition:
    """Load and validate a versioned two-step expansion panel."""

    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, Mapping):
        raise ValueError("precedent-route definition must be an object")
    routes = []
    for raw in value.get("routes") or ():
        if not isinstance(raw, Mapping):
            raise ValueError("precedent-route entries must be objects")
        first = raw.get("first_step") or {}
        second = raw.get("second_step") or {}
        routes.append(
            TwoStepPrecedentRouteDefinition(
                route_id=str(raw["route_id"]),
                reference_id=str(raw["reference_id"]),
                first_step_id=str(first["step_id"]),
                first_reaction_smiles=str(first["reaction_smiles"]),
                second_step_id=str(second["step_id"]),
                second_reaction_smiles=str(second["reaction_smiles"]),
                intermediate_smiles=str(raw["intermediate_smiles"]),
                target_smiles=str(raw["target_smiles"]),
                first_step_inputs=tuple(
                    _input(item) for item in first.get("input_sets") or ()
                ),
                second_step_partners=tuple(
                    _input(item) for item in second.get("partner_sets") or ()
                ),
                source_dataset_id=str(raw.get("source_dataset_id") or ""),
                source_route_id=str(raw.get("source_route_id") or ""),
                patent_id=str(raw.get("patent_id") or ""),
                source_tree_id=str(raw.get("source_tree_id") or ""),
                source_route_core_id=str(
                    raw.get("source_route_core_id") or ""
                ),
                source_lineage_id=str(raw.get("source_lineage_id") or ""),
                first_source_reaction_id=str(
                    first.get("source_reaction_id") or ""
                ),
                second_source_reaction_id=str(
                    second.get("source_reaction_id") or ""
                ),
            )
        )
    return PrecedentRouteExpansionDefinition(
        definition_id=str(value["definition_id"]),
        routes=tuple(routes),
        schema_version=str(value["schema_version"]),
    )


def _admitted_inputs(
    values: Iterable[PrecedentRouteExpansionInput],
    level: ExpansionLevel,
) -> tuple[PrecedentRouteExpansionInput, ...]:
    return tuple(
        sorted(
            (
                item
                for item in values
                if _LEVEL_RANK[item.minimum_level] <= _LEVEL_RANK[level]
            ),
            key=lambda item: (item.smiles, item.minimum_level),
        )
    )


def _input_evidence(
    item: PrecedentRouteExpansionInput,
    role: Literal["first_step", "second_step_partner"],
    stock_lookup: Any | None,
) -> PrecedentRouteInputEvidence:
    if item.source_kind == "curated_stock" and stock_lookup is None:
        raise ValueError("curated-stock expansion inputs require a stock lookup")
    component_evidence = []
    for component in item.stock_components:
        match = stock_lookup.lookup(component) if stock_lookup is not None else None
        component_evidence.append(
            PrecedentRouteStockComponentEvidence(
                smiles=component,
                available=match is not None,
                canonical_smiles=(match.canonical_smiles if match else ""),
                offer_count=(int(match.occurrence_count) if match else 0),
                source_records=(
                    tuple(dict(record) for record in match.source_records)
                    if match
                    else ()
                ),
            )
        )
    complete = (
        all(component.available for component in component_evidence)
        if component_evidence
        else None
    )
    return PrecedentRouteInputEvidence(
        role=role,
        input_smiles=item.smiles,
        minimum_level=item.minimum_level,
        source_kind=item.source_kind,
        label=item.label,
        stock_evidence_complete=complete,
        stock_components=tuple(component_evidence),
    )


def _input_is_eligible(
    item: PrecedentRouteExpansionInput,
    evidence: PrecedentRouteInputEvidence,
) -> bool:
    return not (
        item.source_kind == "curated_stock"
        and evidence.stock_evidence_complete is not True
    )


def _all_query_components_participate(candidate: Any, query: str) -> bool:
    canonical_query = canonical_molecule_collection(query)
    canonical_participating = canonical_molecule_collection(
        candidate.participating_precursor_smiles
    )
    return bool(
        canonical_query is not None
        and canonical_participating is not None
        and canonical_query == canonical_participating
        and not candidate.uses_virtual_copies
    )


def _combine_intermediate_and_partner(
    intermediate: str,
    partner: str,
) -> str:
    return f"{intermediate}.{partner}" if partner else intermediate


def _reaction_without_reagents(reaction_smiles: str) -> str:
    if reaction_smiles.count(">>") == 1:
        return reaction_smiles
    fields = reaction_smiles.split(">")
    if len(fields) == 3 and fields[0] and fields[2]:
        return f"{fields[0]}>>{fields[2]}"
    return ""


def _file_sha256(path: str | Path) -> str:
    value = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def validate_precedent_route_sources(
    definition: PrecedentRouteExpansionDefinition,
    route_core_source: str | Path,
) -> dict[str, Any]:
    """Verify declared observed segments against a route-core artifact."""

    from .route_core_conversion import iter_route_core_projections

    expected = {
        route.source_route_id: route
        for route in definition.routes
        if route.source_route_id
    }
    found: dict[str, Any] = {}
    for projection in iter_route_core_projections(route_core_source):
        if projection.source_route_id in expected:
            found[str(projection.source_route_id)] = projection
        if len(found) == len(expected):
            break
    route_results = []
    for route_id, route in sorted(expected.items()):
        issues = []
        projection = found.get(route_id)
        if projection is None:
            issues.append("source_route_missing")
        else:
            comparisons = {
                "patent_id_mismatch": (projection.patent_id, route.patent_id),
                "source_tree_id_mismatch": (
                    projection.source_tree_id,
                    route.source_tree_id,
                ),
                "source_route_core_id_mismatch": (
                    projection.route_core_id,
                    route.source_route_core_id,
                ),
            }
            for issue, (observed, declared) in comparisons.items():
                if declared and observed != declared:
                    issues.append(issue)
            steps = {step.step_id: step for step in projection.steps}
            first = steps.get(route.first_step_id)
            second = steps.get(route.second_step_id)
            if first is None:
                issues.append("first_source_step_missing")
            if second is None:
                issues.append("second_source_step_missing")
            if first is not None:
                if (
                    route.first_source_reaction_id
                    and first.source_reaction_id
                    != route.first_source_reaction_id
                ):
                    issues.append("first_source_reaction_id_mismatch")
                if _reaction_without_reagents(first.reaction_smiles) != (
                    route.first_reaction_smiles
                ):
                    issues.append("first_source_reaction_mismatch")
            if second is not None:
                if (
                    route.second_source_reaction_id
                    and second.source_reaction_id
                    != route.second_source_reaction_id
                ):
                    issues.append("second_source_reaction_id_mismatch")
                if _reaction_without_reagents(second.reaction_smiles) != (
                    route.second_reaction_smiles
                ):
                    issues.append("second_source_reaction_mismatch")
            link = next(
                (
                    value
                    for value in projection.lineage_links
                    if value.lineage_id == route.source_lineage_id
                ),
                None,
            )
            if link is None:
                issues.append("source_lineage_missing")
            elif (
                first is not None
                and second is not None
                and (
                    link.producer_reaction_node_id != first.reaction_node_id
                    or link.consumer_reaction_node_id != second.reaction_node_id
                )
            ):
                issues.append("source_lineage_step_mismatch")
            elif link.status != "unique":
                issues.append("source_lineage_not_unique")
        route_results.append(
            {
                "route_id": route.route_id,
                "source_route_id": route_id,
                "valid": not issues,
                "issues": sorted(issues),
            }
        )
    report = {
        "source_path": str(route_core_source),
        "source_sha256": _file_sha256(route_core_source),
        "declared_route_count": len(expected),
        "validated_route_count": sum(item["valid"] for item in route_results),
        "valid": all(item["valid"] for item in route_results),
        "routes": route_results,
    }
    if not report["valid"]:
        messages = [
            f"{item['route_id']}:{','.join(item['issues'])}"
            for item in route_results
            if not item["valid"]
        ]
        raise ValueError("observed route-source validation failed: " + "; ".join(messages))
    return report


def _step_operator_ids(library: Any, reaction_id: str) -> tuple[str, ...]:
    return tuple(
        sorted(
            {
                template.operator_id
                for template in library.templates
                if any(
                    precedent.reaction_id == reaction_id
                    for precedent in template.precedents
                )
            }
        )
    )


def expand_two_step_precedent_route(
    route: TwoStepPrecedentRouteDefinition,
    *,
    definition_id: str = "two_step_precedent_route_expansion_poc.v1@1.0",
    top_k_per_application: int = 100,
    stock_lookup: Any | None = None,
) -> PrecedentRouteExpansionResult:
    """Generate cumulative R0-R2 chemical space from one observed route.

    R0 replays only the observed inputs with L2 operators. R1 adds declared
    building blocks without relaxing L2 context. R2 adds the remaining declared
    blocks and L1 operators. Every propagated product is emitted by the existing
    forward predictor after reverse recovery and edit-signature agreement.
    """

    if top_k_per_application < 1:
        raise ValueError("top-k per application must be positive")
    rows = (
        {
            "reaction_id": route.first_step_id,
            "reference_id": route.reference_id,
            "reaction_smiles": route.first_reaction_smiles,
        },
        {
            "reaction_id": route.second_step_id,
            "reference_id": route.reference_id,
            "reaction_smiles": route.second_reaction_smiles,
        },
    )
    generic_library = build_generic_library(
        rows,
        levels=("L1", "L2"),
        admission_mode="data_driven",
    )
    forward_library = build_forward_library(generic_library)
    first_operator_ids = _step_operator_ids(
        generic_library, route.first_step_id
    )
    second_operator_ids = _step_operator_ids(
        generic_library, route.second_step_id
    )
    if not first_operator_ids or not second_operator_ids:
        raise ValueError(
            "both precedent steps must compile into forward-admitted operators"
        )

    target = canonical_molecule_collection(route.target_smiles)
    assert target is not None
    carried_pathways: dict[str, PrecedentRouteExpansionPathway] = {}
    level_results = []
    previous_products: set[str] = set()

    for level in EXPANSION_LEVELS:
        allowed_operator_levels = _OPERATOR_LEVELS[level]
        declared_inputs = _admitted_inputs(route.first_step_inputs, level)
        declared_partners = _admitted_inputs(route.second_step_partners, level)
        rejections: Counter[str] = Counter()
        input_evidence = tuple(
            _input_evidence(item, "first_step", stock_lookup)
            for item in declared_inputs
        )
        partner_evidence = tuple(
            _input_evidence(item, "second_step_partner", stock_lookup)
            for item in declared_partners
        )
        inputs = tuple(
            item
            for item, evidence in zip(declared_inputs, input_evidence)
            if _input_is_eligible(item, evidence)
        )
        partners = tuple(
            item
            for item, evidence in zip(declared_partners, partner_evidence)
            if _input_is_eligible(item, evidence)
        )
        rejections["first_step_stock_evidence_missing"] += (
            len(declared_inputs) - len(inputs)
        )
        rejections["second_step_partner_stock_evidence_missing"] += (
            len(declared_partners) - len(partners)
        )
        intermediate_seeds: dict[
            tuple[str, str], tuple[tuple[str, str, str], Any]
        ] = {}
        for input_set in inputs:
            prediction = predict_products(
                input_set.smiles,
                forward_library,
                operator_ids=first_operator_ids,
                levels=allowed_operator_levels,
                include_self_reactions=False,
                top_k=top_k_per_application,
            )
            participating_candidates = tuple(
                candidate
                for candidate in prediction.candidates
                if _all_query_components_participate(
                    candidate, prediction.canonical_starting_materials
                )
            )
            if level == "R0":
                expected_intermediate = canonical_molecule_collection(
                    route.intermediate_smiles
                )
                participating_candidates = tuple(
                    candidate
                    for candidate in participating_candidates
                    if candidate.product_smiles == expected_intermediate
                )
            if prediction.candidates and not participating_candidates:
                rejections["first_step_incomplete_component_participation"] += 1
            if not participating_candidates:
                rejections["first_step_no_supported_product"] += 1
            for candidate in participating_candidates:
                key = (candidate.product_smiles, input_set.smiles)
                current = intermediate_seeds.get(key)
                candidate_key = (
                    candidate.abstraction_level,
                    candidate.forward_operator_id,
                    candidate.pathway_id,
                )
                if current is None or candidate_key < current[0]:
                    intermediate_seeds[key] = (candidate_key, candidate)

        for (intermediate, first_input), (_, first_candidate) in sorted(
            intermediate_seeds.items()
        ):
            for partner in partners:
                second_start = _combine_intermediate_and_partner(
                    intermediate, partner.smiles
                )
                prediction = predict_products(
                    second_start,
                    forward_library,
                    operator_ids=second_operator_ids,
                    levels=allowed_operator_levels,
                    include_self_reactions=False,
                    top_k=top_k_per_application,
                )
                participating_candidates = tuple(
                    candidate
                    for candidate in prediction.candidates
                    if _all_query_components_participate(
                        candidate, prediction.canonical_starting_materials
                    )
                )
                if level == "R0":
                    participating_candidates = tuple(
                        candidate
                        for candidate in participating_candidates
                        if candidate.product_smiles == target
                    )
                if prediction.candidates and not participating_candidates:
                    rejections[
                        "second_step_incomplete_component_participation"
                    ] += 1
                if not participating_candidates:
                    rejections["second_step_no_supported_product"] += 1
                for second_candidate in participating_candidates:
                    pathway_id = digest(
                        "PRECRTEXP1",
                        route.route_id,
                        first_input,
                        intermediate,
                        partner.smiles,
                        second_candidate.product_smiles,
                        first_candidate.forward_operator_id,
                        second_candidate.forward_operator_id,
                    )
                    if pathway_id in carried_pathways:
                        continue
                    pathway = PrecedentRouteExpansionPathway(
                        pathway_id=pathway_id,
                        first_discovered_level=level,
                        product_smiles=second_candidate.product_smiles,
                        intermediate_smiles=intermediate,
                        first_step_starting_materials=(
                            first_candidate.participating_precursor_smiles
                        ),
                        second_step_partner_smiles=partner.smiles,
                        second_step_starting_materials=(
                            second_candidate.participating_precursor_smiles
                        ),
                        first_operator_id=first_candidate.operator_id,
                        first_forward_operator_id=(
                            first_candidate.forward_operator_id
                        ),
                        first_abstraction_level=(
                            first_candidate.abstraction_level
                        ),
                        first_reverse_round_trip=(
                            first_candidate.reverse_round_trip
                        ),
                        first_operator_edit_agreement=(
                            first_candidate.operator_edit_agreement
                        ),
                        second_operator_id=second_candidate.operator_id,
                        second_forward_operator_id=(
                            second_candidate.forward_operator_id
                        ),
                        second_abstraction_level=(
                            second_candidate.abstraction_level
                        ),
                        second_reverse_round_trip=(
                            second_candidate.reverse_round_trip
                        ),
                        second_operator_edit_agreement=(
                            second_candidate.operator_edit_agreement
                        ),
                        warnings=tuple(
                            sorted(
                                set(
                                    first_candidate.warnings
                                    + second_candidate.warnings
                                )
                            )
                        ),
                    )
                    carried_pathways[pathway_id] = pathway

        product_smiles = tuple(
            sorted({pathway.product_smiles for pathway in carried_pathways.values()})
        )
        level_results.append(
            PrecedentRouteExpansionLevelResult(
                level=level,
                claim=_LEVEL_CLAIMS[level],
                allowed_operator_levels=allowed_operator_levels,
                first_step_input_count=len(inputs),
                second_step_partner_count=len(partners),
                intermediate_count=len(
                    {
                        pathway.intermediate_smiles
                        for pathway in carried_pathways.values()
                    }
                ),
                product_count=len(product_smiles),
                new_product_count=len(set(product_smiles) - previous_products),
                product_smiles=product_smiles,
                pathways=tuple(
                    carried_pathways[key] for key in sorted(carried_pathways)
                ),
                input_evidence=tuple(
                    sorted(
                        (*input_evidence, *partner_evidence),
                        key=lambda evidence: (
                            evidence.role,
                            evidence.input_smiles,
                            evidence.minimum_level,
                        ),
                    )
                ),
                rejection_counts=dict(
                    sorted(
                        (key, count)
                        for key, count in rejections.items()
                        if count
                    )
                ),
            )
        )
        previous_products = set(product_smiles)

    exact_replay_valid = target in set(level_results[0].product_smiles)
    warnings = []
    if not exact_replay_valid:
        warnings.append("EXACT_ROUTE_REPLAY_FAILED")
    if generic_library.rejection_counts:
        warnings.append("SOURCE_STEP_COMPILATION_HAS_REJECTIONS")
    expansion_id = digest(
        "PRECRTSPACE1",
        PRECEDENT_ROUTE_EXPANSION_ALGORITHM_VERSION,
        definition_id,
        route.route_id,
        *(level_results[-1].product_smiles),
    )
    return PrecedentRouteExpansionResult(
        expansion_id=expansion_id,
        route_id=route.route_id,
        reference_id=route.reference_id,
        source_dataset_id=route.source_dataset_id,
        source_route_id=route.source_route_id,
        patent_id=route.patent_id,
        source_tree_id=route.source_tree_id,
        source_route_core_id=route.source_route_core_id,
        source_lineage_id=route.source_lineage_id,
        first_source_reaction_id=route.first_source_reaction_id,
        second_source_reaction_id=route.second_source_reaction_id,
        intermediate_smiles=canonical_molecule_collection(
            route.intermediate_smiles
        )
        or route.intermediate_smiles,
        target_smiles=target,
        exact_replay_valid=exact_replay_valid,
        compiled_template_count=len(generic_library.templates),
        admitted_forward_operator_count=len(forward_library.operators),
        first_step_operator_ids=first_operator_ids,
        second_step_operator_ids=second_operator_ids,
        levels=tuple(level_results),
        warnings=tuple(warnings),
        definition_id=definition_id,
    )


def run_precedent_route_expansion_poc(
    definition: PrecedentRouteExpansionDefinition | str | Path = (
        DEFAULT_PRECEDENT_ROUTE_EXPANSION_DEFINITION
    ),
    *,
    output_path: str | Path | None = None,
    stock_index_path: str | Path | None = None,
    route_core_source: str | Path | None = None,
) -> dict[str, Any]:
    """Run a definition panel and optionally write deterministic JSON."""

    resolved = (
        load_precedent_route_expansion_definition(definition)
        if isinstance(definition, (str, Path))
        else definition
    )
    requires_stock = any(
        item.source_kind == "curated_stock"
        for route in resolved.routes
        for item in (*route.first_step_inputs, *route.second_step_partners)
    )
    if requires_stock and stock_index_path is None:
        raise ValueError("this expansion definition requires a stock index")
    requires_route_core = any(route.source_route_id for route in resolved.routes)
    if requires_route_core and route_core_source is None:
        raise ValueError(
            "observed route definitions require their route-core source"
        )
    source_validation = (
        validate_precedent_route_sources(resolved, route_core_source)
        if route_core_source is not None
        else {}
    )
    stock_metadata: dict[str, Any] = {}
    if stock_index_path is not None:
        with open_stock_lookup(stock_index_path) as stock_lookup:
            if hasattr(stock_lookup, "metadata"):
                stock_metadata = dict(stock_lookup.metadata())
            results = tuple(
                expand_two_step_precedent_route(
                    route,
                    definition_id=resolved.definition_id,
                    stock_lookup=stock_lookup,
                )
                for route in resolved.routes
            )
    else:
        results = tuple(
            expand_two_step_precedent_route(
                route,
                definition_id=resolved.definition_id,
            )
            for route in resolved.routes
        )
    report = {
        "artifact_type": "two_step_precedent_route_expansion_poc",
        "schema_version": PRECEDENT_ROUTE_EXPANSION_SCHEMA_VERSION,
        "algorithm_version": PRECEDENT_ROUTE_EXPANSION_ALGORITHM_VERSION,
        "definition_id": resolved.definition_id,
        "definition_schema_version": resolved.schema_version,
        "stock_index_path": str(stock_index_path or ""),
        "stock_index_metadata": stock_metadata,
        "source_validation": source_validation,
        "route_count": len(results),
        "exact_replay_count": sum(result.exact_replay_valid for result in results),
        "level_product_counts": {
            level: sum(
                next(
                    item.product_count
                    for item in result.levels
                    if item.level == level
                )
                for result in results
            )
            for level in EXPANSION_LEVELS
        },
        "routes": [result.to_dict() for result in results],
    }
    if output_path is not None:
        destination = Path(output_path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(
            json.dumps(
                report,
                indent=2,
                sort_keys=True,
                ensure_ascii=True,
            )
            + "\n",
            encoding="utf-8",
        )
    return report


__all__ = [
    "DEFAULT_PRECEDENT_ROUTE_EXPANSION_DEFINITION",
    "EXPANSION_LEVELS",
    "PRECEDENT_ROUTE_EXPANSION_ALGORITHM_VERSION",
    "PRECEDENT_ROUTE_EXPANSION_SCHEMA_VERSION",
    "PrecedentRouteExpansionDefinition",
    "PrecedentRouteExpansionInput",
    "PrecedentRouteExpansionLevelResult",
    "PrecedentRouteExpansionPathway",
    "PrecedentRouteExpansionResult",
    "PrecedentRouteInputEvidence",
    "PrecedentRouteStockComponentEvidence",
    "TwoStepPrecedentRouteDefinition",
    "expand_two_step_precedent_route",
    "load_precedent_route_expansion_definition",
    "run_precedent_route_expansion_poc",
    "validate_precedent_route_sources",
]

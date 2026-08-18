"""Nested chemical-space expansion from observed two-step routes.

This proof of concept treats a precedent route as executable structural evidence:
each observed step is compiled into source-round-tripped forward operators, and
only products that pass the normal forward/reverse graph checks are propagated.
The expansion levels are cumulative and make increasing extrapolation explicit.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Any, Iterable, Literal, Mapping

from reactive_taxonomy import canonical_molecule_collection

from forward_synthesis import build_forward_library, predict_products

from .chemistry import digest, split_reaction_smiles
from .generic_library import build_generic_library


PRECEDENT_ROUTE_EXPANSION_SCHEMA_VERSION = "1.0"
PRECEDENT_ROUTE_EXPANSION_ALGORITHM_VERSION = (
    "two_step_precedent_route_expansion.v1"
)
DEFAULT_PRECEDENT_ROUTE_EXPANSION_DEFINITION = (
    Path(__file__).resolve().parent
    / "definitions"
    / "two_step_precedent_route_expansion_poc.v1.json"
)

ExpansionLevel = Literal["R0", "R1", "R2"]
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

    def __post_init__(self) -> None:
        if self.minimum_level not in EXPANSION_LEVELS:
            raise ValueError(f"unsupported expansion level: {self.minimum_level}")


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

    def __post_init__(self) -> None:
        if not self.route_id or not self.reference_id:
            raise ValueError("route and reference IDs are required")
        if self.first_step_id == self.second_step_id:
            raise ValueError("two-step route step IDs must be distinct")
        if not any(item.minimum_level == "R0" for item in self.first_step_inputs):
            raise ValueError("the first step requires an R0 input")
        if not any(
            item.minimum_level == "R0" for item in self.second_step_partners
        ):
            raise ValueError("the second step requires an R0 partner")
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
        if self.schema_version != PRECEDENT_ROUTE_EXPANSION_SCHEMA_VERSION:
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
            "rejection_counts": dict(sorted(self.rejection_counts.items())),
        }


@dataclass(frozen=True)
class PrecedentRouteExpansionResult:
    """Auditable nested chemical space generated from one observed route."""

    expansion_id: str
    route_id: str
    reference_id: str
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
    return PrecedentRouteExpansionInput(
        smiles=str(value["smiles"]),
        minimum_level=str(value["minimum_level"]),  # type: ignore[arg-type]
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
        inputs = _admitted_inputs(route.first_step_inputs, level)
        partners = _admitted_inputs(route.second_step_partners, level)
        rejections: Counter[str] = Counter()
        intermediate_seeds: dict[
            tuple[str, str], tuple[tuple[str, str, str], Any]
        ] = {}
        for input_set in inputs:
            prediction = predict_products(
                input_set.smiles,
                forward_library,
                operator_ids=first_operator_ids,
                levels=allowed_operator_levels,
                top_k=top_k_per_application,
            )
            if not prediction.candidates:
                rejections["first_step_no_supported_product"] += 1
            for candidate in prediction.candidates:
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
                second_start = f"{intermediate}.{partner.smiles}"
                prediction = predict_products(
                    second_start,
                    forward_library,
                    operator_ids=second_operator_ids,
                    levels=allowed_operator_levels,
                    top_k=top_k_per_application,
                )
                if not prediction.candidates:
                    rejections["second_step_no_supported_product"] += 1
                for second_candidate in prediction.candidates:
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
                rejection_counts=dict(sorted(rejections.items())),
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
        definition_id,
        route.route_id,
        *(level_results[-1].product_smiles),
    )
    return PrecedentRouteExpansionResult(
        expansion_id=expansion_id,
        route_id=route.route_id,
        reference_id=route.reference_id,
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
) -> dict[str, Any]:
    """Run a definition panel and optionally write deterministic JSON."""

    resolved = (
        load_precedent_route_expansion_definition(definition)
        if isinstance(definition, (str, Path))
        else definition
    )
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
    "TwoStepPrecedentRouteDefinition",
    "expand_two_step_precedent_route",
    "load_precedent_route_expansion_definition",
    "run_precedent_route_expansion_poc",
]

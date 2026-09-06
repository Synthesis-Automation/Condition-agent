"""Bounded multistep evaluation for targets selected by partition review."""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Callable, Mapping, Optional

from .chemistry import canonical_smiles, digest
from .generic_models import GenericTemplateLibrary
from .multistep import MultistepRetrosynthesisResult, plan_multistep_routes
from .multistep_panel_review import (
    MultistepPanelCase,
    MultistepPanelTarget,
    render_multistep_panel_html,
)
from .observed_route_action import build_observed_route_action_label
from .route_contract import ReactionRouteTree, iter_molecule_occurrences


MULTISTEP_DATASET_EVALUATION_SCHEMA_VERSION = "1.0"
MULTISTEP_DATASET_EVALUATION_ALGORITHM_VERSION = (
    "guarded_dataset_multistep_evaluation.v1"
)


@dataclass(frozen=True)
class MultistepDatasetEvaluationConfig:
    """Explicit limits and terminal assumptions for the bounded search."""

    max_depth: int = 3
    molecular_weight_threshold: float = 150.0
    top_k_routes: int = 3
    per_step_top_k: int = 3
    beam_width: int = 8
    max_expansions: int = 15
    max_templates_to_apply: int = 100
    max_candidates_to_validate: int = 25
    use_context: bool = True
    include_l0: bool = True
    diversify: bool = True
    use_hierarchical_ranking: bool = True
    allow_untyped_literature_terminals: bool = False
    route_state_definition_id: str = ""
    route_state_catalog_sha256: str = ""
    route_state_ordering_enabled: bool = False

    def __post_init__(self) -> None:
        for value, label in (
            (self.max_depth, "maximum depth"),
            (self.top_k_routes, "top-k routes"),
            (self.per_step_top_k, "per-step top-k"),
            (self.beam_width, "beam width"),
            (self.max_expansions, "maximum expansions"),
            (self.max_templates_to_apply, "maximum templates"),
            (self.max_candidates_to_validate, "maximum validations"),
        ):
            if value < 1:
                raise ValueError(f"{label} must be positive")
        if self.molecular_weight_threshold <= 0:
            raise ValueError("molecular-weight threshold must be positive")

    @property
    def config_id(self) -> str:
        """Return a deterministic identity for all search controls."""

        return digest(
            "MSDECFG1",
            json.dumps(asdict(self), sort_keys=True, separators=(",", ":")),
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible configuration."""

        value = asdict(self)
        value["config_id"] = self.config_id
        return value


@dataclass(frozen=True)
class MultistepDatasetCaseResult:
    """Search result and observed-route recovery measurements for one target."""

    case_id: str
    selection_rank: int
    split: Optional[str]
    observed_reaction_count: int
    verified_observed_action_count: int
    maximum_observed_action_matches: int
    exact_root_action_recovered: bool
    full_observed_route_recovered: bool
    panel_case: MultistepPanelCase

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible case result."""

        return {
            "case_id": self.case_id,
            "selection_rank": self.selection_rank,
            "split": self.split,
            "observed_reaction_count": self.observed_reaction_count,
            "verified_observed_action_count": self.verified_observed_action_count,
            "maximum_observed_action_matches": (
                self.maximum_observed_action_matches
            ),
            "exact_root_action_recovered": self.exact_root_action_recovered,
            "full_observed_route_recovered": self.full_observed_route_recovered,
            "planner_solved": bool(self.panel_case.baseline.routes),
            "planner_partial_only": (
                not self.panel_case.baseline.routes
                and bool(self.panel_case.baseline.partial_routes)
            ),
            "panel_case": self.panel_case.to_dict(),
        }


@dataclass(frozen=True)
class MultistepDatasetEvaluation:
    """Self-contained evaluation of a fixed partition-review target set."""

    evaluation_id: str
    source_review_path: str
    source_review_sha256: str
    library_path: str
    library_sha256: str
    library_admission_policy: str
    stock_index_path: str
    config: MultistepDatasetEvaluationConfig
    cases: tuple[MultistepDatasetCaseResult, ...]
    warnings: tuple[str, ...]
    schema_version: str = MULTISTEP_DATASET_EVALUATION_SCHEMA_VERSION
    algorithm_version: str = MULTISTEP_DATASET_EVALUATION_ALGORITHM_VERSION

    @property
    def summary(self) -> dict[str, int]:
        """Return aggregate search and recovery counts."""

        return {
            "target_count": len(self.cases),
            "planner_solved_count": sum(
                bool(case.panel_case.baseline.routes) for case in self.cases
            ),
            "planner_partial_only_count": sum(
                not case.panel_case.baseline.routes
                and bool(case.panel_case.baseline.partial_routes)
                for case in self.cases
            ),
            "exact_root_action_recovery_count": sum(
                case.exact_root_action_recovered for case in self.cases
            ),
            "any_observed_action_recovery_count": sum(
                case.maximum_observed_action_matches > 0 for case in self.cases
            ),
            "full_observed_route_recovery_count": sum(
                case.full_observed_route_recovered for case in self.cases
            ),
            "verified_observed_action_count": sum(
                case.verified_observed_action_count for case in self.cases
            ),
            "matched_observed_action_count": sum(
                case.maximum_observed_action_matches for case in self.cases
            ),
        }

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible evaluation."""

        return {
            "evaluation_id": self.evaluation_id,
            "source_review_path": self.source_review_path,
            "source_review_sha256": self.source_review_sha256,
            "library_path": self.library_path,
            "library_sha256": self.library_sha256,
            "library_admission_policy": self.library_admission_policy,
            "stock_index_path": self.stock_index_path,
            "config": self.config.to_dict(),
            "summary": self.summary,
            "cases": [case.to_dict() for case in self.cases],
            "warnings": list(self.warnings),
            "schema_version": self.schema_version,
            "algorithm_version": self.algorithm_version,
        }


MultistepPlanner = Callable[..., MultistepRetrosynthesisResult]
ProgressCallback = Callable[[int, int, str], None]


def _sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _observed_actions(
    tree: ReactionRouteTree,
) -> tuple[tuple[str, str], ...]:
    actions = []
    for molecule in iter_molecule_occurrences(tree):
        reaction = molecule.reaction
        if reaction is None:
            continue
        label = build_observed_route_action_label(
            reaction.reaction_smiles,
            route_product_smiles=molecule.smiles,
            reaction_id=reaction.evidence.source_reaction_id or reaction.step_id,
            reference_id=tree.patent_id or tree.source_route_id or tree.tree_id,
        )
        if not label.exact_precursors_verified:
            continue
        product = canonical_smiles(molecule.smiles)
        precursors = label.expected_precursor_smiles
        if product is not None and precursors is not None:
            actions.append((product, precursors))
    return tuple(actions)


def _generated_action_counter(route: Any) -> Counter[tuple[str, str]]:
    values: Counter[tuple[str, str]] = Counter()
    for step in route.steps:
        product = canonical_smiles(step.product_smiles)
        precursors = canonical_smiles(step.candidate.precursor_smiles)
        if product is not None and precursors is not None:
            values[(product, precursors)] += 1
    return values


def _maximum_matches(
    result: MultistepRetrosynthesisResult,
    observed: tuple[tuple[str, str], ...],
) -> int:
    expected = Counter(observed)
    routes = result.routes or result.partial_routes
    return max(
        (sum((expected & _generated_action_counter(route)).values()) for route in routes),
        default=0,
    )


def evaluate_partition_review_routes(
    review: Mapping[str, Any],
    library: GenericTemplateLibrary,
    literature_index: Any,
    *,
    config: MultistepDatasetEvaluationConfig = MultistepDatasetEvaluationConfig(),
    source_review_path: str = "<in-memory>",
    source_review_sha256: str = "",
    library_path: str = "<in-memory>",
    library_sha256: str = "",
    stock_index_path: str = "<in-memory>",
    planner: MultistepPlanner = plan_multistep_routes,
    progress: Optional[ProgressCallback] = None,
    search_guidance: Any = None,
    route_action_selector: Any = None,
) -> MultistepDatasetEvaluation:
    """Run the same bounded multistep search for every selected route root."""

    raw_cases = review.get("cases")
    if not isinstance(raw_cases, list) or not raw_cases:
        raise ValueError("partition review must contain selected cases")
    results = []
    total = len(raw_cases)
    for index, raw_case in enumerate(raw_cases, start=1):
        raw_tree = raw_case.get("route_tree")
        if not isinstance(raw_tree, Mapping):
            raise ValueError("partition review case omits its route tree")
        tree = ReactionRouteTree.from_dict(dict(raw_tree))
        if progress is not None:
            progress(index, total, tree.root.smiles)
        result = planner(
            tree.root.smiles,
            library,
            literature_index,
            max_depth=config.max_depth,
            molecular_weight_threshold=config.molecular_weight_threshold,
            top_k_routes=config.top_k_routes,
            per_step_top_k=config.per_step_top_k,
            beam_width=config.beam_width,
            max_expansions=config.max_expansions,
            max_templates_to_apply=config.max_templates_to_apply,
            max_candidates_to_validate=config.max_candidates_to_validate,
            use_context=config.use_context,
            include_l0=config.include_l0,
            diversify=config.diversify,
            use_hierarchical_ranking=config.use_hierarchical_ranking,
            allow_untyped_literature_terminals=(
                config.allow_untyped_literature_terminals
            ),
            search_guidance=search_guidance,
            route_action_selector=route_action_selector,
        )
        observed = _observed_actions(tree)
        matched = _maximum_matches(result, observed)
        root_recovered = bool(observed) and any(
            _generated_action_counter(route)[observed[0]] > 0
            for route in (result.routes or result.partial_routes)
        )
        reference_reactions = tuple(
            molecule.reaction.reaction_smiles
            for molecule in iter_molecule_occurrences(tree)
            if molecule.reaction is not None
        )
        panel_case = MultistepPanelCase(
            target=MultistepPanelTarget(
                target_id=str(raw_case.get("case_id") or tree.tree_id),
                name=(
                    f"Case {int(raw_case.get('selection_rank') or index)} · "
                    f"observed {tree.reaction_count}-step route"
                ),
                category=f"{tree.split or 'unknown'}_observed_route",
                smiles=tree.root.smiles,
            ),
            baseline=result,
            reference_route_id=(tree.source_route_id or tree.tree_id),
            reference_reactions=reference_reactions,
            reference_maximum_depth=tree.maximum_depth,
            evaluation_metrics=(
                (
                    "Known root action",
                    "recovered" if root_recovered else "not retained",
                ),
                ("Known actions", f"{matched}/{len(observed)}"),
                (
                    "Search result",
                    "solved by heuristics" if result.routes else "partial only",
                ),
            ),
        )
        results.append(
            MultistepDatasetCaseResult(
                case_id=str(raw_case.get("case_id") or tree.tree_id),
                selection_rank=int(raw_case.get("selection_rank") or index),
                split=tree.split,
                observed_reaction_count=tree.reaction_count,
                verified_observed_action_count=len(observed),
                maximum_observed_action_matches=matched,
                exact_root_action_recovered=root_recovered,
                full_observed_route_recovered=(
                    bool(observed)
                    and len(observed) == tree.reaction_count
                    and matched == len(observed)
                ),
                panel_case=panel_case,
            )
        )
    warnings = [
        "BOUNDED_SEARCH_IS_NOT_A_COMPLETE_RETROSYNTHESIS_PROOF",
        "MOLECULAR_WEIGHT_TERMINALS_ARE_HEURISTIC",
        "CONDITIONS_AND_EXPERIMENTAL_SELECTIVITY_NOT_ASSESSED",
        "KNOWN_ACTION_RECOVERY_IS_NOT_ROUTE_FEASIBILITY",
        "SAMPLE_INCLUDES_TRAIN_ROUTES_AND_IS_NOT_A_BLIND_ACCURACY_ESTIMATE",
    ]
    if config.allow_untyped_literature_terminals:
        warnings.append("UNTYPED_LITERATURE_TERMINALS_ENABLED")
    if config.route_state_definition_id:
        warnings.append(
            "TRAIN_ONLY_ROUTE_STATE_GUIDANCE_ENABLED_"
            f"{config.route_state_definition_id}"
        )
    if config.route_state_ordering_enabled:
        warnings.append("EXPERIMENTAL_ROUTE_STATE_ORDERING_ENABLED")
    policy = str(library.definition.get("core_admission_policy") or "pass_only")
    identity = digest(
        "MSDE1",
        source_review_sha256 or source_review_path,
        library_sha256 or library_path,
        stock_index_path,
        config.config_id,
    )
    return MultistepDatasetEvaluation(
        evaluation_id=identity,
        source_review_path=source_review_path,
        source_review_sha256=source_review_sha256,
        library_path=library_path,
        library_sha256=library_sha256,
        library_admission_policy=policy,
        stock_index_path=stock_index_path,
        config=config,
        cases=tuple(results),
        warnings=tuple(warnings),
    )


def build_multistep_dataset_evaluation(
    source_review: str | Path,
    library_path: str | Path,
    literature_index: Any,
    *,
    stock_index_path: str | Path,
    config: MultistepDatasetEvaluationConfig = MultistepDatasetEvaluationConfig(),
    progress: Optional[ProgressCallback] = None,
    route_state_catalog_path: str | Path | None = None,
    enable_route_state_ordering: bool = False,
) -> MultistepDatasetEvaluation:
    """Load file inputs and run a reproducible ten-target-style evaluation."""

    review_path = Path(source_review).resolve()
    resolved_library_path = Path(library_path).resolve()
    review = json.loads(review_path.read_text(encoding="utf-8"))
    from .generic_library import load_generic_library

    library = load_generic_library(resolved_library_path)
    search_guidance = None
    route_action_selector = None
    if route_state_catalog_path is not None:
        from .route_state_learning import (
            LiteratureRouteActionSelector,
            LiteratureRouteOrderingGuidance,
            load_route_state_learning_catalog,
        )

        route_state_catalog = load_route_state_learning_catalog(
            route_state_catalog_path
        )
        if enable_route_state_ordering:
            search_guidance = LiteratureRouteOrderingGuidance(
                route_state_catalog
            )
        route_action_selector = LiteratureRouteActionSelector(route_state_catalog)
    return evaluate_partition_review_routes(
        review,
        library,
        literature_index,
        config=config,
        source_review_path=str(review_path),
        source_review_sha256=_sha256(review_path),
        library_path=str(resolved_library_path),
        library_sha256=_sha256(resolved_library_path),
        stock_index_path=str(Path(stock_index_path).resolve()),
        progress=progress,
        search_guidance=search_guidance,
        route_action_selector=route_action_selector,
    )


def write_multistep_dataset_evaluation(
    evaluation: MultistepDatasetEvaluation,
    output_json: str | Path,
    output_html: str | Path,
) -> dict[str, Any]:
    """Write a structured result and structure-rich chemist review."""

    summary = evaluation.summary
    metrics = (
        {"label": "Exact root actions", "value": summary["exact_root_action_recovery_count"]},
        {"label": "Any known action", "value": summary["any_observed_action_recovery_count"]},
        {"label": "Full known route", "value": summary["full_observed_route_recovery_count"]},
        {
            "label": "Known actions matched",
            "value": (
                f"{summary['matched_observed_action_count']}/"
                f"{summary['verified_observed_action_count']}"
            ),
        },
    )
    document = render_multistep_panel_html(
        (case.panel_case for case in evaluation.cases),
        title="Dataset-10 bounded multistep retrosynthesis evaluation",
        top_k=evaluation.config.top_k_routes,
        metadata={
            "panel_id": evaluation.evaluation_id,
            "summary_metrics": metrics,
            "warnings": evaluation.warnings,
        },
    )
    json_path = Path(output_json)
    html_path = Path(output_html)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    html_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(
        json.dumps(evaluation.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    html_path.write_text(document, encoding="utf-8", newline="\n")
    return {
        "evaluation_id": evaluation.evaluation_id,
        "summary": summary,
        "output_json": str(json_path.resolve()),
        "output_html": str(html_path.resolve()),
    }


__all__ = [
    "MULTISTEP_DATASET_EVALUATION_ALGORITHM_VERSION",
    "MULTISTEP_DATASET_EVALUATION_SCHEMA_VERSION",
    "MultistepDatasetCaseResult",
    "MultistepDatasetEvaluation",
    "MultistepDatasetEvaluationConfig",
    "build_multistep_dataset_evaluation",
    "evaluate_partition_review_routes",
    "write_multistep_dataset_evaluation",
]

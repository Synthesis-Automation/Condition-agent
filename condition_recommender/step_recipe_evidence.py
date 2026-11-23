"""Inspect recipe-specific experiments alongside the unchanged structural assessment."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping
from dataclasses import asdict, dataclass, replace
from functools import lru_cache
from hashlib import sha256
import json
from math import isclose, isfinite
from pathlib import Path
from typing import Any

from condition_registry.models import ResolvedConditionRecipe
from condition_registry.recipe_evidence import (
    RecipeComparison, compare_evidence_recipes, normalize_evidence_recipe,
)
from reactive_taxonomy.fragment_search import indexed_product
from reactive_taxonomy.shared_reaction_core import LEVELS
from reactive_taxonomy.step_precedent_context import StepContextComparison, build_step_context

from .generic_indexing import GenericReactionIndex, _valid_yield
from .processed_catalog import EXPERIMENT_FIELDS, ProcessedCatalog
from .reaction_precedent_support import (
    ReactionPrecedentEvidence, qualify_reaction_precedent, reaction_evidence_core,
)
from .shared_core_index import SharedCoreIndex
from .step_plausibility import (
    StepPlausibilityAssessment, StepPlausibilityMatch, _assess_step_and_pool,
    _assess_structural_match, _classify_outcome,
)


_COMPARISONS = ("matches_stated_fields", "partially_reported", "different", "unknown")
_TARGET_YIELDS = frozenset(("reported_yield_pct", "reported_product_yield_pct", "uv_area_yield_pct"))


@dataclass(frozen=True)
class RecipeExperiment:
    """One source experiment. The intended target graph can occur in a failed run."""

    precedent: ReactionPrecedentEvidence
    context: StepContextComparison
    eligible_for_measured_transfer: bool
    recipe_comparison: RecipeComparison
    outcome: str
    reported_outcomes: tuple[dict[str, Any], ...]
    reported_conditions: dict[str, Any]
    source_locator: dict[str, Any]
    reference: dict[str, Any]
    record_sha256: str
    index_eligibility: str | None
    admission_reasons: tuple[str, ...]
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class RecipeOutcomeCounts:
    """Outcome buckets cover structurally transferable observations in this stratum."""

    status: str
    experiment_count: int
    transferable_count: int
    above_threshold: int
    low_yield: int
    zero_yield: int
    trace: int
    unreported: int
    unusable: int


@dataclass(frozen=True)
class ExperimentSearchScope:
    """Visited candidates from one lookup before qualification and deduplication."""

    source: str
    scope: str
    candidates_examined: int
    lookup_complete: bool
    truncated: bool


@dataclass(frozen=True)
class StepRecipeEvidence:
    """Source observations partitioned by recorded conditions, alongside the step grade."""

    step_assessment: StepPlausibilityAssessment
    query_recipe: Mapping[str, Any]
    experiments: tuple[RecipeExperiment, ...]
    outcome_counts: tuple[RecipeOutcomeCounts, ...]
    experiment_count: int
    display_truncated: bool
    search_scopes: tuple[ExperimentSearchScope, ...]
    warnings: tuple[str, ...]
    schema_version: str = "step_recipe_evidence.v1"

    def to_dict(self) -> dict[str, Any]:
        """Serialize source locators, experiment conditions and all comparison evidence."""
        return asdict(self)


def _reported_outcome(
    precedent: ReactionPrecedentEvidence, raw: tuple[dict, ...], threshold: float,
) -> tuple[str, float | None, set[str]]:
    """Interpret a yield only for the specified source target and measurement basis."""
    target = indexed_product(precedent.reaction_smiles.split(">")[-1])[0]
    values: set[float] = set()
    trace = invalid = False
    warnings: set[str] = set()
    product_fields = ("product_smiles", "target_smiles", "product_position", "product_cas")
    attributed = any(
        item.get("outcome_type") in _TARGET_YIELDS and
        any((item.get("metadata") or {}).get(key) is not None for key in product_fields)
        for item in raw
    )
    for item in raw:
        if item.get("outcome_type") not in _TARGET_YIELDS:
            continue
        metadata = item.get("metadata") or {}
        if attributed and not any(metadata.get(key) is not None for key in product_fields):
            warnings.add("UNATTRIBUTED_YIELD_WITH_PRODUCT_SPECIFIC_ASSAYS")
            continue
        product = metadata.get("product_smiles") or metadata.get("target_smiles")
        if product:
            try:
                is_target = indexed_product(product)[0] == target
            except (TypeError, ValueError):
                is_target = False
            if not is_target:
                warnings.add("NONMATCHING_MEASUREMENT_PRODUCT")
                continue
        elif ("." in target or metadata.get("product_position") not in (None, 1, "1")
              or (metadata.get("product_cas") and metadata.get("product_position") is None)):
            warnings.add("MEASUREMENT_TARGET_UNRESOLVED")
            continue
        if str(item.get("unit") or "").casefold() not in ("percent", "%", "pct"):
            invalid = True
            continue
        number = item.get("value")
        reported = str(item.get("raw_value") or "").strip().casefold()
        if reported in ("trace", "traces"):
            trace = True
            if number is not None:
                invalid = True
        elif number is not None:
            if type(number) not in (float, int) or not isfinite(number) or not 0 <= number <= 100:
                invalid = True
            else:
                if reported:
                    try:
                        # Published text and a parsed float32 may differ in trailing digits.
                        if not isclose(float(reported.removesuffix("%")), number, rel_tol=1e-6):
                            invalid = True
                    except ValueError:
                        invalid = True
                values.add(number)
        elif reported:
            invalid = True
    if invalid or len(values) > 1 or (values and trace):
        return "unusable", None, warnings | {"TARGET_YIELD_AMBIGUOUS_OR_UNUSABLE"}
    if trace:
        return "trace", None, warnings
    if values:
        value = values.pop()
        status = _classify_outcome(replace(precedent, yield_pct=value, outcome_status="usable"), threshold)
        return status, value, warnings
    return "unreported", None, warnings


def _experiment(match: StepPlausibilityMatch, record: Mapping[str, Any], query: dict,
                threshold: float) -> RecipeExperiment:
    row = match.precedent
    source = record.get("source") or {}
    fields = source.get("raw_fields") or {}
    original = fields.get("source_observation") or {}
    origin = original.get("source") or {}
    raw = tuple(dict(value) for value in original.get("outcomes") or ())
    locator = {key: record[key] for key in ("source_dataset", "source_path", "source_row_number")
               if record.get(key) is not None}
    locator.update(origin)
    raw_fields = original.get("raw_fields") or fields.get("source_raw_fields") or fields
    for key in ("source_location", "source_url", "source_binary_sha256",
                "table", "entry", "well", "plate", "plate_id", "well_id"):
        if key in raw_fields:
            locator[key] = raw_fields[key]
    if source.get("experimental_procedure"):
        locator["experimental_procedure"] = source["experimental_procedure"]
    reference = dict(record.get("reference_identity") or {})
    reference.setdefault("reference_id", row.reference_id)
    if origin.get("reference") or source.get("reference"):
        reference.setdefault("raw_reference", origin.get("reference") or source["reference"])
    warnings: set[str] = set()
    outcome = match.outcome
    if "outcomes" in original:
        outcome, value, warnings = _reported_outcome(row, raw, threshold)
        status = {"unreported": "unreported", "unusable": "invalid", "trace": "reported_trace"}
        row = replace(row, yield_pct=value, outcome_status=status.get(outcome, "usable"))
    elif outcome != "unreported":
        warnings.add("CANONICAL_YIELD_WITHOUT_RAW_ASSAY")
    if row.resolved_recipe is None:
        warnings.add("SOURCE_RECIPE_UNAVAILABLE")
    try:
        comparison = compare_evidence_recipes(query, row.resolved_recipe)
    except (AttributeError, KeyError, TypeError, ValueError):
        comparison = compare_evidence_recipes(query, {})
        warnings.add("SOURCE_RECIPE_UNINTERPRETABLE")
    if comparison.source_only:
        warnings.add("SOURCE_CONTAINS_ADDITIONAL_RECORDED_CONDITIONS")
    reported_conditions = dict(original.get("conditions") or {})
    for stage in reported_conditions.get("stages") or ():
        if (stage.get("pressure_bar") is not None or stage.get("concentration_m") is not None
                or (stage.get("provenance") or {}).get("reactor")):
            warnings.add("INSPECT_ADDITIONAL_REPORTED_STAGE_OPERATIONS")
    payload = (row.observation_id, row.reaction_id, row.reaction_smiles, row.reference_id,
               locator, row.resolved_recipe, row.outcome_status, row.yield_pct, raw, reported_conditions)
    digest = sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return RecipeExperiment(
        row, match.context, match.comparable_for_measured_support, comparison, outcome,
        raw, reported_conditions, locator, reference, digest, record.get("index_eligibility"),
        tuple(record.get("admission_reasons") or ()), tuple(sorted(warnings)),
    )


def _counts(status: str, experiments: list[RecipeExperiment]) -> RecipeOutcomeCounts:
    sub = [e for e in experiments if e.recipe_comparison.status == status]
    counts = Counter(e.outcome for e in sub if e.eligible_for_measured_transfer)
    return RecipeOutcomeCounts(status, len(sub), sum(counts.values()), counts["above_threshold"],
                               counts["low_yield"], counts["zero_yield"], counts["trace"],
                               counts["unreported"], counts["unusable"])


def _projection_key(row: ReactionPrecedentEvidence) -> tuple:
    """Identify exactly the source data available from an unhydrated index row."""
    return (row.observation_id, row.reaction_id, row.reaction_smiles, row.reference_id,
            row.outcome_status, _valid_yield(row.yield_pct),
            json.dumps(row.resolved_recipe, sort_keys=True, separators=(",", ":")))


def _display(experiments: list[RecipeExperiment], limit: int) -> tuple[RecipeExperiment, ...]:
    def priority(e: RecipeExperiment) -> tuple:
        sim = e.context.minimum_partner_similarity
        return (not e.eligible_for_measured_transfer, _COMPARISONS.index(e.recipe_comparison.status),
                ("whole_reaction", *LEVELS).index(e.precedent.level),
                -(sim if sim is not None else -1), len(e.recipe_comparison.differences),
                len(e.recipe_comparison.unknowns), -len(e.recipe_comparison.agreements), e.precedent.reference_id,
                e.precedent.reaction_id, e.record_sha256)

    ordered = sorted(experiments, key=priority)
    # Keep opposed outcomes visible among equally close, recipe-matched experiments.
    if ordered:
        best = priority(ordered[0])[:-3]
        peers = [e for e in ordered if priority(e)[:-3] == best]
        adverse = next((e for e in peers if e.outcome in ("zero_yield", "low_yield", "trace")), None)
        positive = next((e for e in peers if e.outcome == "above_threshold"), None)
        leading = [e for e in (adverse, positive) if e is not None]
        seen = {e.record_sha256 for e in leading}
        ordered = leading + [e for e in ordered if e.record_sha256 not in seen]
    return tuple(ordered[:limit])


def assess_step_recipe_evidence(
    reaction_smiles: str, recipe: ResolvedConditionRecipe | Mapping[str, Any],
    index: GenericReactionIndex, shared_index: SharedCoreIndex, *,
    catalog: ProcessedCatalog | None = None, fragment_index: str | Path | None = None,
    experiment_catalog: ProcessedCatalog | None = None, candidate_limit: int = 128,
    match_limit: int = 20,
) -> StepRecipeEvidence:
    """Compare graph-qualified source experiments to the fields in a proposed recipe.

    Each lookup is bounded independently. Included controls can be absent from the
    condition recommendation index. Counts describe observed target outcomes in
    the visited pool and preserve the measurement basis of the original experiment.
    """
    query = normalize_evidence_recipe(recipe)
    assessment, pool = _assess_step_and_pool(
        reaction_smiles, index, shared_index, candidate_limit=candidate_limit, match_limit=match_limit,
    )
    core_for = lru_cache(maxsize=None)(reaction_evidence_core)
    context_for = lru_cache(maxsize=None)(build_step_context)
    threshold = assessment.policy.positive_yield_strictly_above_pct
    graph_resolved = assessment.grade != "unresolved"
    scopes = [ExperimentSearchScope("canonical_condition_index",
                                    "bounded_shared_cores" if graph_resolved else "query_graph_unresolved",
                                    assessment.candidate_count, graph_resolved and not assessment.candidate_truncated,
                                    assessment.candidate_truncated)]
    warnings = {"MATCHES_APPLY_ONLY_TO_EXPLICIT_QUERY_FIELDS", "COUNTS_DESCRIBE_VISITED_EXPERIMENTS"}
    if not graph_resolved:
        warnings.add("QUERY_GRAPH_UNRESOLVED_EXPERIMENT_SEARCH_NOT_PERFORMED")
    found: dict[str, RecipeExperiment] = {}
    unhydrated: dict[tuple, str] = {}
    for match in pool:
        row = match.precedent
        record = {}
        if catalog is not None and row.observation_id:
            try:
                record = catalog.observation(row.observation_id, EXPERIMENT_FIELDS)
            except KeyError:
                warnings.add("INDEXED_OBSERVATION_MISSING_FROM_CATALOG")
            if record:
                source_row = replace(
                    row, observation_id=record.get("observation_id"), reaction_id=record["reaction_id"],
                    reaction_smiles=record["reaction_smiles"], reference_id=record.get("reference_id") or "",
                    resolved_recipe=record.get("resolved_recipe"), yield_pct=record.get("yield_pct"),
                    outcome_status=record.get("outcome_status") or "",
                )
                if _projection_key(row) != _projection_key(source_row):
                    raise ValueError("Catalog identity, recipe or outcome disagrees with indexed observation")
        experiment = _experiment(match, record, query, threshold)
        found[experiment.record_sha256] = experiment
        if not record and row.observation_id:
            unhydrated[_projection_key(row)] = experiment.record_sha256

    for store, label in ((catalog, "release_experiment_catalog"),
                         (experiment_catalog, "supplemental_experiment_catalog")):
        if store is None:
            continue
        if not graph_resolved:
            scopes.append(ExperimentSearchScope(label, "query_graph_unresolved", 0, False, False))
            continue
        lookup = store.experiment_candidates(
            reaction_smiles, limit=candidate_limit,
            include_analogues=store.supports_experiment_analogues,
            fragment_index=fragment_index if label == "release_experiment_catalog" else None,
        )
        scopes.append(ExperimentSearchScope(label, lookup.scope, lookup.candidates_examined,
                                            lookup.lookup_complete, lookup.truncated))
        query_core = core_for(reaction_smiles)
        for record in lookup.records:
            smiles = record["reaction_smiles"]
            core = core_for(smiles)
            if core is None or query_core is None:
                continue
            row = qualify_reaction_precedent(
                query_core, core, reaction_smiles=smiles, reaction_id=record["reaction_id"],
                observation_id=record["observation_id"], reference_id=record.get("reference_id") or "",
                source=label, yield_pct=record.get("yield_pct"), resolved_recipe=record.get("resolved_recipe"),
                outcome_status=record.get("outcome_status") or "unreported",
            )
            if row is None:
                continue
            match = _assess_structural_match(assessment.query_context, context_for(smiles), row, assessment.policy)
            experiment = _experiment(match, record, query, threshold)
            placeholder = unhydrated.pop(_projection_key(row), None)
            if placeholder:
                del found[placeholder]
            found.setdefault(experiment.record_sha256, experiment)
    observations = Counter(e.precedent.observation_id for e in found.values() if e.precedent.observation_id)
    if any(count > 1 for count in observations.values()):
        warnings.add("OBSERVATION_ID_REUSED_FOR_DISTINCT_EXPERIMENT_RECORDS")
    if any(not scope.lookup_complete for scope in scopes):
        warnings.add("SOME_CANDIDATE_SEARCHES_INCOMPLETE")
    if catalog is None:
        warnings.add("RELEASE_EXPERIMENT_CATALOG_UNAVAILABLE")
    experiments = list(found.values())
    return StepRecipeEvidence(
        assessment, query, _display(experiments, match_limit),
        tuple(_counts(status, experiments) for status in _COMPARISONS), len(experiments),
        len(experiments) > match_limit, tuple(scopes), tuple(sorted(warnings)),
    )

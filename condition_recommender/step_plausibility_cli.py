"""Assess precedent support for a proposed reaction using a bound local index pair."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sqlite3
from textwrap import shorten

from .generic_indexing import load_generic_index
from .processed_catalog import ProcessedCatalog
from .processed_release import DEFAULT_PROCESSED_ROOT, resolve_processed_release
from .reaction_precedent_support import ReactionPrecedentEvidence
from .shared_core_index import load_shared_core_index
from .step_plausibility import StepPlausibilityAssessment, assess_step_plausibility
from .step_recipe_evidence import StepRecipeEvidence, assess_step_recipe_evidence


def _graph_diagnostic(assessment: StepPlausibilityAssessment) -> str | None:
    recovery = assessment.mapping_recovery
    if recovery and recovery.status == "conflicting_mappings":
        return "Exact source reactions disagree on the atom correspondence."
    if recovery and recovery.status == "truncated":
        return "The exact-source lookup reached its limit before agreement could be established."
    if assessment.query_context.product_element_excess:
        deficits = ", ".join(
            f"{symbol} ({count})" for symbol, count in assessment.query_context.product_element_excess
        )
        return f"Product atoms unaccounted for by the submitted reactants: {deficits}"
    if "REACTION_COMPLETENESS_UNRESOLVED" in assessment.query_context.warnings:
        return "Product atom correspondence could not be verified from the submitted structures."
    return None


def _brief_precedent(
    source: ReactionPrecedentEvidence, outcome: str, comparable: bool, reference: dict,
    *, agreement: str = "", location: str = "",
) -> str:
    """Fit the structural match, observed yield and a findable citation on one line."""
    fields = [agreement] if agreement else []
    fields.append("Exact reaction" if source.level == "whole_reaction" else "Reaction-core analogue")
    if not comparable:
        fields.append("outside transfer criteria")
    if outcome == "trace":
        fields.append("trace yield")
    elif outcome == "unusable":
        fields.append("yield unusable")
    elif source.yield_pct is None:
        fields.append("yield unreported")
    else:
        fields.append(f"{source.yield_pct:g}% yield")

    citation = reference.get("doi") or reference.get("patent_number")
    if not citation:
        raw = reference.get("raw_reference") or reference.get("normalized_citation") or ""
        # The catalog's "title | authors | journal" form has a short, locatable last field.
        parts = raw.split(" | ")
        label = parts[-1] if len(parts) == 3 and parts[-1].strip() else raw
        citation = shorten(label, width=76, placeholder="...") or source.reference_id or "reference unavailable"
    fields.append(citation)
    if location:
        fields.append(shorten(str(location), width=50, placeholder="..."))
    return " | ".join(fields)


def _concise_summary(
    step: StepPlausibilityAssessment, citations: dict[str, dict], recipe: StepRecipeEvidence | None = None,
) -> str:
    """Print the step grade and one line per selected source measurement."""
    if step.evidence_rank is None:
        reason = _graph_diagnostic(step)
        if not reason:
            reason = (
                "Supply a forward reaction as reactants>>product."
                if "INVALID_REACTION_FORMAT" in step.warnings else
                "The reaction graph could not be verified from the supplied structures."
            )
        return f"Ungraded  |  Unresolved reaction\n{reason}"

    counts = step.support_counts
    threshold = step.policy.low_yield_below_pct
    lines = [f"{step.evidence_rank}/4  |  {step.grade.replace('_', ' ').capitalize()}"]
    if recipe is None:
        groups = counts.supporting_source_group_count
        summary = (f"{step.qualified_count} qualified {'precedent' if step.qualified_count == 1 else 'precedents'}"
                   f" | {groups} supporting source {'group' if groups == 1 else 'groups'}")
        adverse = counts.low_yield_count + counts.zero_yield_count
        if adverse:
            summary += f" | {adverse} comparable {'yield' if adverse == 1 else 'yields'} below {threshold:g}%"
        if "LOW_YIELD_REPORTED_FOR_SAME_OBSERVATION" in step.warnings:
            summary += " | low yield also applies to duplicate observations"
        rows = [
            _brief_precedent(m.precedent, m.outcome, m.structurally_comparable,
                             citations.get(m.precedent.reference_id, {}))
            for m in step.matches
        ]
        truncated = step.candidate_truncated
    else:
        group = next(g for g in recipe.outcome_counts if g.status == "matches_stated_fields")
        subject = "experiment matches" if group.transferable_count == 1 else "experiments match"
        summary = f"{group.transferable_count} comparable {subject} the stated recipe"
        outcomes = (
            (group.at_or_above_threshold, f"at or above {threshold:g}%"),
            (group.low_yield, f"nonzero below {threshold:g}%"), (group.zero_yield, "zero"), (group.trace, "trace"),
            (group.unreported, "unreported"), (group.unusable, "unusable"),
        )
        if group.transferable_count:
            summary += " | " + ", ".join(f"{n} {label}" for n, label in outcomes if n)
        rows = [
            _brief_precedent(m.precedent, m.outcome, m.eligible_for_measured_transfer, m.reference,
                             agreement=m.recipe_comparison.status.replace("_", " "),
                             location=m.source_locator.get("source_location") or m.source_locator.get("source_record_id", ""))
            for m in recipe.experiments
        ]
        truncated = any(scope.truncated for scope in recipe.search_scopes)
        if not truncated and any(not scope.lookup_complete for scope in recipe.search_scopes):
            summary += " | partial search coverage"
    if truncated:
        summary += " | search limit reached"
    lines.append(summary)
    if step.mapping_recovery and step.mapping_recovery.mapped_reaction_smiles:
        lines.append("Atom correspondence recovered from an identical reported reaction.")
    lines.extend(f"  {number}. {row}" for number, row in enumerate(rows, 1))
    return "\n".join(lines)


def _summary(assessment: StepPlausibilityAssessment, citations: dict[str, dict]) -> str:
    counts = assessment.support_counts
    threshold = assessment.policy.low_yield_below_pct
    rank = "ungraded" if assessment.evidence_rank is None else f"{assessment.evidence_rank}/4"
    lines = [
        f"Precedent support: {assessment.grade} ({rank})",
        f"Strongest graph match: {assessment.strongest_level or 'none'}",
        f"Qualified observations: {assessment.qualified_count} of {assessment.candidate_count} visited",
        f"References among qualified observations: {counts.distinct_reference_count}",
        f"Close observations covering the query context: {counts.context_covered_close_observation_count}",
        f"Comparable observations eligible for source grouping: {counts.supporting_observation_count}"
        f" ({counts.missing_reference_supporting_count} missing a reference)",
        f"Source groups supporting the rank: {counts.supporting_source_group_count}",
        f"Distinct reactions represented by supporting references: {counts.supporting_distinct_reaction_count}",
        f"Yield rule: a reported target yield below {threshold:g}% excludes that observation from rank 3/4 support.",
        f"Outcomes in {counts.comparable_observation_count} comparable observations: "
        f"{counts.at_or_above_threshold_yield_count} at or above {threshold:g}%, "
        f"{counts.low_yield_count} nonzero below {threshold:g}%, {counts.zero_yield_count} zero, "
        f"{counts.unreported_yield_count} unreported, {counts.unusable_yield_count} unusable",
        f"Candidate search truncated: {assessment.candidate_truncated}",
        f"Displayed matches truncated: {assessment.matches_truncated}",
    ]
    diagnostic = _graph_diagnostic(assessment)
    if diagnostic:
        lines.append(diagnostic)
    if assessment.mapping_recovery and assessment.mapping_recovery.mapped_reaction_smiles:
        sources = ", ".join(source.reaction_id for source in assessment.mapping_recovery.sources)
        lines.append(f"Atom correspondence recovered from exact source reactions: {sources}")
    if assessment.query_context.requirements:
        lines.append("Retained context required on the reacting partners:")
        lines.extend(f"  {requirement}" for requirement in assessment.query_context.requirements)
    if assessment.matches:
        lines.append("Source matches (use --json for structures, context differences and recorded recipes):")
    for match in assessment.matches:
        source = match.precedent
        if match.outcome == "unusable":
            outcome = "yield unusable"
            if source.yield_pct is not None:
                outcome += f" (raw {source.yield_pct}%)"
        elif source.yield_pct is None:
            outcome = "yield unreported"
        else:
            outcome = f"source yield {source.yield_pct:g}%"
        similarity = match.context.minimum_partner_similarity
        similarity_text = "unresolved" if similarity is None else f"{similarity:.2f}"
        lines.append(
            f"  {source.level} | context {match.context.status} | "
            f"weakest partner similarity {similarity_text} | {outcome} | "
            f"{source.reference_id or 'reference unknown'} | {source.observation_id or source.reaction_id}"
        )
        if match.excluded_by_low_yield and match.outcome not in ("zero_yield", "low_yield"):
            lines.append(f"    The same observation also has a reported yield below {threshold:g}%.")
        if match.context.missing:
            lines.append("    Missing source context: " + ", ".join(match.context.missing))
        if match.context.unresolved:
            lines.append("    Unresolved source context: " + ", ".join(match.context.unresolved))
        if match.context.source_only:
            lines.append("    Extra source context absent in query: " + ", ".join(match.context.source_only))
        if match.context.source_unresolved:
            lines.append("    Extra source context unresolved: " + ", ".join(match.context.source_unresolved))
    if citations:
        lines.append("Source references:")
        for ref_id, ref in citations.items():
            label = ref.get("raw_reference") or ref.get("normalized_citation") or ref_id
            if ref.get("doi"):
                label = f"{label}  https://doi.org/{ref['doi']}"
            elif ref.get("patent_number"):
                label = f"{label}  https://patents.google.com/patent/{ref['patent_number']}"
            lines.append(f"  {ref_id}  {label}")
    if assessment.warnings:
        lines.append("Warnings: " + ", ".join(assessment.warnings))
    return "\n".join(lines)


def _recipe_summary(view: StepRecipeEvidence) -> str:
    """Show the measured outcomes behind each recipe comparison."""
    step = view.step_assessment
    threshold = step.policy.low_yield_below_pct
    rank = "ungraded" if step.evidence_rank is None else f"{step.evidence_rank}/4"
    lines = [
        f"Step precedent support: {step.grade} ({rank})",
        f"Recipe evidence: showing {len(view.experiments)} of {view.experiment_count} experiments in the bounded search",
    ]
    diagnostic = _graph_diagnostic(step)
    if diagnostic:
        lines.append(diagnostic)
    if step.grade == "unresolved" and step.warnings:
        lines.append("Structural evidence notes: " + ", ".join(step.warnings))
    lines.append("Outcome counts for structurally comparable experiments by recipe agreement:")
    for group in view.outcome_counts:
        lines.append(
            f"  {group.status.replace('_', ' ')} | {group.transferable_count} comparable of "
            f"{group.experiment_count} | {group.at_or_above_threshold} at or above {threshold:g}%, "
            f"{group.low_yield} nonzero below {threshold:g}%, {group.zero_yield} zero, {group.trace} trace, "
            f"{group.unreported} unreported, {group.unusable} unusable"
        )
    for match in view.experiments:
        source = match.precedent
        recipe = match.recipe_comparison
        similarity = match.context.minimum_partner_similarity
        similarity_text = "unresolved" if similarity is None else f"{similarity:.2f}"
        outcome = match.outcome.replace("_", " ")
        if source.yield_pct is not None and match.outcome != "unusable":
            outcome = f"{source.yield_pct:g}% ({outcome})"
        transfer = "comparable substrate" if match.eligible_for_measured_transfer else "outside transfer criteria"
        lines.append(
            f"  {recipe.status.replace('_', ' ')} | {source.level} | weakest partner "
            f"{similarity_text} | {transfer} | {outcome} | {source.observation_id or source.reaction_id}"
        )
        lines.append(f"    Source reaction: {source.reaction_smiles}")
        for measurement in match.reported_outcomes:
            value = measurement.get("raw_value")
            if value in (None, ""):
                value = measurement.get("value")
            if value is not None:
                basis = (measurement.get("metadata") or {}).get("measurement_basis")
                field = measurement.get("source_field") or measurement.get("outcome_type") or "source"
                lines.append(f"    Reported measurement ({field}): {value}" + (f" ({basis})" if basis else ""))
        for label, values in (
            ("Agrees", recipe.agreements), ("Differs", recipe.differences),
            ("Unknown", recipe.unknowns), ("Additional source conditions", recipe.source_only),
            ("Missing substrate context", match.context.missing),
            ("Additional substrate context", match.context.source_only),
            ("Unresolved substrate context", (*match.context.unresolved, *match.context.source_unresolved)),
        ):
            if values:
                lines.append(f"    {label}: {', '.join(values)}")
        ref = match.reference or {}
        citation = (
            ref.get("raw_reference") or ref.get("normalized_citation")
            or match.source_locator.get("reference") or source.reference_id
        )
        if ref.get("doi") and ref["doi"] not in (citation or ""):
            citation = f"{citation or ''}  https://doi.org/{ref['doi']}".strip()
        elif ref.get("patent_number") and ref["patent_number"] not in (citation or ""):
            citation = f"{citation or ''}  https://patents.google.com/patent/{ref['patent_number']}".strip()
        if citation:
            lines.append(f"    Reference: {citation}")
        locator = match.source_locator
        location = locator.get("source_location") or locator.get("source_record_id")
        if not location and locator.get("source_row_number") is not None:
            location = f"{locator.get('source_file') or locator.get('source_dataset') or 'source'}, row {locator['source_row_number']}"
        if location:
            lines.append(f"    Location: {location}")
        if match.warnings:
            lines.append("    Evidence notes: " + ", ".join(match.warnings))
    for scope in view.search_scopes:
        coverage = "complete in scope" if scope.lookup_complete else "incomplete"
        lines.append(f"Search coverage: {scope.source} | {scope.scope} | {scope.candidates_examined} examined | {coverage}")
    if view.warnings:
        lines.append("Coverage: " + ", ".join(view.warnings))
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    """Print an interpretable summary or the complete versioned JSON assessment."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reaction_smiles", help="Proposed reaction, e.g. 'precursor1.precursor2>>product'")
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--release", type=Path, help="Processed dataset root or pinned release manifest (default: active local release)")
    source.add_argument("--index", type=Path, help="Explicit canonical reaction index; also supply --shared-index")
    parser.add_argument("--shared-index", type=Path, help="Shared reaction-core index bound to --index")
    parser.add_argument("--recipe", type=Path, help="Compare the conditions in this recipe JSON with reported experiments")
    parser.add_argument("--experiment-catalog", type=Path, help="Additional evidence SQLite catalog to search with --recipe")
    parser.add_argument("--candidate-limit", type=int, default=128, help="Maximum candidates visited")
    parser.add_argument("--match-limit", type=int, help="Maximum records displayed (default: 3 with --concise, otherwise 20)")
    output = parser.add_mutually_exclusive_group()
    output.add_argument("--concise", action="store_true", help="Show the score and one short line per precedent")
    output.add_argument("--json", action="store_true", help="Print complete machine-readable evidence")
    args = parser.parse_args(argv)
    if args.match_limit is None:
        args.match_limit = 3 if args.concise else 20
    if (args.index is None) != (args.shared_index is None):
        parser.error("--index and --shared-index must be supplied together")
    if args.experiment_catalog and not args.recipe:
        parser.error("--experiment-catalog requires --recipe")
    try:
        catalog = None
        fragment_path = None
        if args.index is None:
            release = resolve_processed_release(args.release or DEFAULT_PROCESSED_ROOT)
            index_path = release.artifact("condition_index")
            shared_path = release.artifact("shared_core_index")
            if (args.recipe or not args.json) and "evidence_catalog" in release.manifest["artifacts"]:
                catalog = ProcessedCatalog(release.artifact("evidence_catalog"))
            if args.recipe and "fragment_index" in release.manifest["artifacts"]:
                fragment_path = release.artifact("fragment_index")
        else:
            index_path, shared_path = args.index, args.shared_index
        index = load_generic_index(index_path)
        shared_index = load_shared_core_index(shared_path, index)
        if args.recipe:
            recipe = json.loads(args.recipe.read_text(encoding="utf-8-sig"))
            if not isinstance(recipe, dict):
                parser.error("--recipe must contain a JSON object")
            view = assess_step_recipe_evidence(
                args.reaction_smiles, recipe, index, shared_index,
                catalog=catalog, fragment_index=fragment_path,
                experiment_catalog=ProcessedCatalog(args.experiment_catalog) if args.experiment_catalog else None,
                candidate_limit=args.candidate_limit, match_limit=args.match_limit,
            )
            if args.concise:
                print(_concise_summary(view.step_assessment, {}, recipe=view))
            else:
                print(json.dumps(view.to_dict(), indent=2, ensure_ascii=True) if args.json else _recipe_summary(view))
            return 0
        result = assess_step_plausibility(
            args.reaction_smiles, index, shared_index,
            candidate_limit=args.candidate_limit, match_limit=args.match_limit,
        )
        citations = catalog.references(
            m.precedent.reference_id for m in result.matches if m.precedent.reference_id
        ) if catalog is not None and result.matches else {}
    except (OSError, ValueError, sqlite3.Error) as exc:
        parser.error(str(exc))
    if args.concise:
        print(_concise_summary(result, citations))
    else:
        print(json.dumps(result.to_dict(), indent=2, ensure_ascii=True) if args.json else _summary(result, citations))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

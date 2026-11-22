"""Assess precedent support for a proposed reaction using a bound local index pair."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sqlite3

from .generic_indexing import load_generic_index
from .processed_catalog import ProcessedCatalog
from .processed_release import DEFAULT_PROCESSED_ROOT, resolve_processed_release
from .shared_core_index import load_shared_core_index
from .step_plausibility import StepPlausibilityAssessment, assess_step_plausibility


def _summary(assessment: StepPlausibilityAssessment, citations: dict[str, dict]) -> str:
    counts = assessment.support_counts
    rank = "ungraded" if assessment.evidence_rank is None else f"{assessment.evidence_rank}/4"
    lines = [
        f"Precedent support: {assessment.grade} ({rank})",
        f"Strongest graph match: {assessment.strongest_level or 'none'}",
        f"Qualified observations: {assessment.qualified_count} of {assessment.candidate_count} visited",
        f"References among qualified observations: {counts.distinct_reference_count}",
        f"Close observations covering the query context: {counts.context_covered_close_observation_count}",
        f"Measured evidence: {assessment.measured_evidence_status}",
        f"Outcomes in {counts.comparable_observation_count} comparable observations: "
        f"{counts.positive_yield_count} above {assessment.policy.positive_yield_strictly_above_pct:g}%, "
        f"{counts.low_yield_count} low but nonzero, {counts.zero_yield_count} zero, "
        f"{counts.unreported_yield_count} unreported, {counts.unusable_yield_count} unusable",
        f"Positive source groups contributing to the rank: {counts.positive_source_group_count}",
        f"Distinct positive source reactions: {counts.positive_distinct_reaction_count}",
        f"Candidate search truncated: {assessment.candidate_truncated}",
        f"Displayed matches truncated: {assessment.matches_truncated}",
    ]
    if assessment.query_context is not None and assessment.query_context.requirements:
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


def main(argv: list[str] | None = None) -> int:
    """Print an interpretable summary or the complete versioned JSON assessment."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reaction_smiles", help="Proposed reaction, e.g. 'precursor1.precursor2>>product'")
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--release", type=Path, help="Processed dataset root or pinned release manifest (default: active local release)")
    source.add_argument("--index", type=Path, help="Explicit canonical reaction index; also supply --shared-index")
    parser.add_argument("--shared-index", type=Path, help="Shared reaction-core index bound to --index")
    parser.add_argument("--candidate-limit", type=int, default=128, help="Maximum candidates visited")
    parser.add_argument("--match-limit", type=int, default=20, help="Maximum matching source records displayed")
    parser.add_argument("--json", action="store_true", help="Print complete machine-readable evidence")
    args = parser.parse_args(argv)
    if (args.index is None) != (args.shared_index is None):
        parser.error("--index and --shared-index must be supplied together")
    try:
        catalog = None
        if args.index is None:
            release = resolve_processed_release(args.release or DEFAULT_PROCESSED_ROOT)
            index_path = release.artifact("condition_index")
            shared_path = release.artifact("shared_core_index")
            if not args.json and "evidence_catalog" in release.manifest["artifacts"]:
                catalog = ProcessedCatalog(release.artifact("evidence_catalog"))
        else:
            index_path, shared_path = args.index, args.shared_index
        index = load_generic_index(index_path)
        shared_index = load_shared_core_index(shared_path, index)
        result = assess_step_plausibility(
            args.reaction_smiles, index, shared_index,
            candidate_limit=args.candidate_limit, match_limit=args.match_limit,
        )
        citations = catalog.references(
            m.precedent.reference_id for m in result.matches if m.precedent.reference_id
        ) if catalog is not None and result.matches else {}
    except (OSError, ValueError, sqlite3.Error) as exc:
        parser.error(str(exc))
    print(json.dumps(result.to_dict(), indent=2, ensure_ascii=True) if args.json else _summary(result, citations))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

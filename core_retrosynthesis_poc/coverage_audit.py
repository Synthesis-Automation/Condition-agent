"""Stage-wise held-out coverage accounting for operator libraries."""

from __future__ import annotations

import gzip
import io
import json
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable

from .generic_compiler import analyze_generic_reaction, compile_generic_templates
from .generic_models import GenericTemplateLibrary
from .generic_search import disconnect_generic_target_detailed
from .retrieval_index import indexed_template_ids


def _rank(values: Iterable[str], expected: str) -> int | None:
    for rank, value in enumerate(values, start=1):
        if expected and value == expected:
            return rank
    return None


def _preferred_template(templates: Iterable[Any]) -> Any:
    level_rank = {"L2": 0, "L1": 1, "L0": 2, "RDCHIRAL": 3}
    return min(
        templates,
        key=lambda template: (
            level_rank.get(template.abstraction_level, 9),
            template.template_id,
        ),
    )


def _failure_stage(
    *,
    operator_present: bool,
    expected_operator_indexed: bool,
    expected_diagnostics: Any,
    operator_rank: int | None,
) -> str:
    if not operator_present:
        return "operator_absent_from_library"
    if not expected_operator_indexed:
        return "expected_operator_not_retrieved_by_index"
    if expected_diagnostics.product_query_match_count == 0:
        return "expected_operator_no_product_query_match"
    if expected_diagnostics.generated_precursor_count == 0:
        return "expected_operator_no_precursor_generated"
    if expected_diagnostics.valid_candidate_count == 0:
        return "expected_operator_all_precursors_failed_validation"
    if operator_rank is None:
        return "correct_operator_lost_by_global_ranking"
    return "operator_recovered"


def audit_operator_library_coverage(
    rows: Iterable[Dict[str, Any]],
    library: GenericTemplateLibrary,
    output_directory: str | Path,
    *,
    top_k: int = 25,
    max_templates_to_apply: int = 500,
    max_candidates_to_validate: int = 100,
) -> Dict[str, Any]:
    """Audit held-out rows and attribute every miss to one pipeline stage."""

    if top_k < 1:
        raise ValueError("top-k must be positive")
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    case_path = output / "coverage_cases.jsonl.gz"
    stage_counts: Counter[str] = Counter()
    rejection_counts: Counter[str] = Counter()
    totals: Counter[str] = Counter()
    recall_counts: Counter[str] = Counter()
    operators = {operator.operator_id for operator in library.operators}
    operator_template_ids: dict[str, set[str]] = {}
    for template in library.templates:
        operator_template_ids.setdefault(template.operator_id, set()).add(
            template.template_id
        )

    with case_path.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8") as handle:
                for row in rows:
                    totals["source_targets"] += 1
                    compilation = compile_generic_templates(
                        row,
                        levels=("L0", "L1", "L2"),
                        admission_mode="data_driven",
                    )
                    if not compilation.templates:
                        stage = f"source_not_compilable:{compilation.rejection_stage}"
                        stage_counts[stage] += 1
                        rejection_counts[
                            str(compilation.rejection_reason or "unknown_rejection")
                        ] += 1
                        case = {
                            "reaction_id": str(row.get("reaction_id") or ""),
                            "status": stage,
                            "rejection_reason": compilation.rejection_reason,
                            "compilation_diagnostics": compilation.diagnostics,
                        }
                        handle.write(
                            json.dumps(case, sort_keys=True, separators=(",", ":"))
                            + "\n"
                        )
                        continue

                    totals["compilable_targets"] += 1
                    template = _preferred_template(compilation.templates)
                    precedent = template.precedents[0]
                    identity = analyze_generic_reaction(
                        precedent.mapped_reaction_smiles
                    )
                    if identity is None:
                        stage_counts["source_identity_unavailable"] += 1
                        continue
                    indexed = (
                        indexed_template_ids(
                            precedent.product_smiles,
                            library.retrieval_index,
                        )
                        if library.retrieval_index is not None
                        else frozenset(
                            item.template_id for item in library.templates
                        )
                    )
                    operator_present = template.operator_id in operators
                    expected_operator_indexed = bool(
                        indexed.intersection(
                            operator_template_ids.get(template.operator_id, set())
                        )
                    )
                    candidates, diagnostics = disconnect_generic_target_detailed(
                        precedent.product_smiles,
                        library,
                        levels=("L1", "L2"),
                        top_k=top_k,
                        max_templates_to_apply=max_templates_to_apply,
                        max_candidates_to_validate=max_candidates_to_validate,
                    )
                    _, expected_diagnostics = disconnect_generic_target_detailed(
                        precedent.product_smiles,
                        library,
                        operator_ids=(template.operator_id,),
                        levels=("L1", "L2"),
                        top_k=top_k,
                        max_templates_to_apply=max_templates_to_apply,
                        max_candidates_to_validate=max_candidates_to_validate,
                    )
                    totals["targets_with_candidates"] += int(bool(candidates))
                    totals["candidate_count"] += len(candidates)
                    for key, value in diagnostics.to_dict().items():
                        totals[f"search_{key}"] += int(value)
                    ranks = {
                        "exact": _rank(
                            (
                                candidate.precursor_smiles
                                for candidate in candidates
                            ),
                            precedent.precursor_smiles,
                        ),
                        "synthon": _rank(
                            (
                                candidate.synthon_signature
                                for candidate in candidates
                            ),
                            identity.synthon_signature,
                        ),
                        "operator": _rank(
                            (candidate.operator_id for candidate in candidates),
                            template.operator_id,
                        ),
                        "site": _rank(
                            (
                                candidate.disconnection_site_key
                                for candidate in candidates
                            ),
                            identity.disconnection_site_key,
                        ),
                    }
                    for label, rank in ranks.items():
                        for limit in (1, 5, 10, 25):
                            recall_counts[f"top{limit}_{label}"] += int(
                                rank is not None and rank <= limit
                            )
                    stage = _failure_stage(
                        operator_present=operator_present,
                        expected_operator_indexed=expected_operator_indexed,
                        expected_diagnostics=expected_diagnostics,
                        operator_rank=ranks["operator"],
                    )
                    stage_counts[stage] += 1
                    case = {
                        "reaction_id": precedent.reaction_id,
                        "reference_id": precedent.reference_id,
                        "target_smiles": precedent.product_smiles,
                        "expected_precursor_smiles": precedent.precursor_smiles,
                        "expected_operator_id": template.operator_id,
                        "expected_synthon_signature": identity.synthon_signature,
                        "expected_site_key": identity.disconnection_site_key,
                        "status": stage,
                        "ranks": ranks,
                        "search_diagnostics": diagnostics.to_dict(),
                        "expected_operator_diagnostics": (
                            expected_diagnostics.to_dict()
                        ),
                    }
                    handle.write(
                        json.dumps(case, sort_keys=True, separators=(",", ":"))
                        + "\n"
                    )

    compilable = max(1, totals["compilable_targets"])
    metrics = {
        "source_compilable_fraction": totals["compilable_targets"]
        / max(1, totals["source_targets"]),
        "target_candidate_coverage": totals["targets_with_candidates"]
        / compilable,
        "mean_candidates_per_compilable_target": totals["candidate_count"]
        / compilable,
    }
    for label in ("exact", "synthon", "operator", "site"):
        for limit in (1, 5, 10, 25):
            metrics[f"top{limit}_{label}_recall"] = (
                recall_counts[f"top{limit}_{label}"] / compilable
            )
    report = {
        "definition_id": "operator_coverage_stage_audit.v1",
        "library_definition": library.definition,
        "top_k": top_k,
        "totals": dict(sorted(totals.items())),
        "metrics": metrics,
        "stage_counts": dict(sorted(stage_counts.items())),
        "source_rejection_counts": dict(sorted(rejection_counts.items())),
        "cases_path": str(case_path.resolve()),
    }
    (output / "coverage_audit.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return report


__all__ = ["audit_operator_library_coverage"]

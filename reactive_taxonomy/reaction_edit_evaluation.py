"""Phase 2 reaction-edit benchmark and blind chemist-review artifacts."""

from __future__ import annotations

import csv
import html
import json
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Sequence

from .reaction_api import featurize_reaction
from .reaction_edits import normalize_reaction_edits

_DEFAULT_BENCHMARK = (
    Path(__file__).parents[1]
    / "benchmarks"
    / "reaction_edits"
    / "benchmark_manifest.v1.json"
)


def load_reaction_edit_benchmark(
    path: str | Path = _DEFAULT_BENCHMARK,
) -> Dict[str, Any]:
    """Load and minimally validate the versioned Phase 2 benchmark."""
    with Path(path).open("r", encoding="utf-8") as handle:
        benchmark = json.load(handle)
    if benchmark.get("schema_version") != "1.0":
        raise ValueError("Unsupported reaction-edit benchmark schema")
    cases = benchmark.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ValueError("Reaction-edit benchmark has no cases")
    case_ids = [str(case.get("case_id") or "") for case in cases]
    if not all(case_ids) or len(case_ids) != len(set(case_ids)):
        raise ValueError("Reaction-edit case IDs must be non-empty and unique")
    partitioned_ids = [
        str(case_id)
        for case_ids_in_partition in (benchmark.get("partitions") or {}).values()
        for case_id in case_ids_in_partition
    ]
    if sorted(partitioned_ids) != sorted(case_ids):
        raise ValueError("Reaction-edit partitions must contain every case exactly once")
    return dict(benchmark)


def _counter(values: Iterable[str]) -> Counter[str]:
    return Counter(str(value) for value in values)


def _expected_counts(case: Mapping[str, Any]) -> Counter[str]:
    if "expected_edit_counts" in case:
        return Counter(
            {
                str(key): int(value)
                for key, value in case["expected_edit_counts"].items()
            }
        )
    return _counter(edit["edit_type"] for edit in case.get("expected_edits") or ())


def _overlap(expected: Counter[str], observed: Counter[str]) -> int:
    return sum(min(count, observed[key]) for key, count in expected.items())


def _edits(result: Any) -> tuple[Any, ...]:
    signature = result.reaction_signature
    return tuple(signature.edits) if signature is not None else ()


def _edit_evidence(result: Any) -> str:
    signature = result.reaction_signature
    return str(signature.evidence_quality if signature is not None else result.evidence_quality)


def _edit_record(edit: Any) -> Dict[str, Any]:
    return {
        "edit_type": edit.edit_type,
        "atom_1_map": edit.atom_1.atom_map_number,
        "atom_2_map": edit.atom_2.atom_map_number if edit.atom_2 else None,
        "atom_1_element": edit.atom_1.element,
        "atom_2_element": edit.atom_2.element if edit.atom_2 else "H",
        "old_order": edit.old_order,
        "new_order": edit.new_order,
    }


def _records_equal(expected: Sequence[Mapping[str, Any]], observed: Sequence[Any]) -> bool:
    expected_tokens = sorted(
        json.dumps(dict(record), sort_keys=True, separators=(",", ":"))
        for record in expected
    )
    observed_tokens = sorted(
        json.dumps(_edit_record(edit), sort_keys=True, separators=(",", ":"))
        for edit in observed
    )
    return expected_tokens == observed_tokens


def _parity_fingerprint(result: Any) -> tuple[tuple[Any, ...], ...]:
    tokens = []
    for edit in _edits(result):
        elements = tuple(
            sorted((edit.atom_1.element, edit.atom_2.element if edit.atom_2 else "H"))
        )
        tokens.append(
            (
                edit.edit_type,
                elements,
                edit.old_order or "NONE",
                edit.new_order or "NONE",
            )
        )
    return tuple(sorted(tokens))


def _result_fingerprint(result: Any) -> str:
    payload = {
        "valid": result.valid,
        "evidence": _edit_evidence(result),
        "edits": [_edit_record(edit) for edit in _edits(result)],
        "warnings": list(result.warnings),
        "exact": bool(
            result.selected_candidate
            and result.selected_candidate.verification == "exact_product_reconstruction"
        ),
    }
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))


def _connectivity_shadow_record(result: Any) -> Dict[str, Any]:
    if not result.valid:
        return {
            "shadow_key": None,
            "edit_component_keys": [],
            "scope_counts": {},
            "bond_transition_count": 0,
            "hydrogen_delta_count": 0,
            "atom_state_transition_count": 0,
            "compatibility_parity": True,
            "warnings": [],
        }
    normalized = normalize_reaction_edits(
        result.reactants,
        result.products,
        result.selected_candidate,
        result.selected_events,
        result.candidates,
    )
    graph = normalized.connectivity_edit_graph
    if graph is None:
        return {
            "shadow_key": None,
            "edit_component_keys": [],
            "scope_counts": {},
            "bond_transition_count": 0,
            "hydrogen_delta_count": 0,
            "atom_state_transition_count": 0,
            "compatibility_parity": not _edits(result),
            "warnings": [],
        }
    scopes = Counter(
        transition.observation_scope for transition in graph.bond_transitions
    )
    scopes.update(delta.observation_scope for delta in graph.hydrogen_deltas)
    scopes.update(
        transition.observation_scope
        for transition in graph.atom_state_transitions
    )
    normalized_fingerprint = tuple(
        sorted(
            (
                edit.edit_type,
                tuple(
                    sorted(
                        (
                            edit.atom_1.element,
                            edit.atom_2.element if edit.atom_2 else "H",
                        )
                    )
                ),
                edit.old_order or "NONE",
                edit.new_order or "NONE",
            )
            for edit in normalized.edits
        )
    )
    return {
        "shadow_key": graph.shadow_key,
        "edit_component_keys": list(graph.edit_component_keys),
        "scope_counts": dict(sorted(scopes.items())),
        "bond_transition_count": len(graph.bond_transitions),
        "hydrogen_delta_count": len(graph.hydrogen_deltas),
        "atom_state_transition_count": len(graph.atom_state_transitions),
        "compatibility_parity": normalized_fingerprint
        == _parity_fingerprint(result),
        "warnings": list(graph.warnings),
    }


def _draw_molecule(
    smiles: str,
    *,
    highlight_indices: Iterable[int] = (),
    highlight_maps: Iterable[int] = (),
    color: tuple[float, float, float],
) -> str:
    from rdkit.Chem.Draw import rdMolDraw2D
    from .chemistry.rdkit_utils import parse_smiles

    molecule = parse_smiles(smiles)
    if molecule is None:
        return "<p class='invalid'>Invalid structure</p>"
    indices = {int(index) for index in highlight_indices}
    maps = {int(value) for value in highlight_maps}
    indices.update(
        atom.GetIdx()
        for atom in molecule.GetAtoms()
        if atom.GetAtomMapNum() in maps
    )
    highlighted = sorted(indices)
    drawer = rdMolDraw2D.MolDraw2DSVG(360, 220)
    drawer.drawOptions().addAtomIndices = True
    drawer.DrawMolecule(
        molecule,
        highlightAtoms=highlighted,
        highlightAtomColors={index: color for index in highlighted},
    )
    drawer.FinishDrawing()
    return drawer.GetDrawingText().replace(
        "<?xml version='1.0' encoding='iso-8859-1'?>", ""
    )


def _reaction_diagram(result: Any) -> str:
    reactant_indices: Dict[int, set[int]] = {}
    product_maps: set[int] = set()
    for edit in _edits(result):
        for atom in (edit.atom_1, edit.atom_2):
            if atom is None:
                continue
            reactant_indices.setdefault(atom.component_index, set()).add(atom.atom_index)
            if atom.atom_map_number is not None:
                product_maps.add(atom.atom_map_number)
    reactants = "".join(
        "<div class='molecule'><code>"
        + html.escape(component.input_smiles)
        + "</code>"
        + _draw_molecule(
            component.input_smiles,
            highlight_indices=reactant_indices.get(component.component_index, ()),
            color=(1.0, 0.55, 0.15),
        )
        + "</div>"
        for component in result.reactants
    )
    products = "".join(
        "<div class='molecule'><code>"
        + html.escape(component.input_smiles)
        + "</code>"
        + _draw_molecule(
            component.input_smiles,
            highlight_maps=product_maps,
            color=(0.35, 0.75, 0.45),
        )
        + "</div>"
        for component in result.products
    )
    return (
        "<div class='reaction'><div><h3>Reactants</h3>"
        + reactants
        + "</div><div class='arrow'>→</div><div><h3>Products</h3>"
        + products
        + "</div></div>"
    )


def _review_html(rows: Sequence[Mapping[str, Any]]) -> str:
    cards = []
    for row in rows:
        cards.append(
            "<section><h2>"
            + html.escape(str(row["case_id"]))
            + "</h2><code>"
            + html.escape(str(row["reaction_smiles"]))
            + "</code>"
            + str(row["diagram"])
            + "<p><b>Orange:</b> reactant edit atoms; <b>green:</b> mapped product atoms.</p>"
            + "<p><b>Reaction label:</b> "
            + html.escape(str(row["reaction_label"]))
            + "</p><p><b>Structured-label status:</b> "
            + html.escape(str(row["display_label_status"]))
            + "</p><p><b>Observed edits:</b> "
            + html.escape(str(row["edit_descriptions"]))
            + "</p><p><b>Edit evidence:</b> "
            + html.escape(str(row["edit_evidence"]))
            + "</p><p><b>Warnings:</b> "
            + html.escape(str(row["warnings"]))
            + "</p></section>"
        )
    return """<!doctype html>
<html><head><meta charset="utf-8"><title>Phase 2 reaction-edit review</title>
<style>
body { font-family: sans-serif; margin: 2rem; color: #17202a; }
section { border: 1px solid #ccd1d1; border-radius: 8px; margin: 1rem 0;
  padding: 1rem; break-inside: avoid; }
.reaction { display: grid; grid-template-columns: minmax(0, 1fr) 3rem minmax(0, 1fr);
  gap: 1rem; align-items: center; }
.arrow { font-size: 2rem; text-align: center; }
.molecule { display: inline-block; vertical-align: top; max-width: 370px; }
svg { max-width: 360px; display: block; margin: 0.5rem 0; }
code { background: #f4f6f7; padding: 0.2rem 0.4rem; overflow-wrap: anywhere; }
.invalid { color: #a93226; font-weight: bold; }
</style></head><body><h1>Phase 2 reaction-edit chemist review</h1>
<p>This report intentionally omits benchmark expectations. Compare highlighted
reaction centers and edit descriptions with the supplied products. Review bond
formation, cleavage, order changes, hydrogen changes, evidence, and warnings.</p>""" + "".join(cards) + "</body></html>\n"


def evaluate_reaction_edits(
    output_dir: str | Path,
    *,
    benchmark_path: str | Path = _DEFAULT_BENCHMARK,
) -> Dict[str, Any]:
    """Run Phase 2 evaluation and write machine and chemist-review artifacts."""
    benchmark = load_reaction_edit_benchmark(benchmark_path)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    partition_by_case = {
        str(case_id): str(partition)
        for partition, case_ids in benchmark["partitions"].items()
        for case_id in case_ids
    }
    expected_total: Counter[str] = Counter()
    observed_total: Counter[str] = Counter()
    true_positive = 0
    detail_passed = detail_total = 0
    evidence_passed = resolution_passed = reconstruction_passed = 0
    valence_passed = valence_total = 0
    conflict_passed = conflict_total = 0
    deterministic = 0
    case_results = []
    shadow_results = []
    review_rows = []
    parity_results: Dict[str, list[tuple[tuple[Any, ...], ...]]] = {}
    partition_counts: Dict[str, Counter[str]] = {
        name: Counter(total=0, passed=0) for name in benchmark["partitions"]
    }
    for case in benchmark["cases"]:
        reaction_smiles = str(case["reaction_smiles"])
        result = featurize_reaction(reaction_smiles)
        repeated = featurize_reaction(reaction_smiles)
        deterministic_case = _result_fingerprint(result) == _result_fingerprint(repeated)
        deterministic += int(deterministic_case)
        edits = _edits(result)
        shadow = _connectivity_shadow_record(result)
        expected_counts = _expected_counts(case)
        observed_counts = _counter(edit.edit_type for edit in edits)
        expected_total.update(expected_counts)
        observed_total.update(observed_counts)
        true_positive += _overlap(expected_counts, observed_counts)
        counts_passed = expected_counts == observed_counts
        if "expected_edits" in case:
            detail_total += 1
            detail_case_passed = _records_equal(case["expected_edits"], edits)
            detail_passed += int(detail_case_passed)
        else:
            detail_case_passed = True
        observed_evidence = _edit_evidence(result)
        evidence_case_passed = observed_evidence == case["expected_evidence"]
        evidence_passed += int(evidence_case_passed)
        resolved = bool(edits)
        resolution_case_passed = resolved == bool(case["expected_edit_resolved"])
        resolution_passed += int(resolution_case_passed)
        exact = bool(
            result.selected_candidate
            and result.selected_candidate.verification == "exact_product_reconstruction"
        )
        reconstruction_case_passed = exact == bool(case["expected_exact_reconstruction"])
        reconstruction_passed += int(reconstruction_case_passed)
        valid_case_passed = result.valid == bool(case["expected_valid"])
        if not case["expected_valid"]:
            valence_total += 1
            valence_passed += int(not result.valid)
        expected_warnings = set(case.get("expected_warning_contains") or ())
        warnings_case_passed = expected_warnings <= set(result.warnings)
        if case["expected_evidence"] == "conflicting_edit_evidence":
            conflict_total += 1
            conflict_passed += int(
                observed_evidence == "conflicting_edit_evidence"
                and "MAPPING_RECONSTRUCTION_CONFLICT" in result.warnings
            )
        parity_group = case.get("parity_group")
        if parity_group:
            parity_results.setdefault(str(parity_group), []).append(
                _parity_fingerprint(result)
            )
        case_passed = all(
            (
                counts_passed,
                detail_case_passed,
                evidence_case_passed,
                resolution_case_passed,
                reconstruction_case_passed,
                valid_case_passed,
                warnings_case_passed,
                deterministic_case,
            )
        )
        partition = partition_by_case[str(case["case_id"])]
        partition_counts[partition]["total"] += 1
        partition_counts[partition]["passed"] += int(case_passed)
        descriptions = (
            "; ".join(
                clause.detailed for clause in result.display_label.clauses
            )
            if result.display_label is not None and result.display_label.clauses
            else "none"
        )
        case_results.append(
            {
                "case_id": case["case_id"],
                "partition": partition,
                "reaction_smiles": reaction_smiles,
                "valid": result.valid,
                "expected_edit_counts": dict(expected_counts),
                "observed_edit_counts": dict(observed_counts),
                "observed_edits": [_edit_record(edit) for edit in edits],
                "expected_evidence": case["expected_evidence"],
                "observed_evidence": observed_evidence,
                "expected_exact_reconstruction": case["expected_exact_reconstruction"],
                "observed_exact_reconstruction": exact,
                "warnings": list(result.warnings),
                "case_passed": case_passed,
                "error": result.error,
            }
        )
        shadow_results.append(
            {
                "case_id": case["case_id"],
                "partition": partition,
                **shadow,
            }
        )
        review_rows.append(
            {
                "case_id": case["case_id"],
                "partition": partition,
                "reaction_smiles": reaction_smiles,
                "reaction_label": result.reaction_label or "",
                "display_label_status": (
                    result.display_label.status if result.display_label else "unavailable"
                ),
                "edit_descriptions": descriptions,
                "edit_evidence": observed_evidence,
                "warnings": "; ".join(result.warnings),
                "reviewer_id": "",
                "reaction_center_correct": "",
                "bond_edits_correct": "",
                "hydrogen_changes_correct": "",
                "evidence_assessment_correct": "",
                "missing_edits": "",
                "incorrect_edits": "",
                "comments": "",
                "diagram": _reaction_diagram(result),
            }
        )
    case_count = len(case_results)

    def ratio(numerator: int, denominator: int) -> float:
        return round(numerator / denominator, 6) if denominator else 1.0

    parity_group_passes = [
        len(fingerprints) >= 2 and len(set(fingerprints)) == 1
        for fingerprints in parity_results.values()
    ]
    shadow_eligible = tuple(
        row
        for row, case in zip(shadow_results, case_results)
        if case["observed_edits"]
    )
    metrics = {
        "case_count": case_count,
        "passed_case_count": sum(case["case_passed"] for case in case_results),
        "critical_case_pass_rate": ratio(
            sum(case["case_passed"] for case in case_results), case_count
        ),
        "edit_precision": ratio(true_positive, sum(observed_total.values())),
        "edit_recall": ratio(true_positive, sum(expected_total.values())),
        "edit_detail_accuracy": ratio(detail_passed, detail_total),
        "evidence_accuracy": ratio(evidence_passed, case_count),
        "edit_resolution_accuracy": ratio(resolution_passed, case_count),
        "product_reconstruction_accuracy": ratio(reconstruction_passed, case_count),
        "valence_rejection_accuracy": ratio(valence_passed, valence_total),
        "mapped_unmapped_parity_rate": ratio(
            sum(parity_group_passes), len(parity_group_passes)
        ),
        "conflict_retention_rate": ratio(conflict_passed, conflict_total),
        "determinism_rate": ratio(deterministic, case_count),
        "connectivity_shadow_coverage": ratio(
            sum(bool(row["shadow_key"]) for row in shadow_eligible),
            len(shadow_eligible),
        ),
        "connectivity_compatibility_parity_rate": ratio(
            sum(bool(row["compatibility_parity"]) for row in shadow_results),
            len(shadow_results),
        ),
    }
    thresholds = dict(benchmark["thresholds"])
    threshold_results = {
        name: metrics[name] >= float(threshold)
        for name, threshold in thresholds.items()
    }
    report: Dict[str, Any] = {
        "schema_version": "1.1",
        "benchmark_id": benchmark["benchmark_id"],
        "benchmark_path": str(Path(benchmark_path)),
        "metrics": metrics,
        "partition_metrics": {
            partition: {
                "case_count": counts["total"],
                "passed_case_count": counts["passed"],
                "case_pass_rate": ratio(counts["passed"], counts["total"]),
            }
            for partition, counts in sorted(partition_counts.items())
        },
        "thresholds": thresholds,
        "threshold_results": threshold_results,
        "machine_gate_passed": all(threshold_results.values()),
        "human_gate_status": "pending_chemist_review",
        "failed_case_ids": [
            case["case_id"] for case in case_results if not case["case_passed"]
        ],
        "connectivity_shadow": {
            "schema_version": "1.0",
            "scope_counts": dict(
                sorted(
                    sum(
                        (
                            Counter(row["scope_counts"])
                            for row in shadow_results
                        ),
                        Counter(),
                    ).items()
                )
            ),
            "unsupported_bond_domain_count": sum(
                "UNSUPPORTED_BOND_DOMAIN" in row["warnings"]
                for row in shadow_results
            ),
            "canonicalization_overflow_count": sum(
                "CONNECTIVITY_CANONICALIZATION_OVERFLOW"
                in row["warnings"]
                for row in shadow_results
            ),
        },
    }
    (destination / "machine_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (destination / "case_results.jsonl").write_text(
        "".join(
            json.dumps(case, ensure_ascii=False, sort_keys=True) + "\n"
            for case in case_results
        ),
        encoding="utf-8",
    )
    (destination / "connectivity_shadow_report.json").write_text(
        json.dumps(
            {
                "schema_version": "1.0",
                "benchmark_id": benchmark["benchmark_id"],
                "summary": report["connectivity_shadow"],
                "cases": shadow_results,
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    review_fields = [key for key in review_rows[0] if key != "diagram"]
    with (destination / "chemist_review.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=review_fields)
        writer.writeheader()
        writer.writerows(
            {key: value for key, value in row.items() if key != "diagram"}
            for row in review_rows
        )
    (destination / "review_structures.html").write_text(
        _review_html(review_rows), encoding="utf-8"
    )
    with (destination / "disagreements.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(
            (
                "case_id",
                "reviewer_id",
                "category",
                "machine_output",
                "chemist_assessment",
                "resolution_status",
                "resolution_notes",
                "regression_test_id",
            )
        )
    return report


__all__ = ["evaluate_reaction_edits", "load_reaction_edit_benchmark"]

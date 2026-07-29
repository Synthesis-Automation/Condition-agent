"""Phase 1 molecular-feature benchmark and chemist-review artifacts."""

from __future__ import annotations

import csv
import html
import json
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Sequence

from .api import featurize_molecule
from .validation import validate_taxonomy

_DEFAULT_BENCHMARK = (
    Path(__file__).parents[1]
    / "benchmarks"
    / "molecular_features"
    / "benchmark_manifest.v1.json"
)


def load_molecular_feature_benchmark(
    path: str | Path = _DEFAULT_BENCHMARK,
) -> Dict[str, Any]:
    """Load and minimally validate the versioned Phase 1 benchmark."""
    with Path(path).open("r", encoding="utf-8") as handle:
        benchmark = json.load(handle)
    if benchmark.get("schema_version") != "1.0":
        raise ValueError("Unsupported molecular-feature benchmark schema")
    cases = benchmark.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ValueError("Molecular-feature benchmark has no cases")
    case_ids = [str(case.get("case_id") or "") for case in cases]
    if not all(case_ids) or len(case_ids) != len(set(case_ids)):
        raise ValueError("Benchmark case IDs must be non-empty and unique")
    partitions = benchmark.get("partitions") or {}
    partitioned_ids = [
        str(case_id)
        for partition in partitions.values()
        for case_id in partition
    ]
    if sorted(partitioned_ids) != sorted(case_ids):
        raise ValueError("Benchmark partitions must contain every case exactly once")
    atom_case_ids = set((benchmark.get("reactive_atom_expectations") or {}).keys())
    if not atom_case_ids <= set(case_ids):
        raise ValueError("Reactive-atom expectations reference unknown cases")
    return dict(benchmark)


def _counter(values: Iterable[str]) -> Counter[str]:
    return Counter(str(value) for value in values)


def _expected_counter(value: Mapping[str, Any]) -> Counter[str]:
    return Counter({str(key): int(count) for key, count in value.items()})


def _overlap(expected: Counter[Any], observed: Counter[Any]) -> int:
    return sum(min(count, observed[key]) for key, count in expected.items())


def _environment_by_signature(result: Any) -> Dict[str, list[Any]]:
    environments = {item.site_id: item for item in result.site_environments}
    grouped: Dict[str, list[Any]] = {}
    for site in result.sites:
        environment = environments.get(site.site_id)
        if environment is not None:
            grouped.setdefault(site.canonical_signature, []).append(environment)
    return grouped


def _benchmark_environment_projection(
    environment: Any,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Project v1 benchmark expectations from the canonical typed profile."""
    profile = environment.reactivity_profile
    if profile is None:
        return {}, {}
    context = profile.context
    if profile.context_kind == "aromatic":
        steric_class = (
            "ortho_hindered"
            if context.ortho_occupancy_count >= 2
            else "ortho_substituted"
            if context.ortho_occupancy_count == 1
            else "open"
        )
        steric = {
            "class": steric_class,
            "ortho_substituent_count": context.ortho_occupancy_count,
        }
    elif profile.context_kind == "alkyl":
        steric = {"class": context.carbon_substitution}
    elif profile.context_kind == "heteroatom":
        steric = {
            "class": profile.steric.accessibility_class,
            "center_substitution_class": context.substitution_class,
            "attached_groups": [
                {
                    "context": group.context,
                    "attachment_carbon_class": group.attachment_carbon_class,
                    "alpha_branched": group.alpha_branched,
                }
                for group in context.attached_groups
            ],
        }
    else:
        steric = {"class": profile.steric.accessibility_class}
    activation = profile.electronic.activation_class
    electronic_class = {
        "electron_rich": "electron_rich",
        "slightly_rich": "electron_rich",
        "balanced": "neutral",
        "slightly_poor": "electron_poor",
        "electron_poor": "electron_poor",
        "high": "neutral",
        "medium": "neutral",
        "low": "neutral",
    }.get(activation, activation)
    return steric, {"class": electronic_class}


def _check_environment(
    expectation: Mapping[str, Any], grouped: Mapping[str, Sequence[Any]]
) -> tuple[bool, list[str]]:
    signature = str(expectation["site_signature"])
    candidates = list(grouped.get(signature) or ())
    if not candidates:
        return False, [f"missing environment for {signature}"]
    failures = []
    for environment in candidates:
        candidate_failures = []
        steric, electronic = _benchmark_environment_projection(environment)
        for key, expected in (expectation.get("steric_equals") or {}).items():
            if steric.get(key) != expected:
                candidate_failures.append(
                    f"steric.{key}: expected {expected!r}, observed {steric.get(key)!r}"
                )
        for key, minimum in (expectation.get("steric_minimum") or {}).items():
            observed = steric.get(key)
            if observed is None or float(observed) < float(minimum):
                candidate_failures.append(
                    f"steric.{key}: expected >= {minimum}, observed {observed!r}"
                )
        expected_attached = sorted(expectation.get("attached_group_classes") or ())
        if expected_attached:
            observed_attached = sorted(
                str(group.get("attachment_carbon_class"))
                for group in steric.get("attached_groups") or ()
                if group.get("attachment_carbon_class")
            )
            if observed_attached != expected_attached:
                candidate_failures.append(
                    f"attached groups: expected {expected_attached}, observed {observed_attached}"
                )
        if expectation.get("attached_alpha_branched") is not None:
            observed_branched = any(
                bool(group.get("alpha_branched"))
                for group in steric.get("attached_groups") or ()
            )
            if observed_branched != bool(expectation["attached_alpha_branched"]):
                candidate_failures.append(
                    "attached alpha branching did not match expectation"
                )
        expected_electronic = expectation.get("electronic_class")
        if expected_electronic and electronic.get("class") != expected_electronic:
            candidate_failures.append(
                "electronic.class: expected "
                f"{expected_electronic!r}, observed {electronic.get('class')!r}"
            )
        nearby = {str(group.get("group_id")) for group in environment.nearby_groups}
        missing_nearby = set(expectation.get("nearby_group_ids_present") or ()) - nearby
        if missing_nearby:
            candidate_failures.append(
                f"missing nearby groups: {sorted(missing_nearby)}"
            )
        if not candidate_failures:
            return True, []
        failures = candidate_failures
    return False, failures


def _feature_fingerprint(result: Any) -> Dict[str, Any]:
    environments = {item.site_id: item for item in result.site_environments}
    environment_tokens = []
    for site in result.sites:
        environment = environments.get(site.site_id)
        if environment is None:
            continue
        profile = environment.reactivity_profile
        attached = (
            tuple(
                sorted(
                    (
                        group.context,
                        str(group.attachment_carbon_class or ""),
                        bool(group.alpha_branched),
                    )
                    for group in profile.context.attached_groups
                )
            )
            if profile is not None and profile.context_kind == "heteroatom"
            else ()
        )
        nearby = tuple(
            sorted(
                (str(group.get("group_id")), int(group.get("distance", -1)))
                for group in environment.nearby_groups
            )
        )
        steric, electronic = _benchmark_environment_projection(environment)
        environment_tokens.append(
            (
                site.canonical_signature,
                str(steric.get("class") or ""),
                str(steric.get("center_substitution_class") or ""),
                attached,
                str(electronic.get("class") or ""),
                nearby,
            )
        )
    return {
        "groups": sorted(_counter(group.group_id for group in result.functional_groups).items()),
        "sites": sorted(_counter(site.canonical_signature for site in result.sites).items()),
        "environments": sorted(environment_tokens),
    }


def _draw_svg(smiles: str, result: Any) -> str:
    if not result.valid:
        return "<p class='invalid'>Invalid structure</p>"
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return "<p class='invalid'>Structure could not be rendered</p>"
    site_atoms = {
        int(index)
        for site in result.sites
        for index in site.atom_indices
    }
    group_atoms = {
        int(index)
        for group in result.functional_groups
        for index in group.atom_indices
    }
    highlighted = sorted(site_atoms | group_atoms)
    colors = {
        index: (1.0, 0.55, 0.15) if index in site_atoms else (0.35, 0.65, 1.0)
        for index in highlighted
    }
    drawer = rdMolDraw2D.MolDraw2DSVG(500, 300)
    drawer.drawOptions().addAtomIndices = True
    drawer.DrawMolecule(
        molecule,
        highlightAtoms=highlighted,
        highlightAtomColors=colors,
    )
    drawer.FinishDrawing()
    return drawer.GetDrawingText().replace("<?xml version='1.0' encoding='iso-8859-1'?>", "")


def _review_html(rows: Sequence[Mapping[str, Any]]) -> str:
    cards = []
    for row in rows:
        cards.append(
            "<section><h2>"
            + html.escape(str(row["case_id"]))
            + "</h2><code>"
            + html.escape(str(row["smiles"]))
            + "</code>"
            + str(row["svg"])
            + "<p><b>Orange:</b> detected reactive-site atoms; "
            + "<b>blue:</b> other functional-group atoms.</p>"
            + "<p><b>Detected groups:</b> "
            + html.escape(str(row["detected_groups"]))
            + "</p><p><b>Functional-group labels:</b> "
            + html.escape(str(row["detected_group_labels"]))
            + "</p><p><b>Reactive-site labels:</b> "
            + html.escape(str(row["detected_site_labels"]))
            + "</p><p><b>Canonical site signatures:</b> "
            + html.escape(str(row["detected_sites"]))
            + "</p><p><b>Environments:</b> "
            + html.escape(str(row["detected_environments"]))
            + "</p></section>"
        )
    return """<!doctype html>
<html><head><meta charset="utf-8"><title>Phase 1 chemist review</title>
<style>
body { font-family: sans-serif; margin: 2rem; color: #17202a; }
section { border: 1px solid #ccd1d1; border-radius: 8px; margin: 1rem 0;
  padding: 1rem; break-inside: avoid; }
svg { max-width: 500px; display: block; margin: 0.5rem 0; }
code { background: #f4f6f7; padding: 0.2rem 0.4rem; }
.invalid { color: #a93226; font-weight: bold; }
</style></head><body><h1>Phase 1 molecular-feature chemist review</h1>
<p>This report intentionally omits benchmark expectations. Review atom highlights,
functional groups, reactive sites, sterics, electronic context, and missing features
without using the automated answer key.</p>""" + "".join(cards) + "</body></html>\n"


def evaluate_molecular_features(
    output_dir: str | Path,
    *,
    benchmark_path: str | Path = _DEFAULT_BENCHMARK,
) -> Dict[str, Any]:
    """Run Phase 1 evaluation and write machine and chemist-review artifacts."""
    benchmark = load_molecular_feature_benchmark(benchmark_path)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    group_expected = group_observed = group_true_positive = 0
    site_expected = site_observed = site_true_positive = 0
    reactive_atom_total = reactive_atom_passed = 0
    environment_total = environment_passed = 0
    deterministic = invariant_total = invariant_passed = 0
    case_results = []
    review_rows = []
    partition_by_case = {
        str(case_id): str(partition)
        for partition, case_ids in benchmark["partitions"].items()
        for case_id in case_ids
    }
    partition_counts: Dict[str, Counter[str]] = {
        str(partition): Counter() for partition in benchmark["partitions"]
    }
    for case in benchmark["cases"]:
        result = featurize_molecule(str(case["smiles"]))
        repeated = featurize_molecule(str(case["smiles"]))
        deterministic += int(result.to_dict() == repeated.to_dict())
        expected_groups = _expected_counter(case.get("expected_group_counts") or {})
        observed_groups = _counter(group.group_id for group in result.functional_groups)
        expected_sites = _expected_counter(case.get("expected_site_counts") or {})
        observed_sites = _counter(site.canonical_signature for site in result.sites)
        group_expected += sum(expected_groups.values())
        group_observed += sum(observed_groups.values())
        group_true_positive += _overlap(expected_groups, observed_groups)
        site_expected += sum(expected_sites.values())
        site_observed += sum(observed_sites.values())
        site_true_positive += _overlap(expected_sites, observed_sites)
        environment_failures = []
        grouped_environments = _environment_by_signature(result)
        for expectation in case.get("environment_expectations") or ():
            environment_total += 1
            passed, failures = _check_environment(expectation, grouped_environments)
            environment_passed += int(passed)
            environment_failures.extend(failures)
        invariance_failures = []
        if result.valid:
            baseline_fingerprint = _feature_fingerprint(result)
            for equivalent in case.get("equivalent_smiles") or ():
                invariant_total += 1
                equivalent_result = featurize_molecule(str(equivalent))
                passed = (
                    equivalent_result.valid
                    and _feature_fingerprint(equivalent_result) == baseline_fingerprint
                )
                invariant_passed += int(passed)
                if not passed:
                    invariance_failures.append(str(equivalent))
        validity_passed = result.valid == bool(case["expected_valid"])
        error_passed = (
            case.get("expected_error") is None
            or result.error == case.get("expected_error")
        )
        groups_passed = observed_groups == expected_groups
        sites_passed = observed_sites == expected_sites
        atom_failures = []
        observed_atoms: Dict[str, Counter[tuple[int, ...]]] = {}
        for site in result.sites:
            observed_atoms.setdefault(site.canonical_signature, Counter())[
                tuple(sorted(int(index) for index in site.atom_indices))
            ] += 1
        for signature, expected_values in (
            (benchmark.get("reactive_atom_expectations") or {})
            .get(case["case_id"], {})
            .items()
        ):
            expected_atoms = Counter(
                tuple(sorted(int(index) for index in indices))
                for indices in expected_values
            )
            reactive_atom_total += sum(expected_atoms.values())
            matched = _overlap(expected_atoms, observed_atoms.get(signature, Counter()))
            reactive_atom_passed += matched
            if matched != sum(expected_atoms.values()):
                atom_failures.append(
                    f"{signature}: expected atoms {list(expected_atoms.elements())}, "
                    f"observed {list(observed_atoms.get(signature, Counter()).elements())}"
                )
        case_passed = all(
            (
                validity_passed,
                error_passed,
                groups_passed,
                sites_passed,
                not atom_failures,
                not environment_failures,
                not invariance_failures,
            )
        )
        environments_summary = [
            {
                "site_signature": site.canonical_signature,
                "reactivity_profile": (
                    environment.reactivity_profile.to_dict()
                    if environment.reactivity_profile is not None
                    else None
                ),
                "nearby_groups": [
                    {
                        "group_id": group.get("group_id"),
                        "distance": group.get("distance"),
                    }
                    for group in environment.nearby_groups
                ],
            }
            for site in result.sites
            for environment in result.site_environments
            if environment.site_id == site.site_id
            and site.site_type != "aromatic_CH"
        ]
        case_result = {
            "case_id": case["case_id"],
            "partition": partition_by_case[str(case["case_id"])],
            "smiles": case["smiles"],
            "valid": result.valid,
            "expected_group_counts": dict(expected_groups),
            "observed_group_counts": dict(observed_groups),
            "expected_site_counts": dict(expected_sites),
            "observed_site_counts": dict(observed_sites),
            "reactive_atom_failures": atom_failures,
            "environment_failures": environment_failures,
            "invariance_failures": invariance_failures,
            "case_passed": case_passed,
            "error": result.error,
        }
        case_results.append(case_result)
        partition_counts[case_result["partition"]]["total"] += 1
        partition_counts[case_result["partition"]]["passed"] += int(case_passed)
        review_rows.append(
            {
                "case_id": case["case_id"],
                "partition": partition_by_case[str(case["case_id"])],
                "smiles": case["smiles"],
                "canonical_smiles": result.canonical_smiles or "",
                "detected_groups": "; ".join(
                    f"{key} ({count})" for key, count in sorted(observed_groups.items())
                ),
                "detected_group_labels": "; ".join(
                    f"{key} ({count})"
                    for key, count in sorted(
                        _counter(
                            group.chemist_label for group in result.functional_groups
                        ).items()
                    )
                ),
                "detected_sites": "; ".join(
                    f"{key} ({count})" for key, count in sorted(observed_sites.items())
                ),
                "detected_site_labels": "; ".join(
                    f"{key} ({count})"
                    for key, count in sorted(
                        _counter(site.chemist_label for site in result.sites).items()
                    )
                ),
                "detected_environments": json.dumps(
                    environments_summary, ensure_ascii=False, sort_keys=True
                ),
                "reviewer_id": "",
                "reactive_sites_correct": "",
                "functional_groups_correct": "",
                "steric_correct": "",
                "electronic_context_correct": "",
                "missing_features": "",
                "incorrect_features": "",
                "comments": "",
                "svg": _draw_svg(str(case["smiles"]), result),
            }
        )
    case_count = len(case_results)
    taxonomy_errors = validate_taxonomy()

    def ratio(numerator: int, denominator: int) -> float:
        return round(numerator / denominator, 6) if denominator else 1.0

    metrics = {
        "case_count": case_count,
        "passed_case_count": sum(case["case_passed"] for case in case_results),
        "critical_case_pass_rate": ratio(
            sum(case["case_passed"] for case in case_results), case_count
        ),
        "functional_group_precision": ratio(group_true_positive, group_observed),
        "functional_group_recall": ratio(group_true_positive, group_expected),
        "reactive_site_precision": ratio(site_true_positive, site_observed),
        "reactive_site_recall": ratio(site_true_positive, site_expected),
        "reactive_atom_accuracy": ratio(
            reactive_atom_passed, reactive_atom_total
        ),
        "reactive_atom_check_count": reactive_atom_total,
        "environment_check_count": environment_total,
        "environment_check_accuracy": ratio(environment_passed, environment_total),
        "determinism_rate": ratio(deterministic, case_count),
        "equivalent_smiles_check_count": invariant_total,
        "equivalent_smiles_invariance_rate": ratio(invariant_passed, invariant_total),
        "taxonomy_error_count": len(taxonomy_errors),
    }
    thresholds = dict(benchmark["thresholds"])
    threshold_results = {
        name: metrics[name] <= threshold
        if name == "taxonomy_error_count"
        else metrics[name] >= threshold
        for name, threshold in thresholds.items()
    }
    report: Dict[str, Any] = {
        "schema_version": "1.0",
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
        "taxonomy_errors": taxonomy_errors,
        "failed_case_ids": [
            case["case_id"] for case in case_results if not case["case_passed"]
        ],
    }
    (destination / "machine_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    (destination / "case_results.jsonl").write_text(
        "".join(
            json.dumps(case, ensure_ascii=False, sort_keys=True) + "\n"
            for case in case_results
        ),
        encoding="utf-8",
    )
    review_fields = [key for key in review_rows[0] if key != "svg"]
    with (destination / "chemist_review.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=review_fields)
        writer.writeheader()
        writer.writerows(
            {key: value for key, value in row.items() if key != "svg"}
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


__all__ = [
    "evaluate_molecular_features",
    "load_molecular_feature_benchmark",
]

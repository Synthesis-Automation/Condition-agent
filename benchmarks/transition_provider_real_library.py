"""Run one transition provider against a real library and target CSV panel."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import time
from collections import Counter
from pathlib import Path
from typing import Any, Sequence

from core_retrosynthesis import (
    ExpandLeafAction,
    ExpansionLeaf,
    ExpansionState,
    OperatorLadderTransitionProvider,
    TransitionProviderOrchestrator,
    load_generic_library,
)


REPORT_SCHEMA_VERSION = "1.0"


def _integer(value: Any) -> int:
    return int(value) if isinstance(value, (int, float)) else 0


def _percentile(values: Sequence[float], fraction: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    index = min(len(ordered) - 1, int((len(ordered) - 1) * fraction))
    return round(ordered[index], 6)


def run_single_provider_benchmark(
    library_path: str | Path,
    target_csv_path: str | Path,
    *,
    proposal_budget: int = 3,
    maximum_targets: int | None = None,
) -> dict[str, Any]:
    """Run the generic operator provider and return an auditable report."""

    library_source = Path(library_path).resolve()
    target_source = Path(target_csv_path).resolve()
    library = load_generic_library(library_source)
    provider = OperatorLadderTransitionProvider(library)
    if proposal_budget > provider.metadata.maximum_proposal_budget:
        raise ValueError("proposal budget exceeds provider capability")
    if maximum_targets is not None and maximum_targets < 1:
        raise ValueError("maximum targets must be positive")
    orchestrator = TransitionProviderOrchestrator((provider,))

    with target_source.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = tuple(csv.DictReader(handle))
    if maximum_targets is not None:
        rows = rows[:maximum_targets]

    cases: list[dict[str, Any]] = []
    for index, row in enumerate(rows, start=1):
        leaf_id = f"target:{index}"
        started = time.perf_counter()
        try:
            state = ExpansionState(
                state_id=f"single-provider-panel:{index}",
                leaves=(ExpansionLeaf(leaf_id, str(row.get("smiles") or "")),),
            )
            outcome = orchestrator.expand(
                state,
                ExpandLeafAction(
                    leaf_id=leaf_id,
                    provider_id=provider.metadata.provider_id,
                    proposal_budget=proposal_budget,
                ),
            )
        except Exception as error:
            cases.append(
                {
                    "case_index": index,
                    "rank": str(row.get("rank") or ""),
                    "drug_name": str(row.get("drug_name") or f"target_{index}"),
                    "input_smiles": str(row.get("smiles") or ""),
                    "execution_status": "error",
                    "error_type": type(error).__name__,
                    "error": str(error),
                    "elapsed_seconds": round(time.perf_counter() - started, 6),
                }
            )
            continue
        cases.append(
            {
                "case_index": index,
                "rank": str(row.get("rank") or ""),
                "drug_name": str(row.get("drug_name") or f"target_{index}"),
                "input_smiles": str(row.get("smiles") or ""),
                "execution_status": "completed",
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "outcome": outcome.to_dict(),
            }
        )

    completed = tuple(
        case for case in cases if case["execution_status"] == "completed"
    )
    elapsed = tuple(float(case["elapsed_seconds"]) for case in completed)
    rejection_reasons = Counter(
        rejected["reason"]
        for case in completed
        for rejected in case["outcome"]["rejected"]
    )
    levels = Counter(
        level
        for case in completed
        for level in case["outcome"]["provider_diagnostics"].get(
            "levels_attempted", ()
        )
    )
    summary = {
        "target_count": len(cases),
        "completed_target_count": len(completed),
        "error_target_count": len(cases) - len(completed),
        "target_with_raw_candidates_count": sum(
            case["outcome"]["raw_candidate_count"] > 0 for case in completed
        ),
        "target_with_admitted_transitions_count": sum(
            bool(case["outcome"]["admitted"]) for case in completed
        ),
        "raw_candidate_count": sum(
            _integer(case["outcome"]["raw_candidate_count"])
            for case in completed
        ),
        "admitted_transition_count": sum(
            len(case["outcome"]["admitted"]) for case in completed
        ),
        "rejected_transition_count": sum(
            len(case["outcome"]["rejected"]) for case in completed
        ),
        "duplicate_transition_count": sum(
            _integer(case["outcome"]["duplicate_candidate_count"])
            for case in completed
        ),
        "provider_proposed_action_count": sum(
            _integer(
                case["outcome"]["provider_diagnostics"].get(
                    "proposed_action_count"
                )
            )
            for case in completed
        ),
        "provider_validation_attempt_count": sum(
            _integer(
                case["outcome"]["provider_diagnostics"].get(
                    "validation_attempt_count"
                )
            )
            for case in completed
        ),
        "provider_valid_action_count": sum(
            _integer(
                case["outcome"]["provider_diagnostics"].get(
                    "valid_action_count"
                )
            )
            for case in completed
        ),
        "rejection_reasons": dict(sorted(rejection_reasons.items())),
        "level_attempt_counts": dict(sorted(levels.items())),
        "total_elapsed_seconds": round(sum(elapsed), 6),
        "median_elapsed_seconds": (
            round(statistics.median(elapsed), 6) if elapsed else None
        ),
        "p90_elapsed_seconds": _percentile(elapsed, 0.90),
    }
    return {
        "schema_version": REPORT_SCHEMA_VERSION,
        "benchmark_kind": "real_library_single_provider_shadow",
        "library": {
            "path": str(library_source),
            "schema_version": library.schema_version,
            "template_count": len(library.templates),
            "operator_count": len(library.operators),
            "accepted_observation_count": library.accepted_observation_count,
        },
        "target_source": str(target_source),
        "provider": {
            "provider_id": provider.metadata.provider_id,
            "display_name": provider.metadata.display_name,
            "capability_tags": list(provider.metadata.capability_tags),
            "proposal_budget": proposal_budget,
        },
        "summary": summary,
        "cases": cases,
    }


def _write_case_summary(report: dict[str, Any], path: Path) -> None:
    fieldnames = (
        "case_index",
        "rank",
        "drug_name",
        "execution_status",
        "raw_candidate_count",
        "admitted_transition_count",
        "rejected_transition_count",
        "duplicate_transition_count",
        "levels_attempted",
        "provider_proposed_action_count",
        "provider_validation_attempt_count",
        "provider_valid_action_count",
        "rejection_reasons",
        "elapsed_seconds",
        "error",
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for case in report["cases"]:
            outcome = case.get("outcome") or {}
            diagnostics = outcome.get("provider_diagnostics") or {}
            reasons = Counter(
                item["reason"] for item in outcome.get("rejected") or ()
            )
            writer.writerow(
                {
                    "case_index": case["case_index"],
                    "rank": case["rank"],
                    "drug_name": case["drug_name"],
                    "execution_status": case["execution_status"],
                    "raw_candidate_count": outcome.get("raw_candidate_count", 0),
                    "admitted_transition_count": len(outcome.get("admitted") or ()),
                    "rejected_transition_count": len(outcome.get("rejected") or ()),
                    "duplicate_transition_count": outcome.get(
                        "duplicate_candidate_count", 0
                    ),
                    "levels_attempted": ";".join(
                        diagnostics.get("levels_attempted") or ()
                    ),
                    "provider_proposed_action_count": diagnostics.get(
                        "proposed_action_count", 0
                    ),
                    "provider_validation_attempt_count": diagnostics.get(
                        "validation_attempt_count", 0
                    ),
                    "provider_valid_action_count": diagnostics.get(
                        "valid_action_count", 0
                    ),
                    "rejection_reasons": ";".join(
                        f"{reason}:{count}"
                        for reason, count in sorted(reasons.items())
                    ),
                    "elapsed_seconds": case["elapsed_seconds"],
                    "error": case.get("error", ""),
                }
            )


def main(argv: Sequence[str] | None = None) -> int:
    """Run the benchmark CLI and write JSON plus concise CSV artifacts."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("library")
    parser.add_argument("targets")
    parser.add_argument("output_json")
    parser.add_argument("--output-csv")
    parser.add_argument("--proposal-budget", type=int, default=3)
    parser.add_argument("--maximum-targets", type=int)
    arguments = parser.parse_args(argv)

    report = run_single_provider_benchmark(
        arguments.library,
        arguments.targets,
        proposal_budget=arguments.proposal_budget,
        maximum_targets=arguments.maximum_targets,
    )
    output_json = Path(arguments.output_json).resolve()
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(
        json.dumps(report, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    if arguments.output_csv:
        output_csv = Path(arguments.output_csv).resolve()
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        _write_case_summary(report, output_csv)
    print(json.dumps(report["summary"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


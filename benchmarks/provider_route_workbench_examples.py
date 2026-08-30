"""Run provider-backed route workbench examples against real local artifacts."""

from __future__ import annotations

import argparse
import csv
import json
import time
from collections import Counter
from pathlib import Path
from typing import Any, Sequence

from cas_tools import open_stock_lookup
from core_retrosynthesis import (
    OperatorLadderTransitionProvider,
    RouteWorkbenchSettings,
    TransitionProviderOrchestrator,
    load_generic_library,
    run_provider_route_workbench,
)


REPORT_SCHEMA_VERSION = "1.0"


def _case_summary(case: dict[str, Any]) -> dict[str, Any]:
    if case["execution_status"] != "completed":
        return {
            "case_id": case["case_id"],
            "target_name": case["target_name"],
            "execution_status": case["execution_status"],
            "error": case.get("error", ""),
            "elapsed_seconds": case["elapsed_seconds"],
        }
    report = case["report"]
    routes = report["routes"]
    first = routes[0] if routes else None
    route = first["route"] if first else None
    issues = first["issues"] if first else []
    weakest_id = first["weakest_issue_id"] if first else None
    weakest = next(
        (item for item in issues if item["issue_id"] == weakest_id),
        None,
    )
    outcomes = report["provider_outcomes"]
    return {
        "case_id": case["case_id"],
        "target_name": case["target_name"],
        "execution_status": case["execution_status"],
        "route_kind": report["route_kind"],
        "retained_route_count": len(routes),
        "solved_route_count": len(report["search_result"]["routes"]),
        "partial_route_count": len(report["search_result"]["partial_routes"]),
        "best_route_id": route["route_id"] if route else "",
        "best_route_reaction_count": route["reaction_count"] if route else 0,
        "best_route_unresolved_leaf_count": (
            sum(not leaf["terminal"] for leaf in route["leaves"]) if route else 0
        ),
        "best_route_verification": (
            first["verification"]["status"] if first else "unavailable"
        ),
        "best_route_issue_count": len(issues),
        "best_route_strong_issue_count": sum(
            item["severity"] == "strong" for item in issues
        ),
        "weakest_issue_kind": weakest["kind"] if weakest else "",
        "weakest_issue_severity": weakest["severity"] if weakest else "",
        "repair_proposal_count": len(first["repair_proposals"]) if first else 0,
        "actionable_repair_count": (
            sum(
                item["status"] == "actionable"
                for item in first["repair_proposals"]
            )
            if first
            else 0
        ),
        "provider_call_count": len(outcomes),
        "provider_raw_candidate_count": sum(
            item["raw_candidate_count"] for item in outcomes
        ),
        "provider_admitted_transition_count": sum(
            len(item["admitted"]) for item in outcomes
        ),
        "provider_rejected_transition_count": sum(
            len(item["rejected"]) for item in outcomes
        ),
        "elapsed_seconds": case["elapsed_seconds"],
        "error": "",
    }


def run_panel(
    library_path: str | Path,
    stock_path: str | Path,
    targets_path: str | Path,
    *,
    settings: RouteWorkbenchSettings,
) -> dict[str, Any]:
    """Run every CSV target through the provider-backed route workbench."""

    library_source = Path(library_path).resolve()
    stock_source = Path(stock_path).resolve()
    target_source = Path(targets_path).resolve()
    library = load_generic_library(library_source)
    provider = OperatorLadderTransitionProvider(library)
    with target_source.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = tuple(csv.DictReader(handle))

    cases = []
    with open_stock_lookup(stock_source) as stock_index:
        for row in rows:
            started = time.perf_counter()
            try:
                orchestrator = TransitionProviderOrchestrator((provider,))
                result = run_provider_route_workbench(
                    str(row.get("target_smiles") or ""),
                    library,
                    stock_index,
                    orchestrator,
                    provider.metadata.provider_id,
                    settings=settings,
                )
            except Exception as error:
                cases.append(
                    {
                        "case_id": str(row.get("case_id") or ""),
                        "target_name": str(row.get("target_name") or ""),
                        "challenge_focus": str(row.get("challenge_focus") or ""),
                        "target_smiles": str(row.get("target_smiles") or ""),
                        "execution_status": "error",
                        "error_type": type(error).__name__,
                        "error": str(error),
                        "elapsed_seconds": round(
                            time.perf_counter() - started, 6
                        ),
                    }
                )
                continue
            cases.append(
                {
                    "case_id": str(row.get("case_id") or ""),
                    "target_name": str(row.get("target_name") or ""),
                    "challenge_focus": str(row.get("challenge_focus") or ""),
                    "target_smiles": str(row.get("target_smiles") or ""),
                    "execution_status": "completed",
                    "elapsed_seconds": round(time.perf_counter() - started, 6),
                    "report": result.to_dict(),
                }
            )

    summaries = tuple(_case_summary(case) for case in cases)
    route_kinds = Counter(
        item.get("route_kind", "error") for item in summaries
    )
    return {
        "schema_version": REPORT_SCHEMA_VERSION,
        "benchmark_kind": "provider_route_workbench_real_examples",
        "library": {
            "path": str(library_source),
            "template_count": len(library.templates),
            "operator_count": len(library.operators),
            "accepted_observation_count": library.accepted_observation_count,
        },
        "stock_path": str(stock_source),
        "target_source": str(target_source),
        "provider_id": provider.metadata.provider_id,
        "settings": {
            "planner": settings.__dict__,
            "provider": {
                "max_templates_to_apply": provider.max_templates_to_apply,
                "max_candidates_to_validate": (
                    provider.max_candidates_to_validate
                ),
            },
        },
        "summary": {
            "target_count": len(cases),
            "completed_count": sum(
                case["execution_status"] == "completed" for case in cases
            ),
            "error_count": sum(
                case["execution_status"] == "error" for case in cases
            ),
            "route_kind_counts": dict(sorted(route_kinds.items())),
            "target_with_solved_route_count": sum(
                item.get("solved_route_count", 0) > 0 for item in summaries
            ),
            "target_with_actionable_repair_count": sum(
                item.get("actionable_repair_count", 0) > 0
                for item in summaries
            ),
            "total_elapsed_seconds": round(
                sum(float(case["elapsed_seconds"]) for case in cases), 6
            ),
        },
        "case_summaries": summaries,
        "cases": cases,
    }


def _write_summary_csv(report: dict[str, Any], path: Path) -> None:
    summaries = tuple(report["case_summaries"])
    fieldnames = tuple(
        dict.fromkeys(
            key for item in summaries for key in item
        )
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summaries)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the example panel and write JSON plus concise CSV artifacts."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("library")
    parser.add_argument("stock")
    parser.add_argument("targets")
    parser.add_argument("output_json")
    parser.add_argument("--output-csv")
    parser.add_argument("--max-depth", type=int, default=2)
    parser.add_argument("--top-k-routes", type=int, default=2)
    parser.add_argument("--per-step-top-k", type=int, default=3)
    parser.add_argument("--beam-width", type=int, default=8)
    parser.add_argument("--max-expansions", type=int, default=8)
    parser.add_argument("--molecular-weight-threshold", type=float, default=150.0)
    arguments = parser.parse_args(argv)

    report = run_panel(
        arguments.library,
        arguments.stock,
        arguments.targets,
        settings=RouteWorkbenchSettings(
            max_depth=arguments.max_depth,
            molecular_weight_threshold=arguments.molecular_weight_threshold,
            top_k_routes=arguments.top_k_routes,
            per_step_top_k=arguments.per_step_top_k,
            beam_width=arguments.beam_width,
            max_expansions=arguments.max_expansions,
        ),
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
        _write_summary_csv(report, output_csv)
    print(json.dumps(report["summary"], indent=2, sort_keys=True))
    print(json.dumps(report["case_summaries"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


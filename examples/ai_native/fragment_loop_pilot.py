"""Opt-in live development pilot of agent-chosen iterative fragment searching.

Uses the production investigator/runtime and recorded workspace operations. No
fragment, route, or automatic search ladder is supplied. This is an integration
pilot, not a paired quality comparison or an untouched chemistry release gate.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
import json
import math
from pathlib import Path
from typing import Any

from chem_coworker.scientific_workspace import ScientificWorkspace
from chem_coworker.scientific_workspace.agent_context.learning import build_learning_context
from chem_coworker.scientific_workspace.core.baseline import capture_baseline
from chem_coworker.scientific_workspace.runtime.agent_runtime import CodexRuntime
from examples.ai_native.fragment_agent_comparison import ComparisonCase, load_cases, run_trial


def pilot_question(case: ComparisonCase, timeout: float) -> str:
    """Supply a target and bounded task, without prescribing a fragment or route."""
    return f"""Development pilot of iterative, agent-chosen fragment searching.
Target SMILES: {case.target_smiles}
{case.question}
Choose your own connected key fragment and explain the chemical question it tests.
Use the persistent scientific console described in the supplied fragment guide.
Read each result before deciding whether to refine, broaden, follow a precursor,
or stop. Inspect actual reaction structures for promising hits and explain what
the source constructs versus carries through. A returned hit is not itself proof
of a feasible route. Record reasons and parent evidence for consequential revisions.
No minimum query count is required. Do not force a positive result.

Budget: {timeout:g} seconds including answer submission; at most six fragment
searches and twelve scientific calls. Inspections and notes do not consume calls.
Do not use automatic discovery, automatic fragment suggestions, web search,
planners, forward predictors, condition search, or build an index in this pilot.
Use only this investigation and its pinned artifacts. Do not read other trials,
shared lessons, case answers, or the pilot's source code. Repository code is read
only; create any scripts under this investigation. Do not ask for clarification.
Keep time for a bounded partial answer and quit the console when finished.

Answer using the normal scientific answer schema and finalization handoff. Report
the chosen bottleneck, inspected evidence, mismatches and unresolved gaps. Propose
at most two steps if supported; omit invented intermediates. Source-supported
observations and synthesis hypotheses must remain distinct.
"""


def loop_metrics(workspace: ScientificWorkspace) -> dict[str, Any]:
    """Expose query choices, timings and opened fields without scoring chemistry."""
    searches, inspections, decisions = [], [], []
    prohibited = []
    for event in workspace.store.events():
        if event.kind not in {"call", "console_inspection", "decision", "hypothesis", "observation"}:
            continue
        value = workspace.store.read_artifact(event.artifact_ref)
        if event.kind == "console_inspection":
            inspections.append({"ref": event.artifact_ref, "evidence_refs": event.evidence_refs, **value})
        elif event.kind != "call":
            decisions.append({"ref": event.artifact_ref, "kind": event.kind,
                              "evidence_refs": event.evidence_refs, **value})
        elif value.get("operation") == "search_fragment_precedents":
            result = value.get("result") or {}
            searches.append({"ref": event.artifact_ref, "arguments": value.get("arguments"),
                             "evidence_refs": event.evidence_refs,
                             "duration_seconds": value.get("duration_seconds"),
                             "timings": value.get("timings"),
                             "status": value.get("execution_status"), "counts": result.get("counts"),
                             "execution": result.get("execution"), "error": value.get("error"),
                             "hit_count": len(result.get("hits", []))})
        elif value.get("operation") in {"find_synthesis_precedents", "suggest_search_fragments",
                                        "plan_synthesis", "predict_forward", "recommend_conditions"}:
            prohibited.append({"ref": event.artifact_ref, "operation": value["operation"]})
    return {"searches": searches, "console_inspections": inspections, "notes": decisions,
            "search_budget_exceeded": len(searches) > 6, "recorded_prohibited_calls": prohibited,
            "manual_chemistry_review": "pending; opened fields do not establish useful support"}


def main() -> int:
    """Run the configured runtime without changing its model or security settings."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-live", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fragment-index", type=Path, required=True)
    parser.add_argument("--cases", type=Path, default=Path(__file__).with_name("fragment_loop_cases.json"))
    parser.add_argument("--timeout", type=float, default=600)
    args = parser.parse_args()
    if not args.run_live:
        parser.error("Add --run-live to run the configured agent")
    if not math.isfinite(args.timeout) or not 1 <= args.timeout <= 1800:
        parser.error("timeout must be finite and in 1..1800")
    if not args.fragment_index.is_file():
        parser.error("Provide an existing fragment index")
    cases = load_cases(args.cases)
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    repository = Path(__file__).resolve().parents[2]
    baseline = capture_baseline(repository, {"fragment_index": args.fragment_index.resolve()})
    baseline["evaluation_partition"] = "development_fragment_loop_not_an_untouched_evaluation"
    baseline["learning_context"] = build_learning_context(baseline)
    runtime = CodexRuntime(timeout_seconds=args.timeout)
    report: dict[str, Any] = {"schema_version": "fragment_loop_pilot.v1", "trials": [],
                             "limitations": ["Unpaired two-case development pilot; no quality uplift claim.",
                                             "Targets newly authored for this pilot; training novelty unknown.",
                                             "Tool restrictions are instructions, not process isolation.",
                                             "No independent chemistry review or untouched evaluation."]}
    for case in cases:
        row = run_trial(args.output / case.case_id, case, "fragment_assisted", deepcopy(baseline),
                        runtime, args.timeout, question_override=pilot_question(case, args.timeout),
                        task_names=("retrosynthesis_fragments",), call_budget=12)
        row["loop"] = loop_metrics(ScientificWorkspace(args.output / case.case_id))
        report["trials"].append(row)
        (args.output / "pilot.json").write_text(json.dumps(report, indent=2), "utf-8")
    return 0 if all(row["status"] == "completed" for row in report["trials"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())

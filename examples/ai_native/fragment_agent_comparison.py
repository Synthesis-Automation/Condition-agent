"""Opt-in paired agent development runs with isolated evidence and no shared learning.

Uses the production investigator prompt, runtime and answer validators. This is a
prompt-based tool ablation, not process isolation or an independent chemistry gate.
One runtime attempt per arm is used; the conversation service's repair loop is not
part of this evaluation. Chemistry usefulness and unsupported steps need review.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
from threading import Event
from time import perf_counter
from typing import Any

from rdkit import Chem

from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.agent_runtime import AgentRuntime, AgentStopped, CodexRuntime
from chem_coworker.scientific_workspace.answer_contracts import ScientificAnswer, validate_answer_evidence
from chem_coworker.scientific_workspace.baseline import capture_baseline, verify_baseline
from chem_coworker.scientific_workspace.conversation import investigation_prompt
from chem_coworker.scientific_workspace.evidence_review import matching_evidence_review
from chem_coworker.scientific_workspace.learning import build_learning_context
from chem_coworker.scientific_workspace.store import canonical_bytes


FRAGMENT_OPERATIONS = frozenset({"suggest_search_fragments", "search_fragment_precedents"})
ARMS = ("agent_only", "fragment_assisted")
ACTIVE_LIMITATIONS = [
    "Development comparison only; targets are not a held-out chemistry evaluation or proven novel structures.",
    "Tool counts and returned construction witnesses do not establish source inspection, transferability or route quality.",
    "Arm restrictions are prompt-based; recorded prohibited calls are flagged, but unrecorded access cannot be excluded.",
    "Requested model/settings do not confirm the effective provider model; inspect saved runtime observations.",
    "New threads and disabled lesson memory prevent intentional history transfer; both arms still share local files.",
    "No automatic chemistry score or improvement claim is made; manual review fields remain pending.",
]


@dataclass(frozen=True)
class ComparisonCase:
    """Structure-only development question, without a supplied route or core answer."""

    case_id: str
    target_smiles: str
    question: str


def load_cases(path: Path) -> tuple[ComparisonCase, ...]:
    """Validate case IDs and connected structures before starting a paid runtime."""
    data = json.loads(path.read_text("utf-8"))
    if data.get("schema_version") != "fragment_agent_cases.v1":
        raise ValueError("Unsupported comparison case schema")
    cases = tuple(ComparisonCase(**item) for item in data["cases"])
    if not cases or len({c.case_id for c in cases}) != len(cases):
        raise ValueError("Supply nonempty cases with unique case IDs")
    for case in cases:
        if not isinstance(case.case_id, str) or not re.fullmatch(r"[a-z][a-z0-9_]{0,63}", case.case_id):
            raise ValueError("Case IDs must be safe lowercase identifiers")
        if not isinstance(case.question, str) or not 1 <= len(case.question.strip()) <= 3000:
            raise ValueError("Each case needs a bounded question")
        if not isinstance(case.target_smiles, str) or len(case.target_smiles) > 5000:
            raise ValueError("Invalid case target")
        mol = Chem.MolFromSmiles(case.target_smiles)
        if mol is None or not 0 < mol.GetNumAtoms() <= 200 or len(Chem.GetMolFrags(mol)) != 1:
            raise ValueError("Each case needs one connected target with at most 200 atoms")
    return cases


def comparison_question(case: ComparisonCase, arm: str, timeout_seconds: float) -> str:
    """Keep the scientific task and budgets identical apart from fragment-tool access."""
    if arm not in ARMS:
        raise ValueError("Unknown comparison arm")
    access = (
        "For this comparison arm, do not use suggest_search_fragments or search_fragment_precedents, "
        "their underlying APIs, the fragment index, or k-way fragmentation. You may identify cores "
        "yourself and use the other configured research tools. These experiment constraints take "
        "precedence over optional fragment advice in the investigator prompt."
        if arm == "agent_only" else
        "For this comparison arm, the optional suggest_search_fragments and search_fragment_precedents "
        "tools are available. Choose whether they resolve a useful uncertainty; you may choose your own core "
        "and skip suggestions. Do not search every candidate merely to demonstrate tool use."
    )
    return f"""Development comparison, not a request to implement code or provide a complete experimental protocol.
Target SMILES: {case.target_smiles}
{case.question}
You have {timeout_seconds:g} seconds for one runtime attempt including final answer submission.
Keep the investigation small enough to finish: at most four recorded scientific
operation calls, two native web searches and two source opens/fetches. Reading
saved fields and recording notes do not consume scientific-operation calls.
These are matched evaluation limits, not production workflow rules.
Preserve failures. If evidence is unavailable, submit a bounded partial answer.
Do not run a forward predictor, multistep planner, library build or condition search.
Use only this investigation's evidence and configured artifacts. Do not read other
investigations, pilot reports, case answers, shared lessons or the evaluation source.
Do not ask the user to clarify this benchmark question.

{access}

In the final answer explain the chosen bottleneck, any inspected construction
precedent (or absence of support), what could transfer and what would differ.
Show at most two proposed steps when structures and evidence permit; do not invent
intermediates to fill a route. Clearly label hypotheses and unresolved gaps.
Record consequential decisions with existing evidence-linked notes. Use the normal
scientific answer schema and runtime handoff, retaining actual evidence references.
"""


def collect_metrics(workspace: ScientificWorkspace, arm: str, answer: dict[str, Any] | None) -> dict[str, Any]:
    """Count observable events; leave scientific quality and actual inspection to review."""
    calls, constructed, excerpts = [], set(), 0
    for event in workspace.store.events():
        if event.kind == "literature_excerpt":
            excerpts += 1
        if event.kind != "call":
            continue
        value = workspace.store.read_artifact(event.artifact_ref)
        calls.append({"artifact_ref": event.artifact_ref, "operation": value.get("operation"),
                      "execution_status": value.get("execution_status"),
                      "duration_seconds": value.get("duration_seconds"), "error": value.get("error")})
        if value.get("operation") == "search_fragment_precedents":
            for hit in (value.get("result") or {}).get("hits", []):
                if "constructed" in hit.get("relationships", []) and hit.get("observation_id"):
                    constructed.add(hit["observation_id"])
    prohibited = [row for row in calls if arm == "agent_only" and row["operation"] in FRAGMENT_OPERATIONS]
    return {
        "calls": calls, "scientific_call_count": len(calls),
        "failed_call_count": sum(row["execution_status"] != "completed" for row in calls),
        "fragment_call_count": sum(row["operation"] in FRAGMENT_OPERATIONS for row in calls),
        "returned_construction_observation_ids": sorted(constructed),
        "recorded_excerpt_count": excerpts, "recorded_arm_violations": prohibited,
        "scientific_call_budget_exceeded": len(calls) > 4,
        "answer_step_count": len((answer or {}).get("steps", [])),
        "manual_review": {"status": "pending", "useful_inspected_construction_precedents": None,
                          "unsupported_route_steps": None, "transfer_argument_supported": None},
    }


def run_trial(
    directory: Path, case: ComparisonCase, arm: str, baseline: dict[str, Any],
    runtime: AgentRuntime, timeout_seconds: float,
) -> dict[str, Any]:
    """Run one fresh agent thread with production validation and persistent failure evidence."""
    question = comparison_question(case, arm, timeout_seconds)
    InvestigationStore.create(directory, objective=question, baseline=deepcopy(baseline),
                              agent_metadata=runtime.describe())
    workspace = ScientificWorkspace(directory)
    turn = directory / "turns" / "attempt-1"
    turn.mkdir(parents=True)
    user = workspace.store.append("user_message", {"text": question, "origin": "development_comparison"})
    started = perf_counter()
    answer, answer_ref, thread_id, error = None, None, None, None
    usage: dict[str, Any] = {}
    status = "failed"
    print(json.dumps({"case": case.case_id, "arm": arm, "status": "running"}), flush=True)
    try:
        verify_baseline(baseline)
        result = runtime.run(prompt=investigation_prompt(workspace, question), workspace=directory,
                             turn_directory=turn, thread_id=None, cancel=Event(), on_event=lambda _: None)
        usage, thread_id = result.usage, result.thread_id
        (turn / "submitted-answer.json").write_text(json.dumps(result.answer, indent=2), "utf-8")
        verify_baseline(workspace.store.manifest["baseline"])
        validated = ScientificAnswer.model_validate(result.answer)
        validated.evidence_refs = validate_answer_evidence(validated, workspace.store)
        review_ref = matching_evidence_review(workspace.store, validated, after_sequence=user.sequence)
        answer = validated.model_dump()
        saved = workspace.store.append("agent_answer", {
            **answer, "thread_id": thread_id, "usage": usage, "review_status": "unreviewed",
            "self_review_ref": review_ref,
        }, evidence_refs=tuple(validated.evidence_refs))
        answer_ref, status = saved.artifact_ref, "completed"
    except Exception as exc:
        status = exc.status if isinstance(exc, AgentStopped) else "failed"
        error = {"type": type(exc).__name__, "message": str(exc)[:8000]}
        workspace.store.append("agent_turn_error", {"status": status, "error": error})
    report = {"case_id": case.case_id, "arm": arm, "directory": str(directory), "status": status,
              "elapsed_seconds": round(perf_counter() - started, 3), "error": error,
              "thread_id": thread_id, "usage": usage, "answer_ref": answer_ref,
              "runtime_requested": runtime.describe(),
              "baseline_sha256": hashlib.sha256(canonical_bytes(baseline)).hexdigest(),
              "metrics": collect_metrics(workspace, arm, answer)}
    (directory / "trial_report.json").write_text(json.dumps(report, indent=2), "utf-8")
    (directory / "answer.md").write_text((answer or {}).get("answer_markdown") or
                                         f"No validated answer. Status: {status}. See trial_report.json.", "utf-8")
    print(json.dumps({"case": case.case_id, "arm": arm, "status": status,
                      "elapsed_seconds": report["elapsed_seconds"], "error": error}), flush=True)
    return report


def run_comparison(
    root: Path, cases: tuple[ComparisonCase, ...], baseline: dict[str, Any],
    runtime: AgentRuntime, timeout_seconds: float,
) -> dict[str, Any]:
    """Run paired arms sequentially, alternating order and disabling shared lesson memory."""
    if not math.isfinite(timeout_seconds) or not 1 <= timeout_seconds <= 1800:
        raise ValueError("Comparison timeout must be between 1 and 1800 seconds")
    if not cases or len({case.case_id for case in cases}) != len(cases):
        raise ValueError("Supply nonempty, unique comparison cases")
    root = root.resolve()
    root.mkdir(parents=True, exist_ok=False)
    baseline = deepcopy(baseline)
    baseline["evaluation_partition"] = "development_fragment_comparison_not_an_untouched_evaluation"
    baseline["learning_context"] = build_learning_context(baseline)
    rows: list[dict[str, Any]] = []
    report: dict[str, Any] = {"schema_version": "fragment_agent_comparison.v1", "trials": rows,
                              "pairs": [], "limitations": ACTIVE_LIMITATIONS, "chemistry_review_status": "pending"}
    for index, case in enumerate(cases):
        for arm in ARMS if index % 2 == 0 else reversed(ARMS):
            rows.append(run_trial(root / case.case_id / arm, case, arm, baseline, runtime, timeout_seconds))
            # Preserve incremental results even if a later run is interrupted.
            (root / "comparison.json").write_text(json.dumps(report, indent=2), "utf-8")
        pair = [row for row in rows if row["case_id"] == case.case_id]
        report["pairs"].append({
            "case_id": case.case_id, "both_completed": all(row["status"] == "completed" for row in pair),
            "same_scientific_baseline": pair[0]["baseline_sha256"] == pair[1]["baseline_sha256"],
            "same_requested_runtime": pair[0]["runtime_requested"] == pair[1]["runtime_requested"],
            "recorded_constraint_violation": any(row["metrics"]["recorded_arm_violations"] or
                                                 row["metrics"]["scientific_call_budget_exceeded"] for row in pair),
            "chemistry_review_status": "pending",
        })
        (root / "comparison.json").write_text(json.dumps(report, indent=2), "utf-8")
    lines = ["# Fragment tools: paired agent development comparison", "",
             "Execution and retrieval counts are observations, not chemistry-quality scores.", "",
             "| Case | Arm | Status | Seconds | Calls | Returned construction observations |",
             "|---|---|---|---:|---:|---:|"]
    for row in rows:
        metrics = row["metrics"]
        lines.append(f"| {row['case_id']} | {row['arm']} | {row['status']} | {row['elapsed_seconds']} | "
                     f"{metrics['scientific_call_count']} | {len(metrics['returned_construction_observation_ids'])} |")
    lines += ["", "## Review each paired answer", "",
              "Read each arm's answer.md, trial_report.json and cited artifacts. Check whether the source",
              "actually constructs the selected core, whether the agent inspected the relevant record/passage,",
              "whether proposed transfer preserves chemistry/stereo, and which route steps lack support.",
              "Record useful inspected precedents, unsupported steps and the rationale separately from latency.",
              "Timed-out or failed trials remain in the comparison; do not score them as zero-quality chemistry.",
              "", *[f"- {item}" for item in ACTIVE_LIMITATIONS]]
    (root / "comparison.md").write_text("\n".join(lines) + "\n", "utf-8")
    return report


def main() -> int:
    """Require explicit live opt-in; never silently change model or security settings."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-live", action="store_true")
    parser.add_argument("--output", required=True, help="New comparison directory")
    parser.add_argument("--cases", type=Path, default=Path(__file__).with_name("fragment_agent_cases.json"))
    parser.add_argument("--case", action="append", help="Choose case IDs; defaults to all cases")
    parser.add_argument("--artifacts", type=Path, default=Path(__file__).with_name("artifacts.local.example.json"))
    parser.add_argument("--timeout", type=float, default=240, help="Runtime seconds per arm, including final answer")
    parser.add_argument("--codex", help="Optional Codex executable")
    args = parser.parse_args()
    if not args.run_live:
        parser.error("Add --run-live to run the configured agent and research tools")
    if not math.isfinite(args.timeout) or not 1 <= args.timeout <= 1800:
        parser.error("--timeout must be finite and between 1 and 1800")
    cases = load_cases(args.cases)
    if args.case:
        if set(args.case) - {case.case_id for case in cases}:
            parser.error("Unknown case ID")
        cases = tuple(case for case in cases if case.case_id in args.case)
    artifacts = json.loads(args.artifacts.read_text("utf-8"))
    if not isinstance(artifacts, dict) or not all(isinstance(k, str) and isinstance(v, str) for k, v in artifacts.items()):
        parser.error("Artifact configuration must map names to paths")
    repository = Path(__file__).resolve().parents[2]
    selected = {key: (repository / value).resolve() for key, value in artifacts.items()
                if key in {"fragment_index", "retro_library", "condition_index", "procedure_catalog"}}
    if "fragment_index" not in selected or any(not p.is_file() for p in selected.values()):
        parser.error("Provide an existing fragment_index and any configured retro_library")
    runtime = CodexRuntime(executable=args.codex, timeout_seconds=args.timeout)
    baseline = capture_baseline(repository, selected)
    report = run_comparison(Path(args.output), cases, baseline, runtime, args.timeout)
    print(json.dumps({"report": str(Path(args.output) / "comparison.md"),
                      "completed": sum(r["status"] == "completed" for r in report["trials"]),
                      "trials": len(report["trials"])}))
    return 0 if all(r["status"] == "completed" for r in report["trials"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())

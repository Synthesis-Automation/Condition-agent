"""Recorded development smoke test or Tools + guidance supplier investigation.

This is not a blind chemistry benchmark or an independent release evaluation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from threading import Event
from time import monotonic

from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.adapters.operations import ScientificOperations
from chem_coworker.scientific_workspace.adapters.retrochimera import RetroChimeraWorker
from chem_coworker.scientific_workspace.agent_context.context import record_application_context
from chem_coworker.scientific_workspace.agent_context.prompts import investigation_prompt
from chem_coworker.scientific_workspace.core.baseline import capture_baseline
from chem_coworker.scientific_workspace.runtime.agent_runtime import CodexRuntime
from chem_coworker.scientific_workspace.runtime.workspace_modes import WorkspaceMode, runtime_configuration

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    """Run a bounded, recorded supplier check with explicit test provenance."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--artifacts", type=Path)
    parser.add_argument("--agent", action="store_true")
    parser.add_argument("--seconds", type=int, default=900)
    parser.add_argument("--target", default="COc1ccc2c(c1)OCCN2C(=O)C")
    args = parser.parse_args()
    paths = json.loads(args.artifacts.read_text("utf-8")) if args.artifacts else {}
    paths["retrochimera_config"] = str(args.config.resolve())
    baseline = capture_baseline(
        ROOT, {k: (ROOT / v).resolve() for k, v in paths.items()},
        fingerprint_cache=ROOT / "results/ai_native/conversations/.fingerprint-cache",
    )
    store = InvestigationStore.create(
        args.output, objective="Evaluate the optional RetroChimera supplier on development targets",
        baseline=baseline, constraints=("Development only; no independent chemistry validation claim",),
        agent_metadata={"pilot": "retrochimera_supplier", "mode": "tools_guidance" if args.agent else "authored_smoke"},
    )
    store.attach_file(__file__, description="Supplier pilot harness")
    if not args.agent:
        with RetroChimeraWorker() as worker:
            w = ScientificWorkspace(store.root, operations=ScientificOperations(store, retrochimera_worker=worker))
            summaries = []
            for smiles, k in [("Oc1ccc(OCc2ccccc2)c(Br)c1", 5), (args.target, 10)]:
                event = w.run("predict_retrochimera", {"target_smiles": smiles, "num_results": k, "timeout_seconds": 300})
                payload = store.read_artifact(event.artifact_ref)
                summaries.append({"ref": event.artifact_ref, "status": payload["execution_status"],
                                  "result": payload.get("result"), "error": payload.get("error")})
                print(json.dumps(w.call_summary(event)), flush=True)
            (store.root / "smoke_results.json").write_text(json.dumps(summaries, indent=2), "utf-8")
        store.set_status("completed" if all(s["status"] == "completed" for s in summaries) else "failed", "Authored live supplier smoke complete")
        return
    runtime = CodexRuntime(profile="research", reasoning_effort="high", timeout_seconds=args.seconds)
    mode = WorkspaceMode.TOOLS_GUIDANCE
    record_application_context(store, turn_id="supplier-pilot", runtime=runtime_configuration(runtime.describe(), mode))
    w = ScientificWorkspace(store.root)
    question = f"""User-authorized development test of an additional single-step supplier.
Target: {args.target}
Use Tools + guidance mode and finish within {args.seconds} seconds. Do not modify source,
install packages, read prior generated answers, or change pinned artifacts. Research only.
Compare disconnect_target and predict_retrochimera on this target. The new supplier is
experimental: its model probability is NOT feasibility confidence. Use the persistent
scientific console so a second RetroChimera call reuses its loaded model. Request 10
candidates initially, timeout_seconds=300. Inspect the returned structures, choose one
strategically relevant precursor to expand with a second RetroChimera call, and review
at least one consequential proposed step with assess_retro_validity and available
precedent tools. You may use literature search when helpful. Keep evidence gaps explicit.
Read operation help for exact schemas. Do not pass external candidate IDs as operator
realization IDs; submit explicit proposed reaction structures for assessment and cite
the saved supplier call in evidence_refs. Do not treat missing operator support as proof
of chemical impossibility. Also do not call parseable structures experimentally verified.
Explain what each supplier actually contributed. Distinguish a new supplier suggestion
from a route idea you already knew, and variation in reagents from a new disconnection.
Use at most two short proposed routes, clearly identifying incomplete branches. Record
calls and review notes in the workspace. Close the console with quit when finished.
Reserve the last 120 seconds for reporting. Final response: concise Markdown with a
supplier comparison, explicit reaction SMILES, critical review findings, call references,
timings where recorded, and a recommendation. This is one development case, not evidence
of general superiority. Do not create or modify application files.
"""
    store.append("user_message", {"text": question, "origin": "supplier_development_test"})
    turn = store.root / "turns" / "supplier-pilot"
    turn.mkdir(parents=True)
    started = monotonic()
    try:
        result = runtime.run(
            prompt=investigation_prompt(w, question, mode=mode), workspace=store.root,
            turn_directory=turn, thread_id=None, cancel=Event(),
            on_event=lambda event: print(json.dumps({"type": event.get("type")}), flush=True)
            if event.get("type") in {"turn.completed", "turn.failed"} else None,
            mode=mode,
        )
        report = {"status": "completed", "elapsed_seconds": monotonic() - started,
                  "mode": mode.value, "usage": result.usage, "answer": result.answer}
        (store.root / "assistant_review.md").write_text(result.answer.get("answer_markdown", ""), "utf-8")
        store.set_status("completed", "Tools + guidance development investigation completed; chemistry review remains provisional")
    except Exception as exc:
        report = {"status": "error", "elapsed_seconds": monotonic() - started,
                  "mode": mode.value, "error": {"type": type(exc).__name__, "message": str(exc)}}
        store.set_status("failed", "Agent runtime did not complete; inspect preserved calls and logs")
    (store.root / "pilot_report.json").write_text(json.dumps(report, indent=2), "utf-8")
    print(json.dumps({k: v for k, v in report.items() if k != "answer"}), flush=True)


if __name__ == "__main__":
    main()

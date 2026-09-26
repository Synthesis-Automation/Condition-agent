"""Explicitly opted-in live research workflow check using an isolated conversation.

This invokes the configured model and public web services, and can incur usage.
It records execution and evidence plumbing, not chemistry correctness, model
superiority, or independent review. It is never part of deterministic pytest.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from time import monotonic, sleep
from typing import Any


TARGET = "CCOC(=O)N1CC[C@H](NC(=O)c2cc(Cl)c(N)cc2OC)[C@H](OC)C1"
PATENT_URL = "https://patents.google.com/patent/EP0076530A2/en"
ACTIVE_STATUSES = {"queued", "preparing", "running"}
PROMPT = f"""Perform a bounded development integration smoke check, not a synthesis recommendation.
Use at most one native web search for EP0076530A2 Example LXVIII. Report if the
search capability is unavailable; do not pretend an HTTP fetch was native search.
Use the recorded workspace w.run('analyze_molecule', {{'smiles': {TARGET!r}}})
and inspect selected identity/stereochemistry fields. Do not run route planning,
condition retrieval, or a broad literature investigation.

Attempt w.fetch_source({PATENT_URL!r}, title='EP0076530A2'). Retain the returned
artifact even if retrieval/extraction fails or the full patent exceeds limits.
Call w.inspect_source(source.artifact_ref, query='LXVIII', limit=1600). If direct
retrieval/extraction fails or lacks a relevant passage, you may make at most one
native web open of that public patent URL, or the native search result for that
same patent, to inspect the relevant passage. Native search and native open are
separate actions; the total budget is one search and one open. Keep the failed
HTTP source artifact and its error even if native open succeeds.

Only when native open actually returns readable source text, save that exact
returned text using w.capture_source(text, url=..., title=..., locator=...).
Preserve its honest acquisition=agent_supplied_excerpt provenance; this does not
convert it into an HTTP-verified source. If native open returns only a truncated
window near a search hit, capture only that returned window, retain the actual
locator and truncation limit, and do not fill any missing words. Search snippets
alone, model memory, and generated summaries are not substitutes for opened text.

Inspect the usable snapshot with query='LXVIII', then save an exact selection
using w.record_source_excerpt with the returned location start/end and an honest
locator. If LXVIII is absent but the opened text contains a different explicit
experiment, you may inspect and capture that actual experiment within the same
bounded passage budget, clearly identifying its different locator. Never force
an irrelevant passage to satisfy the smoke check. If a hit is only a
cross-reference, or no relevant experiment is readable, say so. Missing evidence
is an acceptable recorded outcome; this smoke check does not require proving a
route. Do not reconstruct missing text or change source/network limits.

Return a concise status answer with the actual artifacts and material limits,
not a procedure or inferred experimental result. Record the required evidence
self-review for the exact final draft, explicitly using not_checked or
not_applicable where appropriate. A completed call or an exact excerpt is not
proof of scientific correctness. Do not claim improved intelligence or chemistry
validation from this smoke check. Preserve failures and give a final answer
using the remaining evidence rather than spending the budget retrying.
"""


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text("utf-8"))


def collect_report(root: Path, submitted: dict[str, str], turn: dict[str, Any]) -> dict[str, Any]:
    """Summarize recorded outcomes without inferring chemistry or hidden runtime state."""
    from chem_coworker.scientific_workspace import ScientificWorkspace

    directory = root / "conversations" / submitted["conversation_id"]
    turn_directory = directory / "turns" / submitted["turn_id"]
    observations = [
        {"path": str(path.relative_to(directory)), "record": _read_json(path)}
        for path in sorted(turn_directory.rglob("runtime-observations.json"))
    ]
    requests = [
        {"path": str(path.relative_to(directory)), "record": _read_json(path)}
        for path in sorted(turn_directory.rglob("runtime-request.json"))
    ]
    calls, sources, excerpts, reviews = [], [], [], []
    if (directory / "investigation.json").is_file():
        workspace = ScientificWorkspace(directory)
        for event in workspace.store.events():
            value = workspace.store.read_artifact(event.artifact_ref)
            if event.kind == "call":
                calls.append({"artifact_ref": event.artifact_ref, **{
                    key: value.get(key) for key in (
                        "operation", "execution_status", "duration_seconds", "timings", "result_bytes", "error",
                    )
                }})
            elif event.kind == "literature_source":
                inspection = workspace.inspect_source(event.artifact_ref, query="LXVIII", limit=1600)
                sources.append({
                    "artifact_ref": event.artifact_ref,
                    "source_url": value.get("source_url"), "acquisition": value.get("acquisition"),
                    "retrieval_status": value.get("retrieval_status"),
                    "extraction_status": value.get("extraction", {}).get("status"),
                    "text_characters": value.get("text_characters"), "error": value.get("error"),
                    "report_inspection_query_found": inspection["query_found"],
                    "report_inspection_scope": "Report-time read; does not prove the agent inspected the source",
                })
            elif event.kind == "literature_excerpt":
                excerpts.append({"artifact_ref": event.artifact_ref, **{
                    key: value.get(key) for key in ("source_ref", "verification", "location", "reported_locator")
                }})
            elif event.kind == "evidence_review":
                reviews.append({"artifact_ref": event.artifact_ref, **{
                    key: value.get(key) for key in ("answer_sha256", "review_status", "unresolved_count")
                }})
    answer = turn.get("answer", {})
    return {
        "schema_version": "scientific_research_smoke.v1",
        "development_smoke_only": True,
        "conversation_id": submitted["conversation_id"], "turn_id": submitted["turn_id"],
        "conversation_directory": str(directory), "turn_status": turn.get("status", "unknown"),
        "error": turn.get("error"), "runtime_requests": requests,
        "runtime_observations": observations, "local_calls": calls,
        "sources": sources, "excerpts": excerpts, "self_reviews": reviews,
        "self_review_status": answer.get("self_review_status", "not_reported"),
        "limitations": [
            "This checks one bounded integration workflow, not comparative intelligence or chemical correctness.",
            "Native search events, source extraction, graph analysis, and self-review are distinct observations.",
            "Failures and missing observations remain visible; successful retrieval does not verify source claims.",
            "No independent chemist review, held-out evaluation, or experimental verification was performed.",
        ],
    }


def main() -> int:
    """Run only with explicit opt-in, and cancel owned work before closing."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-live", action="store_true", help="Opt in to model usage and external web access")
    parser.add_argument("--root", default="results/ai_native/research_smoke",
                        help="Isolated output root; conversations live in its conversations subdirectory")
    parser.add_argument("--timeout", type=float, default=300, help="Overall polling/runtime budget, 1..300 seconds")
    parser.add_argument("--artifacts", help="Optional artifact mapping JSON; defaults to no large corpus files")
    parser.add_argument("--codex", help="Optional native Codex executable")
    parser.add_argument("--model", help="Optional model override; no model is substituted by this script")
    args = parser.parse_args()
    if not args.run_live:
        parser.error("Add --run-live to explicitly authorize this live development smoke check")
    if not math.isfinite(args.timeout) or not 1 <= args.timeout <= 300:
        parser.error("--timeout must be a finite number from 1 to 300 seconds")

    from chem_coworker.scientific_workspace.agent_runtime import CodexRuntime
    from chem_coworker.scientific_workspace.conversation import ConversationService

    root = Path(args.root).resolve()
    artifacts = _read_json(Path(args.artifacts)) if args.artifacts else {}
    if not isinstance(artifacts, dict) or not all(
        isinstance(key, str) and isinstance(value, str) for key, value in artifacts.items()
    ):
        parser.error("--artifacts must contain an object mapping artifact names to file paths")
    runtime = CodexRuntime(executable=args.codex, model=args.model,
                           profile="research", timeout_seconds=args.timeout)
    service = ConversationService(root / "conversations", runtime=runtime, artifacts=artifacts)
    submitted = None
    turn: dict[str, Any] = {}
    interrupted = None
    started = monotonic()
    try:
        submitted = service.submit(PROMPT)
        print(json.dumps({"submitted": submitted, "runtime": runtime.describe()}), flush=True)
        last_status = None
        while monotonic() - started < args.timeout:
            conversation = service.get(submitted["conversation_id"])
            turn = next(item for item in conversation["turns"] if item["id"] == submitted["turn_id"])
            if turn["status"] != last_status:
                last_status = turn["status"]
                print(json.dumps({"status": last_status, "elapsed_seconds": round(monotonic() - started, 1)}), flush=True)
            if turn["status"] not in ACTIVE_STATUSES:
                break
            sleep(0.5)
        else:
            interrupted = "smoke_deadline_exceeded"
    except KeyboardInterrupt:
        interrupted = "keyboard_interrupt"
    finally:
        if submitted is not None:
            service.cancel(submitted["conversation_id"])
        service.close()
    if submitted is None:
        return 1
    conversation = service.get(submitted["conversation_id"])
    turn = next(item for item in conversation["turns"] if item["id"] == submitted["turn_id"])
    report = collect_report(root, submitted, turn)
    report["smoke_interruption"] = interrupted
    report["elapsed_seconds"] = round(monotonic() - started, 3)
    report_path = root / "reports" / f"{submitted['conversation_id']}.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), "utf-8")
    print(json.dumps({"report": str(report_path), "turn_status": turn["status"],
                      "self_review_status": report["self_review_status"], "interruption": interrupted}), flush=True)
    # Completion is an execution status only; individual missing/failing tools
    # are reported above rather than promoted into a chemistry pass/fail score.
    return 0 if turn["status"] == "completed" and interrupted is None else 1


if __name__ == "__main__":
    raise SystemExit(main())

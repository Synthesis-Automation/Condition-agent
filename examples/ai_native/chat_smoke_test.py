"""Live, opt-in question-to-evidence check against the local scientific chat API.

This uses the configured model service and can incur usage. It is not a chemistry
quality benchmark and is never run as part of the deterministic pytest suite.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from time import monotonic, sleep
from typing import Any
from urllib.parse import urlsplit
from urllib.request import Request, urlopen


def main() -> int:
    """Submit a real question, inspect its evidence, and optionally verify a follow-up."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", default="http://127.0.0.1:8000")
    parser.add_argument("--follow-up", action="store_true")
    parser.add_argument("--output", help="Optional local JSON report path")
    args = parser.parse_args()
    if urlsplit(args.url).hostname not in {"127.0.0.1", "localhost", "::1"}:
        parser.error("The development smoke test only connects to loopback")
    base = args.url.rstrip("/") + "/api/v1/scientific"
    token = ""

    def request(path: str, payload: dict[str, Any] | None = None) -> Any:
        body = json.dumps(payload).encode("utf-8") if payload is not None else None
        with urlopen(Request(base + path, data=body, headers={
            "Content-Type": "application/json", "X-Scientific-Token": token,
        }), timeout=30) as response:
            return json.load(response)

    def ask(question: str, conversation_id: str | None = None) -> dict[str, Any]:
        submitted = request("/turns", {"question": question, "conversation_id": conversation_id})
        print(json.dumps({"submitted": submitted}), flush=True)
        deadline = monotonic() + float(configuration["timeout_seconds"]) + 180
        while monotonic() < deadline:
            conversation = request("/conversations/" + submitted["conversation_id"])
            turn = next(item for item in conversation["turns"] if item["id"] == submitted["turn_id"])
            if turn["status"] not in {"queued", "preparing", "running"}:
                if turn["status"] != "completed":
                    raise RuntimeError(json.dumps(turn.get("error", turn["status"])))
                for reference in turn["answer"]["evidence_refs"]:
                    request(f"/conversations/{submitted['conversation_id']}/artifacts/{reference}")
                return turn
            sleep(1)
        request("/conversations/" + submitted["conversation_id"] + "/cancel", {})
        raise TimeoutError("Live smoke test deadline exceeded")

    configuration = request("/config")
    token = configuration.pop("token")
    first = ask(
        "Analyze CCBr.N>>CCN using the recorded analyze_reaction operation. "
        "Briefly explain the structural evidence and one limitation. "
        "Inspect selected output fields, not the whole analysis. "
        "Do not run conditions or routes for this small development smoke test."
    )
    if not first["answer"]["evidence_refs"]:
        raise RuntimeError("The structure-specific smoke answer lacks recorded evidence")
    report = {"configuration": configuration, "turns": [first], "development_smoke_only": True}
    if args.follow_up:
        second = ask(
            "Does the analysis above prove experimental feasibility? "
            "Answer in two sentences using the existing evidence without new calculations.",
            first["conversation_id"],
        )
        if second["thread_id"] != first["thread_id"]:
            raise RuntimeError("Follow-up did not resume the original runtime thread")
        report["turns"].append(second)
    if args.output:
        path = Path(args.output)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=2, ensure_ascii=False), "utf-8")
    print(json.dumps({
        "status": "passed", "conversation_id": first["conversation_id"],
        "turns": len(report["turns"]), "answer": first["answer"]["answer_markdown"],
    }, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

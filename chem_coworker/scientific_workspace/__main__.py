"""Run recorded scientific operations with JSON request files."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from .store import canonical_bytes
from .workspace import ScientificWorkspace


def main(argv: Sequence[str] | None = None) -> int:
    """Create, inspect, invoke, annotate, and replay scientific investigations."""
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    create = sub.add_parser("init")
    create.add_argument("workspace")
    create.add_argument("--objective", required=True)
    create.add_argument("--repository", default=str(Path(__file__).resolve().parents[2]))
    create.add_argument("--artifacts", help="JSON file mapping artifact names to local paths")
    for name in ("catalog", "summary", "capabilities", "run", "run-python", "show", "replay", "note", "attach", "status", "fetch-source", "inspect-source"):
        command = sub.add_parser(name)
        command.add_argument("workspace")
        if name == "run":
            command.add_argument("operation")
            command.add_argument("--input", required=True, help="JSON object with operation arguments")
        elif name == "fetch-source":
            command.add_argument("url")
            command.add_argument("--title")
        elif name == "inspect-source":
            command.add_argument("reference")
            command.add_argument("--query")
            command.add_argument("--offset", type=int, default=0)
            command.add_argument("--limit", type=int, default=4000)
        elif name == "run-python":
            command.add_argument("script", help="Python file relative to the investigation")
            command.add_argument("--input", required=True, help="JSON parameters file")
            command.add_argument("--evidence", action="append", default=[])
            command.add_argument("--timeout", type=int, default=60)
        elif name in {"show", "replay"}:
            command.add_argument("reference")
        elif name == "note":
            command.add_argument("kind", choices=("hypothesis", "decision", "question", "limitation", "review"))
            command.add_argument("--text", required=True)
            command.add_argument("--evidence", action="append", default=[])
        elif name == "attach":
            command.add_argument("file")
            command.add_argument("--description", required=True)
            command.add_argument("--evidence", action="append", default=[])
        elif name == "status":
            command.add_argument("status")
            command.add_argument("--reason", required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "init":
            paths = json.loads(Path(args.artifacts).read_text("utf-8")) if args.artifacts else {}
            workspace = ScientificWorkspace.create(
                args.workspace, objective=args.objective, repository=args.repository, artifacts=paths,
            )
            result = workspace.store.summary()
        else:
            workspace = ScientificWorkspace(args.workspace)
            if args.command == "catalog":
                result = workspace.operations.catalog()
            elif args.command == "summary":
                result = workspace.store.summary()
            elif args.command == "capabilities":
                result = workspace.capabilities()
            elif args.command == "fetch-source":
                event = workspace.fetch_source(args.url, title=args.title)
                value = workspace.store.read_artifact(event.artifact_ref)
                result = {"event": event, "retrieval_status": value["retrieval_status"],
                          "extraction_status": value["extraction"]["status"], "error": value.get("error")}
                print(canonical_bytes(result).decode("utf-8"))
                return 0 if value["retrieval_status"] == "completed" else 1
            elif args.command == "inspect-source":
                result = workspace.inspect_source(args.reference, query=args.query, offset=args.offset, limit=args.limit)
            elif args.command == "show":
                result = workspace.store.read_artifact(args.reference)
            elif args.command == "run":
                event = workspace.run(args.operation, json.loads(Path(args.input).read_text("utf-8")))
                result = workspace.call_summary(event)
                print(canonical_bytes(result).decode("utf-8"))
                return 0 if result["execution_status"] == "completed" else 1
            elif args.command == "run-python":
                event = workspace.run_python(args.script, json.loads(Path(args.input).read_text("utf-8")),
                                             evidence_refs=tuple(args.evidence), timeout_seconds=args.timeout)
                result = workspace.call_summary(event)
                print(canonical_bytes(result).decode("utf-8"))
                return 0 if result["execution_status"] == "completed" else 1
            elif args.command == "replay":
                event = workspace.replay(args.reference)
                result = workspace.store.read_artifact(event.artifact_ref)
                print(canonical_bytes({"event": event, "matches": result["matches"]}).decode("utf-8"))
                return 0 if result["matches"] else 1
            elif args.command == "note":
                result = workspace.store.note(args.kind, args.text, evidence_refs=tuple(args.evidence))
            elif args.command == "attach":
                result = workspace.store.attach_file(args.file, description=args.description, evidence_refs=tuple(args.evidence))
            else:
                result = workspace.store.set_status(args.status, args.reason)
        print(canonical_bytes(result).decode("utf-8"))
        return 0
    except (OSError, ValueError, TypeError, RuntimeError) as exc:
        print(canonical_bytes({"error": str(exc), "error_type": type(exc).__name__}).decode("utf-8"))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

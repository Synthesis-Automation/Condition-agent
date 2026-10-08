"""Agent-steered search pilot using existing recorded scientific operations.

CLI: python -m examples.ai_native.iterative_search_poc WORKSPACE REQUEST.json
Requests use actions init, prepare, search, baseline, inspect, or summary.
This development harness adds no chemistry rules and never chooses the next query.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

from chem_coworker.scientific_workspace import ScientificWorkspace


SCHEMA = "iterative_search_pilot.v1"


def search_card(hit: dict[str, Any]) -> dict[str, Any]:
    """Show structures and evidence distinctions without loading procedure prose."""
    matches = hit.get("matches", [])
    witnesses = [w for match in matches for w in match.get("witnesses", [])]
    return {
        "observation_id": hit["observation_id"],
        "reference_id": hit.get("reference_id"),
        "reaction_smiles": hit.get("record", {}).get("reaction_smiles"),
        "product_smiles": hit.get("product_smiles"),
        "relationships": hit.get("relationships", []),
        "discovery": hit.get("discovery"),
        "witnesses": witnesses[:3], "witnesses_omitted": max(0, len(witnesses) - 3),
        "procedure_availability": hit.get("procedure_availability"),
        "warnings": hit.get("warnings", []),
    }


class SearchPilot:
    """Thin experiment harness; the caller supplies each query and revision reason."""

    def __init__(self, root: str | Path) -> None:
        self.workspace = ScientificWorkspace(root)

    @classmethod
    def create(cls, root: str | Path, request: dict[str, Any]) -> "SearchPilot":
        """Freeze inputs and equal per-arm search budgets before the first search."""
        budget = request.get("budget_seconds", 90)
        if type(budget) is not int or not 1 <= budget <= 120:
            raise ValueError("budget_seconds must be an integer in 1..120")
        workspace = ScientificWorkspace.create(
            root, objective=request["objective"], repository=Path(__file__).resolve().parents[2],
            artifacts=request["artifacts"],
            constraints=("Known development target; agent-authored review, not blind evaluation.",
                         "No web retrieval or prior pilot results during search.",
                         "Equal cumulative operation-time budgets; reasoning and inspection measured separately."),
        )
        workspace.store.append("pilot_config", {
            "schema_version": SCHEMA, "target_smiles": request["target_smiles"],
            "budget_seconds": budget,
            "budget_scope": "search operation time, including worker startup; preparation and review excluded",
        })
        workspace.store.attach_file(__file__, description="Exact pilot harness at initialization")
        return cls(root)

    def records(self, kind: str) -> list[dict[str, Any]]:
        """Read immutable experiment records in event order."""
        return [self.workspace.store.read_artifact(e.artifact_ref)
                for e in self.workspace.store.events() if e.kind == kind]

    def _result(self, reference: str, operation: str) -> dict[str, Any]:
        payload = self.workspace.store.read_artifact(reference)
        if payload.get("operation") != operation or payload.get("execution_status") != "completed":
            raise ValueError(f"Expected a completed {operation} artifact")
        return payload["result"]

    def dispatch(self, request: dict[str, Any]) -> dict[str, Any]:
        """Execute exactly one caller-selected action and preserve its rationale."""
        action = request["action"]
        config = self.records("pilot_config")[0]
        target = config["target_smiles"]
        w = self.workspace
        if action == "summary":
            return {"config": config, "runs": self.records("pilot_run"),
                    "inspections": self.records("pilot_inspection")}
        reason = request.get("reason", "").strip()
        if not reason:
            raise ValueError("Provide the search question or reason for this action")
        if action == "prepare":
            parent = request.get("parent_search_ref")
            if parent:
                previous = w.store.read_artifact(parent)
                if previous.get("operation") != "search_fragment_precedents":
                    raise ValueError("A revision must reference an actual fragment search")
            refs = (parent,) if parent else ()
            w.store.note("decision", reason, evidence_refs=refs)
            event = w.run("propose_fragment_queries", {
                "target_smiles": target, "query": request["query"],
                "query_format": request.get("query_format", "smiles"),
                "topology": request.get("topology", "preserve_rings"),
                **({"aromatic_atom_ids": request["aromatic_atom_ids"]}
                   if "aromatic_atom_ids" in request else {}),
            }, evidence_refs=refs)
            return {"query_ref": event.artifact_ref,
                    "payload": w.store.read_artifact(event.artifact_ref)}
        if action == "inspect":
            ref = request["search_ref"]
            payload = w.store.read_artifact(ref)
            if payload.get("operation") not in {"search_fragment_precedents", "find_synthesis_precedents"}:
                raise ValueError("Inspection requires a search result")
            hits = payload.get("result", {}).get("hits", [])
            hit = next(h for h in hits if h["observation_id"] == request["observation_id"])
            event = w.store.append("pilot_inspection", {
                "schema_version": SCHEMA, "search_ref": ref,
                "observation_id": hit["observation_id"], "reason": reason,
                "status": "opened_for_agent_review", "usefulness": "not_yet_adjudicated",
            }, evidence_refs=(ref,))
            return {"inspection_ref": event.artifact_ref, "hit": hit}
        if action not in {"search", "baseline"}:
            raise ValueError("Unknown action")
        arm = "iterative" if action == "search" else "automatic"
        used = sum(r["operation_seconds"] for r in self.records("pilot_run") if r["arm"] == arm)
        remaining = math.floor(config["budget_seconds"] - used)
        if remaining < 1:
            raise ValueError("Search budget exhausted")
        maximum = 30 if action == "search" else 120
        timeout = min(request.get("timeout_seconds", maximum), maximum, remaining)
        refs: tuple[str, ...] = ()
        if action == "search":
            ref = request["query_ref"]
            preview = self._result(ref, "propose_fragment_queries")
            variant_id = request.get("variant_id")
            if variant_id:
                selected = next(v for v in preview["variants"] if v["variant_id"] == variant_id)
            else:
                selected = {**preview["parent_query"], "query": preview["parent_query"]["expression"]}
            arguments = {k: selected[k] for k in ("query", "query_format", "topology")}
            arguments.update(target_smiles=target, limit=10, timeout_seconds=timeout)
            if variant_id:
                arguments.update(query_variant_ref=ref, query_variant_id=variant_id)
            operation = "search_fragment_precedents"
            refs = (ref,)
        else:
            arguments = {"target_smiles": target, "limit": 10, "timeout_seconds": timeout}
            operation = "find_synthesis_precedents"
        w.store.note("decision", reason, evidence_refs=refs)
        event = w.run(operation, arguments, evidence_refs=refs)
        payload = w.store.read_artifact(event.artifact_ref)
        result = payload.get("result", {})
        cards = [search_card(hit) for hit in result.get("hits", [])]
        run = {
            "schema_version": SCHEMA, "arm": arm, "reason": reason,
            "search_ref": event.artifact_ref, "query_ref": request.get("query_ref"),
            "variant_id": request.get("variant_id"), "timeout_seconds": timeout,
            "operation_seconds": payload.get("timings", {}).get("operation_seconds", payload["duration_seconds"]),
            "workspace_seconds": payload["duration_seconds"],
            "execution_status": payload["execution_status"],
            "worker_status": result.get("execution_status"),
            "search_status": result.get("search_status"), "stop_reason": result.get("stop_reason"),
            "error": payload.get("error") or result.get("error"),
            "counts": result.get("counts"), "relationship_groups": result.get("relationship_groups"),
            "execution": result.get("execution"), "source_scope": result.get("source_scope"),
            "cards": cards,
            "full_result_bytes": len(json.dumps(result).encode()),
            "card_bytes": len(json.dumps(cards).encode()),
            "limitations": result.get("limitations", []),
        }
        w.store.append("pilot_run", run, evidence_refs=(event.artifact_ref,))
        return run


def main() -> None:
    """Read one JSON request; emit one result, leaving strategy to the caller."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workspace")
    parser.add_argument("request")
    parser.add_argument("--output", help="Save full response; stdout then prints its path")
    args = parser.parse_args()
    request = json.loads(Path(args.request).read_text("utf-8"))
    if request["action"] == "init":
        pilot = SearchPilot.create(args.workspace, request)
        result = pilot.dispatch({"action": "summary"})
    else:
        result = SearchPilot(args.workspace).dispatch(request)
    text = json.dumps(result, indent=2)
    if args.output:
        Path(args.output).write_text(text, "utf-8")
        print(json.dumps({"output": args.output}))
    else:
        print(text)


if __name__ == "__main__":
    main()

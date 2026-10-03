"""Deterministic fragment-guided single-step development investigation.

Use existing graph query, source compiler and reverse/forward validators. This
experiment changes search guidance, never chemistry admission. Run as a module
with --output results/ai_native/<new-directory>. No LLM or network is invoked.
"""

from __future__ import annotations

import argparse
from html import escape
import json
from pathlib import Path
import sys
from typing import Any

from core_retrosynthesis.fragment_guidance import (
    DEFAULT_POLICY, evaluate_transfers, load_policy, project_construction_guidance,
    restore_source_smiles, scientific_projection,
)


DEFAULT_TARGET = "Fc(cn1)cc2c1c(c3ccccc3OC)n[nH]2"
DEFAULT_INDEX = "results/ai_native/indexes/fragment_precedents.sqlite"
DEFAULT_LIBRARY = (
    "results/operator_retrosynthesis_poc/full_scale_v3/compact/"
    "operator_library_v3.json.gz"
)


def worker(inputs: dict[str, Any]) -> dict[str, Any]:
    """Consume saved evidence or run canonical chemistry inside a recorded subprocess."""
    parameters = inputs["parameters"]
    if parameters["action"] == "project":
        evidence = inputs["evidence"]
        suggestions = evidence[parameters["suggestions_ref"]]["result"]
        searches = [{**item, **evidence[item["artifact_ref"]]}
                    for item in parameters["searches"]]
        return project_construction_guidance(suggestions, searches, parameters["policy"])
    if parameters["action"] == "evaluate":
        parameters["guidance"] = inputs["evidence"][parameters["guidance_ref"]]["result"]
        return evaluate_transfers(parameters)
    raise ValueError("Unsupported POC worker action")


def render_report(report: dict[str, Any], directory: Path) -> None:
    """Write an inspectable local HTML report using the shared molecular renderer."""
    from visualization import render_route_scheme_svg

    def scheme(reaction: str) -> str:
        parts = reaction.split(">")
        return render_route_scheme_svg([f"{parts[0]}>>{parts[-1]}"]).decode("utf-8")

    def card(candidate: dict[str, Any]) -> str:
        return ("<article>" + scheme(candidate["proposed_reaction_smiles"])
                + "<p>Graph validation: " + escape(candidate["forward_validation_status"])
                + "; context level: " + escape(candidate["abstraction_level"]) + "</p>"
                + "<p>Supporting operator reactions: "
                + escape(", ".join(candidate["precedent_reaction_ids"])) + "</p>"
                + "<details><summary>SMILES and validation evidence</summary><pre>"
                + escape(json.dumps(candidate, indent=2)) + "</pre></details></article>")

    transfer = report["transfers"]
    sections = ["<h2>Unrestricted baseline</h2>",
                *(card(c) for c in transfer["baseline"]["candidates"])]
    for branch in transfer["guided"]:
        sections.append("<h2>Construction-site hypothesis: target atoms "
                        + escape(str(branch["target_atom_ids"])) + "</h2>")
        for arm, title in (("direct_source_transfer", "Transfer of compiled search-hit operators"),
                           ("witness_directed_library", "Witness-directed prepared-library search")):
            sections.append("<h3>" + title + "</h3>")
            candidates = branch[arm]["candidates"]
            sections.extend(card(c) for c in candidates)
            if not candidates:
                sections.append("<p>No supported candidate in this bounded arm.</p>")
    sections.append("<h2>Selected source reactions</h2>")
    for record in report["guidance"]["source_records"]:
        sections.extend(["<h3>" + escape(record["reaction_id"]) + "</h3>",
                         scheme(record["reaction_smiles"])])
    sections.append("<h2>Compiler admission and rejection evidence</h2><pre>"
                    + escape(json.dumps(transfer["source_admissions"], indent=2)) + "</pre>")
    sections.append("<h2>Limits</h2><ul>" + "".join(
        "<li>" + escape(value) + "</li>" for value in report["limitations"]) + "</ul>")
    page = ("<!doctype html><html lang='en'><meta charset='utf-8'>"
            "<title>Deterministic fragment-guided retrosynthesis POC</title>"
            "<style>body{font:16px system-ui;max-width:1100px;margin:32px auto;padding:0 20px}"
            "article{border:1px solid #ccd7df;padding:18px;margin:16px 0;border-radius:8px}"
            "svg{max-width:100%;height:auto}pre{white-space:pre-wrap;overflow-wrap:anywhere}"
            "code{overflow-wrap:anywhere}</style><h1>Fragment-guided retrosynthesis POC</h1>"
            "<p>Target: <code>" + escape(report["target_smiles"]) + "</code></p>"
            "<p>Development proposals. Graph reconstruction does not establish experimental feasibility.</p>"
            + "".join(sections) + "</html>")
    (directory / "report.html").write_text(page, "utf-8")


def main() -> int:
    """Freeze inputs, record fragment calls and custom computations, and save results."""
    from chem_coworker.scientific_workspace import ScientificWorkspace

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="New investigation directory")
    parser.add_argument("--target", default=DEFAULT_TARGET)
    parser.add_argument("--index", default=DEFAULT_INDEX)
    parser.add_argument("--library", default=DEFAULT_LIBRARY)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    args = parser.parse_args()
    repository = Path(__file__).resolve().parents[2]
    policy = load_policy(args.policy)
    workspace = ScientificWorkspace.create(
        args.output, repository=repository,
        objective="Deterministic fragment-guided single-step retrosynthesis POC",
        artifacts={"fragment_index": args.index, "retro_library": args.library,
                   "poc_policy": args.policy.resolve()},
        constraints=("Development only; no independent chemistry or untouched evaluation claim",),
    )
    directory = workspace.store.root
    (directory / "poc.py").write_text(Path(__file__).read_text("utf-8"), "utf-8")

    def completed(event: Any) -> dict[str, Any]:
        payload = workspace.store.read_artifact(event.artifact_ref)
        if payload["execution_status"] != "completed":
            raise RuntimeError(f"Recorded failure {event.artifact_ref}: {payload.get('error')}")
        return payload["result"]

    event = workspace.run("suggest_search_fragments", {"target_smiles": args.target})
    suggestions_ref = event.artifact_ref
    suggestions = completed(event)
    search_rows, refs = [], [suggestions_ref]
    search_repeats = []
    for candidate in suggestions["candidates"][:policy["query_limit"]]:
        arguments = {"query": candidate["query"], "query_format": candidate["query_format"],
                     "topology": candidate["topology"], "target_smiles": suggestions["target_smiles"],
                     "limit": policy["hit_limit"], "timeout_seconds": policy["search_timeout_seconds"]}
        event = workspace.run("search_fragment_precedents", arguments,
                              evidence_refs=(suggestions_ref,))
        result = workspace.store.read_artifact(event.artifact_ref)
        refs.append(event.artifact_ref)
        search_rows.append({"candidate_id": candidate["candidate_id"], "artifact_ref": event.artifact_ref})
        q = result.get("result") or {}
        print(json.dumps({"query": candidate["query"], "status": result["execution_status"],
                          "search_status": q.get("search_status"),
                          "constructed": q.get("relationship_groups", {}).get("constructed")}), flush=True)
        if result["execution_status"] == "completed" and q.get("search_status") == "complete":
            repeat_event = workspace.run("search_fragment_precedents", arguments,
                                         evidence_refs=(suggestions_ref,))
            repeat_payload = workspace.store.read_artifact(repeat_event.artifact_ref)
            repeat_result = repeat_payload.get("result") or {}
            search_repeats.append({"candidate_id": candidate["candidate_id"],
                                   "artifact_ref": repeat_event.artifact_ref,
                                   "matches": repeat_payload["execution_status"] == "completed"
                                   and scientific_projection(q) == scientific_projection(repeat_result)})
            refs.append(repeat_event.artifact_ref)
    event = workspace.run_python("poc.py", {
        "action": "project", "suggestions_ref": suggestions_ref,
        "searches": search_rows, "policy": policy,
    }, evidence_refs=tuple(refs), timeout_seconds=30)
    guidance_ref, guidance = event.artifact_ref, completed(event)
    print(json.dumps({"construction_sites": [b["target_atom_ids"] for b in guidance["focus_bonds"]],
                      "eligible_sources": guidance["eligible_source_count"]}), flush=True)
    event = workspace.run_python("poc.py", {
        "action": "evaluate", "guidance_ref": guidance_ref, "policy": policy,
        "library": workspace.store.manifest["baseline"]["artifacts"]["retro_library"]["path"],
    }, evidence_refs=(guidance_ref,), timeout_seconds=120)
    transfer_ref, transfers = event.artifact_ref, completed(event)
    report = {
        "schema_version": "fragment_guided_retrosynthesis_poc.v1", "development_only": True,
        "input_target_smiles": args.target, "target_smiles": suggestions["target_smiles"],
        "policy": policy, "suggestions_ref": suggestions_ref, "searches": search_rows,
        "guidance_ref": guidance_ref, "transfer_ref": transfer_ref,
        "guidance": guidance, "transfers": transfers, "search_repeat_checks": search_repeats,
        "limitations": [
            "Authored development target; no independent chemistry review or untouched evaluation.",
            "Single-step proposals only; precursor preparation, stock and complete route feasibility are unassessed.",
            "A local construction witness does not verify its complete source reaction; compiler admission is separate.",
            "A witness-directed library candidate may use a different source operator and different other edits.",
            "Only direct-source candidates transfer a compiler-admitted search-hit operator.",
            "Query alignments are analogue hypotheses, never observed target atom correspondence.",
            "Returned hits and representative query embeddings are bounded; absence is not synthetic impossibility.",
            "Search has wall-clock deadlines; incomplete searches cannot seed this POC.",
            "Repeated complete-result equality is observed here, not a guarantee across machines or versions.",
            "Arms have different total work budgets; candidate counts do not establish an improvement.",
            "Conditions, selectivity and experimental feasibility remain unvalidated.",
        ],
    }
    path = directory / "report.json"
    path.write_text(json.dumps(report, indent=2), "utf-8")
    render_report(report, directory)
    for artifact in (path, directory / "report.html"):
        workspace.store.attach_file(artifact, description="Deterministic fragment-guided development POC",
                                    evidence_refs=(guidance_ref, transfer_ref))
    workspace.store.note("limitation", "Direct source admission and witness-directed search are distinct; "
                         "these development results do not validate experimental synthesis.",
                         evidence_refs=(transfer_ref,))
    print(json.dumps({"report": str(path), "html": str(directory / "report.html"),
                      "source_templates": transfers["compiled_source_template_count"],
                      "search_repeat_matches": all(c["matches"] for c in search_repeats),
                      "transfer_repeat_matches": transfers["repeat_scientific_results_identical"]}), flush=True)
    return 0 if transfers["repeat_scientific_results_identical"] and all(
        c["matches"] for c in search_repeats) else 1


if __name__ == "__main__":
    if len(sys.argv) == 3 and Path(sys.argv[1]).name == "input.json":
        inputs = json.loads(Path(sys.argv[1]).read_text("utf-8"))
        Path(sys.argv[2]).write_text(json.dumps(worker(inputs)), "utf-8")
    else:
        raise SystemExit(main())

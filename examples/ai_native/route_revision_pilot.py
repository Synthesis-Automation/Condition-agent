"""Authored development route revision using recorded canonical assessments.

This deterministic example exercises the contracts, not agent intelligence or
experimental chemistry. It never edits the operator library or source records.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from chem_coworker.scientific_workspace import ScientificWorkspace


def main() -> int:
    """Assess a sketch with an unavailable leaf and investigate making that leaf."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="New investigation directory")
    parser.add_argument("--library", required=True, help="Canonical generic operator library")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    workspace = ScientificWorkspace.create(
        args.output, objective="Development example: revise an unavailable ethylamine branch",
        repository=root, artifacts={"retro_library": args.library},
        constraints=("Authored development case; no independent chemistry review or experimental validation",),
    )

    def call(operation: str, arguments: dict[str, Any]) -> str:
        event = workspace.run(operation, arguments)
        summary = workspace.call_summary(event)
        print(json.dumps(summary, ensure_ascii=False), flush=True)
        if summary["execution_status"] != "completed":
            raise RuntimeError(f"Recorded {operation} failed; inspect {event.artifact_ref}")
        return event.artifact_ref

    original = call("assess_route_proposal", {
        "proposal": {"target_smiles": "O=C(NCC)c1ccccc1", "steps": [{
            "external_step_id": "amide", "target_smiles": "O=C(NCC)c1ccccc1",
            "precursor_smiles": "O=C(O)c1ccccc1.NCC",
        }]}, "unavailable_starting_materials": ["NCC"],
    })
    revised = call("revise_route_branch", {
        "source_ref": original, "remove_step_ids": [], "replacement_steps": [{
            "external_step_id": "amine", "target_smiles": "NCC", "precursor_smiles": "CC=O.N",
        }], "reason": "Investigate preparing the user-declared unavailable intermediate",
        "assumptions": ["The proposed amine preparation is a development hypothesis"],
        "risks": ["Selectivity, operating conditions and actual starting-material stock remain unverified"],
    })
    inspection = call("inspect_route_step", {"source_ref": revised, "step_id": "amine"})
    comparison = call("compare_route_proposals", {"source_refs": [original, revised]})
    replay = ScientificWorkspace(workspace.store.root).replay(revised)
    replay_matches = workspace.store.read_artifact(replay.artifact_ref)["matches"]
    report = {
        "development_example_only": True, "original_ref": original, "revision_ref": revised,
        "inspection_ref": inspection, "comparison_ref": comparison,
        "replay_ref": replay.artifact_ref, "replay_matches": replay_matches,
        "comparison": workspace.store.read_artifact(comparison)["result"],
    }
    report_path = workspace.store.root / "route_revision_report.json"
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), "utf-8")
    workspace.store.attach_file(report_path, description="Authored route revision development report",
                                evidence_refs=(original, revised, inspection, comparison, replay.artifact_ref))
    print(json.dumps({"report": str(report_path), "replay_matches": replay_matches}), flush=True)
    return 0 if replay_matches else 1


if __name__ == "__main__":
    raise SystemExit(main())

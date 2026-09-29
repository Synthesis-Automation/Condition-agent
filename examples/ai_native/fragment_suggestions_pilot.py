"""Recorded development panel for optional fragment queries, not agent evaluation.

Run from the repository root with ``python -m examples.ai_native.fragment_suggestions_pilot
--output results/ai_native/fragment_suggestions_pilot``. An optional ``--index``
enables one explicit end-to-end search for the cyclic-ether core, not every query.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from statistics import median

from rdkit import Chem

from chem_coworker.scientific_workspace import ScientificWorkspace
from reactive_taxonomy.fragment_search import compile_fragment_query, fragment_embeddings


CASES = {
    "acetyl_cyclic_ether": "CC(=O)c1ccc2c(c1)COc1ccccc1-2",
    "simple_biphenyl": "c1ccccc1-c1ccccc1",
    "substituted_biphenyl": "CCOc1ccc(-c2ccc(N)cc2)cc1",
    "fused_aromatic": "CCc1ccc2ccccc2c1",
    "fused_heterocycle": "CCc1ccc2ncccc2c1",
    "spiro": "CC1CCC2(CC1)CCCC2",
    "bridged": "CC1CC2CCC1C2",
    "lactam": "CCC1CCNC1=O",
    "lactone": "CCC1CCOC1=O",
    "macrocycle": "CC1CCCCCCCCCCC1",
    "linked_rings": "c1ccccc1CCOc1ccncc1",
    "morpholine_substituent": "CCOc1ccc(CN2CCOCC2)cc1",
    "amide_chain": "CC(=O)NCCO",
    "acyclic_chiral": "CC[C@H](O)CCCCCCN",
    "alkene_stereo": "CC/C=C/CCO",
    "ammonium": "C[N+](C)(C)CCO",
    "zwitterion": "[NH3+]CC(=O)[O-]",
    "carbon_isotope": "[13CH3]COCCO",
    "hydrogen_isotope": "[2H]C(O)CCN",
    "aromatic_cation": "CC[n+]1ccccc1",
    "sulfonamide": "CCS(=O)(=O)NCCO",
    "nitrile": "N#CCCc1ccccc1",
    "urea": "CCNC(=O)NCCOc1ccccc1",
    "small_ether": "COC",
}


def main() -> int:
    """Record candidate validity, timings, replay and an optional chosen-core search."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="New investigation directory")
    parser.add_argument("--index", help="Existing fragment precedent index; never built here")
    args = parser.parse_args()
    workspace = ScientificWorkspace.create(
        args.output, objective="Development validation of target-derived search fragments",
        repository=Path(__file__).resolve().parents[2],
        artifacts={"fragment_index": args.index} if args.index else {},
        constraints=("Development panel, not independent chemistry review or an agent comparison",),
    )
    rows, first_ref, first_query = [], None, None
    for name, smiles in CASES.items():
        event = workspace.run("suggest_search_fragments", {"target_smiles": smiles})
        payload = workspace.store.read_artifact(event.artifact_ref)
        if payload["execution_status"] != "completed":
            raise RuntimeError(f"{name}: inspect {event.artifact_ref}")
        result = payload["result"]
        target = Chem.MolFromSmiles(result["target_smiles"])
        valid = []
        for candidate in result["candidates"]:
            query = compile_fragment_query(candidate["query"])
            matches, _ = fragment_embeddings(query, target, maximum=10000)
            valid.append(tuple(candidate["query_atom_target_ids"]) in matches)
        if not valid or not all(valid):
            raise RuntimeError(f"Invalid or empty panel result: {name}")
        rows.append({"name": name, "input": smiles, "artifact_ref": event.artifact_ref,
                     "candidate_count": len(valid), "all_queries_match_selected_atoms": all(valid),
                     "rejected_count": result["rejected_count"], "timings": payload["timings"],
                     "workspace_seconds": payload["duration_seconds"]})
        if first_ref is None:
            first_ref, first_query = event.artifact_ref, result["candidates"][0]["query"]
    replay = workspace.replay(first_ref)
    replay_matches = workspace.store.read_artifact(replay.artifact_ref)["matches"]
    report = {"development_panel_only": True, "cases": rows, "replay_matches": replay_matches,
              "median_operation_seconds": median(row["timings"]["operation_seconds"] for row in rows),
              "max_operation_seconds": max(row["timings"]["operation_seconds"] for row in rows)}
    if args.index:
        event = workspace.run("search_fragment_precedents", {
            "query": first_query, "limit": 5, "timeout_seconds": 30,
        }, evidence_refs=(first_ref,))
        payload = workspace.store.read_artifact(event.artifact_ref)
        report["chosen_core_search"] = {"artifact_ref": event.artifact_ref,
                                       "summary": workspace.call_summary(event)}
        if payload["execution_status"] != "completed":
            raise RuntimeError(f"Chosen-core search failed: inspect {event.artifact_ref}")
    path = workspace.store.root / "fragment_suggestions_report.json"
    path.write_text(json.dumps(report, indent=2), "utf-8")
    workspace.store.attach_file(path, description="Optional fragment suggestions development panel",
                                evidence_refs=tuple(row["artifact_ref"] for row in rows))
    print(json.dumps({"report": str(path), "cases": len(rows), "replay_matches": replay_matches,
                      "median_operation_seconds": report["median_operation_seconds"],
                      "max_operation_seconds": report["max_operation_seconds"]}))
    return 0 if replay_matches else 1


if __name__ == "__main__":
    raise SystemExit(main())

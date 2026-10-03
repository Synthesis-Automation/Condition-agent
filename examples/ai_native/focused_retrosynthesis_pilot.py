"""Recorded paired development evaluation of focused single-step search.

Run with --output results/ai_native/<new-name> --library <existing-library>.
Cases are authored development probes, not a held-out or independent chemistry
evaluation. A fixed single search pool gives both arms the same total budgets;
the agent-facing specificity ladder has per-level budgets and is checked separately.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
from statistics import median
import sys
from time import perf_counter
from typing import Any

from rdkit import Chem
from rdchiral.main import rdchiralReactants

from reactive_taxonomy.chemistry.smarts_cache import compile_smarts
from reactive_taxonomy.fragment_search import indexed_product


# Atom pairs are selected by graph queries before either search runs.
CASES = (
    {"name": "ethyl_methyl_amine_ethyl", "target": "CCNC", "bond": [1, 2]},
    {"name": "ethyl_methyl_amine_methyl", "target": "CCNC", "bond": [2, 3]},
    {"name": "symmetric_diethylamine_left", "target": "CCNCC", "bond": [1, 2]},
    {"name": "symmetric_diethylamine_right", "target": "CCNCC", "bond": [2, 3]},
    {"name": "biaryl_amide_amide", "target": "CC(=O)Nc1ccc(-c2ccccc2)cc1",
     "query": "[C](=[O])-[N]", "endpoints": [0, 2]},
    {"name": "biaryl_amide_biaryl", "target": "CC(=O)Nc1ccc(-c2ccccc2)cc1",
     "query": "[c]-[c]", "endpoints": [0, 1]},
    {"name": "ether_amide_ether", "target": "COc1ccc(NC(C)=O)cc1",
     "query": "[CH3]-[O]-[c]", "endpoints": [0, 1]},
    {"name": "ether_amide_amide", "target": "COc1ccc(NC(C)=O)cc1",
     "query": "[C](=[O])-[N]", "endpoints": [0, 2]},
    {"name": "lactam_ring", "target": "O=C1CCCN1",
     "query": "[C;R](=[O])-[N;R]", "endpoints": [0, 2]},
    {"name": "ide_target_amide", "target": "CCN(CC)CCNC(C1=COC(NC(C)=O)=N1)=O",
     "query": "[C](=[O])-[N]", "endpoints": [0, 2]},
)


def selected_target(case: dict[str, Any]) -> tuple[str, tuple[int, int]]:
    """Resolve a predefined canonical target bond without looking at search results."""
    target, _ = indexed_product(case["target"])
    if "bond" in case:
        return target, tuple(case["bond"])
    molecule = Chem.MolFromSmiles(target)
    query = compile_smarts(case["query"], validate=True)
    matches = sorted(molecule.GetSubstructMatches(query))
    if not matches:
        raise ValueError(f"Development case has no selected bond: {case['name']}")
    match = matches[0]
    return target, tuple(sorted(match[i] for i in case["endpoints"]))


def independently_matches(candidate: Any, target: str, bond: tuple[int, int]) -> bool:
    """Audit generator-labelled structures directly, independently of focus helpers."""
    generated = rdchiralReactants(target)
    selected = tuple(generated.idx_to_mapnum(i) for i in bond)
    before, after = candidate.condition_query_reaction_smiles.split(">>")
    reactants, product = Chem.MolFromSmiles(before), Chem.MolFromSmiles(after)
    if Chem.MolToSmiles(product) != Chem.MolToSmiles(generated.reactants):
        return False
    lookup = {atom.GetAtomMapNum(): atom.GetIdx() for atom in reactants.GetAtoms()
              if atom.GetAtomMapNum() > 0}
    if not set(selected) <= lookup.keys():
        return False
    return (candidate.forward_validation_status == "verified_signature"
            and reactants.GetBondBetweenAtoms(*(lookup[number] for number in selected)) is None)


def evaluate_pair(parameters: dict[str, Any]) -> dict[str, Any]:
    """Run both arms against one pinned library with equal total numerical budgets."""
    from core_retrosynthesis import load_generic_library, disconnect_generic_target_detailed
    from core_retrosynthesis.generic_compiler import clear_generic_compilation_cache
    from core_retrosynthesis.generic_search import _reaction, _strategic_complexity
    from core_retrosynthesis.search import _forward_analysis

    loaded = perf_counter()
    library = load_generic_library(parameters["library"])
    load_seconds = perf_counter() - loaded
    target, bond = selected_target(parameters["case"])
    result = {"case": parameters["case"]["name"], "target": target, "bond": bond,
              "library_load_seconds": load_seconds, "budgets": parameters["budgets"],
              "arm_order": parameters["arm_order"], "arms": {}}
    for arm in parameters["arm_order"]:
        clear_generic_compilation_cache()
        _forward_analysis.cache_clear()
        _strategic_complexity.cache_clear()
        _reaction.cache_clear()
        options = ({"required_disconnection_bond": bond, "focus_target_smiles": target}
                   if arm == "focused" else {})
        started = perf_counter()
        candidates, diagnostics = disconnect_generic_target_detailed(
            target, library, **parameters["budgets"], **options,
        )
        seconds = perf_counter() - started
        matches = [independently_matches(candidate, target, bond) for candidate in candidates]
        result["arms"][arm] = {
            "operation_seconds": seconds, "returned_count": len(candidates),
            "matching_count": sum(matches), "has_matching_candidate": any(matches),
            "diagnostics": diagnostics.to_dict(),
            "candidates": [{"precursor_smiles": candidate.precursor_smiles,
                            "matches_selected_bond": match,
                            "strategy_id": candidate.strategy_id,
                            "template_id": candidate.template_id,
                            "bond_focus_check": (asdict(candidate.bond_focus_check)
                                                 if candidate.bond_focus_check else None)}
                           for candidate, match in zip(candidates, matches)],
            "focus_rejection_examples": [asdict(check) for check in diagnostics.focus_rejection_examples],
        }
    return result


def main() -> int:
    """Save scientific calls, an independent graph audit and the paired report locally."""
    from chem_coworker.scientific_workspace import ScientificWorkspace

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--library", required=True)
    parser.add_argument("--templates", type=int, default=40)
    parser.add_argument("--validations", type=int, default=10)
    parser.add_argument("--top-k", type=int, default=3)
    args = parser.parse_args()
    library = str(Path(args.library).resolve())
    workspace = ScientificWorkspace.create(
        args.output, repository=Path(__file__).resolve().parents[2],
        objective="Paired development evaluation of a required retrosynthetic bond",
        artifacts={"retro_library": library},
        constraints=("Authored development probes; no feasibility or held-out evaluation claim",),
    )
    # Snapshot the complete evaluator, not a wrapper importing unpinned example code.
    (workspace.store.root / "evaluate_pair.py").write_text(Path(__file__).read_text("utf-8"), "utf-8")
    budgets = {"max_templates_to_apply": args.templates,
               "max_candidates_to_validate": args.validations, "top_k": args.top_k}
    rows, failures, refs = [], [], []
    for ordinal, case in enumerate(CASES):
        event = workspace.run_python("evaluate_pair.py", {
            "case": case, "library": library, "budgets": budgets,
            "arm_order": ["focused", "unrestricted"] if ordinal % 2 else ["unrestricted", "focused"],
        }, timeout_seconds=120)
        refs.append(event.artifact_ref)
        payload = workspace.store.read_artifact(event.artifact_ref)
        if payload["execution_status"] != "completed":
            failures.append({"case": case["name"], "artifact_ref": event.artifact_ref,
                             "error": payload.get("error")})
            print(json.dumps(failures[-1]), flush=True)
            continue
        row = payload["result"]
        row["artifact_ref"] = event.artifact_ref
        rows.append(row)
        print(json.dumps({"case": case["name"], "matching": {
            arm: value["matching_count"] for arm, value in row["arms"].items()}}), flush=True)
    summary = {}
    for arm in ("unrestricted", "focused"):
        values = [row["arms"][arm] for row in rows]
        summary[arm] = {
            "cases_with_matching_candidate": sum(value["has_matching_candidate"] for value in values),
            "matching_candidates": sum(value["matching_count"] for value in values),
            "returned_candidates": sum(value["returned_count"] for value in values),
            "validation_attempts": sum(value["diagnostics"]["validation_attempt_count"] for value in values),
            "median_operation_seconds": median(value["operation_seconds"] for value in values) if values else None,
        }
    report = {"schema_version": "focused_retrosynthesis_development.v1", "development_only": True,
              "scope": "single_pool_equal_total_budgets_not_ladder_or_agent_benchmark",
              "limitations": ["Authored probes, not independent chemistry review or an untouched evaluation.",
                              "Library coverage and experimental feasibility are not established.",
                              "Timing is observational; arm order is counterbalanced and major search caches cleared."],
              "budgets": budgets, "summary": summary, "cases": rows, "failures": failures}
    path = workspace.store.root / "report.json"
    path.write_text(json.dumps(report, indent=2), "utf-8")
    workspace.store.attach_file(path, description="Focused retrosynthesis development comparison",
                                evidence_refs=tuple(refs))
    print(json.dumps({"report": str(path), "summary": summary, "failures": len(failures)}), flush=True)
    return 1 if failures or any(row["arms"]["focused"]["returned_count"] !=
                               row["arms"]["focused"]["matching_count"] for row in rows) else 0


if __name__ == "__main__":
    if len(sys.argv) == 3 and Path(sys.argv[1]).name == "input.json":
        parameters = json.loads(Path(sys.argv[1]).read_text("utf-8"))["parameters"]
        Path(sys.argv[2]).write_text(json.dumps(evaluate_pair(parameters)), encoding="utf-8")
    else:
        raise SystemExit(main())

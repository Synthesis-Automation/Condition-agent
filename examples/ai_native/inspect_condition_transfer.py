"""Derived pilot analysis of substrate differences and evidence coverage.

This script compares saved graph descriptors and source completeness. It does not
estimate chemical reactivity, propose a new recipe, or change recommendation ranks.
"""

from __future__ import annotations

import argparse
import json
from typing import Any

from rdkit import Chem

from chem_coworker.scientific_workspace import ScientificWorkspace


def reactant_summary(reaction: str) -> list[dict[str, Any]]:
    """Describe supplied reactant graphs without assigning mechanism or reactive roles."""
    summaries = []
    for component in reaction.split(">")[0].split("."):
        molecule = Chem.MolFromSmiles(component)
        if molecule is None:
            summaries.append({"smiles": component, "parse_status": "invalid"})
            continue
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(0)
        symbols = [atom.GetSymbol() for atom in molecule.GetAtoms()]
        summaries.append({
            "smiles": Chem.MolToSmiles(molecule),
            "element_counts": {symbol: symbols.count(symbol) for symbol in sorted(set(symbols))},
            "ring_count": molecule.GetRingInfo().NumRings(),
        })
    return sorted(summaries, key=lambda value: value["smiles"])


def main() -> None:
    """Inspect a recorded recommendation and save the exact derived comparison."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workspace")
    parser.add_argument("recommendation_ref")
    args = parser.parse_args()
    workspace = ScientificWorkspace(args.workspace)
    call = workspace.store.read_artifact(args.recommendation_ref)
    if call.get("operation") != "recommend_conditions" or call.get("execution_status") != "completed":
        raise ValueError("Expected a completed recommendation call")
    result = call["result"]
    query = reactant_summary(call["arguments"]["reaction_smiles"])
    comparisons = []
    for recipe in result["recommendations"]:
        contexts = recipe.get("precedent_reaction_contexts", [])
        comparisons.append({
            "recipe_id": recipe["recipe_id"],
            "reported_support": recipe.get("support"),
            "independent_support": recipe.get("score_trace", {}).get("independent_evidence_count"),
            "missing_recipe_fields": [key for key in ("temperature_c", "time_h", "concentration_m", "atmosphere")
                                      if recipe["resolved_recipe"].get(key) is None],
            "inspected_contexts": [{
                "reaction_id": item["reaction_id"], "reference_id": item.get("reference_id"),
                "reactants": reactant_summary(item["reaction_smiles"]),
                "same_reactants_as_query": reactant_summary(item["reaction_smiles"]) == query,
            } for item in contexts],
            "cautions": recipe.get("cautions", []),
        })
    script = workspace.store.attach_file(
        __file__, description="Derived graph-description and evidence-coverage comparison",
        evidence_refs=(args.recommendation_ref,),
    )
    event = workspace.store.append("derived_analysis", {
        "origin": "derived_analysis", "review_status": "unreviewed",
        "query_reactants": query, "comparisons": comparisons,
        "limitations": ["Only stored representative precedent contexts were inspected",
                        "Graph differences do not establish reaction failure or a superior recipe"],
    }, evidence_refs=(args.recommendation_ref, script.artifact_ref))
    print(json.dumps({"artifact_ref": event.artifact_ref,
                      "recipes": [{k:row[k] for k in ("recipe_id", "reported_support", "independent_support", "missing_recipe_fields")}
                                  for row in comparisons]}))


if __name__ == "__main__":
    main()

"""Benchmark helper for taxonomy-first general reaction inference."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List

from .detection import detect_reaction_type
from .reaction_inference import analyze_reaction_general


@dataclass(frozen=True)
class BenchmarkCase:
    name: str
    reaction_smiles: str
    expected_reaction_type: str


_DEFAULT_CASES: List[BenchmarkCase] = [
    BenchmarkCase(
        name="amide_formation",
        reaction_smiles="O=C(O)c1ccccc1.NCc1ccccc1>>O=C(NCc1ccccc1)c1ccccc1",
        expected_reaction_type="Amide_formation",
    ),
    BenchmarkCase(
        name="oxidation_primary_alcohol",
        reaction_smiles="CCO>>CC=O",
        expected_reaction_type="Oxidation_primary_alcohol_to_aldehyde",
    ),
    BenchmarkCase(
        name="alkyl_nucleophilic_substitution",
        reaction_smiles="CN.CBr>>CNC",
        expected_reaction_type="Alkyl_Nucleophilic_Substitution",
    ),
    BenchmarkCase(
        name="suzuki_miyaura",
        reaction_smiles="Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
        expected_reaction_type="Suzuki_miyaura",
    ),
    BenchmarkCase(
        name="taxonomy_gap_snar_like",
        reaction_smiles="Clc1ncc(-c2ccccc2)cn1.NN>>NN=c1ncc(-c2ccccc2)c[nH]1",
        expected_reaction_type="unknown",
    ),
]


def _baseline_predict(reaction_smiles: str) -> str:
    detection = detect_reaction_type(reaction_smiles)
    top = detection.top_match
    return str(top.reaction_type) if top else "unknown"


def run_benchmark() -> Dict[str, Any]:
    """
    Compare baseline taxonomy detection with the general inference workflow.

    Baseline ("v1_mapped") is direct `detect_reaction_type`.
    V2 is `analyze_reaction_general`.
    """
    rows: List[Dict[str, Any]] = []
    baseline_correct = 0
    v2_correct = 0
    for case in _DEFAULT_CASES:
        expected = case.expected_reaction_type
        v1_pred = _baseline_predict(case.reaction_smiles)
        v2_pred = analyze_reaction_general(case.reaction_smiles).decision.reaction_type
        v1_exact = v1_pred == expected
        v2_exact = v2_pred == expected
        baseline_correct += int(v1_exact)
        v2_correct += int(v2_exact)
        rows.append(
            {
                "name": case.name,
                "reaction_smiles": case.reaction_smiles,
                "expected_reaction_type": expected,
                "v1_mapped_prediction": v1_pred,
                "v2_prediction": v2_pred,
                "v1_exact": v1_exact,
                "v2_exact": v2_exact,
            }
        )

    n_cases = max(1, len(_DEFAULT_CASES))
    return {
        "summary": {
            "n_cases": len(_DEFAULT_CASES),
            "v1_mapped_exact_accuracy": baseline_correct / float(n_cases),
            "v2_exact_accuracy": v2_correct / float(n_cases),
        },
        "cases": rows,
    }


__all__ = ["run_benchmark"]


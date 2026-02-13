"""Small deterministic benchmark: v1 family-specific PoC vs v2 general PoC."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from statistics import mean
from typing import Any, Dict, List, Optional

from poc_gpt52_reaction import analyze_reaction as analyze_v1
from poc_gpt52_reaction_v2 import analyze_reaction_general as analyze_v2


@dataclass(frozen=True)
class BenchmarkCase:
    """One curated benchmark entry."""

    case_id: str
    reaction_smiles: str
    expected_taxonomy: str
    note: str


@dataclass
class CaseResult:
    """Combined v1 + v2 result for one case."""

    case_id: str
    expected_taxonomy: str
    v1_raw: str
    v1_confidence: float
    v1_mapped_taxonomy: str
    v2_decision: str
    v2_confidence: float
    v2_candidate_count: int
    v2_detection_error: Optional[str]
    v1_correct: bool
    v2_correct: bool


CASES: List[BenchmarkCase] = [
    BenchmarkCase(
        case_id="snar_hydrazinolysis",
        reaction_smiles="Clc1ncc(-c2ccccc2)cn1.NN>>NN=c1ncc(-c2ccccc2)c[nH]1",
        expected_taxonomy="SNAr_CN",
        note="Taxonomy has SNAr_CN, but detector currently returns no candidate.",
    ),
    BenchmarkCase(
        case_id="amide_formation",
        reaction_smiles="O=C(O)c1ccccc1.NCc1ccccc1>>O=C(NCc1ccccc1)c1ccccc1",
        expected_taxonomy="Amide_formation",
        note="Typical acid + amine coupling.",
    ),
    BenchmarkCase(
        case_id="oxidation_primary_alcohol",
        reaction_smiles="CCO>>CC=O",
        expected_taxonomy="Oxidation_primary_alcohol_to_aldehyde",
        note="Primary alcohol oxidation.",
    ),
    BenchmarkCase(
        case_id="reduction_carbonyl",
        reaction_smiles="O=CC1=CC=CC=C1>>OCC1=CC=CC=C1",
        expected_taxonomy="Reduction_carbonyl_to_alcohol",
        note="Carbonyl reduction.",
    ),
    BenchmarkCase(
        case_id="alkyl_substitution",
        reaction_smiles="CCO.CBr>>CCOC",
        expected_taxonomy="Alkyl_Nucleophilic_Substitution",
        note="Simple alkyl substitution.",
    ),
    BenchmarkCase(
        case_id="acyl_halide_formation",
        reaction_smiles="O=C(O)CN1C(=O)c2ccccc2C1=O>>O=C(Cl)CN1C(=O)c2ccccc2C1=O",
        expected_taxonomy="Acyl_Halides_formation",
        note="Acid to acyl chloride conversion.",
    ),
    BenchmarkCase(
        case_id="cn_coupling",
        reaction_smiles="Brc1ccccc1.Nc1ccccc1>>Nc1ccccc1c1ccccc1",
        expected_taxonomy="C_N_Coupling",
        note="Representative aryl-halide amination.",
    ),
]


def _map_v1_to_taxonomy(v1_label: str) -> str:
    if v1_label == "SNAr_hydrazinolysis":
        return "SNAr_CN"
    return "unknown"


def run_benchmark() -> Dict[str, Any]:
    case_results: List[CaseResult] = []
    v1_fallback_confidences: List[float] = []

    for case in CASES:
        v1 = analyze_v1(case.reaction_smiles, use_llm=False)
        v1_raw = str(v1.best_hypothesis.get("reaction_class", "unknown"))
        v1_conf = float(v1.confidence)
        v1_mapped = _map_v1_to_taxonomy(v1_raw)
        if v1_raw == "other_amination_or_annotation_error":
            v1_fallback_confidences.append(v1_conf)

        v2 = analyze_v2(case.reaction_smiles, use_llm=False)
        v2_label = str(v2.decision.reaction_type)
        v2_conf = float(v2.decision.confidence)

        case_results.append(
            CaseResult(
                case_id=case.case_id,
                expected_taxonomy=case.expected_taxonomy,
                v1_raw=v1_raw,
                v1_confidence=round(v1_conf, 2),
                v1_mapped_taxonomy=v1_mapped,
                v2_decision=v2_label,
                v2_confidence=round(v2_conf, 2),
                v2_candidate_count=len(v2.taxonomy_candidates),
                v2_detection_error=v2.detection_error,
                v1_correct=(v1_mapped == case.expected_taxonomy),
                v2_correct=(v2_label == case.expected_taxonomy),
            )
        )

    total = len(case_results)
    v1_correct = sum(1 for row in case_results if row.v1_correct)
    v2_correct = sum(1 for row in case_results if row.v2_correct)
    v1_specific = sum(1 for row in case_results if row.v1_raw != "other_amination_or_annotation_error")
    v2_unknown = sum(1 for row in case_results if row.v2_decision == "unknown")
    v2_with_candidates = sum(1 for row in case_results if row.v2_candidate_count > 0)

    summary = {
        "n_cases": total,
        "v1_mapped_exact_accuracy": round(v1_correct / total, 3),
        "v2_exact_accuracy": round(v2_correct / total, 3),
        "v1_specific_call_rate": round(v1_specific / total, 3),
        "v2_unknown_rate": round(v2_unknown / total, 3),
        "v2_candidate_coverage": round(v2_with_candidates / total, 3),
        "v1_fallback_mean_confidence": round(mean(v1_fallback_confidences), 3)
        if v1_fallback_confidences
        else None,
        "v2_failure_cases": [
            row.case_id for row in case_results if not row.v2_correct
        ],
        "v1_failure_cases": [
            row.case_id for row in case_results if not row.v1_correct
        ],
    }

    return {
        "summary": summary,
        "cases": [asdict(row) for row in case_results],
        "case_definitions": [asdict(case) for case in CASES],
    }


def _print_human(report: Dict[str, Any]) -> None:
    summary = report["summary"]
    print("Benchmark summary")
    print(f"  Cases: {summary['n_cases']}")
    print(f"  v1 mapped exact accuracy: {summary['v1_mapped_exact_accuracy']:.3f}")
    print(f"  v2 exact accuracy: {summary['v2_exact_accuracy']:.3f}")
    print(f"  v1 specific-call rate: {summary['v1_specific_call_rate']:.3f}")
    print(f"  v2 unknown rate: {summary['v2_unknown_rate']:.3f}")
    print(f"  v2 candidate coverage: {summary['v2_candidate_coverage']:.3f}")
    if summary["v1_fallback_mean_confidence"] is not None:
        print(f"  v1 fallback mean confidence: {summary['v1_fallback_mean_confidence']:.3f}")

    print("\nPer-case results")
    header = (
        "case_id | expected | v1_raw | v1_mapped | v2 | "
        "v1_ok | v2_ok | v2_candidates"
    )
    print(header)
    print("-" * len(header))
    for row in report["cases"]:
        print(
            f"{row['case_id']} | {row['expected_taxonomy']} | {row['v1_raw']} | "
            f"{row['v1_mapped_taxonomy']} | {row['v2_decision']} | "
            f"{row['v1_correct']} | {row['v2_correct']} | {row['v2_candidate_count']}"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description="Run v1 vs v2 reaction analysis benchmark.")
    parser.add_argument("--json", action="store_true", help="Print JSON report.")
    args = parser.parse_args()

    report = run_benchmark()
    if args.json:
        print(json.dumps(report, ensure_ascii=True, indent=2, sort_keys=True))
    else:
        _print_human(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

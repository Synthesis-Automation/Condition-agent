"""
Reaction evaluation ToolPlugins for ChemCoworker.

Two tools:

  evaluate_reaction        — full 8-check suite on any reaction SMILES
  check_retro_consistency  — retrosynthesis-specific shortcut (product + two precursors)

Both return a structured EvalReport dict with:
  - verdict: "PASS" | "PASS_WITH_WARNINGS" | "FAIL"
  - per-check results (name, severity, score, message, details)
  - critical_failures / warnings lists
  - llm_eval_prompt: pre-built brief for LLM expert follow-up

The LLM is expected to supplement these RDKit checks with its own assessment
of regiochemistry, mechanism, stereochemistry, and practical feasibility.
"""
from __future__ import annotations

from typing import Any, Dict

from ._helpers import _clean_rxn_smiles, _error, _success
from ._base import ToolPlugin


# ---------------------------------------------------------------------------
# Tool A: evaluate_reaction
# ---------------------------------------------------------------------------

def _evaluate_reaction(reaction_smiles: str, reaction_type: str = "") -> Dict[str, Any]:
    """Run all 8 RDKit-based sanity checks on a reaction SMILES string.

    Evaluates the forward reaction (reactants >> product) for:
      1. Component validity    — all SMILES parse and sanitize
      2. Atom balance          — no elements invented in product
      3. Charge balance        — formal charge conserved
      4. HDI consistency       — no unexplained gain of unsaturation
      5. MW direction          — product lighter than combined reactants
      6. Ring change           — Δ(rings) matches claimed reaction type
      7. Reaction-type pattern — correct FGs present (Suzuki: Ar-X + Ar-B)
      8. Fragment reattachment — precursor skeletons MCS-found in product

    Checks 6 and 7 are skipped when reaction_type is not provided or not in
    the known-type table. This is intentional — the LLM should assess unknown
    types via its own expert evaluation.

    Args:
        reaction_smiles: Reaction SMILES in 'reactant1.reactant2>>product' format.
                         Multiple reactants are dot-separated. Example:
                         "N#Cc1ccc(Br)cc1.OBOc1ccoc1>>N#Cc1ccc(-c2ccoc2)cc1"
        reaction_type:   Optional taxonomy reaction ID for type-specific checks.
                         Examples: "suzuki_miyaura", "buchwald_hartwig",
                         "reductive_amination", "amide_coupling", "diels_alder".
                         Leave empty for generic checks only.

    Returns:
        dict with verdict, overall_score, checks (list), critical_failures,
        warnings, and llm_eval_prompt for LLM follow-up expert assessment.
    """
    reaction_smiles = _clean_rxn_smiles(reaction_smiles)
    if not reaction_smiles or not reaction_smiles.strip():
        return _error("reaction_smiles cannot be empty")

    if ">>" not in reaction_smiles:
        return _error(
            "reaction_smiles must be in 'reactants>>product' format. "
            "For retrosynthesis evaluation, use the FORWARD direction: "
            "'precursor1.precursor2>>product'."
        )

    try:
        from chemtools.eval.reaction_evaluator import evaluate_reaction
        report = evaluate_reaction(
            reaction_smiles.strip(),
            reaction_type=reaction_type.strip() if reaction_type else None,
        )

        # Serialize the report dataclass to a plain dict
        checks_dicts = [
            {
                "name":     c.name,
                "passed":   c.passed,
                "severity": c.severity,
                "score":    c.score,
                "message":  c.message,
                "details":  c.details,
            }
            for c in report.checks
        ]

        return _success({
            "verdict":           report.verdict,
            "overall_score":     report.overall_score,
            "overall_passed":    report.overall_passed,
            "reaction_smiles":   report.reaction_smiles,
            "reaction_type":     report.reaction_type or "",
            "checks":            checks_dicts,
            "critical_failures": report.critical_failures,
            "warnings":          report.warnings,
            "llm_eval_prompt":   report.llm_eval_prompt,
            "summary": (
                f"{report.verdict} (score={report.overall_score:.2f}). "
                + (f"Failures: {'; '.join(report.critical_failures)}. " if report.critical_failures else "")
                + (f"Warnings: {'; '.join(report.warnings[:2])}." if report.warnings else "")
            ),
        })

    except Exception as exc:
        return _error(f"Reaction evaluation failed: {exc}")


evaluate_reaction_tool = ToolPlugin(
    name="evaluate_reaction",
    category="reaction_eval",
    description=(
        "Run 8 RDKit-based sanity checks on a reaction SMILES (reactants>>product): "
        "SMILES validity, atom balance (no invented elements), charge conservation, "
        "HDI/unsaturation consistency, MW direction, ring-count change, "
        "reaction-type FG patterns (Suzuki: aryl halide + boronate), and "
        "MCS fragment reattachment. Returns verdict (PASS/PASS_WITH_WARNINGS/FAIL), "
        "per-check details, and an llm_eval_prompt for follow-up expert assessment. "
        "Use after generate_disconnections to validate retrosynthetic proposals."
    ),
    prerequisites=[],
    fn=_evaluate_reaction,
)


# ---------------------------------------------------------------------------
# Tool B: check_retro_consistency
# ---------------------------------------------------------------------------

def _check_retro_consistency(
    product_smiles: str,
    precursor_1: str,
    precursor_2: str,
    reaction_name: str = "",
) -> Dict[str, Any]:
    """Validate a retrosynthetic disconnection: do these two precursors give this product?

    Assembles the forward reaction 'precursor_1.precursor_2>>product_smiles' and runs
    the full evaluate_reaction suite. Also checks that both precursors are simpler
    than the product (complexity reduction — the point of retrosynthesis).

    Args:
        product_smiles: The target molecule SMILES.
                        Example: "N#Cc1ccc(-c2ccoc2)cc1"
        precursor_1:    First precursor SMILES (electrophile side for coupling).
                        Example: "N#Cc1ccc(Br)cc1"
        precursor_2:    Second precursor SMILES (nucleophile side for coupling).
                        Example: "OBOc1ccoc1"
        reaction_name:  Optional reaction type ID for type-specific checks.
                        Example: "suzuki_miyaura"

    Returns:
        Same structure as evaluate_reaction, with an additional
        'complexity_check' field: {target_complexity, p1_complexity, p2_complexity,
        max_precursor_complexity, simplification_achieved}.
    """
    for label, val in [("product_smiles", product_smiles),
                        ("precursor_1", precursor_1),
                        ("precursor_2", precursor_2)]:
        if not val or not val.strip():
            return _error(f"'{label}' cannot be empty")

    try:
        # Build forward reaction SMILES
        rxn = f"{precursor_1.strip()}.{precursor_2.strip()}>>{product_smiles.strip()}"
        result = _evaluate_reaction(rxn, reaction_type=reaction_name)

        if not result.get("success"):
            return result

        # Additional complexity check (BertzCT)
        complexity_check: Dict[str, Any] = {}
        try:
            from rdkit import Chem
            from rdkit.Chem.GraphDescriptors import BertzCT

            def bertz(smi: str) -> float:
                mol = Chem.MolFromSmiles(smi)
                return round(BertzCT(mol), 1) if mol else 0.0

            t_cx  = bertz(product_smiles)
            p1_cx = bertz(precursor_1)
            p2_cx = bertz(precursor_2)
            max_p = max(p1_cx, p2_cx)
            delta = t_cx - max_p

            complexity_check = {
                "target_complexity":      t_cx,
                "precursor_1_complexity": p1_cx,
                "precursor_2_complexity": p2_cx,
                "max_precursor_complexity": max_p,
                "complexity_delta":       round(delta, 1),
                "simplification_achieved": delta > 0,
                "note": (
                    f"Target ({t_cx}) > max(precursor)({max_p}): good simplification."
                    if delta > 0
                    else (
                        f"Target ({t_cx}) ≤ max(precursor)({max_p}): "
                        "precursors not simpler than target — check disconnection."
                    )
                ),
            }
        except Exception as cx_err:
            complexity_check = {"error": str(cx_err)}

        result["complexity_check"] = complexity_check
        result["reaction_smiles_used"] = rxn
        return result

    except Exception as exc:
        return _error(f"Retro consistency check failed: {exc}")


check_retro_consistency_tool = ToolPlugin(
    name="check_retro_consistency",
    category="reaction_eval",
    description=(
        "Validate a retrosynthetic disconnection: given a product and two precursor SMILES, "
        "runs the full 8-check evaluation suite on 'precursor_1.precursor_2>>product' plus "
        "a BertzCT complexity check to confirm precursors are simpler than the target. "
        "Use immediately after generate_disconnections to catch atom balance errors, "
        "wrong FG assignments, and structural mismatches before writing the final answer."
    ),
    prerequisites=[],
    fn=_check_retro_consistency,
)


# ---------------------------------------------------------------------------
# Module export
# ---------------------------------------------------------------------------

EVAL_TOOLS = [
    evaluate_reaction_tool,
    check_retro_consistency_tool,
]

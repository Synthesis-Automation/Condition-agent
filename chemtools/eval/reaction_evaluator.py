"""
RDKit-based reaction evaluation engine.

Runs up to 8 computable sanity checks on a reaction SMILES string
(format: "reactant1.reactant2>>product"). Each check returns a CheckResult
with a severity level, score, human-readable message, and raw details dict
for LLM reasoning.

Design philosophy
-----------------
- Checks verify HARD facts RDKit can prove; they do NOT replace LLM expertise
- Leniency by default: reactants may have extra atoms (omitted reagents OK)
- Severity tiers: "error" (definite violation) vs "warning" (suspicious)
- The EvalReport.llm_eval_prompt field gives the LLM a structured brief
  for its own mechanistic/regiochemical assessment

Usage
-----
    from chemtools.eval.reaction_evaluator import evaluate_reaction

    report = evaluate_reaction(
        "N#Cc1ccc(Br)cc1.OBOc1ccoc1>>N#Cc1ccc(-c2ccoc2)cc1",
        reaction_type="suzuki_miyaura",
    )
    print(report.verdict)          # "PASS" / "PASS_WITH_WARNINGS" / "FAIL"
    print(report.llm_eval_prompt)  # ready-made brief for LLM expert assessment
"""
from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Reaction-type knowledge tables
# ---------------------------------------------------------------------------

# Expected ring-count Δ = rings(product) − sum(rings(reactants))
_EXPECTED_RING_DELTA: Dict[str, int] = {
    "suzuki_miyaura":           0,
    "buchwald_hartwig":         0,
    "ullmann_coupling":         0,
    "negishi_coupling":         0,
    "heck_reaction":            0,
    "sonogashira_coupling":     0,
    "reductive_amination":      0,
    "amide_coupling":           0,
    "mitsunobu":                0,
    "williamson_ether":         0,
    "chan_lam_coupling":        0,
    "diels_alder":              1,
    "ring_closing_metathesis":  1,
    "staudinger_ketene":        1,
    "pauson_khand":             1,
}

# SMARTS patterns that MUST appear in the reactants for the claimed reaction type.
# Format: {rxn_type: {"r1": [smarts_list], "r2": [smarts_list]}}
# r1 / r2 are searched across ALL reactants (not positionally fixed).
# A reaction passes if at least one reactant matches r1 AND at least one matches r2.
_REACTION_PATTERNS: Dict[str, Dict[str, List[str]]] = {
    "suzuki_miyaura": {
        "r1": ["[c,C][Br,I]", "[c,C][Cl]", "[c,C]OC(F)(F)F"],  # aryl/vinyl halide or triflate
        "r2": ["[c,C]B", "[B](O)O"],                              # boronate / boronic acid
    },
    "buchwald_hartwig": {
        "r1": ["[c][Br,I,Cl]"],              # aryl halide
        "r2": ["[NH1,NH2]", "[nH]"],          # primary/secondary amine or NH heterocycle
    },
    "ullmann_coupling": {
        "r1": ["[c][Br,I,Cl]"],              # aryl halide
        "r2": ["[NH1,NH2]", "[OH1]", "[SH1]"],
    },
    "negishi_coupling": {
        "r1": ["[c,C][Br,I,Cl]"],
        "r2": ["[Zn]"],
    },
    "heck_reaction": {
        "r1": ["[c,C][Br,I,Cl]"],
        "r2": ["C=C"],
    },
    "sonogashira_coupling": {
        "r1": ["[c,C][Br,I,Cl]"],
        "r2": ["[C]#[C]", "C#C"],
    },
    "reductive_amination": {
        "r1": ["[CX3](=O)", "C=O"],          # aldehyde or ketone
        "r2": ["[NH1,NH2]"],                  # amine
    },
    "amide_coupling": {
        "r1": ["C(=O)[OH1]", "C(=O)Cl", "C(=O)OC(=O)"],   # carboxylic acid / acyl chloride / anhydride
        "r2": ["[NH1,NH2]"],
    },
    "mitsunobu": {
        "r1": ["[OH1][CX4]"],                # alcohol (aliphatic)
        "r2": ["OC(=O)", "[NH1,NH2]"],        # carboxylic acid or amine pronucleophile
    },
    "williamson_ether": {
        "r1": ["[OH1,O-]"],                  # alcohol or alkoxide
        "r2": ["[C][Br,I,Cl]", "[C]OS(=O)"],  # alkyl halide or sulfonate
    },
    "wittig": {
        "r1": ["C=O"],                        # aldehyde or ketone
        "r2": ["[P+]", "P=C"],               # phosphonium ylide
    },
    "chan_lam_coupling": {
        "r1": ["[NH1,NH2]", "OH"],
        "r2": ["[c,C]B"],                    # aryl boronate
    },
    "diels_alder": {
        "r1": ["C=CC=C", "c"],              # diene (conjugated)
        "r2": ["C=C"],                       # dienophile
    },
}

# HDI tolerance: allow product HDI to exceed reactant HDI sum by at most this value
# (Small positive tolerance handles imprecision from omitted H atoms in SMILES)
_HDI_TOLERANCE = 1.5


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class CheckResult:
    """Result of a single evaluation check."""
    name: str
    passed: bool
    severity: str        # "error" | "warning" | "info"
    score: float         # 0.0 (fail) – 1.0 (perfect)
    message: str
    details: Dict = field(default_factory=dict)


@dataclass
class EvalReport:
    """Aggregated evaluation report for a reaction SMILES."""
    reaction_smiles: str
    reaction_type: Optional[str]
    checks: List[CheckResult]
    overall_passed: bool    # True if zero error-severity failures
    overall_score: float    # weighted mean of check scores
    critical_failures: List[str]
    warnings: List[str]
    verdict: str            # "PASS" | "PASS_WITH_WARNINGS" | "FAIL"
    llm_eval_prompt: str    # structured brief for LLM expert follow-up


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _parse_reaction(reaction_smiles: str) -> Tuple[List[str], List[str]]:
    """Split 'A.B.C>>D.E' into ([reactants], [products]) SMILES lists."""
    if ">>" not in reaction_smiles:
        raise ValueError("Reaction SMILES must contain '>>'")
    parts = reaction_smiles.split(">>")
    reactants_str = parts[0]
    products_str  = parts[-1]

    def split_smiles(s: str) -> List[str]:
        # Split on '.' but not inside brackets or ring closures
        result, buf, depth = [], [], 0
        for ch in s:
            if ch == '.' and depth == 0:
                if buf:
                    result.append("".join(buf))
                buf = []
            else:
                if ch in "([":
                    depth += 1
                elif ch in ")]":
                    depth -= 1
                buf.append(ch)
        if buf:
            result.append("".join(buf))
        return [x for x in result if x.strip()]

    return split_smiles(reactants_str), split_smiles(products_str)


def _mol(smiles: str):
    """Parse SMILES, return mol or None."""
    try:
        from rdkit import Chem
        return Chem.MolFromSmiles(smiles)
    except Exception:
        return None


def _mol_with_h(smiles: str):
    """Parse SMILES and add explicit H."""
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        return Chem.AddHs(mol) if mol else None
    except Exception:
        return None


def _atom_counts(smiles_list: List[str], include_h: bool = True) -> Counter:
    """Count all atoms across a list of SMILES."""
    total: Counter = Counter()
    for smi in smiles_list:
        mol = _mol_with_h(smi) if include_h else _mol(smi)
        if mol:
            for atom in mol.GetAtoms():
                total[atom.GetSymbol()] += 1
    return total


def _charge_sum(smiles_list: List[str]) -> int:
    """Sum formal charges across all atoms in a SMILES list."""
    total = 0
    for smi in smiles_list:
        mol = _mol(smi)
        if mol:
            total += sum(a.GetFormalCharge() for a in mol.GetAtoms())
    return total


def _hdi(atom_counts: Counter) -> float:
    """Hydrogen Deficiency Index = (2C + 2 + N - H - halogens) / 2."""
    C = atom_counts.get("C", 0)
    H = atom_counts.get("H", 0)
    N = atom_counts.get("N", 0)
    X = sum(atom_counts.get(e, 0) for e in ("F", "Cl", "Br", "I"))
    return (2 * C + 2 + N - H - X) / 2.0


def _exact_mw(smiles: str) -> float:
    try:
        from rdkit.Chem import Descriptors
        mol = _mol(smiles)
        return Descriptors.ExactMolWt(mol) if mol else 0.0
    except Exception:
        return 0.0


def _ring_count(smiles: str) -> int:
    try:
        from rdkit.Chem import rdMolDescriptors
        mol = _mol(smiles)
        return rdMolDescriptors.CalcNumRings(mol) if mol else 0
    except Exception:
        return 0


def _smarts_match(smiles: str, smarts_list: List[str]) -> bool:
    """Return True if the molecule matches ANY SMARTS in the list."""
    try:
        from rdkit import Chem
        mol = _mol(smiles)
        if mol is None:
            return False
        for sma in smarts_list:
            try:
                patt = Chem.MolFromSmarts(sma)
                if patt and mol.HasSubstructMatch(patt):
                    return True
            except Exception:
                continue
    except Exception:
        pass
    return False


def _mcs_coverage(query_smiles: str, target_smiles: str) -> float:
    """
    Fraction of heavy atoms in query_smiles covered by MCS with target_smiles.
    Returns 0.0 on failure.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdFMCS
        mol_q = _mol(query_smiles)
        mol_t = _mol(target_smiles)
        if mol_q is None or mol_t is None:
            return 0.0
        q_atoms = mol_q.GetNumHeavyAtoms()
        if q_atoms == 0:
            return 1.0
        result = rdFMCS.FindMCS(
            [mol_q, mol_t],
            timeout=3,
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            ringMatchesRingOnly=False,
            completeRingsOnly=False,
        )
        return result.numAtoms / q_atoms if result.numAtoms > 0 else 0.0
    except Exception:
        return 0.0


# ---------------------------------------------------------------------------
# The 8 checks
# ---------------------------------------------------------------------------

def _check_validate_components(reactants: List[str], products: List[str]) -> CheckResult:
    """Check 1: All SMILES components parse and sanitize correctly."""
    bad = []
    for smi in reactants + products:
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                bad.append(smi)
            else:
                Chem.SanitizeMol(mol)
        except Exception:
            bad.append(smi)

    passed = len(bad) == 0
    return CheckResult(
        name="validate_components",
        passed=passed,
        severity="error",
        score=1.0 if passed else 0.0,
        message=(
            "All SMILES components are valid."
            if passed
            else f"Invalid/unparseable SMILES: {bad}"
        ),
        details={"invalid_smiles": bad, "total_components": len(reactants + products)},
    )


def _check_atom_balance(reactants: List[str], products: List[str]) -> CheckResult:
    """
    Check 2: Atom balance (lenient).
    - ERROR if product contains an element absent from all reactants.
    - WARNING if a reactant element appears 0× in product (lost atoms, but could be OK).
    - Reactants are allowed to have MORE atoms than product (omitted reagents, byproducts).
    """
    r_counts = _atom_counts(reactants)
    p_counts = _atom_counts(products)

    invented: List[str] = []   # in product but not reactants → hard error
    lost: List[str] = []       # in reactants but not product → usually OK

    for elem, count in p_counts.items():
        if elem == "H":
            continue  # H is often implicit; skip
        if r_counts.get(elem, 0) == 0:
            invented.append(elem)

    for elem, count in r_counts.items():
        if elem == "H":
            continue
        if p_counts.get(elem, 0) == 0:
            lost.append(elem)

    if invented:
        return CheckResult(
            name="atom_balance",
            passed=False,
            severity="error",
            score=0.0,
            message=f"Elements appear in product but NOT in reactants (atoms invented): {invented}",
            details={"invented_elements": invented, "lost_elements": lost,
                     "reactant_counts": dict(r_counts), "product_counts": dict(p_counts)},
        )

    # Compute per-element ratios for scoring
    score = 1.0
    if lost:
        score = 0.8  # warn: lost atoms, but OK if reagents are omitted
    return CheckResult(
        name="atom_balance",
        passed=True,
        severity="warning" if lost else "info",
        score=score,
        message=(
            f"Atom balance OK. Elements lost to byproducts (expected): {lost}"
            if lost else "Atom balance: all product atoms accounted for in reactants."
        ),
        details={"invented_elements": [], "lost_elements": lost,
                 "reactant_counts": dict(r_counts), "product_counts": dict(p_counts)},
    )


def _check_charge_balance(reactants: List[str], products: List[str]) -> CheckResult:
    """Check 3: Formal charge must be conserved exactly."""
    r_charge = _charge_sum(reactants)
    p_charge = _charge_sum(products)
    delta = p_charge - r_charge
    passed = delta == 0
    return CheckResult(
        name="charge_balance",
        passed=passed,
        severity="error" if not passed else "info",
        score=1.0 if passed else 0.0,
        message=(
            f"Charge balanced (both sides = {r_charge})."
            if passed
            else f"Charge imbalance: reactants={r_charge}, products={p_charge}, Δ={delta:+d}"
        ),
        details={"reactant_charge": r_charge, "product_charge": p_charge, "delta": delta},
    )


def _check_hdi_consistency(reactants: List[str], products: List[str]) -> CheckResult:
    """
    Check 4: Hydrogen Deficiency Index consistency.
    Product HDI cannot exceed sum of reactant HDIs by more than _HDI_TOLERANCE,
    because you cannot gain degrees of unsaturation (rings, pi bonds) in a
    simple bond-forming step without a ring-forming or aromatization event.
    """
    r_hdi = sum(_hdi(_atom_counts([smi])) for smi in reactants)
    p_hdi = sum(_hdi(_atom_counts([smi])) for smi in products)
    delta = p_hdi - r_hdi  # positive = product MORE unsaturated than expected

    passed = delta <= _HDI_TOLERANCE
    score = max(0.0, 1.0 - max(0.0, delta - _HDI_TOLERANCE) / 5.0)

    return CheckResult(
        name="hdi_consistency",
        passed=passed,
        severity="warning" if not passed else "info",
        score=score,
        message=(
            f"HDI consistent (reactants={r_hdi:.1f}, products={p_hdi:.1f}, Δ={delta:+.1f})."
            if passed
            else (
                f"HDI anomaly: products have {delta:.1f} more degrees of unsaturation than reactants "
                f"(reactants={r_hdi:.1f}, products={p_hdi:.1f}). "
                f"Unexpected ring formation or aromatization? Tolerance=±{_HDI_TOLERANCE}."
            )
        ),
        details={"reactant_hdi": r_hdi, "product_hdi": p_hdi, "delta": delta,
                 "tolerance": _HDI_TOLERANCE},
    )


def _check_mw_direction(reactants: List[str], products: List[str]) -> CheckResult:
    """
    Check 5: Molecular weight direction.
    Product MW must be ≤ sum of reactant MWs. Bond formation always eliminates
    atoms (HX, H2O, etc.) so the product must be lighter than combined reactants.
    """
    r_mw = sum(_exact_mw(smi) for smi in reactants)
    p_mw = sum(_exact_mw(smi) for smi in products)
    delta = p_mw - r_mw  # positive = product heavier than reactants (violation)
    # Allow small tolerance for floating point
    passed = delta <= 0.5
    score = 1.0 if passed else max(0.0, 1.0 - delta / 50.0)

    return CheckResult(
        name="mw_direction",
        passed=passed,
        severity="warning" if not passed else "info",
        score=score,
        message=(
            f"MW direction OK (reactants={r_mw:.1f}, products={p_mw:.1f}, Δ={delta:+.1f} Da)."
            if passed
            else (
                f"MW violation: products ({p_mw:.1f} Da) heavier than reactants ({r_mw:.1f} Da) "
                f"by {delta:.1f} Da. Missing byproducts or incorrect SMILES?"
            )
        ),
        details={"reactant_mw": r_mw, "product_mw": p_mw, "delta_da": delta},
    )


def _check_ring_change(
    reactants: List[str], products: List[str], reaction_type: Optional[str]
) -> CheckResult:
    """
    Check 6: Ring count change vs expected for the claimed reaction type.
    Skipped if reaction_type is unknown or not in the expected-Δ table.
    """
    expected = _EXPECTED_RING_DELTA.get(reaction_type)
    if expected is None:
        return CheckResult(
            name="ring_change",
            passed=True,
            severity="info",
            score=1.0,
            message=f"Ring change check skipped (reaction type '{reaction_type}' not in table).",
            details={"reaction_type": reaction_type, "skipped": True},
        )

    r_rings = sum(_ring_count(smi) for smi in reactants)
    p_rings = sum(_ring_count(smi) for smi in products)
    actual_delta = p_rings - r_rings
    passed = actual_delta == expected

    return CheckResult(
        name="ring_change",
        passed=passed,
        severity="warning" if not passed else "info",
        score=1.0 if passed else 0.5,
        message=(
            f"Ring count change as expected for {reaction_type} "
            f"(reactants={r_rings}, products={p_rings}, Δ={actual_delta:+d})."
            if passed
            else (
                f"Unexpected ring change for {reaction_type}: "
                f"expected Δ={expected:+d}, got Δ={actual_delta:+d} "
                f"(reactants={r_rings} rings, products={p_rings} rings)."
            )
        ),
        details={
            "reaction_type": reaction_type,
            "reactant_rings": r_rings, "product_rings": p_rings,
            "actual_delta": actual_delta, "expected_delta": expected,
        },
    )


def _check_reaction_type_pattern(
    reactants: List[str], reaction_type: Optional[str]
) -> CheckResult:
    """
    Check 7: Required functional group patterns for the claimed reaction type.
    Most powerful check for catching regiochemistry and FG assignment errors.
    Skipped if reaction_type has no defined patterns.
    """
    patterns = _REACTION_PATTERNS.get(reaction_type)
    if not patterns:
        return CheckResult(
            name="reaction_type_pattern",
            passed=True,
            severity="info",
            score=1.0,
            message=f"Pattern check skipped (no patterns defined for '{reaction_type}').",
            details={"reaction_type": reaction_type, "skipped": True},
        )

    failures: List[str] = []
    matched: Dict[str, bool] = {}

    for role, smarts_list in patterns.items():
        # At least one reactant must match any SMARTS for this role
        role_matched = any(_smarts_match(smi, smarts_list) for smi in reactants)
        matched[role] = role_matched
        if not role_matched:
            failures.append(
                f"No reactant matches {role} patterns {smarts_list} "
                f"(required for {reaction_type})"
            )

    passed = len(failures) == 0
    return CheckResult(
        name="reaction_type_pattern",
        passed=passed,
        severity="error" if not passed else "info",
        score=1.0 if passed else (0.5 if len(failures) == 1 else 0.0),
        message=(
            f"All required FG patterns present for {reaction_type}."
            if passed
            else f"FG pattern mismatch for {reaction_type}: {'; '.join(failures)}"
        ),
        details={
            "reaction_type": reaction_type,
            "roles_matched": matched,
            "failures": failures,
            "patterns_checked": {r: p for r, p in patterns.items()},
        },
    )


def _check_fragment_reattachment(reactants: List[str], products: List[str]) -> CheckResult:
    """
    Check 8: MCS structural sanity.
    For each reactant (with ≥ 3 heavy atoms), its carbon skeleton should appear
    in the product at ≥ 60% coverage via MCS. Catches cases where the LLM
    suggests precursors that are structurally unrelated to the product.
    Skipped for catalysts / small reagents (< 3 heavy atoms).

    Threshold is 60% (not higher) to tolerate leaving-group atoms (B(OH)2,
    halides, TMS, etc.) that are present in the precursor but absent from the
    product — these correctly disappear during cross-coupling reactions.
    """
    MIN_ATOMS = 3
    MIN_COVERAGE = 0.60

    if not products:
        return CheckResult(
            name="fragment_reattachment",
            passed=False, severity="error", score=0.0,
            message="No product SMILES found.",
            details={},
        )

    product = products[0]  # evaluate against first product
    results: List[Dict] = []
    low_coverage: List[str] = []

    for smi in reactants:
        mol = _mol(smi)
        if mol is None or mol.GetNumHeavyAtoms() < MIN_ATOMS:
            continue
        coverage = _mcs_coverage(smi, product)
        results.append({"smiles": smi, "mcs_coverage": round(coverage, 3)})
        if coverage < MIN_COVERAGE:
            low_coverage.append(f"{smi} ({coverage:.0%} of atoms found in product)")

    passed = len(low_coverage) == 0
    avg_cov = (sum(r["mcs_coverage"] for r in results) / len(results)) if results else 1.0

    return CheckResult(
        name="fragment_reattachment",
        passed=passed,
        severity="warning" if not passed else "info",
        score=avg_cov,
        message=(
            f"All precursor skeletons found in product (avg MCS coverage {avg_cov:.0%})."
            if passed
            else (
                f"Low structural overlap between precursors and product: {low_coverage}. "
                f"Check regiochemistry or disconnection site."
            )
        ),
        details={"per_reactant": results, "low_coverage": low_coverage,
                 "threshold": MIN_COVERAGE, "product": product},
    )


# ---------------------------------------------------------------------------
# LLM eval prompt builder
# ---------------------------------------------------------------------------

def _build_llm_eval_prompt(report: "EvalReport") -> str:
    """Auto-generate a structured brief for LLM expert follow-up evaluation."""
    check_lines = []
    for c in report.checks:
        status = "✓" if c.passed else ("✗" if c.severity == "error" else "⚠")
        check_lines.append(f"  {status} [{c.severity.upper():<7}] {c.name}: {c.message}")

    return (
        f"REACTION TO EVALUATE: {report.reaction_smiles}\n"
        f"CLAIMED TYPE: {report.reaction_type or 'unspecified'}\n"
        f"\nRDKit CHECKS ({report.verdict}, score={report.overall_score:.2f}):\n"
        + "\n".join(check_lines)
        + "\n\nNow provide your EXPERT CHEMICAL ASSESSMENT:\n"
        "1. REGIOCHEMISTRY — Is the disconnection at the correct/most reactive position?\n"
        "2. MECHANISM — Is this transformation mechanistically sound for the claimed conditions?\n"
        "3. STEREOCHEMISTRY — Any stereocontrol issues? Is the proposed outcome achievable?\n"
        "4. PRACTICALITY — Commercial availability of precursors, scalability, known precedent?\n"
        "5. ALTERNATIVES — Are there better disconnections for this target?\n"
        "Flag any RDKit warnings above that need chemical context to interpret correctly."
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def evaluate_reaction(
    reaction_smiles: str,
    reaction_type: Optional[str] = None,
) -> "EvalReport":
    """
    Run all 8 evaluation checks on a reaction SMILES string.

    Args:
        reaction_smiles : Reaction in 'reactant1.reactant2>>product' format.
                          Multiple reactants separated by '.'.
        reaction_type   : Optional taxonomy ID (e.g. 'suzuki_miyaura') to
                          enable type-specific pattern and ring-change checks.

    Returns:
        EvalReport dataclass with per-check results, overall verdict,
        and a pre-built LLM evaluation prompt.
    """
    reaction_smiles = reaction_smiles.strip()

    # Parse
    try:
        reactants, products = _parse_reaction(reaction_smiles)
    except ValueError as e:
        dummy = CheckResult(
            name="parse", passed=False, severity="error", score=0.0,
            message=str(e), details={},
        )
        return EvalReport(
            reaction_smiles=reaction_smiles,
            reaction_type=reaction_type,
            checks=[dummy],
            overall_passed=False,
            overall_score=0.0,
            critical_failures=[str(e)],
            warnings=[],
            verdict="FAIL",
            llm_eval_prompt=f"REACTION: {reaction_smiles}\nPARSE ERROR: {e}",
        )

    # Run all 8 checks
    checks: List[CheckResult] = [
        _check_validate_components(reactants, products),
        _check_atom_balance(reactants, products),
        _check_charge_balance(reactants, products),
        _check_hdi_consistency(reactants, products),
        _check_mw_direction(reactants, products),
        _check_ring_change(reactants, products, reaction_type),
        _check_reaction_type_pattern(reactants, reaction_type),
        _check_fragment_reattachment(reactants, products),
    ]

    # Aggregate
    critical_failures = [c.message for c in checks if not c.passed and c.severity == "error"]
    warnings          = [c.message for c in checks if not c.passed and c.severity == "warning"]
    overall_passed    = len(critical_failures) == 0
    overall_score     = sum(c.score for c in checks) / len(checks) if checks else 0.0

    if not overall_passed:
        verdict = "FAIL"
    elif warnings:
        verdict = "PASS_WITH_WARNINGS"
    else:
        verdict = "PASS"

    report = EvalReport(
        reaction_smiles=reaction_smiles,
        reaction_type=reaction_type,
        checks=checks,
        overall_passed=overall_passed,
        overall_score=round(overall_score, 3),
        critical_failures=critical_failures,
        warnings=warnings,
        verdict=verdict,
        llm_eval_prompt="",  # filled below
    )
    report.llm_eval_prompt = _build_llm_eval_prompt(report)
    return report

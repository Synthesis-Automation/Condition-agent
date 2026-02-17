"""
Validation layer (Tier 4) for reaction analysis.

Provides RDKit-based structural validation, atom balance checks,
and cross-tier consensus testing.
"""

from typing import Dict, Any, List, Optional
import logging

# Reuse chemtools utilities
from chemtools.util.rdkit_helpers import parse_smiles, rdkit_available
from chemtools.featurizers.formatters.reaction_precheck import _count_elements
from chemtools.featurizers.analysis.smiles import normalize_reaction

logger = logging.getLogger(__name__)


def validate_with_rdkit(
    rxn_smiles_clean: str,
    tier2_result: Optional[Dict[str, Any]] = None,
    tier3_result: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    RDKit-based structural validation of reaction analysis.

    Checks:
    1. Structural validity (parseable molecules)
    2. Atom balance (reactants vs products)
    3. Functional group consistency
    4. Plausibility checks based on reported reaction type

    Args:
        rxn_smiles_clean: Clean reaction SMILES (reactants>>products)
        tier2_result: Tier 2 (DeepSeek) analysis results
        tier3_result: Tier 3 (gpt-4o-mini) interpretation results

    Returns:
        {
            "valid": True/False,
            "issues": [list of critical problems],
            "warnings": [list of concerns],
            "atom_balance": {...},
            "confidence_adjustment": -0.1,
            "checks_performed": [...]
        }
    """
    issues = []
    warnings = []
    checks_performed = []

    if not rdkit_available():
        return {
            "valid": True,  # Don't fail if RDKit unavailable
            "issues": [],
            "warnings": ["RDKit not available - structural validation skipped"],
            "checks_performed": ["rdkit_availability"]
        }

    # Parse reaction SMILES
    parts = rxn_smiles_clean.split(">>")
    if len(parts) != 2:
        issues.append("Invalid reaction SMILES format (expected reactants>>products)")
        return {
            "valid": False,
            "issues": issues,
            "warnings": warnings,
            "checks_performed": checks_performed
        }

    reactants_smiles = parts[0]
    products_smiles = parts[1]

    # 1. Validate all structures are parseable
    reactant_mols = []
    for smiles in reactants_smiles.split("."):
        if not smiles.strip():
            continue
        mol = parse_smiles(smiles)
        if mol is None:
            issues.append(f"Invalid reactant structure: {smiles[:50]}")
        else:
            reactant_mols.append(mol)

    product_mols = []
    for smiles in products_smiles.split("."):
        if not smiles.strip():
            continue
        mol = parse_smiles(smiles)
        if mol is None:
            issues.append(f"Invalid product structure: {smiles[:50]}")
        else:
            product_mols.append(mol)

    if issues:
        return {
            "valid": False,
            "issues": issues,
            "warnings": warnings,
            "checks_performed": checks_performed + ["structure_parsing"]
        }

    checks_performed.append("structure_parsing")

    # 2. Count atoms (reuse chemtools function)
    reactant_counts = _count_elements(reactants_smiles.split("."))
    product_counts = _count_elements(products_smiles.split("."))

    reactant_atoms = sum(reactant_counts.values())
    product_atoms = sum(product_counts.values())
    atoms_lost = reactant_atoms - product_atoms

    atom_balance = {
        "reactants": reactant_atoms,
        "products": product_atoms,
        "lost": atoms_lost,
        "reactant_formula": dict(reactant_counts),
        "product_formula": dict(product_counts)
    }
    checks_performed.append("atom_counting")

    # 3. Check atom loss matches reported changes
    if tier2_result:
        all_changes = tier2_result.get("all_changes", [])
        has_deprotection = any(
            "deprotection" in str(c).lower() or
            "removal" in str(c).lower() or
            "cleav" in str(c).lower()
            for c in all_changes
        )

        if atoms_lost > 5 and not has_deprotection:
            warnings.append(
                f"{atoms_lost} atoms lost but no deprotection/cleavage reported in Tier 2"
            )

        checks_performed.append("atom_loss_consistency")

    # 4. Validate reaction type plausibility
    if tier2_result and tier3_result:
        reaction_types = [rt.lower() for rt in tier2_result.get("reaction_types", [])]
        overall_class = tier3_result.get("overall_class", "").lower()

        # Check for Suzuki/cross-coupling
        is_coupling = any(
            "suzuki" in rt or "coupling" in rt or "heck" in rt
            for rt in reaction_types
        ) or "coupling" in overall_class

        if is_coupling:
            # Should have aryl/vinyl halide in reactants
            has_halide = any(
                check_for_aryl_halide(mol) for mol in reactant_mols
            )
            if not has_halide:
                warnings.append(
                    "Cross-coupling reported but no clear aryl/vinyl halide found in reactants"
                )
            checks_performed.append("coupling_plausibility")

        # Check for oxidation
        is_oxidation = any("oxidation" in rt for rt in reaction_types)
        if is_oxidation:
            # Check if oxidation state increased (rough heuristic: O count increased)
            o_increase = product_counts.get('O', 0) - reactant_counts.get('O', 0)
            if o_increase <= 0:
                warnings.append(
                    f"Oxidation reported but O count changed by {o_increase}"
                )
            checks_performed.append("oxidation_plausibility")

    # 5. Calculate confidence adjustment
    confidence_adjustment = 0.0
    if issues:
        confidence_adjustment = -0.3  # Major structural issues
    elif len(warnings) > 2:
        confidence_adjustment = -0.1  # Multiple concerns
    elif len(warnings) == 1:
        confidence_adjustment = -0.05  # Minor concern

    return {
        "valid": len(issues) == 0,
        "issues": issues,
        "warnings": warnings,
        "atom_balance": atom_balance,
        "confidence_adjustment": confidence_adjustment,
        "checks_performed": checks_performed
    }


def check_for_aryl_halide(mol) -> bool:
    """
    Check if molecule contains aryl or vinyl halide.

    Returns True if found, False otherwise.
    """
    if mol is None:
        return False

    try:
        # Look for aromatic carbon or sp2 carbon bonded to Br, I, Cl
        for atom in mol.GetAtoms():
            if atom.GetSymbol() not in ['Br', 'I', 'Cl']:
                continue

            # Check neighbors
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    # Check if aromatic or sp2
                    if neighbor.GetIsAromatic():
                        return True
                    # Check hybridization (2 = SP2)
                    try:
                        from rdkit.Chem import rdchem
                        if neighbor.GetHybridization() == rdchem.HybridizationType.SP2:
                            return True
                    except:
                        pass

        return False
    except Exception as e:
        logger.debug(f"Error checking for aryl halide: {e}")
        return False


def check_consensus(
    tier1_result: Optional[Dict[str, Any]],
    tier2_result: Optional[Dict[str, Any]],
    tier3_result: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Check consensus across all three tiers.

    Returns quality score and issues.
    """
    issues = []
    warnings = []

    # 1. Check Tier 1 → Tier 2 agreement
    if tier1_result and tier2_result:
        tier1_patterns = tier1_result.get('interpretation', {})
        tier1_types = [p.lower() for p in tier1_patterns.get('likely_types', [])]
        tier2_types = [rt.lower() for rt in tier2_result.get('reaction_types', [])]

        # If Tier 1 detected something specific, Tier 2 should confirm
        for t1_type in tier1_types:
            if t1_type and not any(t1_type in t2 for t2 in tier2_types):
                warnings.append(
                    f"Tier 1 detected '{t1_type}' but Tier 2 didn't confirm"
                )

    # 2. Check Tier 2 → Tier 3 agreement
    if tier2_result:
        tier2_types = [rt.lower() for rt in tier2_result.get('reaction_types', [])]
        tier3_class = tier3_result.get('overall_class', '').lower()
        tier3_tags = [str(t).lower() for t in tier3_result.get('tags', [])]

        # Major reaction types that should agree
        suzuki_in_t2 = any('suzuki' in rt or 'coupling' in rt for rt in tier2_types)
        coupling_in_t3 = 'coupling' in tier3_class or any('suzuki' in tag for tag in tier3_tags)

        if suzuki_in_t2 and not coupling_in_t3:
            issues.append(
                f"Major disagreement: Tier 2 says coupling, Tier 3 says {tier3_class}"
            )

        # Check for substitution disagreement
        substitution_in_t2 = any('substitution' in rt for rt in tier2_types)
        substitution_in_t3 = 'substitution' in tier3_class

        if suzuki_in_t2 and substitution_in_t3:
            issues.append(
                "Major disagreement: Tier 2 says Suzuki, Tier 3 says substitution"
            )

    # 3. Check confidence scores
    tier2_conf = tier2_result.get('confidence', 1.0) if tier2_result else 1.0
    tier3_conf = tier3_result.get('confidence', 0.0)

    if tier2_conf < 0.7:
        warnings.append(f"Tier 2 low confidence ({tier2_conf:.2f})")
    if tier3_conf < 0.6:
        warnings.append(f"Tier 3 low confidence ({tier3_conf:.2f})")

    # 4. Calculate overall quality score
    quality_score = 1.0
    quality_score -= len(issues) * 0.2  # -0.2 per major issue
    quality_score -= len(warnings) * 0.05  # -0.05 per warning
    quality_score = max(0.0, min(1.0, quality_score))

    # 5. Determine recommendation
    if quality_score > 0.8:
        recommendation = "accept"
    elif quality_score > 0.6:
        recommendation = "review"
    else:
        recommendation = "re_analyze"

    return {
        "quality_score": quality_score,
        "issues": issues,
        "warnings": warnings,
        "recommendation": recommendation,
        "confidence_scores": {
            "tier2": tier2_conf,
            "tier3": tier3_conf
        }
    }


def quality_gate(
    result: Dict[str, Any],
    rdkit_validation: Dict[str, Any],
    consensus_check: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Quality gate: Decide if analysis is acceptable or needs retry.

    Returns:
        {
            "status": "pass"|"warning"|"fail",
            "action": "accept"|"accept_with_warnings"|"re_analyze",
            "retry_config": {...}  # If re-analysis needed
        }
    """
    # Collect all issues
    all_issues = []
    all_issues.extend(rdkit_validation.get('issues', []))
    all_issues.extend(consensus_check.get('issues', []))

    all_warnings = []
    all_warnings.extend(rdkit_validation.get('warnings', []))
    all_warnings.extend(consensus_check.get('warnings', []))

    # Get quality score
    quality_score = consensus_check.get('quality_score', 0.0)

    # Apply RDKit confidence adjustment
    quality_score += rdkit_validation.get('confidence_adjustment', 0.0)
    quality_score = max(0.0, min(1.0, quality_score))

    # Decision logic
    if len(all_issues) > 0:
        # Critical issues found → Re-analyze with stronger model
        return {
            "status": "fail",
            "action": "re_analyze",
            "issues": all_issues,
            "warnings": all_warnings,
            "quality_score": quality_score,
            "retry_config": {
                "use_stronger_model": True,
                "model": "deepseek-v3.2",  # Use best model for Tier 3
                "include_context": all_issues[:2],  # Pass top issues to model
                "mode": "expert"
            }
        }
    elif quality_score < 0.7 or len(all_warnings) > 3:
        # Multiple warnings → Accept but flag for review
        return {
            "status": "warning",
            "action": "accept_with_warnings",
            "warnings": all_warnings,
            "quality_score": quality_score,
            "suggestion": "Manual review recommended due to low confidence or warnings"
        }
    else:
        # Looks good → Accept
        return {
            "status": "pass",
            "action": "accept",
            "warnings": all_warnings,
            "quality_score": quality_score
        }


# Export public API
__all__ = [
    "validate_with_rdkit",
    "check_consensus",
    "quality_gate",
]

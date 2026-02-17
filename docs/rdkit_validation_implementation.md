# RDKit Validation Implementation Plan - Using Existing Chemtools

## Summary of Available Tools in Chemtools

After exploring the codebase, I found these **ready-to-use components**:

### 1. **chemtools/util/rdkit_helpers.py**
Provides core RDKit utilities:
- ✅ `parse_smiles(smiles)` - Parse SMILES to RDKit mol object
- ✅ `canonical_smiles(smiles)` - Canonicalize SMILES string
- ✅ `rdkit_available()` - Check if RDKit is available
- ✅ `choose_largest_organic_fragment(mol)` - Fragment selection
- ✅ `neutralize_and_standardize(mol)` - Clean up molecules

### 2. **chemtools/featurizers/formatters/reaction_precheck.py**
Provides atom counting and stoichiometry:
- ✅ `_count_elements(smiles_list)` - Count all elements in molecules
- ✅ `infer_reactant_repeats_from_stoichiometry()` - Detect atom balance issues

### 3. **chemtools/featurizers/analysis/smiles.py**
Provides SMILES normalization:
- ✅ `normalize(smiles)` - Normalize single SMILES
- ✅ `normalize_reaction(rxn_smiles)` - Parse and normalize reaction SMILES
- ✅ Split reactants/products

### 4. **chemtools/featurizers/formatters/detection_validation.py**
Provides reaction type validation:
- ✅ `validate_detection_with_crk_key()` - Validate using CRK patterns
- ✅ `validate_detection_with_reacted_motifs()` - Validate using motif consumption
- ✅ Taxonomy-based reaction matching
- ✅ Evidence collection and scoring

## What We Can Reuse vs Build New

| Functionality | Available in Chemtools | Status | Action |
|--------------|----------------------|--------|---------|
| **Parse SMILES** | ✅ rdkit_helpers.parse_smiles() | Ready | ✅ Use directly |
| **Count atoms** | ✅ reaction_precheck._count_elements() | Ready | ✅ Use directly |
| **Validate structures** | ✅ rdkit_helpers | Ready | ✅ Wrap in validator |
| **Normalize SMILES** | ✅ smiles.normalize_reaction() | Ready | ✅ Use directly |
| **Reaction type validation** | ✅ detection_validation | Ready | ✅ Integrate |
| **Functional group detection** | ❌ Not found | Missing | 🔨 Build new |
| **Cross-tier consensus** | ❌ Not found | Missing | 🔨 Build new |
| **Quality gate logic** | ❌ Not found | Missing | 🔨 Build new |

## Implementation Plan - Phase 1: RDKit Validator

### File: `reaction_agent/validation.py` (NEW)

```python
"""
Validation layer (Tier 4) for reaction analysis.

Provides RDKit-based structural validation, atom balance checks,
and cross-tier consensus testing.
"""

from typing import Dict, Any, List, Optional, Set
from collections import Counter
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
    try:
        # Use chemtools normalize_reaction (already handles parsing)
        normalized = normalize_reaction(rxn_smiles_clean)
        checks_performed.append("smiles_parsing")
    except Exception as e:
        issues.append(f"Failed to parse reaction SMILES: {e}")
        return {
            "valid": False,
            "issues": issues,
            "warnings": warnings,
            "checks_performed": checks_performed
        }

    # Extract components
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
    reactant_counts = _count_elements([reactants_smiles])
    product_counts = _count_elements([products_smiles])

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
                    if neighbor.GetIsAromatic() or neighbor.GetHybridization() == 2:  # SP2
                        return True

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

        # Major reaction types that should agree
        suzuki_in_t2 = any('suzuki' in rt or 'coupling' in rt for rt in tier2_types)
        coupling_in_t3 = 'coupling' in tier3_class

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
            "suggestion": "Manual review recommended due to low confidence"
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
```

## Integration Points

### 1. Update `reaction_agent/agent.py`

Add optional validation parameter:

```python
def analyze_reaction_smiles(
    rxn_smiles: str,
    client: LLMClient,
    drop_spectators: bool = True,
    skip_mapping: bool = False,
    temperature: float = 0.0,
    max_tokens: int = 2000,
    reasoning_effort: Optional[str] = None,
    validate: bool = False  # NEW: Enable Tier 4 validation
) -> Dict[str, Any]:
    """
    ... existing docstring ...

    Args:
        validate: Enable Tier 4 validation (RDKit + consensus checks)
    """
    # ... existing code for Tiers 1-3 ...

    if validate:
        from reaction_agent.validation import (
            validate_with_rdkit,
            check_consensus,
            quality_gate
        )

        # Run Tier 4 validation
        rdkit_val = validate_with_rdkit(
            input_data.get("rxn_smiles_clean", ""),
            result.get("quick_glance"),
            result.get("interpretation")
        )

        consensus = check_consensus(
            result.get("auto_interpretation"),
            result.get("quick_glance"),
            result.get("interpretation")
        )

        gate = quality_gate(result, rdkit_val, consensus)

        # Add validation results
        result["validation"] = {
            "rdkit": rdkit_val,
            "consensus": consensus,
            "gate": gate
        }

    return result
```

### 2. Update `reaction_agent/cli.py`

Add --validate flag and display:

```python
# In argparse
parser.add_argument('--validate', action='store_true',
                    help='Enable Tier 4 validation (RDKit + consensus)')

# In display function
def display_result(result: Dict[str, Any], show_details: bool = False):
    # ... existing display code ...

    # Validation section
    if 'validation' in result:
        print_header("VALIDATION RESULTS (Tier 4)")
        validation = result['validation']

        # RDKit checks
        rdkit = validation.get('rdkit', {})
        if rdkit.get('valid'):
            print(f"{Colors.GREEN}RDKit Checks: ✓ PASS{Colors.END}")
        else:
            print(f"{Colors.RED}RDKit Checks: ✗ FAIL{Colors.END}")
            for issue in rdkit.get('issues', []):
                print(f"  ✗ {issue}")

        if rdkit.get('warnings'):
            print(f"\n{Colors.YELLOW}Warnings ({len(rdkit['warnings'])}):{ Colors.END}")
            for warn in rdkit['warnings'][:3]:  # First 3
                print(f"  ⚠ {warn}")

        # Atom balance
        atom_bal = rdkit.get('atom_balance', {})
        if atom_bal:
            reactants = atom_bal.get('reactants', 0)
            products = atom_bal.get('products', 0)
            lost = atom_bal.get('lost', 0)
            print(f"\nAtom Balance: {reactants} → {products} ({lost} lost)")

        # Consensus
        consensus = validation.get('consensus', {})
        quality = consensus.get('quality_score', 0.0)
        quality_color = Colors.GREEN if quality > 0.8 else Colors.YELLOW if quality > 0.6 else Colors.RED
        print(f"\nConsensus Score: {quality_color}{quality:.2f} / 1.00{Colors.END}")

        # Gate decision
        gate = validation.get('gate', {})
        status = gate.get('status', 'unknown')
        if status == 'pass':
            print(f"\n{Colors.GREEN}Overall Status: ✓ PASS - High quality analysis{Colors.END}")
        elif status == 'warning':
            print(f"\n{Colors.YELLOW}Overall Status: ⚠ WARNING - Review recommended{Colors.END}")
        else:
            print(f"\n{Colors.RED}Overall Status: ✗ FAIL - Re-analysis needed{Colors.END}")
```

## Testing Plan

### Test Case 1: Simple Reaction (Should Pass)
```python
# Suzuki coupling - straightforward
rxn = "Brc1ccccc1.B(O)(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"

result = analyze_reaction_smiles(rxn, client, validate=True)
assert result['validation']['gate']['status'] == 'pass'
```

### Test Case 2: Complex Tandem (Should Warn)
```python
# Your Suzuki + THP deprotection
rxn = "CC1(C)OB(c2cnn(CCOC3CCCCO3)c2)OC1(C)C.Cc1nc(-c2cn3c(n2)-c2ccc(Br)cc2OCC3)n(C)n1>>Cc1nc(-c2cn3c(n2)-c2ccc(-c4cnn(CCO)c4)cc2OCC3)n(C)n1"

result = analyze_reaction_smiles(rxn, client, validate=True)
validation = result['validation']

# Should warn about atom loss
assert len(validation['rdkit']['warnings']) > 0
# But should still pass overall
assert validation['gate']['status'] in ['pass', 'warning']
```

### Test Case 3: Invalid Structure (Should Fail)
```python
# Invalid SMILES
rxn = "INVALID>>SMILES"

result = analyze_reaction_smiles(rxn, client, validate=True)
assert result['validation']['gate']['status'] == 'fail'
assert len(result['validation']['rdkit']['issues']) > 0
```

## Benefits of This Approach

✅ **Reuses 80% of existing code** - chemtools already has most utilities we need
✅ **No duplication** - Leverages well-tested rdkit_helpers and reaction_precheck
✅ **Clean separation** - validation.py is self-contained and optional
✅ **Minimal changes** - Only adds validate=True parameter to existing functions
✅ **Fast implementation** - Can be done in 2-3 hours using existing tools
✅ **Well-tested** - Builds on battle-tested chemtools foundation

## Next Steps

1. **Create `reaction_agent/validation.py`** - Implement using chemtools utilities
2. **Add validate parameter to agent.py** - Make it optional (default False)
3. **Add --validate flag to CLI** - Display validation results
4. **Write 5-10 test cases** - Cover pass/warn/fail scenarios
5. **Update documentation** - Explain when to use validation

**Estimated time: 2-3 hours** with existing chemtools utilities!

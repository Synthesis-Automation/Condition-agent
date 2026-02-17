"""
Lightweight reaction interpretation module.

Provides fast, deterministic analysis of reaction patterns to help interpret
mapping results and guide users when methods disagree.
"""

from typing import Dict, List, Tuple, Any
import re


def interpret_reaction_pattern(
    rxn_smiles: str,
    hybrid_result: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Quick interpretation of reaction pattern to help understand mapping results.

    This is a lightweight, deterministic analysis that runs in milliseconds.
    Useful when mapping methods disagree or confidence is borderline.

    Args:
        rxn_smiles: Reaction SMILES
        hybrid_result: Result from analyze_bond_changes_hybrid()

    Returns:
        Dictionary with interpretation insights
    """

    interpretation = {
        'reaction_complexity': 'unknown',
        'likely_reaction_types': [],
        'tandem_reaction_suspected': False,
        'explanation': '',
        'recommendation': '',
        'patterns_detected': []
    }

    # Parse reaction
    if '>>' not in rxn_smiles:
        return interpretation

    reactants, products = rxn_smiles.split('>>')
    reactant_parts = [p.strip() for p in reactants.split('.') if p.strip()]
    product_parts = [p.strip() for p in products.split('.') if p.strip()]

    # Count molecules and estimate atom counts
    n_reactants = len(reactant_parts)
    n_products = len(product_parts)

    # Rough atom count (count uppercase letters as atoms)
    def rough_atom_count(smiles: str) -> int:
        # Remove brackets, numbers, charges
        clean = re.sub(r'\[|\]|\d+|[+-]', '', smiles)
        return len(re.findall(r'[A-Z]', clean))

    total_r_atoms = sum(rough_atom_count(p) for p in reactant_parts)
    total_p_atoms = sum(rough_atom_count(p) for p in product_parts)

    atom_loss = total_r_atoms - total_p_atoms

    # Extract mapping results
    rxn_ok = hybrid_result.get('rxnmapper_result', {}).get('success', False)
    local_ok = hybrid_result.get('local_env_result', {}).get('success', False)
    mcs_ok = hybrid_result.get('mcs_result', {}).get('success', False)

    rxn_broken = len(hybrid_result.get('rxnmapper_result', {}).get('broken_bonds', []))
    local_broken = len(hybrid_result.get('local_env_result', {}).get('broken_bonds', []))

    all_disagree = (
        hybrid_result.get('agreement', {}).get('rxnmapper_vs_local_env') == False and
        hybrid_result.get('agreement', {}).get('rxnmapper_vs_mcs') == False and
        hybrid_result.get('agreement', {}).get('local_env_vs_mcs') == False
    )

    # Pattern detection

    # 1. Check for leaving groups (suggests coupling or substitution)
    leaving_groups = []
    if 'Br' in reactants or '[Br]' in reactants:
        leaving_groups.append('Br')
        interpretation['patterns_detected'].append('aryl/alkyl bromide')
    if 'I' in reactants or '[I]' in reactants:
        leaving_groups.append('I')
        interpretation['patterns_detected'].append('aryl/alkyl iodide')
    if 'Cl' in reactants or '[Cl]' in reactants:
        leaving_groups.append('Cl')
        interpretation['patterns_detected'].append('aryl/alkyl chloride')

    # 2. Check for boron reagent (Suzuki coupling)
    if 'B(' in reactants or 'B1' in reactants:
        interpretation['patterns_detected'].append('boronic ester/acid')
        interpretation['likely_reaction_types'].append('Suzuki coupling')

    # 3. Check for acetal/ketal - multiple patterns
    has_acetal = False
    # Pattern 1: C(OC)(OC) or C(O)(O)
    if re.search(r'C\([^)]*O[^)]*\)\([^)]*O', reactants):
        has_acetal = True
    # Pattern 2: Simple check - carbon with two separate methoxy groups
    if 'COC(OC)' in reactants or 'C(OC)(OC)' in reactants or 'OC(OC)' in reactants:
        has_acetal = True

    if has_acetal:
        interpretation['patterns_detected'].append('acetal/ketal')
        # Check if it becomes carbonyl in product
        if 'C=O' in products or 'C(=O)' in products:
            # Make sure methoxy groups are lost (indicates hydrolysis)
            if reactants.count('OC') > products.count('OC'):
                interpretation['likely_reaction_types'].append('acetal hydrolysis')

    # 4. Check for protected groups
    if 'Si(' in reactants or 'Si[' in reactants:
        interpretation['patterns_detected'].append('silyl protection')
    if 'C(=O)O' in reactants and 'OH' in products:
        interpretation['likely_reaction_types'].append('ester hydrolysis')

    # 5. Atom count changes
    if atom_loss > 10:
        interpretation['patterns_detected'].append(f'significant atom loss ({atom_loss} atoms)')
        interpretation['tandem_reaction_suspected'] = True
    elif atom_loss > 5:
        interpretation['patterns_detected'].append(f'moderate atom loss ({atom_loss} atoms)')

    # Determine complexity
    if len(interpretation['likely_reaction_types']) > 1:
        interpretation['reaction_complexity'] = 'tandem/multi-step'
        interpretation['tandem_reaction_suspected'] = True
    elif rxn_broken > 6 or local_broken > 8:
        interpretation['reaction_complexity'] = 'complex'
    elif rxn_broken <= 3 and local_broken <= 4:
        interpretation['reaction_complexity'] = 'simple'
    else:
        interpretation['reaction_complexity'] = 'moderate'

    # Generate explanation
    explanations = []

    if interpretation['tandem_reaction_suspected']:
        types = ' + '.join(interpretation['likely_reaction_types'])
        explanations.append(f"Tandem reaction detected: {types}")
        explanations.append("Multiple transformations occurring in sequence make atom mapping challenging")

    if all_disagree:
        explanations.append("All mapping methods disagree, indicating:")
        if interpretation['reaction_complexity'] == 'tandem/multi-step':
            explanations.append("  • Multiple reaction steps confuse automated mappers")
            explanations.append("  • Local environment mapper may be most reliable for tracking spectators")
        else:
            explanations.append("  • Complex rearrangement or unusual mechanism")
            explanations.append("  • Leaving groups may be incorrectly mapped")

    if atom_loss > 5:
        explanations.append(f"Significant byproducts generated ({atom_loss} atoms lost)")
        if leaving_groups:
            explanations.append(f"  • Likely leaving groups: {', '.join(leaving_groups)}")

    if not explanations:
        if interpretation['reaction_complexity'] == 'simple':
            explanations.append("Simple transformation - methods should agree")
        else:
            explanations.append("Standard single-step reaction")

    interpretation['explanation'] = '\n'.join(explanations)

    # Generate recommendation
    recommendations = []

    if interpretation['tandem_reaction_suspected']:
        recommendations.append("✓ Trust local environment mapping for spectator identification")
        recommendations.append("✓ Manually verify each transformation step")
        recommendations.append("✓ Consider LLM-assisted analysis for mechanistic insights")
    elif all_disagree:
        recommendations.append("⚠ Manual review required - all methods disagree")
        if local_ok and rxn_ok:
            local_conf = hybrid_result.get('local_env_result', {}).get('confidence', 0)
            rxn_conf = hybrid_result.get('rxnmapper_result', {}).get('mapping_confidence', 0)
            if local_conf > rxn_conf + 0.1:
                recommendations.append("✓ Local environment shows higher confidence - may be more reliable")
            else:
                recommendations.append("✓ RXNMapper typically more accurate for single-step reactions")
    else:
        if hybrid_result.get('combined_confidence', 0) > 0.7:
            recommendations.append("✓ Good confidence - mapping appears reliable")
        else:
            recommendations.append("⚠ Moderate confidence - consider validation")

    interpretation['recommendation'] = '\n'.join(recommendations)

    return interpretation


def format_interpretation_report(
    rxn_smiles: str,
    hybrid_result: Dict[str, Any],
    interpretation: Dict[str, Any]
) -> str:
    """
    Format a human-readable interpretation report.

    Args:
        rxn_smiles: Reaction SMILES
        hybrid_result: Hybrid mapping result
        interpretation: Interpretation from interpret_reaction_pattern()

    Returns:
        Formatted report string
    """

    lines = []
    lines.append("=" * 80)
    lines.append("AUTOMATIC REACTION INTERPRETATION")
    lines.append("=" * 80)
    lines.append("")

    # Reaction complexity
    complexity_emoji = {
        'simple': '🟢',
        'moderate': '🟡',
        'complex': '🟠',
        'tandem/multi-step': '🔴'
    }
    emoji = complexity_emoji.get(interpretation['reaction_complexity'], '⚪')
    lines.append(f"{emoji} Complexity: {interpretation['reaction_complexity'].upper()}")
    lines.append("")

    # Patterns detected
    if interpretation['patterns_detected']:
        lines.append("Patterns Detected:")
        for pattern in interpretation['patterns_detected']:
            lines.append(f"  • {pattern}")
        lines.append("")

    # Likely reaction types
    if interpretation['likely_reaction_types']:
        lines.append("Likely Reaction Type(s):")
        for rxn_type in interpretation['likely_reaction_types']:
            lines.append(f"  • {rxn_type}")
        lines.append("")

    # Tandem reaction flag
    if interpretation['tandem_reaction_suspected']:
        lines.append("⚠️  TANDEM/MULTI-STEP REACTION SUSPECTED")
        lines.append("")

    # Explanation
    if interpretation['explanation']:
        lines.append("Explanation:")
        for line in interpretation['explanation'].split('\n'):
            lines.append(f"  {line}")
        lines.append("")

    # Mapping results summary
    lines.append("Mapping Methods Performance:")
    if hybrid_result.get('rxnmapper_result', {}).get('success'):
        rxn_conf = hybrid_result['rxnmapper_result'].get('mapping_confidence', 0)
        rxn_bonds = len(hybrid_result['rxnmapper_result'].get('broken_bonds', []))
        lines.append(f"  • RXNMapper: {rxn_conf:.3f} confidence, {rxn_bonds} broken bonds")
    if hybrid_result.get('local_env_result', {}).get('success'):
        local_conf = hybrid_result['local_env_result'].get('confidence', 0)
        local_bonds = len(hybrid_result['local_env_result'].get('broken_bonds', []))
        lines.append(f"  • Local Environment: {local_conf:.3f} confidence, {local_bonds} broken bonds")
    if hybrid_result.get('mcs_result', {}).get('success'):
        mcs_conf = hybrid_result['mcs_result'].get('confidence', 0)
        lines.append(f"  • MCS: {mcs_conf:.3f} confidence (approximate)")
    lines.append("")

    # Agreement status
    agreement = hybrid_result.get('agreement', {})
    agrees = sum(1 for v in agreement.values() if v == True)
    total = sum(1 for v in agreement.values() if v is not None)
    if total > 0:
        lines.append(f"Method Agreement: {agrees}/{total} pairs agree")
        lines.append("")

    # Recommendation
    if interpretation['recommendation']:
        lines.append("RECOMMENDATION:")
        for line in interpretation['recommendation'].split('\n'):
            lines.append(f"  {line}")
        lines.append("")

    lines.append("=" * 80)

    return '\n'.join(lines)

#!/usr/bin/env python3
"""
Simple CLI for cross-family reaction condition recommendations.

Usage:
    python app/cross_family_recommendation_cli.py "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    python app/cross_family_recommendation_cli.py --rxn "CCBr.CCN>>CCNCC"
    python app/cross_family_recommendation_cli.py --rxn "..." --k 50
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from chemtools import chem
from chemtools.util.drfp_storage import get_unified_drfp_path


def print_separator(char="=", length=80):
    """Print a separator line."""
    print(char * length)


def extract_chemical_by_role(chemicals: list, role: str) -> str:
    """Extract chemical by role from chemicals list."""
    for chem in chemicals:
        if chem.get("role") == role:
            name = chem.get("name") or chem.get("abbreviation") or chem.get("smiles")
            return name if name else "N/A"
    return "N/A"


def extract_chemicals_by_roles(chemicals: list, roles: list) -> dict:
    """Extract multiple chemicals by roles from chemicals list."""
    result = {}
    for role in roles:
        for chem in chemicals:
            if chem.get("role") == role:
                name = chem.get("name") or chem.get("abbreviation") or chem.get("smiles")
                if name and name != "N/A":
                    result[role] = name
                break
    return result


def find_matching_precedents(rec: dict, all_precedents: list, max_count: int = 5) -> list:
    """
    Find matching precedents for a recommendation from the full precedents list.
    
    Args:
        rec: Recommendation dict with 'summary' and 'combo'
        all_precedents: List of all precedent dicts from result['precedents_used']['top_precedents']
        max_count: Maximum number of precedents to return
    
    Returns:
        List of matching precedent dicts
    """
    # First, check if summary already has precedents
    summary = rec.get('summary', {})
    if summary.get('precedents'):
        return summary['precedents'][:max_count]
    
    # Try to match based on combo (base + solvent)
    combo = rec.get('combo', {})
    base_uid = combo.get('base_uid')
    solvent_uid = combo.get('solvent_uid')
    
    if not base_uid and not solvent_uid:
        return []
    
    # Find precedents that match the combo
    matched = []
    for prec in all_precedents:
        # Check if base matches
        base_match = False
        if base_uid:
            for reagent in prec.get('reagents', []):
                if reagent.get('cas') == base_uid:
                    base_match = True
                    break
        else:
            base_match = True  # No base requirement
        
        # Check if solvent matches
        solvent_match = False
        if solvent_uid:
            for solvent in prec.get('solvents', []):
                if solvent.get('cas') == solvent_uid:
                    solvent_match = True
                    break
        else:
            solvent_match = True  # No solvent requirement
        
        if base_match and solvent_match:
            # Convert to simplified format for display
            matched.append({
                'reaction_id': prec.get('reaction_id'),
                'core': prec.get('core'),
                'yield_pct': prec.get('yield'),
                'reference': prec.get('reference', '')
            })
            
            if len(matched) >= max_count:
                break
    
    return matched


def print_recommendation(result: dict, reaction_smiles: str, max_precedents: int = 3):
    """Print recommendation results in a readable format."""
    print_separator()
    print("CROSS-FAMILY REACTION CONDITION RECOMMENDATION")
    print_separator()
    print(f"\nReaction: {reaction_smiles}\n")
    
    # Check if unified index was used
    unified_path = get_unified_drfp_path()
    if Path(unified_path).exists():
        print("Using unified DRFP index for accurate cross-family search\n")
    else:
        print("Unified DRFP index not found - using feature-based similarity")
        print(f"To enable DRFP: python scripts/build_unified_drfp_index.py\n")
    
    # Extract key information
    strategy = result.get("strategy", "unknown")
    confidence = result.get("confidence", 0.0)
    
    print(f"Strategy: {strategy.upper()}")
    print(f"Confidence: {confidence:.2f}")
    
    # Detection information - show both detected family and analysis results
    if "detection" in result:
        detection = result["detection"]
        detected_family = detection.get("detected_family") or detection.get("family", "Unknown")
        detection_confidence = detection.get("confidence", 0.0)
        analysis_used = detection.get("analysis_module_used", False)
        
        print(f"Detected Family: {detected_family} (confidence: {detection_confidence:.2f})")
        
        # Show reactant classification if available
        if analysis_used and detection.get("reactant_classification"):
            rc = detection["reactant_classification"]
            rxn_type = rc.get("reaction_type", "Unknown")
            num_reactants = rc.get("num_reactants", 0)
            print(f"  - Analysis Module: {rxn_type} ({num_reactants} reactants classified)")
    elif "detected_family" in result:
        # Fallback to top-level detected_family
        detected_family = result.get("detected_family", "Unknown")
        print(f"Detected Family: {detected_family}")
    
    # Recommendations
    recommendations = result.get("recommendations", [])
    if not recommendations:
        print("\nNo recommendations found")
        print_separator()
        return
    
    print(f"\nFound {len(recommendations)} recommendation(s):\n")
    
    # Get all precedents from result for matching
    all_precedents = result.get('precedents_used', {}).get('top_precedents', [])
    
    for idx, rec in enumerate(recommendations, 1):
        print(f"{'─' * 80}")
        print(f"Recommendation #{idx}")
        print(f"{'─' * 80}")
        
        # Extract from summary if available (preferred)
        summary = rec.get("summary", {})
        chemicals = rec.get("chemicals", [])
        conditions = rec.get("conditions", {})
        
        # Catalyst (from summary.core or extract from chemicals)
        catalyst = summary.get("core", "N/A")
        if catalyst == "N/A":
            metal = extract_chemical_by_role(chemicals, "metal_precursor")
            ligand = extract_chemical_by_role(chemicals, "ligand")
            if metal != "N/A" or ligand != "N/A":
                parts = [p for p in [metal, ligand] if p != "N/A"]
                catalyst = "/".join(parts)
        print(f"  Catalyst: {catalyst}")
        
        # Base
        base = summary.get("base", {})
        if isinstance(base, dict):
            base_name = base.get("name") or base.get("abbreviation", "N/A")
        else:
            base_name = extract_chemical_by_role(chemicals, "base")
        if base_name != "N/A":
            print(f"  Base: {base_name}")
        
        # Solvent
        solvent = summary.get("solvent", {})
        if isinstance(solvent, dict):
            solvent_name = solvent.get("name") or solvent.get("abbreviation", "N/A")
        else:
            solvent_name = extract_chemical_by_role(chemicals, "solvent")
        print(f"  Solvent: {solvent_name}")
        
        # Additional reagents
        other_reagents = extract_chemicals_by_roles(
            chemicals, 
            ["additive", "reagent", "coupling_agent", "activating_agent"]
        )
        if other_reagents:
            print(f"  Other Reagents:")
            for role, name in other_reagents.items():
                print(f"    - {role.replace('_', ' ').title()}: {name}")
        
        # Temperature
        temp = conditions.get("temperature")
        if temp:
            print(f"  Temperature: {temp}")
        
        # Time
        time = conditions.get("time")
        if time:
            print(f"  Time: {time}")
        
        # Confidence/support from summary
        if summary:
            rec_confidence = summary.get("confidence")
            if rec_confidence is not None:
                print(f"  Confidence: {rec_confidence:.2f}")
            
            support = summary.get("support", {})
            if support:
                count = support.get("count", 0)
                if count > 0:
                    print(f"  Precedent Support: {count} similar reaction(s)")
        
        # Precedent information - use helper to find matching precedents
        precedents = find_matching_precedents(rec, all_precedents, max_precedents)
        if precedents:
            num_to_show = min(len(precedents), max_precedents)
            support = summary.get('support', {})
            total_support = support.get('count', 0)
            # Use the actual support count if available, otherwise use precedent count
            if total_support > 0:
                print(f"  Top Precedents ({num_to_show} of {total_support}):")
            else:
                print(f"  Top Precedents ({num_to_show}):")
            for i, prec in enumerate(precedents[:max_precedents], 1):
                reaction_id = prec.get("reaction_id", "N/A")
                ref = prec.get("reference", "")
                core = prec.get("core", "")
                yield_pct = prec.get("yield_pct")
                
                # Extract just the title from reference
                ref_title = ref.split("|")[0].strip() if "|" in ref else ref[:80]
                
                print(f"    {i}. {reaction_id}")
                if core:
                    print(f"       Catalyst: {core}")
                if yield_pct is not None:
                    print(f"       Yield: {yield_pct}%")
                if ref_title:
                    print(f"       Ref: {ref_title}...")
        else:
            # If no precedents found
            print(f"  Top Precedents: None available")
        
        print()
    
    print_separator()


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Cross-family reaction condition recommendation CLI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Simple usage (positional argument)
  python app/cross_family_recommendation_cli.py "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
  
  # Using --rxn flag
  python app/cross_family_recommendation_cli.py --rxn "CCBr.CCN>>CCNCC"
  
  # Get more precedents
  python app/cross_family_recommendation_cli.py --rxn "..." --k 100
  
  # Show more precedents per recommendation
  python app/cross_family_recommendation_cli.py --rxn "..." --max-precedents 5
  
  # Disable DRFP (use feature-based similarity only)
  python app/cross_family_recommendation_cli.py --rxn "..." --no-drfp
        """
    )
    
    parser.add_argument(
        "reaction",
        nargs="?",
        help="Reaction SMILES (reactants>>product)"
    )
    parser.add_argument(
        "--rxn",
        dest="reaction_flag",
        help="Reaction SMILES (alternative to positional argument)"
    )
    parser.add_argument(
        "--k",
        type=int,
        default=50,
        help="Number of precedents to retrieve (default: 50)"
    )
    parser.add_argument(
        "--max-precedents",
        type=int,
        default=3,
        help="Maximum number of precedents to show per recommendation (default: 3)"
    )
    parser.add_argument(
        "--no-drfp",
        action="store_true",
        help="Disable DRFP and use feature-based similarity only"
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print debug information"
    )
    
    args = parser.parse_args()
    
    # Get reaction SMILES from either positional or flag argument
    reaction_smiles = args.reaction or args.reaction_flag
    
    if not reaction_smiles:
        parser.print_help()
        print("\nError: No reaction SMILES provided")
        print("Use: python app/cross_family_recommendation_cli.py \"reaction>>product\"")
        sys.exit(1)
    
    # Validate reaction SMILES format
    if ">>" not in reaction_smiles:
        print(f"Error: Invalid reaction SMILES format")
        print(f"Expected format: reactants>>product")
        print(f"Got: {reaction_smiles}")
        sys.exit(1)
    
    try:
        # Build relax options
        relax = {}
        if args.no_drfp:
            relax["use_drfp"] = False
        
        if args.debug:
            print(f"Debug: Calling chem.recommend.conditions()")
            print(f"  reaction: {reaction_smiles}")
            print(f"  k: {args.k}")
            print(f"  search_all_families: True")
            print(f"  relax: {relax}")
            print()
        
        # Get cross-family recommendation
        result = chem.recommend.conditions(
            reaction=reaction_smiles,
            k=args.k,
            search_all_families=True,
            relax=relax if relax else None
        )
        
        if args.debug:
            print(f"Debug: Result keys: {list(result.keys())}")
            print()
        
        # Print results
        print_recommendation(result, reaction_smiles, args.max_precedents)
        
        # Exit with appropriate code
        recommendations = result.get("recommendations", [])
        if recommendations:
            sys.exit(0)
        else:
            sys.exit(1)
    
    except Exception as e:
        print(f"\nError: {str(e)}")
        if args.debug:
            import traceback
            print("\nDebug: Full traceback:")
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

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


def parse_catalyst_from_precedents(precedents: list) -> str:
    """
    Extract catalyst information from precedent list.
    
    Looks at the first few precedents to find catalyst info from:
    1. condition_core field if it contains catalyst info (e.g., "Pd/XPhos")
    2. Parsed catalytic_system (metal + ligand) - but skip if just CAS numbers
    """
    catalyst_names = []
    
    for prec in precedents[:3]:  # Check first 3 precedents
        # Try condition_core first - it often has good catalyst info
        core = prec.get('core', '')
        
        # Skip if it's just a base/acid indicator
        if core and not any(core.startswith(prefix) for prefix in ['Base:', 'Acid:', 'base:', 'acid:']):
            # Check if it looks like a real catalyst (has "/" or chemical name, not just CAS)
            if '/' in core or any(char.isalpha() for char in core):
                # Looks like "Pd/XPhos" or a chemical name
                if core not in catalyst_names:
                    catalyst_names.append(core)
        
        # Also check catalysts array (structured)
        catalysts = prec.get('catalysts', [])
        for cat in catalysts:
            name = cat.get('name') or cat.get('abbreviation')
            # Skip if it's just a CAS number (all digits and hyphens)
            if name and not all(c.isdigit() or c == '-' for c in name.replace('-', '')):
                if name not in catalyst_names:
                    catalyst_names.append(name)
    
    # Return the most common or first found
    if catalyst_names:
        return catalyst_names[0]  # Use first/most common
    
    return "N/A"


def extract_base_from_precedents(precedents: list) -> str:
    """Extract base information from precedent list."""
    for prec in precedents[:3]:
        # Try structured reagents array
        reagents = prec.get('reagents', [])
        for reagent in reagents:
            if reagent.get('role', '').upper() in ['BASE', 'base']:
                name = reagent.get('name') or reagent.get('abbreviation')
                if name and not name.startswith('[Unknown'):
                    return name
        
        # Fallback: check if condition_core starts with "Base:"
        core = prec.get('core', '')
        if core.startswith('Base:'):
            return core.replace('Base:', '').strip()
    
    return "N/A"


def extract_solvent_from_precedents(precedents: list) -> str:
    """Extract solvent information from precedent list."""
    for prec in precedents[:3]:
        # Try structured solvents array
        solvents = prec.get('solvents', [])
        if solvents:
            solvent_names = []
            for solv in solvents:
                name = solv.get('name') or solv.get('abbreviation')
                if name and not name.startswith('[Unknown'):
                    solvent_names.append(name)
            if solvent_names:
                return " + ".join(solvent_names)
    
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
        List of matching precedent dicts with full details
    """
    # First, check if summary already has precedents
    summary = rec.get('summary', {})
    if summary.get('precedents'):
        # Return with full details if available
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
            # Return full precedent details for enhanced display
            matched.append(prec)
            
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
        
        # Get FULL precedents for this recommendation (not just summary.precedents)
        matched_precedents = find_matching_precedents(rec, all_precedents, max_count=10)
        
        # Catalyst - parse from precedents for more accurate info
        catalyst = parse_catalyst_from_precedents(matched_precedents)
        if catalyst == "N/A":
            # Fallback to summary.core but parse it
            core = summary.get("core", "N/A")
            if core != "N/A" and not any(core.startswith(prefix) for prefix in ['Base:', 'Acid:', 'base:', 'acid:']):
                catalyst = core
            else:
                # Last resort: extract from chemicals (but skip CAS-only entries)
                metal = extract_chemical_by_role(chemicals, "metal_precursor")
                ligand = extract_chemical_by_role(chemicals, "ligand")
                # Check if these are just CAS numbers
                if metal != "N/A" and all(c.isdigit() or c == '-' for c in metal.replace('-', '')):
                    metal = "N/A"  # Skip CAS-only metal
                if ligand != "N/A" and all(c.isdigit() or c == '-' for c in ligand.replace('-', '')):
                    ligand = "N/A"  # Skip CAS-only ligand
                if metal != "N/A" or ligand != "N/A":
                    parts = [p for p in [metal, ligand] if p != "N/A"]
                    catalyst = "/".join(parts)
        
        # Only print catalyst if we found a real one (not N/A and not just CAS)
        if catalyst != "N/A":
            print(f"  Catalyst: {catalyst}")
        
        # Base - parse from precedents
        base_name = extract_base_from_precedents(matched_precedents)
        if base_name == "N/A":
            # Fallback to summary
            base = summary.get("base", {})
            if isinstance(base, dict):
                base_name = base.get("name") or base.get("abbreviation", "N/A")
                # Skip if it's an "Unknown" placeholder
                if base_name.startswith("[Unknown"):
                    base_name = extract_chemical_by_role(chemicals, "base")
            else:
                base_name = extract_chemical_by_role(chemicals, "base")
        if base_name != "N/A" and not base_name.startswith("[Unknown"):
            print(f"  Base: {base_name}")
        
        # Solvent - parse from precedents
        solvent_name = extract_solvent_from_precedents(matched_precedents)
        if solvent_name == "N/A":
            # Fallback to summary
            solvent = summary.get("solvent", {})
            if isinstance(solvent, dict):
                solvent_name = solvent.get("name") or solvent.get("abbreviation", "N/A")
                # Skip if it's an "Unknown" placeholder
                if solvent_name.startswith("[Unknown"):
                    solvent_name = extract_chemical_by_role(chemicals, "solvent")
            else:
                solvent_name = extract_chemical_by_role(chemicals, "solvent")
        if solvent_name != "N/A" and not solvent_name.startswith("[Unknown"):
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
                yield_pct = prec.get("yield") or prec.get("yield_pct")
                
                # Extract dataset name from reaction_id (e.g., "31-001-CAS-16722741" -> dataset "31")
                dataset_name = "Unknown"
                if reaction_id and "-" in str(reaction_id):
                    parts = str(reaction_id).split("-")
                    if parts:
                        dataset_name = f"Dataset-{parts[0]}"
                
                # Get reaction SMILES
                reaction_smiles = prec.get("reaction_smiles", "N/A")
                if reaction_smiles == "N/A":
                    # Try to construct from reactants and products
                    smiles_data = prec.get("smiles", {})
                    if isinstance(smiles_data, dict):
                        reactants = smiles_data.get("reactants", "")
                        products = smiles_data.get("products", "")
                        if reactants and products:
                            reaction_smiles = f"{reactants}>>{products}"
                
                # Get full catalytic system
                catalytic_system = prec.get("catalytic_system", [])
                catalyst_display = []
                if catalytic_system:
                    for cat in catalytic_system:
                        if isinstance(cat, dict):
                            cat_name = cat.get("name") or cat.get("abbreviation") or cat.get("cas", "")
                        else:
                            cat_name = str(cat)
                        if cat_name and not cat_name.startswith("[Unknown"):
                            catalyst_display.append(cat_name)
                
                # Fallback to condition_core if no catalytic_system
                if not catalyst_display:
                    core = prec.get("core", "") or prec.get("condition_core", "")
                    if core and not any(core.startswith(prefix) for prefix in ['Base:', 'Acid:', 'base:', 'acid:']):
                        catalyst_display.append(core)
                
                # Get full solvent info
                solvents = prec.get("solvents", [])
                solvent_display = []
                if solvents:
                    for solv in solvents:
                        if isinstance(solv, dict):
                            solv_name = solv.get("name") or solv.get("abbreviation", "")
                        else:
                            solv_name = str(solv)
                        if solv_name and not solv_name.startswith("[Unknown"):
                            solvent_display.append(solv_name)
                
                # Get reagents (base, additives, etc.)
                reagents = prec.get("reagents", [])
                base_display = []
                additive_display = []
                if reagents:
                    for reagent in reagents:
                        if isinstance(reagent, dict):
                            reagent_name = reagent.get("name") or reagent.get("abbreviation", "")
                            reagent_role = reagent.get("role", "").upper()
                            if reagent_name and not reagent_name.startswith("[Unknown"):
                                if reagent_role == "BASE":
                                    base_display.append(reagent_name)
                                elif reagent_role in ["ADDITIVE", "REAGENT", "ACTIVATING_AGENT"]:
                                    additive_display.append(f"{reagent_name} ({reagent_role})")
                
                # Get temperature and time
                conditions = prec.get("conditions", {})
                temperature = conditions.get("temperature_c") if isinstance(conditions, dict) else None
                time_h = conditions.get("time_h") if isinstance(conditions, dict) else None
                
                # Extract just the title from reference
                ref_title = ref.split("|")[0].strip() if "|" in ref else (ref[:80] if ref else "N/A")
                
                # Print precedent details
                print(f"    {i}. [{dataset_name}] {reaction_id}")
                
                # Show reaction SMILES
                if reaction_smiles != "N/A":
                    # Truncate if too long
                    if len(reaction_smiles) > 100:
                        reaction_smiles_display = reaction_smiles[:97] + "..."
                    else:
                        reaction_smiles_display = reaction_smiles
                    print(f"       Reaction: {reaction_smiles_display}")
                
                # Show catalytic system
                if catalyst_display:
                    print(f"       Catalytic System: {', '.join(catalyst_display)}")
                
                # Show base
                if base_display:
                    print(f"       Base: {', '.join(base_display)}")
                
                # Show solvent
                if solvent_display:
                    print(f"       Solvent: {', '.join(solvent_display)}")
                
                # Show additives
                if additive_display:
                    print(f"       Additives: {', '.join(additive_display)}")
                
                # Show temperature
                if temperature is not None:
                    print(f"       Temperature: {temperature}°C")
                
                # Show time
                if time_h is not None:
                    print(f"       Time: {time_h}h")
                
                # Show yield
                if yield_pct is not None:
                    print(f"       Yield: {yield_pct}%")
                
                # Show reference
                if ref_title and ref_title != "N/A":
                    print(f"       Reference: {ref_title}...")
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

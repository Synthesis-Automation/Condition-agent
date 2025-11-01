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


def run_recommendation(
    reaction_smiles: str,
    *,
    k: int,
    max_precedents: int,
    use_drfp: bool,
    debug: bool,
) -> bool:
    """Execute a single recommendation request."""
    relax = {}
    if not use_drfp:
        relax["use_drfp"] = False

    if debug:
        print_separator("-")
        print("Debug: Calling chem.recommend.conditions()")
        print(f"  reaction: {reaction_smiles}")
        print(f"  k: {k}")
        print(f"  search_all_families: True")
        print(f"  relax: {relax if relax else None}")
        print_separator("-")

    result = chem.recommend.conditions(
        reaction=reaction_smiles,
        k=k,
        search_all_families=True,
        relax=relax if relax else None,
    )

    if debug:
        print(f"Debug: Result keys: {list(result.keys())}")
        print_separator("-")

    print_recommendation(result, reaction_smiles, max_precedents)
    recommendations = result.get("recommendations", [])
    return bool(recommendations)


def interactive_session(
    *,
    k: int,
    max_precedents: int,
    use_drfp: bool,
    debug: bool,
) -> None:
    """Interactive REPL for entering reactions and tweaking options."""
    print_separator()
    print("Cross-family Recommendation Interactive Mode")
    print("Type a reaction SMILES (reactants>>product) to get recommendations.")
    print("Commands:")
    print("  /k <int>              -> set number of precedents to search (current: %d)" % k)
    print(
        "  /precedents <int>     -> set precedents to show per recommendation (current: %d)"
        % max_precedents
    )
    print("  /drfp on|off          -> enable or disable DRFP similarity (currently: %s)" %
          ("on" if use_drfp else "off"))
    print("  /debug on|off         -> toggle debug logging (currently: %s)" %
          ("on" if debug else "off"))
    print("  /show                 -> display current settings")
    print("  /help                 -> show this message again")
    print("  /exit or /quit        -> leave interactive mode")
    print_separator()

    current_k = k
    current_max_precedents = max_precedents
    current_use_drfp = use_drfp
    current_debug = debug

    while True:
        try:
            raw = input("reaction> ").strip()
        except KeyboardInterrupt:
            print("\nExiting.")
            break
        except EOFError:
            print("\n")
            break

        if not raw:
            continue

        lowered = raw.lower()
        if lowered in {"/exit", "/quit"}:
            break
        if lowered == "/help":
            print_separator()
            print("Commands:")
            print("  /k <int>")
            print("  /precedents <int>")
            print("  /drfp on|off")
            print("  /debug on|off")
            print("  /show")
            print("  /help")
            print("  /exit or /quit")
            print_separator()
            continue
        if lowered == "/show":
            print_separator("-")
            print(f"k: {current_k}")
            print(f"max_precedents: {current_max_precedents}")
            print(f"use_drfp: {current_use_drfp}")
            print(f"debug: {current_debug}")
            print_separator("-")
            continue
        if lowered.startswith("/k "):
            value = raw.split(maxsplit=1)[1]
            if value.isdigit():
                current_k = max(1, int(value))
                print(f"Updated k to {current_k}")
            else:
                print("Invalid value for /k. Example: /k 50")
            continue
        if lowered.startswith("/precedents "):
            value = raw.split(maxsplit=1)[1]
            if value.isdigit():
                current_max_precedents = max(1, int(value))
                print(f"Updated max precedents to {current_max_precedents}")
            else:
                print("Invalid value for /precedents. Example: /precedents 3")
            continue
        if lowered.startswith("/drfp "):
            value = raw.split(maxsplit=1)[1]
            if value in {"on", "off"}:
                current_use_drfp = value == "on"
                print(f"DRFP similarity {'enabled' if current_use_drfp else 'disabled'}.")
            else:
                print("Invalid value for /drfp. Use /drfp on or /drfp off.")
            continue
        if lowered.startswith("/debug "):
            value = raw.split(maxsplit=1)[1]
            if value in {"on", "off"}:
                current_debug = value == "on"
                print(f"Debug logging {'enabled' if current_debug else 'disabled'}.")
            else:
                print("Invalid value for /debug. Use /debug on or /debug off.")
            continue
        if raw.startswith("/"):
            print("Unknown command. Type /help for options.")
            continue

        if ">>" not in raw:
            print("Invalid reaction format. Expected reactants>>product.")
            continue

        try:
            success = run_recommendation(
                raw,
                k=current_k,
                max_precedents=current_max_precedents,
                use_drfp=current_use_drfp,
                debug=current_debug,
            )
            if not success:
                print("No recommendations returned.")
        except Exception as exc:
            print(f"Error: {exc}")
            if current_debug:
                import traceback

                traceback.print_exc()
        print_separator()

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


ROLE_GROUPS = [
    ("Catalyst(s)", {"metal_catalyst", "catalyst"}),
    ("Ligand(s)", {"ligand"}),
    ("Base(s)", {"base"}),
    ("Acid(s)", {"acid"}),
    ("Oxidant(s)", {"oxidant"}),
    ("Reductant(s)", {"reductant"}),
    ("Solvent(s)", {"solvent"}),
    (
        "Additive(s)",
        {"additive", "activating_agent", "coupling_agent", "reagent", "buffer", "promoter"},
    ),
]


DETAIL_LABELS = [
    ("mol_percent", "mol%"),
    ("mol%", "mol%"),
    ("mol_percent_range", "mol% range"),
    ("wt_percent", "wt%"),
    ("weight_percent", "wt%"),
    ("volume_percent", "vol%"),
    ("v_v_percent", "vol%"),
    ("mass_percent", "mass%"),
]


def _format_number(value):
    """Format numeric values with sensible precision."""
    try:
        if isinstance(value, int):
            return str(value)
        if isinstance(value, float):
            if value.is_integer():
                return str(int(value))
            return f"{value:.3g}"
    except Exception:
        pass
    return str(value)


def _format_loading(loading):
    """Format catalyst loading dictionaries or values."""
    if not loading:
        return None
    if isinstance(loading, dict):
        range_val = loading.get("range")
        if isinstance(range_val, (list, tuple)) and len(range_val) == 2:
            low, high = range_val
            text = f"{_format_number(low)}-{_format_number(high)}"
        else:
            value = loading.get("value")
            if value is not None:
                text = _format_number(value)
            else:
                text = None
        unit = loading.get("unit") or loading.get("units")
        if text and unit:
            return f"{text} {unit}"
        return text
    if isinstance(loading, (int, float, str)):
        return str(loading)
    return None


def _extract_detail_value(chemical, extras, key):
    """Fetch detail value from chemical or extras block."""
    if isinstance(chemical, dict) and key in chemical and chemical[key] is not None:
        return chemical[key]
    if extras and key in extras and extras[key] is not None:
        return extras[key]
    return None


def _format_condition_value(value):
    """Format condition dictionaries/numbers into readable strings."""
    if value is None:
        return None
    if isinstance(value, dict):
        range_val = value.get("range")
        if isinstance(range_val, (list, tuple)) and len(range_val) == 2:
            low, high = range_val
            text = f"{_format_number(low)}-{_format_number(high)}"
        elif "value" in value:
            text = _format_number(value["value"])
        else:
            text = value.get("text") or value.get("note") or None
        unit = value.get("unit") or value.get("units")
        if text and unit:
            return f"{text} {unit}"
        return text or unit
    if isinstance(value, (int, float)):
        return _format_number(value)
    if isinstance(value, (list, tuple)):
        return ", ".join(_format_number(v) for v in value)
    if isinstance(value, str):
        return value
    return str(value)


def build_full_condition_lines(chemicals, conditions) -> list:
    """Build human-readable full condition lines for a recommendation."""
    groups = {label: [] for label, _ in ROLE_GROUPS}
    other_entries = []

    if not isinstance(chemicals, list):
        chemicals = []

    for chemical in chemicals:
        if not isinstance(chemical, dict):
            continue
        role = (chemical.get("role") or "").lower()
        name = (
            chemical.get("name")
            or chemical.get("abbreviation")
            or chemical.get("smiles")
            or chemical.get("cas")
        )
        if not name or name.startswith("[Unknown"):
            continue

        details = []

        equivalents = chemical.get("equivalents")
        if equivalents not in (None, "", "N/A"):
            details.append(f"{_format_number(equivalents)} eq")

        loading_detail = _format_loading(chemical.get("loading"))
        if loading_detail:
            details.append(f"loading {loading_detail}")

        extras = chemical.get("extras") or {}
        for key, label in DETAIL_LABELS:
            value = _extract_detail_value(chemical, extras, key)
            if value is not None:
                if isinstance(value, (list, tuple)) and len(value) == 2:
                    formatted = f"{_format_number(value[0])}-{_format_number(value[1])}"
                elif isinstance(value, (int, float)):
                    formatted = _format_number(value)
                else:
                    formatted = str(value)
                details.append(f"{formatted} {label}")

        entry = name
        if details:
            entry += f" ({', '.join(details)})"

        matched_group = False
        for label, roles in ROLE_GROUPS:
            if role in roles:
                groups[label].append(entry)
                matched_group = True
                break
        if not matched_group:
            descriptor = role.replace("_", " ") if role else "reagent"
            other_entries.append(f"{descriptor.title()}: {entry}")

    lines = []
    for label, _ in ROLE_GROUPS:
        entries = groups[label]
        if entries:
            lines.append(f"{label}: " + "; ".join(entries))

    if other_entries:
        lines.append("Other reagents: " + "; ".join(other_entries))

    condition_data = conditions if isinstance(conditions, dict) else {}
    temp_text = _format_condition_value(
        condition_data.get("temperature") or condition_data.get("temperature_c")
    )
    if temp_text:
        lines.append(f"Temperature: {temp_text}")

    time_text = _format_condition_value(
        condition_data.get("time") or condition_data.get("time_h") or condition_data.get("time_hr")
    )
    if time_text:
        lines.append(f"Time: {time_text}")

    atmosphere = condition_data.get("atmosphere") or condition_data.get("environment")
    if atmosphere:
        lines.append(f"Atmosphere: {atmosphere}")

    extras = condition_data.get("extras") if isinstance(condition_data.get("extras"), dict) else {}
    for key, value in extras.items():
        formatted = _format_condition_value(value)
        label = key.replace("_", " ").title()
        if formatted:
            lines.append(f"{label}: {formatted}")
        else:
            lines.append(f"{label}: {value}")

    note = condition_data.get("note")
    if note:
        lines.append(f"Note: {note}")

    return lines


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
        raw_conditions = rec.get("conditions") or {}
        conditions = raw_conditions if isinstance(raw_conditions, dict) else {"note": raw_conditions}
        
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
                metal = extract_chemical_by_role(chemicals, "metal_catalyst")
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
        temp_text = _format_condition_value(
            conditions.get("temperature") or conditions.get("temperature_c")
        )
        if temp_text:
            print(f"  Temperature: {temp_text}")
    
        # Time
        time_text = _format_condition_value(
            conditions.get("time") or conditions.get("time_h") or conditions.get("time_hr")
        )
        if time_text:
            print(f"  Time: {time_text}")
        
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

        full_condition_lines = build_full_condition_lines(chemicals, conditions)
        if full_condition_lines:
            print("  Full Conditions:")
            for line in full_condition_lines:
                print(f"    - {line}")
        
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
                
                # Extract dataset name from rxn_type field (used internally by precedent loader)
                # The reaction_id format is "31-XXX-CAS-YYYYYYY" where 31 is the SciFinder dataset ID
                # But we want to show the actual reaction family/type from rxn_type field
                dataset_name = prec.get("rxn_type") or prec.get("reaction_type", "Unknown")
                if dataset_name != "Unknown":
                    # Map canonical family names back to dataset filenames
                    # The loader normalizes names like "C_N_Coupling" -> "C_N_Coupling_Pd" or "C_N_Coupling_Cu"
                    # But the actual JSONL files use simpler names like "C_N_Coupling.jsonl"
                    # So we need to reverse-map these
                    family_to_file = {
                        "C_N_Coupling_Pd": "C_N_Coupling",
                        "C_N_Coupling_Cu": "C_N_Coupling",
                        "C_N_Coupling_Ni": "C_N_Coupling",
                        "Suzuki": "Suzuki",
                        "Amide_formation": "Amide_formation",
                        "C_O_Coupling": "C_O_Coupling",
                        "C_S_Coupling": "C_S_Coupling",
                        "SNAr-CN": "SNAr-CN",
                        "SNAr-CO": "SNAr-CO",
                        "SNAr-CS": "SNAr-CS",
                    }
                    # Use mapping if available, otherwise use the raw value
                    file_base = family_to_file.get(dataset_name, dataset_name)
                    dataset_name = f"{file_base}.jsonl"
                
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
                            # Try name, abbreviation, then CAS as fallback
                            cat_name = cat.get("name") or cat.get("abbreviation")
                            if not cat_name:
                                # If no name/abbrev, show CAS with label
                                cas_num = cat.get("cas", "")
                                if cas_num and not cas_num.startswith("[Unknown"):
                                    cat_name = f"CAS:{cas_num}"
                                else:
                                    cat_name = ""
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
                            solv_name = solv.get("name") or solv.get("abbreviation")
                            if not solv_name:
                                # If no name/abbrev, show CAS with label
                                cas_num = solv.get("cas", "")
                                if cas_num and not cas_num.startswith("[Unknown"):
                                    solv_name = f"CAS:{cas_num}"
                                else:
                                    solv_name = ""
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
                            reagent_name = reagent.get("name") or reagent.get("abbreviation")
                            if not reagent_name:
                                # If no name/abbrev, show CAS with label
                                cas_num = reagent.get("cas", "")
                                if cas_num and not cas_num.startswith("[Unknown"):
                                    reagent_name = f"CAS:{cas_num}"
                                else:
                                    reagent_name = ""
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
    
    parser.add_argument("reaction", nargs="?", help="Reaction SMILES (reactants>>product)")
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
    
    use_drfp = not args.no_drfp
    reaction_smiles = args.reaction or args.reaction_flag

    if not reaction_smiles:
        interactive_session(
            k=args.k,
            max_precedents=args.max_precedents,
            use_drfp=use_drfp,
            debug=args.debug,
        )
        return

    if ">>" not in reaction_smiles:
        print("Error: Invalid reaction SMILES format")
        print("Expected format: reactants>>product")
        print(f"Got: {reaction_smiles}")
        sys.exit(1)

    try:
        success = run_recommendation(
            reaction_smiles,
            k=args.k,
            max_precedents=args.max_precedents,
            use_drfp=use_drfp,
            debug=args.debug,
        )
        sys.exit(0 if success else 1)
    except Exception as exc:
        print(f"\nError: {exc}")
        if args.debug:
            import traceback

            print("\nDebug: Full traceback:")
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

"""
Rule-to-Protocol Format Converter.

Converts rule-based condition dictionaries to protocol-like reaction_setup format
for automation-ready output. Reuses the proven protocol structure from protocol_db_v2.

Key Features:
- Converts simple rule conditions → structured reaction_setup
- Maintains protocol-compatible format
- Handles ranges (picks midpoint)
- Handles options (picks first)
- Standard addition order based on chemistry best practices

Usage:
    >>> from chemtools.formatters.rule_to_protocol import rule_conditions_to_reaction_setup
    >>> conditions = {
    ...     "catalyst": "PdCl2(PPh3)2",
    ...     "catalyst_loading_molpct": "1.0-3.0",
    ...     "base": "Et3N",
    ...     "base_equiv": "2.0",
    ...     "solvent": "THF",
    ...     "temperature_C": "60-80"
    ... }
    >>> result = rule_conditions_to_reaction_setup(conditions, scale_mmol=1.0)
    >>> print(result["reaction_setup"][0]["chemicals"])
"""

from __future__ import annotations
from typing import Dict, Any, List, Optional, Tuple
import re


# Standard addition order based on chemistry best practices
# Matches protocol format: chemicals array is in addition order
ADDITION_ORDER_PRIORITY = {
    "solvent": 1,
    "base": 2,
    "ligand": 3,
    "metal_catalyst": 4,
    "catalyst": 5,
    "additive": 6,
    "starting_material": 7,
    "reagent": 8,
}

# Mapping from rule field names to protocol roles
RULE_FIELD_TO_ROLE = {
    # Direct mappings
    "solvent": "solvent",
    "base": "base",
    "ligand": "ligand",
    
    # Metal catalysts (family-specific)
    "catalyst": "metal_catalyst",
    "pd_precatalyst": "metal_catalyst",
    "pd_source": "metal_catalyst",
    "ru_catalyst": "metal_catalyst",
    "cu_source": "metal_catalyst",
    "metal_catalyst": "metal_catalyst",
    
    # Other
    "additives": "additive",
    "additive": "additive",
    "coupling_system": "reagent",
    "reducing_agent": "reagent",
}


def rule_conditions_to_reaction_setup(
    conditions: Dict[str, Any],
    user_substrates: Optional[List[Dict[str, Any]]] = None,
    scale_mmol: float = 1.0,
    reaction_family: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Convert rule conditions to protocol-like reaction_setup format.
    
    Args:
        conditions: Rule conditions dictionary (from rule JSON)
        user_substrates: Optional list of user-provided substrates with:
            - name: Chemical name
            - smiles: SMILES string
            - role: "starting_material" (default)
            - mmol: Amount in mmol (defaults to scale_mmol)
            - equivalents: Stoichiometric equivalents (defaults to 1.0)
        scale_mmol: Reaction scale in mmol (default 1.0)
        reaction_family: Optional reaction family name for context
    
    Returns:
        Dict with "reaction_setup" key containing protocol-compatible structure:
        {
            "reaction_setup": [
                {
                    "chemicals": [
                        {"name": "...", "role": "...", "amount": {...}},
                        ...
                    ],
                    "conditions": [
                        {"temperature_C": 60, "time_h": 4, ...}
                    ]
                }
            ],
            "metadata": {
                "generated_from": "rule",
                "format": "protocol-compatible",
                "scale_mmol": 1.0
            }
        }
    
    Example:
        >>> conditions = {
        ...     "catalyst": "PdCl2(PPh3)2",
        ...     "catalyst_loading_molpct": "1.0",
        ...     "base": "Et3N",
        ...     "base_equiv": "2.0",
        ...     "solvent": "THF"
        ... }
        >>> result = rule_conditions_to_reaction_setup(conditions, scale_mmol=1.0)
        >>> len(result["reaction_setup"][0]["chemicals"])
        4  # solvent, base, catalyst, substrate placeholder
    """
    chemicals = []
    
    # 1. Solvent (priority 1)
    if "solvent" in conditions:
        solvent_name = _pick_first_option(conditions["solvent"])
        chemicals.append({
            "name": solvent_name,
            "abbreviation": _get_abbreviation(solvent_name),
            "role": "solvent",
            "amount": {
                "volume_ml": None,  # To be calculated based on concentration
                "note": "Volume determined by target concentration (typically 0.1-0.5 M)"
            },
            "_addition_order": ADDITION_ORDER_PRIORITY["solvent"]
        })
    
    # 2. Base (priority 2)
    if "base" in conditions:
        base_equiv = _parse_range_midpoint(conditions.get("base_equiv", "2.0"))
        base_name = _pick_first_option(conditions["base"])
        chemicals.append({
            "name": base_name,
            "abbreviation": _get_abbreviation(base_name),
            "role": "base",
            "amount": {
                "mmol": base_equiv * scale_mmol,
                "equivalents": base_equiv,
                "note": f"Based on rule conditions"
            },
            "_addition_order": ADDITION_ORDER_PRIORITY["base"]
        })
    
    # 3. Ligand (priority 3) - if separate from catalyst
    if "ligand" in conditions:
        ligand_str = conditions["ligand"]
        # Skip if it's built-in to catalyst
        if not any(skip in ligand_str.lower() for skip in ["built-in", "built in", "integrated"]):
            ligand_loading = _parse_range_midpoint(conditions.get("ligand_loading_molpct", "2.0"))
            ligand_name = _pick_first_option(ligand_str)
            chemicals.append({
                "name": ligand_name,
                "abbreviation": _get_abbreviation(ligand_name),
                "role": "ligand",
                "amount": {
                    "mmol": (ligand_loading / 100.0) * scale_mmol,
                    "equivalents": ligand_loading / 100.0,
                    "note": f"{ligand_loading} mol%"
                },
                "_addition_order": ADDITION_ORDER_PRIORITY["ligand"]
            })
    
    # 4. Metal catalyst (priority 4)
    # Check multiple possible field names
    catalyst_fields = ["catalyst", "pd_precatalyst", "pd_source", "ru_catalyst", "cu_source", "metal_catalyst", "precatalyst"]
    for field in catalyst_fields:
        if field in conditions:
            cat_loading = _parse_range_midpoint(conditions.get("catalyst_loading_molpct", "1.0"))
            # Also check family-specific loading fields
            if "pd_loading_molpct" in conditions:
                cat_loading = _parse_range_midpoint(conditions["pd_loading_molpct"])
            elif "cu_loading_molpct" in conditions:
                cat_loading = _parse_range_midpoint(conditions["cu_loading_molpct"])
            
            cat_name = _pick_first_option(conditions[field])
            chemicals.append({
                "name": cat_name,
                "abbreviation": _get_abbreviation(cat_name),
                "role": "metal_catalyst",
                "amount": {
                    "mmol": (cat_loading / 100.0) * scale_mmol,
                    "equivalents": cat_loading / 100.0,
                    "note": f"{cat_loading} mol%"
                },
                "_addition_order": ADDITION_ORDER_PRIORITY["metal_catalyst"]
            })
            break  # Only add one catalyst
    
    # 5. Additives (priority 6)
    if "additives" in conditions:
        additives = conditions["additives"]
        if isinstance(additives, list):
            for idx, additive in enumerate(additives):
                chemicals.append({
                    "name": additive,
                    "role": "additive",
                    "amount": {
                        "note": "As specified in rule conditions"
                    },
                    "_addition_order": ADDITION_ORDER_PRIORITY["additive"] + idx * 0.1
                })
        elif isinstance(additives, str):
            chemicals.append({
                "name": additives,
                "role": "additive",
                "amount": {
                    "note": "As specified in rule conditions"
                },
                "_addition_order": ADDITION_ORDER_PRIORITY["additive"]
            })
    
    # 6. User substrates or placeholder (priority 7)
    if user_substrates:
        for idx, substrate in enumerate(user_substrates):
            sub_mmol = substrate.get("mmol", scale_mmol)
            sub_equiv = substrate.get("equivalents", 1.0)
            
            chem = {
                "name": substrate.get("name", f"Substrate {idx + 1}"),
                "role": substrate.get("role", "starting_material"),
                "amount": {
                    "mmol": sub_mmol,
                    "equivalents": sub_equiv
                },
                "_addition_order": ADDITION_ORDER_PRIORITY["starting_material"] + idx * 0.1
            }
            
            if "smiles" in substrate:
                chem["smiles"] = substrate["smiles"]
            if "cas" in substrate:
                chem["cas"] = substrate["cas"]
            
            chemicals.append(chem)
    else:
        # Placeholder for user-provided substrate
        chemicals.append({
            "name": "Substrate (user-provided)",
            "role": "starting_material",
            "amount": {
                "mmol": scale_mmol,
                "equivalents": 1.0,
                "note": "Replace with actual substrate information"
            },
            "_addition_order": ADDITION_ORDER_PRIORITY["starting_material"]
        })
    
    # Sort by addition order
    chemicals.sort(key=lambda x: x["_addition_order"])
    
    # Remove internal _addition_order field
    for chem in chemicals:
        del chem["_addition_order"]
    
    # Build reaction conditions (protocol format)
    reaction_conditions = {}
    
    if "temperature_C" in conditions:
        reaction_conditions["temperature_C"] = _parse_range_midpoint(conditions["temperature_C"])
    
    if "time_h" in conditions:
        reaction_conditions["time_h"] = _parse_range_midpoint(conditions["time_h"])
    
    if "atmosphere" in conditions:
        reaction_conditions["atmosphere"] = conditions["atmosphere"]
    
    if "pressure_bar" in conditions:
        reaction_conditions["pressure_bar"] = _parse_range_midpoint(conditions["pressure_bar"])
    
    # Build final structure (protocol-compatible)
    return {
        "reaction_setup": [
            {
                "chemicals": chemicals,
                "conditions": [reaction_conditions] if reaction_conditions else []
            }
        ],
        "metadata": {
            "generated_from": "rule",
            "format": "protocol-compatible",
            "scale_mmol": scale_mmol,
            "reaction_family": reaction_family
        }
    }


def _pick_first_option(value: str) -> str:
    """
    Pick first option from string with multiple alternatives.
    
    Examples:
        "THF or toluene" → "THF"
        "PdCl2(PPh3)2 or Pd(dppf)Cl2·DCM" → "PdCl2(PPh3)2"
        "THF, toluene, or DMF" → "THF"
    """
    value = str(value).strip()
    
    # Handle comma-separated with "or" (e.g., "A, B, or C")
    if ", or " in value:
        return value.split(",")[0].strip()
    
    # Handle " or " separator
    if " or " in value:
        return value.split(" or ")[0].strip()
    
    # Handle "/ " separator (sometimes used)
    if " / " in value or "/ " in value:
        return value.split("/")[0].strip()
    
    return value


def _parse_range_midpoint(value: str | float | int) -> float:
    """
    Parse range or single value and return midpoint.
    
    Examples:
        "0.5-2.0" → 1.25
        "60-80" → 70.0
        "80–100" → 90.0  (en-dash)
        "1.5" → 1.5
        1.5 → 1.5
    """
    # Already a number
    if isinstance(value, (int, float)):
        return float(value)
    
    value_str = str(value).strip()
    
    # Normalize various dash types to regular hyphen
    value_str = value_str.replace('–', '-')  # en-dash
    value_str = value_str.replace('—', '-')  # em-dash
    
    # Handle range with hyphen
    if "-" in value_str:
        # Split on hyphen (but be careful of negative numbers)
        parts = value_str.split("-")
        # Filter out empty parts (from leading/trailing hyphens)
        parts = [p.strip() for p in parts if p.strip()]
        
        if len(parts) >= 2:
            try:
                low = float(parts[0])
                high = float(parts[1])
                return (low + high) / 2.0
            except ValueError:
                pass  # Fall through to try parsing as single value
    
    # Try parsing as single numeric value
    try:
        return float(value_str)
    except ValueError:
        # Return default if can't parse
        return 1.0


def _get_abbreviation(chemical_name: str) -> Optional[str]:
    """
    Try to extract common abbreviation from chemical name.
    
    Examples:
        "Triethylamine" → "Et3N"
        "N,N-Dimethylformamide" → "DMF"
        "Palladium(II) acetate" → "Pd(OAc)2"
    """
    # Common abbreviations (case-insensitive)
    abbreviations = {
        "triethylamine": "Et3N",
        "diisopropylethylamine": "DIPEA",
        "dimethylformamide": "DMF",
        "tetrahydrofuran": "THF",
        "dimethyl sulfoxide": "DMSO",
        "acetonitrile": "MeCN",
        "dichloromethane": "DCM",
        "ethyl acetate": "EtOAc",
        "methanol": "MeOH",
        "ethanol": "EtOH",
        "tert-butanol": "tBuOH",
        "acetic acid": "AcOH",
    }
    
    name_lower = chemical_name.lower().strip()
    
    for full_name, abbrev in abbreviations.items():
        if full_name in name_lower:
            return abbrev
    
    # If already looks like abbreviation (short, has uppercase), return it
    if len(chemical_name) <= 10 and any(c.isupper() for c in chemical_name):
        # Extract abbreviation pattern
        match = re.search(r'([A-Z][a-z]*\d*[A-Z]*\d*[A-Z]*)', chemical_name)
        if match:
            return match.group(1)
    
    return None

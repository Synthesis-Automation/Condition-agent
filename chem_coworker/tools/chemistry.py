"""
Core chemistry tools for ChemCoworker.

Provides lightweight wrappers over chemtools and reaction_agent internals:
  - normalize_reaction     : parse and validate reaction / molecule SMILES
  - detect_reaction_type   : deterministic taxonomy classification
  - analyze_bond_changes   : which bonds formed/broke (anchor tool)
  - inspect_functional_groups : what FGs are present
  - get_molecular_descriptors : RDKit molecular properties

All functions follow the _success() / _error() contract and are registered
as ToolPlugins so the ToolRegistry can build execution groups automatically.
"""
from __future__ import annotations

from typing import Any, Dict

from ._helpers import _error, _success, _to_jsonable
from ._base import ToolPlugin


# ---------------------------------------------------------------------------
# Tool 1: normalize_reaction
# ---------------------------------------------------------------------------

def _normalize_reaction(smiles: str) -> Dict[str, Any]:
    """Parse, validate, and canonicalize a reaction or molecule SMILES.

    Works for both reaction SMILES (with '>>' or '>') and single molecule SMILES.
    Returns the normalized components and any detected issues.

    Args:
        smiles: Reaction SMILES ('reactants>>products') or molecule SMILES.

    Returns:
        dict with success, is_reaction, normalized_smiles, reactants, products, warnings.
    """
    try:
        from chemtools.featurizers.analysis.smiles import normalize, normalize_reaction

        is_reaction = ">>" in smiles or (">" in smiles and smiles.count(">") >= 1)

        if is_reaction:
            result = normalize_reaction(smiles)
            if "error" in result:
                return _error(f"Invalid reaction SMILES: {result['error']}")
            return _success({
                "is_reaction": True,
                "input_smiles": smiles,
                "reactants": result.get("reactants", []),
                "agents": result.get("agents", []),
                "products": result.get("products", []),
                "warnings": result.get("warnings", []),
            })
        else:
            result = normalize(smiles)
            if "error" in result:
                return _error(f"Invalid SMILES: {result['error']}")
            return _success({
                "is_reaction": False,
                "input_smiles": smiles,
                "normalized_smiles": result.get("smiles_norm", result.get("largest_smiles", smiles)),
                "fragments": result.get("fragments", [smiles]),
                "warnings": [],
            })
    except Exception as exc:
        return _error(f"SMILES normalization failed: {exc}")


normalize_reaction_tool = ToolPlugin(
    name="normalize_reaction",
    category="chemistry",
    description="Parse and canonicalize a reaction SMILES or molecule SMILES. Always run this first for any SMILES input.",
    prerequisites=[],
    fn=_normalize_reaction,
)


# ---------------------------------------------------------------------------
# Tool 2: detect_reaction_type
# ---------------------------------------------------------------------------

def _detect_reaction_type(reaction_smiles: str) -> Dict[str, Any]:
    """Deterministic taxonomy classification of a reaction.

    Uses the chemtools reaction family detector to match the reaction to a
    known reaction type (e.g., Suzuki_miyaura, Buchwald_Hartwig, C_N_Coupling).
    This is fast and does not require atom mapping.

    Args:
        reaction_smiles: Reaction SMILES in 'reactants>>products' format.

    Returns:
        dict with success, reaction_type, family_label, confidence.
    """
    try:
        from chemtools.featurizers.analysis.reactions import analyze_reaction

        result = analyze_reaction(reaction_smiles)
        if not result:
            return _error("Reaction analysis returned no result")

        family = result.get("family") or {}
        return _success({
            "reaction_smiles": reaction_smiles,
            "reaction_type": family.get("canonical_id") or result.get("reaction_type"),
            "family_label": family.get("label") or "",
            "confidence": result.get("confidence", 0.0),
            "matched_slots": result.get("matched_slots", []),
        })
    except Exception as exc:
        # Fallback: try simpler chemtools normalize + detect
        try:
            from chemtools.smiles import normalize_reaction as nr
            norm = nr(reaction_smiles)
            return _success({
                "reaction_smiles": reaction_smiles,
                "reaction_type": norm.get("reaction_type"),
                "family_label": norm.get("family_label", ""),
                "confidence": norm.get("confidence", 0.0),
                "matched_slots": [],
            })
        except Exception:
            return _error(f"Reaction type detection failed: {exc}")


detect_reaction_type_tool = ToolPlugin(
    name="detect_reaction_type",
    category="chemistry",
    description="Deterministic taxonomy classification of a reaction (Suzuki, Buchwald-Hartwig, etc.). Fast, no atom mapping needed.",
    prerequisites=[],
    fn=_detect_reaction_type,
)


# ---------------------------------------------------------------------------
# Tool 3: analyze_bond_changes
# ---------------------------------------------------------------------------

def _analyze_bond_changes(reaction_smiles: str) -> Dict[str, Any]:
    """Identify exactly which bonds were broken and formed in a reaction.

    This is the ANCHOR tool for reaction analysis. Always run this before
    searching the taxonomy. Bond changes reveal:
      - Which bonds formed → reaction family (C-C, C-N, C-O, C-S, ring closure...)
      - Which leaving groups departed
      - Whether it's a substitution, addition, or rearrangement

    Args:
        reaction_smiles: Reaction SMILES in 'reactants>>products' format.

    Returns:
        dict with bonds_broken, bonds_formed, key_bond_type, leaving_groups.
    """
    try:
        from chemtools.util.reaction_center_detector import analyze_bond_changes_hybrid

        result = analyze_bond_changes_hybrid(reaction_smiles)
        if not result:
            return _error("Bond change analysis returned no result")

        # Use recommended_result if available (hybrid analyzer)
        if isinstance(result, dict) and "recommended_result" in result:
            rec = result["recommended_result"]
        else:
            rec = result

        broken = rec.get("broken_bonds") or rec.get("bonds_broken", [])
        formed = rec.get("formed_bonds") or rec.get("bonds_formed", [])

        # Derive key bond type from formed bonds
        key_bond = _infer_key_bond_type(formed)

        return _success({
            "reaction_smiles": reaction_smiles,
            "bonds_broken": _to_jsonable(broken),
            "bonds_formed": _to_jsonable(formed),
            "key_bond_type": key_bond,
            "leaving_groups": _to_jsonable(rec.get("leaving_groups", [])),
            "mapping_confidence": rec.get("confidence", rec.get("mapping_confidence", "")),
        })
    except Exception as exc:
        return _error(f"Bond change analysis failed: {exc}")


def _infer_key_bond_type(formed_bonds: list) -> str:
    """Derive a human-readable key bond type from bond change data."""
    if not formed_bonds:
        return "unknown"
    bond_str = " ".join(str(b) for b in formed_bonds).upper()
    if "C-C" in bond_str or "C–C" in bond_str:
        return "C-C"
    if "C-N" in bond_str or "C–N" in bond_str:
        return "C-N"
    if "C-O" in bond_str or "C–O" in bond_str:
        return "C-O"
    if "C-S" in bond_str or "C–S" in bond_str:
        return "C-S"
    if "C-B" in bond_str or "C–B" in bond_str:
        return "C-B"
    if "RING" in bond_str or "CYCL" in bond_str:
        return "ring_closure"
    return str(formed_bonds[0]) if formed_bonds else "unknown"


analyze_bond_changes_tool = ToolPlugin(
    name="analyze_bond_changes",
    category="chemistry",
    description="Identify bonds broken and formed in a reaction. ALWAYS run this before searching the reaction taxonomy.",
    prerequisites=["normalize_reaction"],
    fn=_analyze_bond_changes,
)


# ---------------------------------------------------------------------------
# Tool 4: inspect_functional_groups
# ---------------------------------------------------------------------------

def _inspect_functional_groups(smiles: str) -> Dict[str, Any]:
    """Detect all functional groups in a molecule or reaction component.

    Returns groups organized by category (oxygen, nitrogen, halides, aromatic,
    leaving_groups, etc.) plus a flat list of all detected group names.
    Use this to understand what reactive sites exist in a molecule.

    Args:
        smiles: Molecule SMILES (single molecule, not a full reaction SMILES).

    Returns:
        dict with detected_groups, categories, total_groups_detected.
    """
    try:
        from chemtools.util.functional_groups import (
            get_functional_groups,
            get_group_categories,
        )
        # Strip reaction arrow if accidentally passed full reaction
        mol_smiles = smiles.split(">>")[0].split(">")[0].strip() if ">" in smiles else smiles

        detected = get_functional_groups(mol_smiles)
        categories = get_group_categories(mol_smiles)
        categories_filtered = {cat: grps for cat, grps in categories.items() if grps}

        return _success({
            "smiles": mol_smiles,
            "detected_groups": detected,
            "categories": categories_filtered,
            "total_groups_detected": len(detected),
        })
    except Exception as exc:
        return _error(f"Functional group detection failed: {exc}")


inspect_functional_groups_tool = ToolPlugin(
    name="inspect_functional_groups",
    category="chemistry",
    description="Detect all functional groups in a molecule (halides, amines, carbonyls, etc.). Use on individual reactants or products.",
    prerequisites=[],
    fn=_inspect_functional_groups,
)


# ---------------------------------------------------------------------------
# Tool 5: get_molecular_descriptors
# ---------------------------------------------------------------------------

def _get_molecular_descriptors(smiles: str) -> Dict[str, Any]:
    """Compute RDKit molecular descriptors for a molecule.

    Returns key physicochemical properties: molecular weight, logP, TPSA,
    HBA, HBD, rotatable bonds, sp3 fraction, ring counts, heavy atoms.
    Useful for drug-likeness assessment (Lipinski's Ro5), ADMET prediction,
    or characterizing reactant properties.

    Args:
        smiles: Molecule SMILES.

    Returns:
        dict with MW, logP, TPSA, HBA, HBD, RotatableBonds, FractionCSP3,
        RingCount, AromaticRings, HeavyAtomCount.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors

        mol_smiles = smiles.split(">>")[0].strip() if ">" in smiles else smiles
        mol = Chem.MolFromSmiles(mol_smiles)
        if mol is None:
            return _error(f"Could not parse SMILES: {mol_smiles!r}")

        descriptors = {
            "MW": round(Descriptors.MolWt(mol), 2),
            "ExactMW": round(Descriptors.ExactMolWt(mol), 4),
            "logP": round(Descriptors.MolLogP(mol), 2),
            "TPSA": round(Descriptors.TPSA(mol), 2),
            "HBA": rdMolDescriptors.CalcNumHBA(mol),
            "HBD": rdMolDescriptors.CalcNumHBD(mol),
            "RotatableBonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "FractionCSP3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
            "RingCount": rdMolDescriptors.CalcNumRings(mol),
            "AromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
            "HeavyAtomCount": mol.GetNumHeavyAtoms(),
            "NumStereocenters": len(rdMolDescriptors.FindPotentialStereo(mol)),
        }

        # Lipinski Ro5 check
        ro5_violations = sum([
            descriptors["MW"] > 500,
            descriptors["logP"] > 5,
            descriptors["HBA"] > 10,
            descriptors["HBD"] > 5,
        ])

        return _success({
            "smiles": mol_smiles,
            "descriptors": descriptors,
            "lipinski_ro5_violations": ro5_violations,
            "is_drug_like": ro5_violations <= 1,
        })
    except Exception as exc:
        return _error(f"Descriptor calculation failed: {exc}")


get_molecular_descriptors_tool = ToolPlugin(
    name="get_molecular_descriptors",
    category="chemistry",
    description="Compute RDKit physicochemical descriptors (MW, logP, TPSA, HBA, HBD, rings). Useful for drug-likeness and property assessment.",
    prerequisites=[],
    fn=_get_molecular_descriptors,
)


# ---------------------------------------------------------------------------
# All tools in this module
# ---------------------------------------------------------------------------

CHEMISTRY_TOOLS = [
    normalize_reaction_tool,
    detect_reaction_type_tool,
    analyze_bond_changes_tool,
    inspect_functional_groups_tool,
    get_molecular_descriptors_tool,
]

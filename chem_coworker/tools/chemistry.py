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

import re
from typing import Any, Dict

from ._helpers import _clean_rxn_smiles, _error, _success, _to_jsonable
from ._base import ToolPlugin

# Pre-warm chemtools to avoid a threading race where parallel tool calls
# simultaneously trigger the first import of chemtools._atom_mapping (a slow
# submodule) from different threads, causing "cannot import name '...' from
# partially initialized module 'chemtools._atom_mapping'".
# Importing here forces full single-threaded initialization at startup.
try:
    import chemtools as _chemtools_warm  # noqa: F401
except Exception:
    pass


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
        try:
            from chem_coworker.tool_runtime import get_current_tool_runtime_context
            _rtx = get_current_tool_runtime_context()
        except Exception:
            _rtx = None
        if _rtx is not None and hasattr(_rtx, "normalize_reaction"):
            cached = _rtx.normalize_reaction(smiles)
            if isinstance(cached, dict):
                return cached

        from chemtools.featurizers.analysis.smiles import normalize, normalize_reaction

        smiles = _clean_rxn_smiles(smiles)
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
        try:
            from chem_coworker.tool_runtime import get_current_tool_runtime_context
            _rtx = get_current_tool_runtime_context()
        except Exception:
            _rtx = None
        if _rtx is not None and hasattr(_rtx, "detect_reaction_type"):
            cached = _rtx.detect_reaction_type(reaction_smiles)
            if isinstance(cached, dict):
                return cached

        from chemtools.featurizers.unified import featurize_reaction
        from chemtools.taxonomy.reaction_catalog import get_reaction_type, resolve_reaction_type

        reaction_smiles = _clean_rxn_smiles(reaction_smiles)
        result = featurize_reaction(reaction_smiles)
        if not result:
            return _error("Reaction featurization returned no result")

        # Extract detection details
        detection = result.get("detection", {})
        validation = detection.get("validation", {})
        evidence = detection.get("evidence", {})

        reaction_type_raw = result.get("reaction_type") or validation.get("validated_detection")
        reaction_type_id = resolve_reaction_type(str(reaction_type_raw)) if reaction_type_raw else None
        reaction_type = reaction_type_id or reaction_type_raw
        confidence = (
            validation.get("validation_confidence")
            or result.get("confidence")
            or 0.0
        )
        # Taxonomy-canonical metadata (fallbacks preserve legacy behavior)
        rt_def = get_reaction_type(str(reaction_type)) if reaction_type else None
        family_label = (rt_def.name if rt_def else (reaction_type.replace("_", " ") if reaction_type else ""))
        taxonomy_metadata = {
            "id": getattr(rt_def, "id", reaction_type_id or reaction_type or ""),
            "name": getattr(rt_def, "name", family_label or ""),
            "category": getattr(rt_def, "category", ""),
            "aliases": list(getattr(rt_def, "aliases", []) or []),
            "has_constraints": bool(getattr(rt_def, "constraints", None)),
        } if reaction_type else {}

        # Collect matched motifs if available
        reacted_motifs = evidence.get("reacted_motifs", [])
        formed_motifs = evidence.get("formed_motifs", [])

        return _success({
            "reaction_smiles": reaction_smiles,
            "reaction_type": reaction_type,
            "reaction_type_id": reaction_type_id or reaction_type,
            "family_label": family_label,
            "confidence": float(confidence),
            "reacted_motifs": reacted_motifs,
            "formed_motifs": formed_motifs,
            "reaction_key": result.get("reaction_key"),
            "reaction_type_metadata": taxonomy_metadata,
        })
    except Exception as exc:
        return _error(f"Reaction type detection failed: {exc}")


detect_reaction_type_tool = ToolPlugin(
    name="detect_reaction_type",
    category="chemistry",
    description="Deterministic taxonomy classification of a reaction (Suzuki, Buchwald-Hartwig, etc.). Fast, no atom mapping needed.",
    prerequisites=[],
    fn=_detect_reaction_type,
    provides=["reaction_type", "reaction_type_id", "reaction_family", "reaction_type_metadata"],
    validators=[],
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
        try:
            from chem_coworker.tool_runtime import get_current_tool_runtime_context
            _rtx = get_current_tool_runtime_context()
        except Exception:
            _rtx = None
        if _rtx is not None and hasattr(_rtx, "analyze_bond_changes"):
            cached = _rtx.analyze_bond_changes(reaction_smiles)
            if isinstance(cached, dict):
                return cached

        from chemtools._atom_mapping import analyze_bond_changes_hybrid

        reaction_smiles = _clean_rxn_smiles(reaction_smiles)
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
        leaving = rec.get("leaving_groups", [])

        # Derive key bond type from leaving groups (more reliable than index pairs)
        key_bond = _infer_key_bond_type(broken, leaving)

        return _success({
            "reaction_smiles": reaction_smiles,
            "bonds_broken": _to_jsonable(broken),
            "bonds_formed": _to_jsonable(formed),
            "key_bond_type": key_bond,
            "leaving_groups": _to_jsonable(leaving),
            "mapping_confidence": rec.get("confidence", rec.get("mapping_confidence", "")),
        })
    except Exception as exc:
        return _error(f"Bond change analysis failed: {exc}")


def _infer_key_bond_type(broken_bonds: list, leaving_groups: list) -> str:
    """Derive a human-readable key bond type from leaving groups.

    Bond change data uses atom indices for formed bonds (no atom type info),
    but leaving groups include the departing atom symbol — use those instead.
    """
    # Collect leaving atom symbols from both broken_bonds and leaving_groups
    leaving_syms: set[str] = set()
    for entry in leaving_groups:
        # entry: [atom_idx, symbol, bond_type]
        if isinstance(entry, (list, tuple)) and len(entry) >= 2:
            sym = str(entry[1]).upper().split()[0]  # "B (leaving group)" → "B"
            leaving_syms.add(sym)
    # Also scan broken_bonds for labelled entries
    for entry in broken_bonds:
        if isinstance(entry, (list, tuple)) and len(entry) >= 2:
            second = str(entry[1]).upper()
            for sym in ("BR", "CL", "I ", "OTF", "OMS", "OTS", " N", " O", " S", " B"):
                if sym.strip() in second:
                    leaving_syms.add(sym.strip())

    # Map leaving group pattern → key bond formed
    # Suzuki: B + Br/Cl/I/OTf → C-C
    if "B" in leaving_syms and leaving_syms & {"BR", "CL", "I", "OTF", "OMS"}:
        return "C-C (Suzuki-type)"
    # Buchwald / C-N: Br/Cl/I leaves, N participates
    if "N" in leaving_syms or (leaving_syms & {"BR", "CL", "I"} and "B" not in leaving_syms):
        # crude heuristic — check if N is in leaving (amination) or just Br/Cl
        if "N" in leaving_syms:
            return "C-N"
        return "C-C or C-heteroatom"
    if "O" in leaving_syms:
        return "C-O"
    if "S" in leaving_syms:
        return "C-S"
    # Fallback
    return "C-C" if leaving_syms else "unknown"


analyze_bond_changes_tool = ToolPlugin(
    name="analyze_bond_changes",
    category="chemistry",
    description="Identify bonds broken and formed in a reaction. ALWAYS run this before searching the reaction taxonomy.",
    prerequisites=["normalize_reaction"],
    fn=_analyze_bond_changes,
    provides=["bonds_broken", "bonds_formed", "key_bond_type", "mapping_confidence"],
    validators=[],
)


def _coerce_confidence_number(value: Any) -> float | None:
    """Best-effort conversion of confidence-like values to float."""
    if value is None:
        return None
    if isinstance(value, (int, float)):
        try:
            return float(value)
        except Exception:
            return None
    text = str(value).strip()
    if not text:
        return None
    if text.endswith("%"):
        try:
            return float(text[:-1]) / 100.0
        except Exception:
            return None
    m = re.search(r"([01](?:\.\d+)?|\d+\.\d+)", text)
    if not m:
        return None
    try:
        num = float(m.group(1))
    except Exception:
        return None
    if num > 1.0 and num <= 100.0:
        return num / 100.0
    return num


def _validate_detect_reaction_type(result: dict) -> object:
    """Surface chemistry-specific caveats for deterministic reaction typing."""
    if not isinstance(result, dict) or not result.get("success"):
        return None

    rt = str(result.get("reaction_type") or result.get("reaction_type_id") or "").strip()
    conf = _coerce_confidence_number(result.get("confidence"))
    meta = result.get("reaction_type_metadata") or {}

    if not rt:
        return (
            "detect_reaction_type did not return a reaction family. "
            "State uncertainty and avoid family-specific claims or conditions."
        )
    if rt.lower() == "unknown":
        return (
            "Reaction type classified as Unknown. "
            "Use bond changes/functional groups and state uncertainty explicitly."
        )
    if conf is not None and conf < 0.40:
        return (
            f"Low reaction-type confidence ({conf:.2f}) for {rt}. "
            "Treat family-specific conclusions as tentative."
        )
    if conf is not None and conf < 0.60 and not meta.get("category"):
        return (
            f"Reaction type {rt} lacks clear taxonomy metadata and confidence is moderate/low ({conf:.2f}). "
            "Cross-check with bond-change evidence before committing."
        )
    return None


def _validate_analyze_bond_changes(result: dict) -> object:
    """Surface chemistry caveats from bond-change/mapping quality."""
    if not isinstance(result, dict) or not result.get("success"):
        return None

    formed = result.get("bonds_formed") or []
    key_bond = str(result.get("key_bond_type") or "").strip().lower()
    map_conf = _coerce_confidence_number(result.get("mapping_confidence"))

    if not formed:
        return (
            "Bond-change analysis returned no formed bonds. "
            "Reaction family inference may be unreliable; verify mapping and reaction SMILES."
        )
    if key_bond in ("", "unknown"):
        return (
            "Bond-change analysis could not infer a clear key bond type. "
            "Taxonomy search and mechanism claims may be ambiguous."
        )
    if map_conf is not None and map_conf < 0.50:
        return (
            f"Low bond-mapping confidence ({map_conf:.2f}). "
            "Treat bond-change-derived conclusions as tentative."
        )
    return None


def _project_detect_reaction_type(result: dict) -> Dict[str, Any]:
    """Project reaction typing result into ChemResponse.structured fields."""
    if not isinstance(result, dict) or not result.get("success"):
        return {}
    out: Dict[str, Any] = {
        "reaction_type": result.get("reaction_type_id") or result.get("reaction_type"),
        "reaction_family": result.get("family_label"),
    }
    if result.get("reaction_type_metadata"):
        out["reaction_type_metadata"] = result.get("reaction_type_metadata")
    return {k: v for k, v in out.items() if v is not None}


def _project_analyze_bond_changes(result: dict) -> Dict[str, Any]:
    """Project bond-change analysis into structured output fields."""
    if not isinstance(result, dict) or not result.get("success"):
        return {}
    return {
        "bonds_formed": result.get("bonds_formed", []),
        "bonds_broken": result.get("bonds_broken", []),
        "key_bond_type": result.get("key_bond_type"),
    }


def _project_get_molecular_descriptors(result: dict) -> Dict[str, Any]:
    """Project descriptor tool outputs into structured fields."""
    if not isinstance(result, dict) or not result.get("success"):
        return {}
    return {
        "descriptors": result.get("descriptors", {}),
        "is_drug_like": result.get("is_drug_like"),
    }


# Register chemistry validators after function definitions to keep tool declarations readable.
detect_reaction_type_tool.validators = [_validate_detect_reaction_type]
analyze_bond_changes_tool.validators = [_validate_analyze_bond_changes]
detect_reaction_type_tool.structured_projection = _project_detect_reaction_type
analyze_bond_changes_tool.structured_projection = _project_analyze_bond_changes


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
        try:
            from chem_coworker.tool_runtime import get_current_tool_runtime_context
            _rtx = get_current_tool_runtime_context()
        except Exception:
            _rtx = None
        from chemtools.util.functional_groups import (
            get_functional_groups,
            get_group_categories,
        )
        # Strip reaction arrow if accidentally passed full reaction
        mol_smiles = smiles.split(">>")[0].split(">")[0].strip() if ">" in smiles else smiles
        if _rtx is not None and hasattr(_rtx, "get_cached_molecule_result"):
            cached = _rtx.get_cached_molecule_result("functional_groups", mol_smiles)
            if isinstance(cached, dict):
                return cached

        detected = get_functional_groups(mol_smiles)
        categories = get_group_categories(mol_smiles)
        categories_filtered = {cat: grps for cat, grps in categories.items() if grps}
        result = _success({
            "smiles": mol_smiles,
            "detected_groups": detected,
            "categories": categories_filtered,
            "total_groups_detected": len(detected),
        })
        if _rtx is not None and hasattr(_rtx, "set_cached_molecule_result"):
            _rtx.set_cached_molecule_result("functional_groups", mol_smiles, result)
        return result
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
        try:
            from chem_coworker.tool_runtime import get_current_tool_runtime_context
            _rtx = get_current_tool_runtime_context()
        except Exception:
            _rtx = None
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors

        mol_smiles = smiles.split(">>")[0].strip() if ">" in smiles else smiles
        if _rtx is not None and hasattr(_rtx, "get_cached_molecule_result"):
            cached = _rtx.get_cached_molecule_result("descriptors", mol_smiles)
            if isinstance(cached, dict):
                return cached
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
            "NumStereocenters": (
                rdMolDescriptors.CalcNumAtomStereoCenters(mol)
                if hasattr(rdMolDescriptors, "CalcNumAtomStereoCenters")
                else len(rdMolDescriptors.FindPotentialStereo(mol))
                if hasattr(rdMolDescriptors, "FindPotentialStereo")
                else len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
            ),
        }

        # Lipinski Ro5 check
        ro5_violations = sum([
            descriptors["MW"] > 500,
            descriptors["logP"] > 5,
            descriptors["HBA"] > 10,
            descriptors["HBD"] > 5,
        ])

        result = _success({
            "smiles": mol_smiles,
            "descriptors": descriptors,
            "lipinski_ro5_violations": ro5_violations,
            "is_drug_like": ro5_violations <= 1,
        })
        if _rtx is not None and hasattr(_rtx, "set_cached_molecule_result"):
            _rtx.set_cached_molecule_result("descriptors", mol_smiles, result)
        return result
    except Exception as exc:
        return _error(f"Descriptor calculation failed: {exc}")


get_molecular_descriptors_tool = ToolPlugin(
    name="get_molecular_descriptors",
    category="chemistry",
    description="Compute RDKit physicochemical descriptors (MW, logP, TPSA, HBA, HBD, rings). Useful for drug-likeness and property assessment.",
    prerequisites=[],
    fn=_get_molecular_descriptors,
)

get_molecular_descriptors_tool.structured_projection = _project_get_molecular_descriptors


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

"""
LangChain tool wrappers for the ReactionReasoningAgent.

These wrap existing chemtools functions into focused tools that the reasoning
agent uses to inspect molecules, analyze reactivity, and map to taxonomy.

Tool categories:
  - Molecular inspection: functional groups, descriptors
  - Electronic/steric analysis: EWG/EDG scoring, steric bulk
  - Bond change analysis: identify broken/formed bonds
  - Mechanism hypothesis testing: SNAr feasibility
  - Structural comparison: reactant vs product diff
  - Taxonomy search: reaction types, motif labels, hierarchy
"""

from __future__ import annotations

import json
import logging
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional

from langchain_core.tools import tool
from pydantic import BaseModel, Field

logger = logging.getLogger(__name__)


# ============================================================================
# Helpers (reuse pattern from chemtools_wrapper.py)
# ============================================================================

def _to_jsonable(value: Any) -> Any:
    from dataclasses import asdict, is_dataclass
    if is_dataclass(value):
        return asdict(value)
    if isinstance(value, dict):
        return {k: _to_jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_to_jsonable(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    to_dict = getattr(value, "to_dict", None)
    if callable(to_dict):
        return _to_jsonable(to_dict())
    to_list = getattr(value, "tolist", None)
    if callable(to_list):
        return _to_jsonable(to_list())
    return value


def _success(data: Any) -> Dict[str, Any]:
    payload: Dict[str, Any] = {"success": True}
    if isinstance(data, dict):
        payload.update(_to_jsonable(data))
    else:
        payload["data"] = _to_jsonable(data)
    return payload


def _error(message: str, extra: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    payload: Dict[str, Any] = {"success": False, "error": message}
    if extra:
        payload.update(_to_jsonable(extra))
    return payload


# ============================================================================
# Input schemas
# ============================================================================

class SmilesInput(BaseModel):
    """Single molecule SMILES."""
    smiles: str = Field(..., description="SMILES string of a molecule.")


class ReactionSmilesInput(BaseModel):
    """Reaction SMILES (reactants>>products)."""
    reaction_smiles: str = Field(
        ..., description="Reaction SMILES in 'reactants>>products' format."
    )


class SearchQueryInput(BaseModel):
    """Free-text search query."""
    query: str = Field(
        ..., description="Search keywords (e.g. 'C-C coupling Pd', 'amide formation')."
    )


class MotifSearchInput(BaseModel):
    """Search for motif labels by scaffold and/or substituent."""
    scaffold: str = Field(
        default="",
        description="Scaffold part (e.g. 'Ar', 'Alkyl', 'Alkenyl'). Empty to search all.",
    )
    substituent: str = Field(
        default="",
        description="Substituent part (e.g. 'Br', 'B(OH)2', 'NH2'). Empty to search all.",
    )


class MotifHierarchyInput(BaseModel):
    """Get sub-motifs for a broad motif."""
    broad_motif: str = Field(
        ..., description="Broad motif label (e.g. 'Alkyl-Br', 'Ar-OH')."
    )


# ============================================================================
# Tool 1: Functional group detection
# ============================================================================

@tool(args_schema=SmilesInput)
def inspect_functional_groups(smiles: str) -> Dict[str, Any]:
    """Detect all functional groups in a molecule, categorized by type.

    Returns groups organized by category (oxygen, nitrogen, halides, aromatic,
    unsaturated, protecting_groups, leaving_groups, etc.) plus a flat list of
    all detected group names. Use this to understand what reactive sites exist.
    """
    try:
        from chemtools.util.functional_groups import (
            get_functional_groups,
            get_group_categories,
        )
        detected = get_functional_groups(smiles)
        categories = get_group_categories(smiles)
        # Filter to only non-empty categories
        categories_filtered = {
            cat: groups for cat, groups in categories.items() if groups
        }
        return _success({
            "smiles": smiles,
            "detected_groups": detected,
            "categories": categories_filtered,
            "total_groups_detected": len(detected),
        })
    except Exception as exc:
        return _error(f"Functional group detection failed: {exc}")


# ============================================================================
# Tool 2: Electronic analysis
# ============================================================================

@tool(args_schema=SmilesInput)
def analyze_electronics(smiles: str) -> Dict[str, Any]:
    """Analyze electronic effects (EWG/EDG) around reactive centers in a molecule.

    Uses the unified featurizer to get motif-level electronic scoring (0-10 scale).
    Score < 5 = electron-poor (EWG-activated), score > 5 = electron-rich (EDG).
    This helps determine: SNAr feasibility, oxidative addition rate, nucleophilicity.
    """
    try:
        from chemtools.featurizers.unified import featurize_molecule
        result = featurize_molecule(smiles)
        if not isinstance(result, dict):
            return _error("featurize_molecule returned non-dict")

        # Extract electronic data from motifs
        electronics_data = []
        reactants = result.get("reactants", [result]) if "reactants" not in result else result["reactants"]
        # featurize_molecule returns a single molecule record
        motifs = result.get("motifs", [])
        for motif in motifs:
            if not isinstance(motif, dict):
                continue
            elec = motif.get("electronics", {})
            if elec:
                electronics_data.append({
                    "compound_id": motif.get("compound_id", ""),
                    "electronics": _to_jsonable(elec),
                })

        return _success({
            "smiles": smiles,
            "motif_electronics": electronics_data,
            "num_motifs": len(motifs),
            "motif_labels": [m.get("compound_id", "") for m in motifs if isinstance(m, dict)],
        })
    except Exception as exc:
        return _error(f"Electronic analysis failed: {exc}")


# ============================================================================
# Tool 3: Steric analysis
# ============================================================================

@tool(args_schema=SmilesInput)
def analyze_steric(smiles: str) -> Dict[str, Any]:
    """Analyze steric environment around reactive centers in a molecule.

    Returns steric scores (0-10 scale) for each motif. High score = more hindered.
    Includes ortho-substituent info for aryl motifs and alpha-branching for alkyl.
    This helps determine: SN2 feasibility, steric selectivity, catalyst choice.
    """
    try:
        from chemtools.featurizers.unified import featurize_molecule
        result = featurize_molecule(smiles)
        if not isinstance(result, dict):
            return _error("featurize_molecule returned non-dict")

        steric_data = []
        motifs = result.get("motifs", [])
        for motif in motifs:
            if not isinstance(motif, dict):
                continue
            steric = motif.get("steric", {})
            if steric:
                steric_data.append({
                    "compound_id": motif.get("compound_id", ""),
                    "steric": _to_jsonable(steric),
                })

        return _success({
            "smiles": smiles,
            "motif_sterics": steric_data,
            "num_motifs": len(motifs),
            "motif_labels": [m.get("compound_id", "") for m in motifs if isinstance(m, dict)],
        })
    except Exception as exc:
        return _error(f"Steric analysis failed: {exc}")


# ============================================================================
# Tool 4: SNAr feasibility check
# ============================================================================

@tool(args_schema=ReactionSmilesInput)
def check_snar_feasibility(reaction_smiles: str) -> Dict[str, Any]:
    """Check if SNAr (nucleophilic aromatic substitution) is feasible for this reaction.

    Analyzes the aryl electrophile's electronic activation. Returns a feasibility
    score and classification. Only meaningful for reactions involving aryl halides/
    pseudohalides with nucleophilic displacement.
    """
    try:
        from chemtools.featurizers.unified import featurize_reaction
        payload = featurize_reaction(reaction_smiles)
        if not isinstance(payload, dict):
            return _error("featurize_reaction returned non-dict")

        reaction = payload.get("reaction") if isinstance(payload.get("reaction"), dict) else payload

        from chemtools.featurizers.analysis.feasibility import analyze_snar_feasibility
        feasibility = analyze_snar_feasibility(reaction)
        return _success({
            "reaction_smiles": reaction_smiles,
            "snar_feasibility": _to_jsonable(feasibility),
        })
    except Exception as exc:
        return _error(f"SNAr feasibility check failed: {exc}")


# ============================================================================
# Tool 5: Bond change analysis
# ============================================================================

@tool(args_schema=ReactionSmilesInput)
def analyze_bond_changes(reaction_smiles: str) -> Dict[str, Any]:
    """Identify which bonds broke and formed in the reaction.

    Uses multiple mapping methods (RXNMapper, local environment, MCS) for best
    accuracy. Returns: bonds broken, bonds formed, mapping confidence, and the
    recommended result. This is essential for determining the reaction mechanism.
    """
    try:
        from chemtools._atom_mapping import analyze_bond_changes_hybrid
        result = analyze_bond_changes_hybrid(reaction_smiles)

        # Extract the recommended result for the agent
        recommended = result.get("recommended_result", {})
        summary = {
            "reaction_smiles": reaction_smiles,
            "method": result.get("method", "unknown"),
            "combined_confidence": result.get("combined_confidence", 0.0),
        }
        if recommended:
            summary["bonds_broken"] = recommended.get("bonds_broken", [])
            summary["bonds_formed"] = recommended.get("bonds_formed", [])
            summary["bond_changes_summary"] = recommended.get("bond_changes_summary", "")
            summary["changed_atoms"] = recommended.get("changed_atoms", [])
            summary["mapping_confidence"] = recommended.get("mapping_confidence", 0.0)

        return _success(summary)
    except Exception as exc:
        return _error(f"Bond change analysis failed: {exc}")


# ============================================================================
# Tool 6: Molecular descriptors
# ============================================================================

@tool(args_schema=SmilesInput)
def compute_molecular_descriptors(smiles: str) -> Dict[str, Any]:
    """Compute RDKit molecular descriptors: MW, logP, TPSA, HBA, HBD, rotatable bonds,
    sp3 fraction, ring counts, and aromatic ring count.

    Useful for understanding molecular complexity, drug-likeness, and reactivity.
    """
    try:
        from chemtools.util.rdkit_helpers import parse_smiles
        mol = parse_smiles(smiles)
        if mol is None:
            return _error(f"Could not parse SMILES: {smiles}")

        from chemtools.featurizers.calculable import _get_property_values
        props = _get_property_values(mol)

        # Add extra descriptors
        try:
            from rdkit.Chem import Descriptors, rdMolDescriptors
            props["sp3_fraction"] = float(Descriptors.FractionCSP3(mol))
            props["num_rings"] = int(rdMolDescriptors.CalcNumRings(mol))
            props["num_aromatic_rings"] = int(rdMolDescriptors.CalcNumAromaticRings(mol))
            props["num_heavy_atoms"] = int(mol.GetNumHeavyAtoms())
        except Exception:
            pass

        return _success({
            "smiles": smiles,
            "descriptors": props,
        })
    except Exception as exc:
        return _error(f"Molecular descriptor computation failed: {exc}")


# ============================================================================
# Tool 7: Reactant vs Product comparison
# ============================================================================

@tool(args_schema=ReactionSmilesInput)
def compare_reactant_product(reaction_smiles: str) -> Dict[str, Any]:
    """Compare reactant and product structures to identify what changed.

    Analyzes: atom count changes, functional groups that appeared/disappeared,
    ring count changes, and overall structural transformation. Use this as a
    first step to understand the net transformation.
    """
    try:
        if ">>" not in reaction_smiles:
            return _error("Invalid reaction SMILES: missing '>>'")

        reactants_str, products_str = reaction_smiles.split(">>", 1)

        from chemtools.util.functional_groups import get_functional_groups
        from chemtools.util.rdkit_helpers import parse_smiles

        # Analyze each reactant
        reactant_smiles_list = [s.strip() for s in reactants_str.split(".") if s.strip()]
        product_smiles_list = [s.strip() for s in products_str.split(".") if s.strip()]

        # Collect FGs for all reactants and products
        reactant_fgs: Dict[str, List[str]] = {}
        product_fgs: Dict[str, List[str]] = {}

        all_reactant_groups = set()
        all_product_groups = set()

        for smi in reactant_smiles_list:
            groups = get_functional_groups(smi)
            reactant_fgs[smi] = groups
            all_reactant_groups.update(groups)

        for smi in product_smiles_list:
            groups = get_functional_groups(smi)
            product_fgs[smi] = groups
            all_product_groups.update(groups)

        # Diff
        groups_removed = sorted(all_reactant_groups - all_product_groups)
        groups_added = sorted(all_product_groups - all_reactant_groups)
        groups_preserved = sorted(all_reactant_groups & all_product_groups)

        # Atom counts
        def count_atoms(smiles_list):
            counts = {}
            for smi in smiles_list:
                mol = parse_smiles(smi)
                if mol is None:
                    continue
                for atom in mol.GetAtoms():
                    sym = atom.GetSymbol()
                    counts[sym] = counts.get(sym, 0) + 1
            return counts

        reactant_atoms = count_atoms(reactant_smiles_list)
        product_atoms = count_atoms(product_smiles_list)

        return _success({
            "reaction_smiles": reaction_smiles,
            "reactants": reactant_smiles_list,
            "products": product_smiles_list,
            "reactant_functional_groups": reactant_fgs,
            "product_functional_groups": product_fgs,
            "groups_removed": groups_removed,
            "groups_added": groups_added,
            "groups_preserved": groups_preserved,
            "reactant_atom_counts": reactant_atoms,
            "product_atom_counts": product_atoms,
        })
    except Exception as exc:
        return _error(f"Reactant-product comparison failed: {exc}")


# ============================================================================
# Tool 8: Search reaction types
# ============================================================================

@lru_cache(maxsize=1)
def _load_reaction_types_data() -> List[Dict[str, Any]]:
    """Load full reaction types data from taxonomy."""
    try:
        from chemtools.taxonomy.loader import load_reaction_types_list
        return load_reaction_types_list()
    except Exception as e:
        logger.warning(f"Could not load reaction types: {e}")
        return []


@tool(args_schema=SearchQueryInput)
def search_reaction_types(query: str) -> Dict[str, Any]:
    """Search the taxonomy for matching reaction types by keyword.

    Searches across: reaction type IDs, aliases, descriptions, catalysts,
    bond constraints, and reactant/product types. Returns matching reaction
    type entries with their full metadata.

    Example queries: "C-C coupling Pd", "amide formation", "SNAr", "borylation"
    """
    try:
        rxn_types = _load_reaction_types_data()
        if not rxn_types:
            return _error("Reaction types taxonomy not available")

        keywords = [kw.lower().strip() for kw in query.split() if kw.strip()]
        matches = []

        for entry in rxn_types:
            if not isinstance(entry, dict):
                continue

            # Build searchable text
            searchable_parts = []
            searchable_parts.append(entry.get("id", ""))
            searchable_parts.extend(entry.get("aliases", []))
            searchable_parts.append(entry.get("description", ""))
            searchable_parts.extend(entry.get("catalysts", []))

            constraints = entry.get("constraints", {})
            if isinstance(constraints, dict):
                for v in constraints.values():
                    if isinstance(v, list):
                        searchable_parts.extend(str(x) for x in v)
                    else:
                        searchable_parts.append(str(v))

            reactants = entry.get("reactants", {})
            if isinstance(reactants, dict):
                searchable_parts.extend(str(v) for v in reactants.values())

            products = entry.get("products", {})
            if isinstance(products, dict):
                for v in products.values():
                    if isinstance(v, list):
                        searchable_parts.extend(str(x) for x in v)
                    else:
                        searchable_parts.append(str(v))

            searchable = " ".join(searchable_parts).lower()

            # Score: count how many keywords match
            score = sum(1 for kw in keywords if kw in searchable)
            if score > 0:
                matches.append({
                    "id": entry.get("id"),
                    "aliases": entry.get("aliases", []),
                    "description": entry.get("description", ""),
                    "catalysts": entry.get("catalysts", []),
                    "constraints": entry.get("constraints", {}),
                    "reactants": entry.get("reactants", {}),
                    "products": entry.get("products", {}),
                    "match_score": score,
                })

        # Sort by match score descending
        matches.sort(key=lambda x: x["match_score"], reverse=True)

        return _success({
            "query": query,
            "num_matches": len(matches),
            "matches": matches[:10],  # Top 10
        })
    except Exception as exc:
        return _error(f"Reaction type search failed: {exc}")


# ============================================================================
# Tool 9: Search motif labels
# ============================================================================

@lru_cache(maxsize=1)
def _load_motif_data() -> Dict[str, Any]:
    """Load organic compounds and motif scope index."""
    data = {"compounds": [], "scope_map": {}, "scaffold_parent_map": {}}
    try:
        from chemtools.taxonomy.loader import load_organic_compounds
        oc = load_organic_compounds()
        data["compounds"] = oc.get("compounds", [])
    except Exception as e:
        logger.warning(f"Could not load organic compounds: {e}")

    try:
        scope_path = Path(__file__).parent.parent / "chemtools" / "taxonomy" / "data" / "motif_scope_index.v1.json"
        if scope_path.exists():
            with open(scope_path, encoding="utf-8") as f:
                scope_data = json.load(f)
            data["scope_map"] = scope_data.get("scope_map", {})
            data["scaffold_parent_map"] = scope_data.get("scaffold_parent_map", {})
    except Exception as e:
        logger.warning(f"Could not load motif scope index: {e}")

    return data


@tool(args_schema=MotifSearchInput)
def search_motifs(scaffold: str = "", substituent: str = "") -> Dict[str, Any]:
    """Search the taxonomy for matching motif labels (scaffold-substituent pairs).

    Motif labels follow the format 'Scaffold-Substituent' (e.g. 'Ar-Br', 'Alkyl-NH2').
    Search by scaffold (e.g. 'Ar'), substituent (e.g. 'Br'), or both.

    Returns matching motif labels with their scaffold/substituent breakdown.
    """
    try:
        motif_data = _load_motif_data()
        compounds = motif_data.get("compounds", [])
        if not compounds:
            return _error("Motif taxonomy not available")

        scaffold_lower = scaffold.lower().strip()
        sub_lower = substituent.lower().strip()
        # Normalize: substituent in taxonomy uses "-" prefix (e.g. "-Br")
        if sub_lower and not sub_lower.startswith("-"):
            sub_lower_search = sub_lower
        else:
            sub_lower_search = sub_lower.lstrip("-")

        matches = []
        for entry in compounds:
            if not isinstance(entry, dict):
                continue
            a = entry.get("A", "")
            b = entry.get("B", "")
            label = f"{a}{b}"

            # Match scaffold
            if scaffold_lower and scaffold_lower not in a.lower():
                continue
            # Match substituent
            if sub_lower_search and sub_lower_search not in b.lower():
                continue

            matches.append({
                "label": label,
                "scaffold": a,
                "substituent": b,
            })

        return _success({
            "scaffold_query": scaffold,
            "substituent_query": substituent,
            "num_matches": len(matches),
            "matches": matches[:30],  # Cap at 30
        })
    except Exception as exc:
        return _error(f"Motif search failed: {exc}")


# ============================================================================
# Tool 10: Motif hierarchy
# ============================================================================

@tool(args_schema=MotifHierarchyInput)
def get_motif_hierarchy(broad_motif: str) -> Dict[str, Any]:
    """Get specific sub-motifs for a broad motif from the scope index.

    E.g. 'Alkyl-Br' expands to ['Propargyl-Br', 'R2CH-Br', 'R3C-Br', 'RCH2-Br'].
    Also shows scaffold parent relationships (e.g. HeteroAr -> Ar, RCH2 -> Alkyl).

    Use this to determine which specific motif label to assign when you know
    the broad category.
    """
    try:
        motif_data = _load_motif_data()
        scope_map = motif_data.get("scope_map", {})
        parent_map = motif_data.get("scaffold_parent_map", {})

        result = {
            "broad_motif": broad_motif,
        }

        # Direct lookup in scope_map
        if broad_motif in scope_map:
            result["specific_motifs"] = scope_map[broad_motif]
        else:
            # Try case-insensitive
            for key, val in scope_map.items():
                if key.lower() == broad_motif.lower():
                    result["specific_motifs"] = val
                    result["canonical_key"] = key
                    break
            else:
                result["specific_motifs"] = []
                result["note"] = f"'{broad_motif}' not found in scope_map. It may already be a specific motif."

        # Show relevant parent relationships
        relevant_parents = {}
        for child, parents in parent_map.items():
            if broad_motif.split("-")[0] in parents or child in broad_motif:
                relevant_parents[child] = parents

        if relevant_parents:
            result["scaffold_parents"] = relevant_parents

        return _success(result)
    except Exception as exc:
        return _error(f"Motif hierarchy lookup failed: {exc}")


# ============================================================================
# Exported tool list
# ============================================================================

REASONING_TOOLS = [
    inspect_functional_groups,
    analyze_electronics,
    analyze_steric,
    check_snar_feasibility,
    analyze_bond_changes,
    compute_molecular_descriptors,
    compare_reactant_product,
    search_reaction_types,
    search_motifs,
    get_motif_hierarchy,
]

__all__ = [
    "REASONING_TOOLS",
    "inspect_functional_groups",
    "analyze_electronics",
    "analyze_steric",
    "check_snar_feasibility",
    "analyze_bond_changes",
    "compute_molecular_descriptors",
    "compare_reactant_product",
    "search_reaction_types",
    "search_motifs",
    "get_motif_hierarchy",
]

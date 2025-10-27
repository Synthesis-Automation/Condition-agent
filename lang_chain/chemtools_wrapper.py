"""
LangChain tool wrappers for ChemTools functions.

This module exposes existing chemtools functionality as LangChain tools
without modifying the original chemtools codebase.

Available Tools:
    - normalize_smiles_tool: Canonicalize SMILES strings
    - normalize_reaction_tool: Canonicalize reaction SMILES
    - detect_reaction_family_tool: Detect reaction family/type
    - recommend_conditions_tool: Get ML-based condition recommendations
    - search_precedents_tool: Search for similar precedent reactions
    - find_reagent_tool: Look up reagent information from database
    - classify_reactant_tool: Classify reactant type (aryl halide, amine, etc.)
    - get_functional_groups_tool: Detect functional groups in a molecule

Usage:
    from lang_chain.chemtools_wrapper import CHEMTOOLS_TOOLS
    from langgraph.prebuilt import create_react_agent
    
    agent = create_react_agent(llm, CHEMTOOLS_TOOLS)
"""

from typing import Dict, Any, List, Optional
import json
from langchain_core.tools import tool

# Import chemtools functions
from chemtools.smiles import normalize, normalize_reaction
from chemtools.router import detect_family_from_reaction
from chemtools.recommend.modules import recommend_from_reaction
from chemtools.precedent import knn as precedent_knn
from chemtools.reagent import find_reagent, classify_reactant_smiles
from chemtools.util.functional_groups import detect_all as detect_functional_groups


# ============================================================================
# SMILES Normalization Tools
# ============================================================================

@tool
def normalize_smiles_tool(smiles: str) -> str:
    """
    Normalize (canonicalize) a SMILES string.
    
    This standardizes molecular representations for consistent comparisons.
    Returns a dict with normalization details.
    
    Args:
        smiles: SMILES string to normalize (e.g., "CCO", "c1ccccc1")
    
    Returns:
        JSON string with normalized SMILES and metadata
    
    Example:
        Input: "c1ccccc1"
        Output: {"smiles_norm": "c1ccccc1", ...}
    """
    try:
        result = normalize(smiles)
        if isinstance(result, dict):
            # Return the normalized SMILES if it's a dict
            return json.dumps(result, indent=2)
        elif result:
            return json.dumps({"smiles_norm": result}, indent=2)
        else:
            return json.dumps({"error": f"Invalid SMILES '{smiles}'"})
    except Exception as e:
        return json.dumps({"error": str(e)})


@tool
def normalize_reaction_tool(reaction_smiles: str) -> str:
    """
    Normalize a reaction SMILES string.
    
    Standardizes both reactants and products, sorting components for consistency.
    Accepts formats: "reactants>>products" or "reactants>agents>products"
    
    Args:
        reaction_smiles: Reaction SMILES (e.g., "CCBr.CCO>>CCOCC")
    
    Returns:
        Canonicalized reaction SMILES string
    
    Example:
        Input: "c1ccccc1Br.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
        Output: Normalized reaction SMILES
    """
    try:
        result = normalize_reaction(reaction_smiles)
        return result if result else f"Error: Invalid reaction SMILES '{reaction_smiles}'"
    except Exception as e:
        return f"Error normalizing reaction: {str(e)}"


# ============================================================================
# Reaction Analysis Tools
# ============================================================================

@tool
def detect_reaction_family_tool(reaction_smiles: str) -> str:
    """
    Detect the reaction family/type from a reaction SMILES.
    
    Analyzes reactants and products to identify the reaction type
    (e.g., Suzuki, Buchwald_CN, Ullmann_CN, etc.)
    
    Args:
        reaction_smiles: Reaction SMILES string
    
    Returns:
        JSON string with detected family and metadata
    
    Example:
        Input: "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
        Output: {"family": "Suzuki", "confidence": "high"}
    """
    try:
        result = detect_family_from_reaction(reaction_smiles)
        return json.dumps(result, indent=2)
    except Exception as e:
        return json.dumps({"error": str(e), "family": "Unknown"})


@tool
def classify_reactant_tool(smiles: str) -> str:
    """
    Classify a reactant molecule type.
    
    Identifies reactant categories like aryl_halide, amine, boronic_acid,
    alkyne, etc. based on structural patterns.
    
    Args:
        smiles: SMILES string of the reactant
    
    Returns:
        JSON string with classification results
    
    Example:
        Input: "Brc1ccccc1"
        Output: {"category": "aryl_halide", "matches": [...]}
    """
    try:
        result = classify_reactant_smiles(smiles)
        return json.dumps(result, indent=2, default=str)
    except Exception as e:
        return json.dumps({"error": str(e), "category": "unknown"})


@tool
def get_functional_groups_tool(smiles: str) -> str:
    """
    Detect functional groups in a molecule.
    
    Identifies 80+ functional groups using SMARTS patterns.
    Uses RDKit when available, falls back to text patterns.
    
    Args:
        smiles: SMILES string to analyze
    
    Returns:
        JSON string with detected functional groups
    
    Example:
        Input: "CCO"
        Output: {"alcohol": true, "primary_alcohol": true}
    """
    try:
        result = detect_functional_groups(smiles)
        return json.dumps(result, indent=2)
    except Exception as e:
        return json.dumps({"error": str(e)})


# ============================================================================
# Recommendation Tools
# ============================================================================

@tool
def recommend_conditions_tool(
    reaction_smiles: str,
    k: int = 25,
    max_variants: int = 3,
    rerank_strategy: str = "rule"
) -> str:
    """
    Recommend reaction conditions using ML-based DRFP similarity search.
    
    Analyzes the reaction and finds the most similar precedent reactions,
    then recommends optimal conditions (catalyst, base, solvent, temperature, time).
    
    Args:
        reaction_smiles: Reaction SMILES string (reactants>>products)
        k: Number of precedents to retrieve (default: 25)
        max_variants: Maximum condition variants to generate (default: 3)
        rerank_strategy: Ranking strategy - "rule", "analytics", or "none" (default: "rule")
    
    Returns:
        JSON string with recommended conditions and alternatives
    
    Example:
        Input: "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
        Output: {
            "family": "Buchwald_CN",
            "recommendation": {
                "core": "Pd/XPhos",
                "base": "Cs2CO3",
                "solvent": "1,4-dioxane",
                "T_C": 100,
                "time_h": 12
            },
            "alternatives": {...},
            "reasons": [...]
        }
    """
    try:
        result = recommend_from_reaction(
            reaction=reaction_smiles,
            k=k,
            max_variants=max_variants,
            rerank_strategy=rerank_strategy
        )
        # Simplify output for LLM consumption
        simplified = {
            "family": result.get("family", "Unknown"),
            "recommendation": result.get("recommendation", {}),
            "alternatives": result.get("alternatives", {}),
            "reasons": result.get("reasons", [])[:5],  # Top 5 reasons
            "precedent_count": len(result.get("precedent_pack", {}).get("precedents", []))
        }
        return json.dumps(simplified, indent=2)
    except Exception as e:
        return json.dumps({"error": str(e), "family": "Unknown"})


@tool
def search_precedents_tool(
    reaction_smiles: str,
    k: int = 10,
    family: Optional[str] = None
) -> str:
    """
    Search for similar precedent reactions using DRFP fingerprints.
    
    Finds the k most similar reactions from the precedent database
    based on reaction fingerprint similarity.
    
    Args:
        reaction_smiles: Query reaction SMILES
        k: Number of similar reactions to find (default: 10)
        family: Optional reaction family filter (e.g., "Suzuki", "Buchwald_CN")
    
    Returns:
        JSON string with similar precedent reactions
    
    Example:
        Input: "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
        Output: {
            "precedents": [
                {
                    "reaction": "...",
                    "similarity": 0.95,
                    "conditions": {...},
                    "yield": 85
                },
                ...
            ]
        }
    """
    try:
        # Auto-detect family if not provided
        if family is None:
            family_result = detect_family_from_reaction(reaction_smiles)
            family = family_result.get("family", "Unknown")
        
        # Use relax parameter for DRFP-based search
        result = precedent_knn(
            family=family,
            features={},  # DRFP uses reaction fingerprint, not substrate features
            k=k,
            relax={"use_drfp": True, "reaction_smiles": reaction_smiles}
        )
        
        # Simplify precedent data
        precedents = result.get("precedents", [])
        simplified_precedents = []
        for p in precedents[:k]:
            simplified_precedents.append({
                "similarity": round(p.get("similarity", 0), 3),
                "conditions": {
                    "catalyst": p.get("catalyst", ""),
                    "base": p.get("base", ""),
                    "solvent": p.get("solvent", ""),
                    "temperature_C": p.get("T_C", ""),
                    "time_h": p.get("time_h", "")
                },
                "yield": p.get("yield", ""),
                "reaction": p.get("rxn_smiles", "")[:100]  # Truncate long SMILES
            })
        
        return json.dumps({
            "family": family,
            "precedent_count": len(simplified_precedents),
            "precedents": simplified_precedents
        }, indent=2)
    except Exception as e:
        return json.dumps({"error": str(e), "precedents": []})


# ============================================================================
# Reagent Database Tools
# ============================================================================

@tool
def find_reagent_tool(query: str, reagent_type: str = "base") -> str:
    """
    Look up reagent information from the reagent database.
    
    Searches by name or abbreviation to find reagent properties,
    roles, and metadata. Supports multiple reagent types.
    
    Args:
        query: Reagent name or abbreviation (e.g., "Cs2CO3", "XPhos")
        reagent_type: Type of reagent database to search 
                     (base, solvent, ligand, metal, additive, etc.)
                     Default: "base"
    
    Returns:
        JSON string with reagent information
    
    Example:
        Input: query="Cs2CO3", reagent_type="base"
        Output: {
            "name": "Cesium carbonate",
            "cas": "534-17-8",
            "role": "base",
            "smiles": "[Cs+].[Cs+].[O-]C([O-])=O"
        }
    """
    try:
        # Try common reagent types if not found
        types_to_try = [reagent_type, "base", "solvent", "ligand", "metal", "additive"]
        
        for r_type in types_to_try:
            result = find_reagent(query, r_type)
            if result:
                # Simplify the output
                simplified = {
                    "name": result.get("name", ""),
                    "cas": result.get("cas", ""),
                    "role": result.get("role", r_type),
                    "family": result.get("family", ""),
                    "smiles": result.get("smiles", ""),
                    "source": result.get("data_source", ""),
                    "reagent_type": r_type
                }
                return json.dumps(simplified, indent=2)
        
        return json.dumps({"error": f"Reagent not found: {query}"})
    except Exception as e:
        return json.dumps({"error": str(e)})


# ============================================================================
# Tool Collection
# ============================================================================

CHEMTOOLS_TOOLS = [
    # SMILES tools
    normalize_smiles_tool,
    normalize_reaction_tool,
    
    # Analysis tools
    detect_reaction_family_tool,
    classify_reactant_tool,
    get_functional_groups_tool,
    
    # Recommendation tools
    recommend_conditions_tool,
    search_precedents_tool,
    
    # Database tools
    find_reagent_tool,
]


# ============================================================================
# Helper Functions
# ============================================================================

def get_tool_descriptions() -> List[Dict[str, str]]:
    """Get descriptions of all available tools."""
    descriptions = []
    for tool_obj in CHEMTOOLS_TOOLS:
        descriptions.append({
            "name": tool_obj.name,
            "description": tool_obj.description or "",
        })
    return descriptions


def print_tool_summary():
    """Print a summary of available tools."""
    print("=" * 70)
    print("ChemTools LangChain Tools")
    print("=" * 70)
    for i, tool_obj in enumerate(CHEMTOOLS_TOOLS, 1):
        print(f"\n{i}. {tool_obj.name}")
        print(f"   {tool_obj.description.split('.')[0] if tool_obj.description else 'No description'}")
    print("\n" + "=" * 70)


if __name__ == "__main__":
    # Print tool summary when run directly
    print_tool_summary()

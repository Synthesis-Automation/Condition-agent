"""
Motif-based steric and electronic analysis using organic compound motifs.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from chemtools.util.rdkit_helpers import parse_smiles, rdkit_available

from .alkyl_steric import analyze_alkyl_steric
from .aryl_electronics import analyze_aryl_electronics
from .aryl_steric import analyze_aryl_steric
from .nearby_groups import analyze_nearby_groups
from .motif_detect import detect_motifs
from .motif_registry import build_compound_registry


def featurize_molecule(
    smiles: str,
    registry_paths: Optional[Dict[str, str | Path]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Analyze motifs, sterics, and electronics for a SMILES string.
    """
    meta = {"rdkit_available": rdkit_available(), "error": None}
    if not meta["rdkit_available"]:
        meta["error"] = "rdkit_unavailable"
        return {
            "schema_version": "v2",
            "smiles": smiles,
            "motifs": [],
            "steric": {"aryl": [], "alkyl": []},
            "electronics": {"aryl": []},
            "analyses": [],
            "meta": meta,
        }

    mol = parse_smiles(smiles)
    if mol is None:
        meta["error"] = "invalid_smiles"
        return {
            "schema_version": "v2",
            "smiles": smiles,
            "motifs": [],
            "steric": {"aryl": [], "alkyl": []},
            "electronics": {"aryl": []},
            "analyses": [],
            "meta": meta,
        }

    registry_paths = registry_paths or _default_registry_paths()
    registry = build_compound_registry(registry_paths)
    compiled = registry["compiled_compounds"]
    groups = registry["groups"]

    options = options or {}
    include_gasteiger = bool(options.get("include_gasteiger", False))
    include_ipso_group = options.get("electronics_include_ipso_group", True)
    max_hits = options.get("max_hits_per_compound")

    motifs = detect_motifs(mol, compiled, max_hits_per_compound=max_hits)
    analyses = []
    for hit in motifs:
        compound_id = hit["compound_id"]
        if compound_id.startswith(("Ar-", "AromN-")):
            steric = analyze_aryl_steric(mol, hit)
            if include_ipso_group == "both":
                electronic = [
                    analyze_aryl_electronics(
                        mol,
                        hit,
                        groups,
                        include_ipso_group=True,
                        include_gasteiger=include_gasteiger,
                    ),
                    analyze_aryl_electronics(
                        mol,
                        hit,
                        groups,
                        include_ipso_group=False,
                        include_gasteiger=include_gasteiger,
                    ),
                ]
            else:
                electronic = analyze_aryl_electronics(
                    mol,
                    hit,
                    groups,
                    include_ipso_group=bool(include_ipso_group),
                    include_gasteiger=include_gasteiger,
                )
            nearby = analyze_nearby_groups(mol, hit, motifs, groups)
            analyses.append(
                {
                    "compound_id": compound_id,
                    "center": {"ipso": hit["a_atom_idx"], "bond": hit["bond"]},
                    "steric": steric,
                    "electronic": electronic,
                    "nearby_groups": nearby,
                }
            )
        elif compound_id.startswith(("R-", "Bn-", "Allyl-")):
            steric = analyze_alkyl_steric(mol, hit)
            nearby = analyze_nearby_groups(mol, hit, motifs, groups)
            analyses.append(
                {
                    "compound_id": compound_id,
                    "center": {"alpha_c": hit["a_atom_idx"], "bond": hit["bond"]},
                    "steric": steric,
                    "nearby_groups": nearby,
                }
            )

    steric_payload = {"aryl": [], "alkyl": []}
    electronic_payload = {"aryl": []}
    nearby_payload = []
    for analysis in analyses:
        steric_entry = {
            "compound_id": analysis.get("compound_id"),
            "center": analysis.get("center"),
            "result": analysis.get("steric"),
        }
        if analysis.get("compound_id", "").startswith(("Ar-", "AromN-")):
            steric_payload["aryl"].append(steric_entry)
        else:
            steric_payload["alkyl"].append(steric_entry)
        if "electronic" in analysis:
            electronic_payload["aryl"].append(
                {
                    "compound_id": analysis.get("compound_id"),
                    "center": analysis.get("center"),
                    "result": analysis.get("electronic"),
                }
            )
        if "nearby_groups" in analysis:
            nearby_payload.append(
                {
                    "compound_id": analysis.get("compound_id"),
                    "center": analysis.get("center"),
                    "result": analysis.get("nearby_groups"),
                }
            )

    return {
        "schema_version": "v2",
        "smiles": smiles,
        "motifs": motifs,
        "steric": steric_payload,
        "electronics": electronic_payload,
        "nearby": nearby_payload,
        "analyses": analyses,
        "meta": meta,
    }


def featurize_reaction(
    reaction_smiles: str,
    registry_paths: Optional[Dict[str, str | Path]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Systematically analyze a reaction: detect type and featurize all reactants.
    """
    from chemtools.smiles import normalize_reaction

    from .reaction_detection import detect_reaction_types

    # 1. Detect reaction types
    max_hits = (options or {}).get("max_hits_per_compound")
    detection = detect_reaction_types(reaction_smiles, max_hits_per_compound=max_hits)

    # 2. Normalize and extract reactants
    normalized = normalize_reaction(reaction_smiles)
    reactants = normalized.get("reactants") or []

    # 3. Featurize each reactant
    reactant_analyses = []
    
    # Aggregated payloads for the whole reaction
    steric_payload = {"aryl": [], "alkyl": []}
    electronic_payload = {"aryl": []}
    nearby_payload = []
    all_motifs = []

    for r in reactants:
        smiles = r.get("smiles_norm") or r.get("input")
        if smiles:
            analysis = featurize_molecule(
                smiles, registry_paths=registry_paths, options=options
            )
            reactant_analyses.append({"smiles": smiles, "analysis": analysis})
            
            # Aggregate features with smiles context
            for entry in analysis.get("steric", {}).get("aryl", []):
                entry_copy = entry.copy()
                entry_copy["smiles"] = smiles
                steric_payload["aryl"].append(entry_copy)
            
            for entry in analysis.get("steric", {}).get("alkyl", []):
                entry_copy = entry.copy()
                entry_copy["smiles"] = smiles
                steric_payload["alkyl"].append(entry_copy)
                
            for entry in analysis.get("electronics", {}).get("aryl", []):
                entry_copy = entry.copy()
                entry_copy["smiles"] = smiles
                electronic_payload["aryl"].append(entry_copy)
                
            for entry in analysis.get("nearby", []):
                entry_copy = entry.copy()
                entry_copy["smiles"] = smiles
                nearby_payload.append(entry_copy)
                
            for motif in analysis.get("motifs", []):
                motif_copy = motif.copy()
                motif_copy["smiles"] = smiles
                all_motifs.append(motif_copy)

    return {
        "reaction_smiles": reaction_smiles,
        "detection": detection.to_dict(),
        "reactants": reactant_analyses,
        "motifs": all_motifs,
        "steric": steric_payload,
        "electronics": electronic_payload,
        "nearby": nearby_payload,
        "meta": {"normalized_reaction": normalized},
    }


def analyze_smiles(
    smiles: str,
    registry_paths: Optional[Dict[str, str | Path]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Compatibility wrapper for motif/steric/electronic analysis."""
    return featurize_molecule(smiles, registry_paths=registry_paths, options=options)


def _default_registry_paths() -> Dict[str, Path]:
    base = Path(__file__).resolve().parent.parent / "taxonomy" / "v2_data"
    return {
        "groups": base / "organic_groups.v1.3.json",
        "compounds": base / "organic_compounds.v1.3.json",
        "templates": base / "smarts_templates.v1.json",
    }


__all__ = ["featurize_molecule", "featurize_reaction", "analyze_smiles"]

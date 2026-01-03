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

    from rdkit import Chem
    mol = Chem.AddHs(mol)

    registry_paths = registry_paths or _default_registry_paths()
    registry = build_compound_registry(registry_paths)
    compiled = registry["compiled_compounds"]
    groups = registry["groups"]
    compound_map = registry.get("compound_map", {})

    options = options or {}
    include_gasteiger = bool(options.get("include_gasteiger", False))
    include_ipso_group = options.get("electronics_include_ipso_group", True)
    max_hits = options.get("max_hits_per_compound")
    discovery_mode = bool(options.get("discovery_mode", False))

    all_motifs = detect_motifs(
        mol,
        compiled,
        max_hits_per_compound=max_hits,
        registry=registry,
        discovery_mode=discovery_mode,
    )
    motifs = list(all_motifs)

    # Filter by target groups if provided
    target_groups = options.get("target_groups")
    if target_groups:
        if isinstance(target_groups, str):
            target_groups = [target_groups]
        
        filtered_motifs = []
        for m in motifs:
            cid = m.get("compound_id", "")
            # Match if cid matches target exactly, or ends with "-target", 
            # or ends with "target" if target already starts with "-"
            for tg in target_groups:
                if cid == tg:
                    filtered_motifs.append(m)
                    break
                if cid.endswith("-" + tg):
                    filtered_motifs.append(m)
                    break
                if tg.startswith("-") and cid.endswith(tg):
                    filtered_motifs.append(m)
                    break
        
        # Only apply filter if we actually found matches for the target groups
        if filtered_motifs:
            motifs = filtered_motifs

    # Identify background motifs (H-motifs)
    background_ids = {"Ar-H", "R-H", "Any-H", "Alkenyl-H", "Alkynyl-H"}
    
    # Determine which motifs were explicitly requested via target_groups
    requested_ids = set()
    if target_groups:
        for m in motifs:
            cid = m.get("compound_id", "")
            for tg in target_groups:
                if cid == tg or cid.endswith("-" + tg) or (tg.startswith("-") and cid.endswith(tg)):
                    requested_ids.add(cid)
                    break

    # Filter background motifs if other motifs exist, unless explicitly requested or include_h_motifs is True
    if not options.get("include_h_motifs", False):
        non_bg_motifs = [m for m in motifs if m.get("compound_id") not in background_ids]
        if non_bg_motifs:
            # Keep non-background motifs + any background motifs that were explicitly requested
            motifs = [m for m in motifs if m.get("compound_id") not in background_ids or m.get("compound_id") in requested_ids]

    analyses = []
    for hit in motifs:
        compound_id = hit["compound_id"]
        # Aryl/Heteroaryl motifs
        if compound_id.startswith((
            "Ar-", "AromN-", "Pyridine-", "Pyrimidine-", "Pyrrole-", "Indole-", 
            "Thiophene-", "Furan-", "Imidazole-", "Quinoline-", "Isoquinoline-"
        )):
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
            nearby = analyze_nearby_groups(mol, hit, all_motifs, groups, compound_map)
            analyses.append(
                {
                    "compound_id": compound_id,
                    "center": {"ipso": hit["a_atom_idx"], "bond": hit["bond"]},
                    "steric": steric,
                    "electronic": electronic,
                    "nearby_groups": nearby,
                    "undocumented": hit.get("undocumented", False),
                }
            )
        # Alkyl/Generic motifs
        elif compound_id.startswith((
            "R-", "Bn-", "Allyl-", "Any-", "RCH2-", "R2CH-", "R3C-", 
            "Vinyl-", "Alkynyl-", "Acyl-", "Propargyl-"
        )):
            steric = analyze_alkyl_steric(mol, hit)
            nearby = analyze_nearby_groups(mol, hit, all_motifs, groups, compound_map)
            analyses.append(
                {
                    "compound_id": compound_id,
                    "center": {"alpha_c": hit["a_atom_idx"], "bond": hit["bond"]},
                    "steric": steric,
                    "nearby_groups": nearby,
                    "undocumented": hit.get("undocumented", False),
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
            "undocumented": analysis.get("undocumented", False),
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
                    "undocumented": analysis.get("undocumented", False),
                }
            )
        if "nearby_groups" in analysis:
            nearby_payload.append(
                {
                    "compound_id": analysis.get("compound_id"),
                    "center": analysis.get("center"),
                    "result": analysis.get("nearby_groups"),
                    "undocumented": analysis.get("undocumented", False),
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


def analyze_smiles(
    smiles: str,
    registry_paths: Optional[Dict[str, str | Path]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Compatibility wrapper for motif/steric/electronic analysis."""
    return featurize_molecule(smiles, registry_paths=registry_paths, options=options)


def _default_registry_paths() -> Dict[str, Path]:
    base = Path(__file__).resolve().parent.parent / "taxonomy" / "data"
    return {
        "groups": base / "organic_groups.v1.3.json",
        "compounds": base / "organic_compounds.v1.3.json",
    }


__all__ = ["featurize_molecule", "analyze_smiles"]

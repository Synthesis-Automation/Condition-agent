"""
Reaction-level structural analysis using motif-based molecule featurization.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from chemtools.smiles import normalize_reaction

from .molecule import featurize_molecule
from .reaction_detection import detect_reaction_types


def featurize_reaction(
    reaction_smiles: str,
    registry_paths: Optional[Dict[str, str | Path]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Systematically analyze a reaction: detect type and featurize all reactants.
    """
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


__all__ = ["featurize_reaction"]

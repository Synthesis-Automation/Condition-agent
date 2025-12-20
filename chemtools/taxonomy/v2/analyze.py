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
from .motif_detect import detect_motifs
from .motif_registry import build_compound_registry


def analyze_smiles(
    smiles: str,
    registry_paths: Optional[Dict[str, str | Path]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Analyze motifs, sterics, and electronics for a SMILES string.
    """
    if not rdkit_available():
        return {"smiles": smiles, "motifs": [], "analyses": [], "error": "rdkit_unavailable"}

    mol = parse_smiles(smiles)
    if mol is None:
        return {"smiles": smiles, "motifs": [], "analyses": [], "error": "invalid_smiles"}

    registry_paths = registry_paths or _default_registry_paths()
    registry = build_compound_registry(registry_paths)
    compiled = registry["compiled_compounds"]
    groups = registry["groups"]

    options = options or {}
    include_gasteiger = bool(options.get("include_gasteiger", False))
    include_ipso_group = options.get("electronics_include_ipso_group", True)
    max_hits = options.get("max_hits_per_compound")

    motifs = detect_motifs(mol, compiled, max_hits_per_compound=max_hits)
    motifs = _filter_arom_duplicates(motifs)
    analyses = []
    for hit in motifs:
        compound_id = hit["compound_id"]
        if compound_id.startswith(("Ar-", "Arom-")):
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
            analyses.append(
                {
                    "compound_id": compound_id,
                    "center": {"ipso": hit["a_atom_idx"], "bond": hit["bond"]},
                    "steric": steric,
                    "electronic": electronic,
                }
            )
        elif compound_id.startswith(("R-", "Bn-", "Allyl-")):
            steric = analyze_alkyl_steric(mol, hit)
            analyses.append(
                {
                    "compound_id": compound_id,
                    "center": {"alpha_c": hit["a_atom_idx"], "bond": hit["bond"]},
                    "steric": steric,
                }
            )

    return {"smiles": smiles, "motifs": motifs, "analyses": analyses}


def _filter_arom_duplicates(motifs: list[Dict[str, Any]]) -> list[Dict[str, Any]]:
    ar_hits = set()
    for hit in motifs:
        compound_id = hit.get("compound_id", "")
        if not compound_id.startswith("Ar-"):
            continue
        suffix = compound_id[3:]
        a_idx = hit.get("a_atom_idx")
        b_idx = hit.get("b_atom_idx")
        if a_idx is None or b_idx is None:
            continue
        ar_hits.add((suffix, a_idx, b_idx))

    filtered: list[Dict[str, Any]] = []
    for hit in motifs:
        compound_id = hit.get("compound_id", "")
        if compound_id.startswith("Arom-"):
            suffix = compound_id[5:]
            a_idx = hit.get("a_atom_idx")
            b_idx = hit.get("b_atom_idx")
            if a_idx is not None and b_idx is not None and (suffix, a_idx, b_idx) in ar_hits:
                continue
        filtered.append(hit)
    return filtered


def _default_registry_paths() -> Dict[str, Path]:
    base = Path(__file__).resolve().parent
    return {
        "groups": base / "organic_groups.v1.2.json",
        "compounds": base / "organic_compounds.v1.2.json",
        "templates": base / "smarts_templates.v1.json",
    }

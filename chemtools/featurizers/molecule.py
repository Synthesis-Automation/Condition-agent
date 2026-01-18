"""
Motif-based steric and electronic analysis using organic compound motifs.
"""

from __future__ import annotations

from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Dict, Optional, List

from chemtools.util.rdkit_helpers import parse_smiles, rdkit_available
from chemtools.util.smarts_cache import compile_smarts

from .alkyl_steric import analyze_alkyl_steric
from .aryl_electronics import analyze_aryl_electronics
from .aryl_steric import analyze_aryl_steric
from .nearby_groups import analyze_nearby_groups
from .motif_detect import detect_motifs
from .motif_registry import build_compound_registry, _default_registry_paths

_INORGANIC_SMARTS = (
    "O=C=O",
    "[CX3](=O)([O-])[O-]",
    "[CX3](=O)([O-])O",
    "[CX3](=O)(O)O",
    "[OX2H]C(=O)[O-]",
)


def _is_inorganic_molecule(mol: Any) -> bool:
    has_carbon = any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if not has_carbon:
        return True
    for smarts in _INORGANIC_SMARTS:
        pattern = compile_smarts(smarts, validate=False)
        if pattern and mol.HasSubstructMatch(pattern):
            return True
    return False


_ORGANOMETAL_B_GROUPS = {
    "B_Any",
    "B(OH)2",
    "B(OR)2",
    "Bpin",
    "BF3K",
    "Zn",
    "Mg",
    "Sn",
    "Si",
}
_ORGANOMETAL_REACTIVITY_BONUS = 4.0
_REACTANT_MOTIF_BONUS = 2.0


@lru_cache(maxsize=1)
def _load_reaction_reactant_motifs() -> set[str]:
    base = Path(__file__).resolve().parents[1] / "taxonomy" / "data"
    rt_path = base / "reaction_types.v4.0.json"
    logic_path = base / "compound_logic.json"
    if not rt_path.exists() or not logic_path.exists():
        return set()
    try:
        with rt_path.open("r", encoding="utf-8") as handle:
            reaction_types = json.load(handle).get("reaction_types", []) or []
        with logic_path.open("r", encoding="utf-8") as handle:
            logic_data = json.load(handle)
    except Exception:
        return set()
    motif_sets = logic_data.get("motif_sets", {}) or {}
    out: set[str] = set()

    def add_value(value: Any) -> None:
        if not value:
            return
        if isinstance(value, str):
            if value.startswith("@"):
                set_name = value[1:]
                members = motif_sets.get(set_name, {}).get("members", []) or []
                for member in members:
                    if member:
                        out.add(str(member))
            else:
                out.add(value)
            return
        if isinstance(value, list):
            for item in value:
                add_value(item)
            return
        if isinstance(value, dict):
            for item in value.values():
                add_value(item)

    for entry in reaction_types:
        if not isinstance(entry, dict):
            continue
        reactants = entry.get("reactants")
        if isinstance(reactants, dict):
            for item in reactants.values():
                add_value(item)
        else:
            add_value(reactants)

    return out


def _motif_rank_score(hit: Dict[str, Any]) -> float:
    reactivity = float(hit.get("reactivity_weight") or 0.0)
    group_b = hit.get("group_b")
    if group_b in _ORGANOMETAL_B_GROUPS:
        reactivity += _ORGANOMETAL_REACTIVITY_BONUS
    compound_id = hit.get("compound_id")
    if compound_id and compound_id in _load_reaction_reactant_motifs():
        reactivity += _REACTANT_MOTIF_BONUS
    priority = int(hit.get("priority") or 0)
    complexity = int(hit.get("complexity") or 0)
    return (reactivity * 100.0) + (priority * 10.0) + complexity


@lru_cache(maxsize=1)
def _load_heterocycle_scaffold_ids() -> set[str]:
    path = Path(__file__).resolve().parents[1] / "taxonomy" / "data" / "scaffold_motifs.v1.3.json"
    if not path.exists():
        return set()
    try:
        with path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except Exception:
        return set()
    heterocycles: set[str] = set()
    for entry in payload.get("compounds", []) or []:
        if not isinstance(entry, dict):
            continue
        desc = str(entry.get("description") or "").lower()
        if "heteroaromatic" not in desc:
            continue
        motif_id = str(entry.get("id") or "").strip()
        if motif_id:
            heterocycles.add(motif_id)
    return heterocycles


def _iter_motif_ids(motif: Dict[str, Any]) -> List[str]:
    ids: List[str] = []
    if not isinstance(motif, dict):
        return ids
    cid = str(motif.get("compound_id") or "").strip()
    if cid:
        ids.append(cid)
    alt_ids = motif.get("alt_compound_ids") or []
    if isinstance(alt_ids, set):
        alt_ids = list(alt_ids)
    for alt_id in alt_ids:
        alt_text = str(alt_id).strip()
        if alt_text:
            ids.append(alt_text)
    return ids


def _collect_heterocycle_types(motifs: List[Dict[str, Any]]) -> List[str]:
    scaffold_ids = _load_heterocycle_scaffold_ids()
    if not scaffold_ids:
        return []
    hits: set[str] = set()
    for motif in motifs:
        for motif_id in _iter_motif_ids(motif):
            if motif_id in scaffold_ids:
                hits.add(motif_id)
    return sorted(hits)


def _aryl_ring_metrics(mol: Any) -> Dict[str, Any]:
    ring_info = mol.GetRingInfo()
    ring_sizes: List[int] = []
    ring_size_counts: Dict[int, int] = {}
    hetero_counts = {"N": 0, "O": 0, "S": 0}
    hetero_ring_sizes: List[int] = []
    aromatic_rings = 0

    for ring in ring_info.AtomRings():
        if not ring:
            continue
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if not atoms or not all(atom.GetIsAromatic() for atom in atoms):
            continue
        aromatic_rings += 1
        size = len(ring)
        ring_sizes.append(size)
        ring_size_counts[size] = ring_size_counts.get(size, 0) + 1
        ring_heteros = {"N": 0, "O": 0, "S": 0}
        for atom in atoms:
            atomic_num = atom.GetAtomicNum()
            if atomic_num == 7:
                ring_heteros["N"] += 1
            elif atomic_num == 8:
                ring_heteros["O"] += 1
            elif atomic_num == 16:
                ring_heteros["S"] += 1
        if sum(ring_heteros.values()) > 0:
            hetero_ring_sizes.append(size)
        for key, value in ring_heteros.items():
            hetero_counts[key] += value

    return {
        "aromatic_ring_count": aromatic_rings,
        "ring_sizes": sorted(set(ring_sizes)),
        "ring_size_counts": {str(size): count for size, count in sorted(ring_size_counts.items())},
        "heteroaromatic": bool(hetero_ring_sizes),
        "hetero_ring_sizes": sorted(set(hetero_ring_sizes)),
        "hetero_counts": hetero_counts,
    }


def _empty_aryl_analysis() -> Dict[str, Any]:
    return {
        "is_heterocycle": False,
        "heterocycle_types": [],
        "aromatic_ring_count": 0,
        "ring_sizes": [],
        "ring_size_counts": {},
        "heteroaromatic": False,
        "hetero_ring_sizes": [],
        "hetero_counts": {"N": 0, "O": 0, "S": 0},
    }


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
            "ranked_motifs": [],
            "steric": {"aryl": [], "alkyl": []},
            "electronics": {"aryl": []},
            "aryl_analysis": _empty_aryl_analysis(),
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
            "ranked_motifs": [],
            "steric": {"aryl": [], "alkyl": []},
            "electronics": {"aryl": []},
            "aryl_analysis": _empty_aryl_analysis(),
            "analyses": [],
            "meta": meta,
        }

    if _is_inorganic_molecule(mol):
        return {
            "schema_version": "v2",
            "smiles": smiles,
            "motifs": [{"compound_id": "Inorganic"}],
            "ranked_motifs": ["Inorganic"],
            "steric": {"aryl": [], "alkyl": []},
            "electronics": {"aryl": []},
            "aryl_analysis": _empty_aryl_analysis(),
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
    include_analysis_details = bool(options.get("include_analysis_details", False))
    include_steric_details = bool(
        options.get("include_steric_details", include_analysis_details)
    )
    include_electronic_details = bool(
        options.get(
            "include_electronic_details",
            include_analysis_details or include_gasteiger,
        )
    )
    max_hits = options.get("max_hits_per_compound")
    discovery_mode = bool(options.get("discovery_mode", False))
    site_filter = options.get("motif_site_filter", "bond")

    all_motifs = detect_motifs(
        mol,
        compiled,
        max_hits_per_compound=max_hits,
        registry=registry,
        discovery_mode=discovery_mode,
        site_filter=site_filter,
    )
    motifs = list(all_motifs)
    heterocycle_types = _collect_heterocycle_types(motifs)
    ring_metrics = _aryl_ring_metrics(mol)
    aryl_analysis = {
        "is_heterocycle": bool(heterocycle_types) or ring_metrics.get("heteroaromatic", False),
        "heterocycle_types": heterocycle_types,
        **ring_metrics,
    }

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
    background_ids = {"Ar-H", "R-H", "Any-H", "Alkyl-H", "Alkenyl-H", "Alkynyl-H"}
    
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
    # Collapse duplicate background motifs (e.g., Alkyl-H per C-H bond).
    if motifs:
        seen_background = set()
        deduped = []
        for m in motifs:
            cid = m.get("compound_id")
            if cid in background_ids:
                if cid in seen_background:
                    continue
                seen_background.add(cid)
            deduped.append(m)
        motifs = deduped

    ranked_motifs = []
    if motifs:
        for m in motifs:
            m["rank_score"] = _motif_rank_score(m)
        motifs.sort(key=lambda m: (m.get("rank_score", 0.0), m.get("compound_id", "")), reverse=True)
        ranked_motifs = [m.get("compound_id") for m in motifs if m.get("compound_id")]

    analyses = []
    for hit in motifs:
        compound_id = hit["compound_id"]
        # Aryl/Heteroaryl motifs
        if compound_id.startswith((
            "Ar-", "AromN-", "Pyridine-", "Pyrimidine-", "Pyrrole-", "Indole-", 
            "Thiophene-", "Furan-", "Imidazole-", "Quinoline-", "Isoquinoline-"
        )):
            steric = analyze_aryl_steric(mol, hit, include_details=include_steric_details)
            if include_ipso_group == "both":
                electronic = [
                    analyze_aryl_electronics(
                        mol,
                        hit,
                        groups,
                        include_ipso_group=True,
                        include_gasteiger=include_gasteiger,
                        include_details=include_electronic_details,
                    ),
                    analyze_aryl_electronics(
                        mol,
                        hit,
                        groups,
                        include_ipso_group=False,
                        include_gasteiger=include_gasteiger,
                        include_details=include_electronic_details,
                    ),
                ]
            else:
                electronic = analyze_aryl_electronics(
                    mol,
                    hit,
                    groups,
                    include_ipso_group=bool(include_ipso_group),
                    include_gasteiger=include_gasteiger,
                    include_details=include_electronic_details,
                )
            nearby = analyze_nearby_groups(mol, hit, all_motifs, groups, compound_map)
            analyses.append(
                {
                    "compound_id": compound_id,
                    "alt_compound_ids": hit.get("alt_compound_ids", set()),
                    "center": {"ipso": hit["a_atom_idx"], "bond": hit["bond"]},
                    "steric": steric,
                    "electronic": electronic,
                    "nearby_groups": nearby,
                    "undocumented": hit.get("undocumented", False),
                }
            )
        # Alkyl/Generic motifs
        elif compound_id.startswith((
            "R-", "Bn-", "Allyl-", "Alkyl-", "Alkenyl-", "Any-", "RCH2-", "R2CH-", "R3C-",
            "Vinyl-", "Alkynyl-", "Acyl-", "Propargyl-", "H-"
        )):
            steric = analyze_alkyl_steric(mol, hit, include_details=include_steric_details)
            nearby = analyze_nearby_groups(mol, hit, all_motifs, groups, compound_map)
            analyses.append(
                {
                    "compound_id": compound_id,
                    "alt_compound_ids": hit.get("alt_compound_ids", set()),
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
        alt_ids = analysis.get("alt_compound_ids", set())
        steric_entry = {
            "compound_id": analysis.get("compound_id"),
            "alt_compound_ids": list(alt_ids) if isinstance(alt_ids, set) else alt_ids,
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
                    "alt_compound_ids": list(alt_ids) if isinstance(alt_ids, set) else alt_ids,
                    "center": analysis.get("center"),
                    "result": analysis.get("electronic"),
                    "undocumented": analysis.get("undocumented", False),
                }
            )
        if "nearby_groups" in analysis:
            nearby_payload.append(
                {
                    "compound_id": analysis.get("compound_id"),
                    "alt_compound_ids": list(alt_ids) if isinstance(alt_ids, set) else alt_ids,
                    "center": analysis.get("center"),
                    "result": analysis.get("nearby_groups"),
                    "undocumented": analysis.get("undocumented", False),
                }
            )

    return {
        "schema_version": "v2",
        "smiles": smiles,
        "motifs": motifs,
        "ranked_motifs": ranked_motifs,
        "steric": steric_payload,
        "electronics": electronic_payload,
        "nearby": nearby_payload,
        "aryl_analysis": aryl_analysis,
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


__all__ = ["featurize_molecule", "analyze_smiles"]

"""
Motif-based steric and electronic analysis using organic compound motifs.
"""

from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, List

from chemtools.core.rdkit import parse_smiles, rdkit_available
from chemtools.core.smarts import compile_smarts

from .electronics import analyze_aryl_electronics
from .sterics_alkyl import analyze_alkyl_steric
from .sterics_aryl import analyze_aryl_steric
from .nearby_groups import analyze_nearby_groups
from .motifs.detection import detect_motifs
from .motifs.registry import build_compound_registry, _default_registry_paths

_DEFAULT_INORGANIC_SMARTS = (
    "O=C=O",
    "[CX3](=O)([O-])[O-]",
    "[CX3](=O)([O-])O",
    "[CX3](=O)(O)O",
    "[OX2H]C(=O)[O-]",
)
_DEFAULT_BACKGROUND_MOTIF_IDS = {"Ar-H", "R-H", "Any-H", "Alkyl-H", "Alkenyl-H", "Alkynyl-H"}
_DEFAULT_ARYL_PREFIXES = (
    "Ar-",
    "AromN-",
    "Pyridine-",
    "Pyrimidine-",
    "Pyrrole-",
    "Indole-",
    "Thiophene-",
    "Furan-",
    "Imidazole-",
    "Quinoline-",
    "Isoquinoline-",
)
_DEFAULT_ALKYL_PREFIXES = (
    "R-",
    "Bn-",
    "Allyl-",
    "Alkyl-",
    "Alkenyl-",
    "Any-",
    "RCH2-",
    "R2CH-",
    "R3C-",
    "Vinyl-",
    "Alkynyl-",
    "Acyl-",
    "Propargyl-",
    "H-",
)

_DEFAULT_ORGANOMETAL_B_GROUPS = {
    "B_Any",
    "B(OH)2",
    "B(OR)2",
    "BF3K",
    "Zn",
    "Mg",
    "Sn",
    "Si",
}
_DEFAULT_ORGANOMETAL_REACTIVITY_BONUS = 4.0
_DEFAULT_REACTANT_MOTIF_BONUS = 2.0
_ORGANOMETAL_HALIDE_REWRITE_RE = re.compile(r"^\[(Mg|Zn)\](Cl|Br|I)(.+)$")


@lru_cache(maxsize=1)
def _load_molecule_logic() -> Dict[str, Any]:
    try:
        from ..taxonomy import loader as taxonomy_loader
    except Exception:
        return {}
    payload = taxonomy_loader.load_featurizer_logic()
    if not payload:
        return {}
    section = payload.get("molecule", {}) or {}
    return section if isinstance(section, dict) else {}


def _inorganic_smarts() -> List[str]:
    logic = _load_molecule_logic()
    values = logic.get("inorganic_smarts", [])
    if isinstance(values, list):
        items = [str(v) for v in values if isinstance(v, str) and v.strip()]
        if items:
            return items
    return list(_DEFAULT_INORGANIC_SMARTS)


def _background_motif_ids() -> set[str]:
    logic = _load_molecule_logic()
    values = logic.get("background_motif_ids", [])
    if isinstance(values, list):
        items = {str(v) for v in values if isinstance(v, str) and v.strip()}
        if items:
            return items
    return set(_DEFAULT_BACKGROUND_MOTIF_IDS)


def _motif_family_prefixes() -> tuple[tuple[str, ...], tuple[str, ...]]:
    logic = _load_molecule_logic()
    family = logic.get("motif_family_prefixes", {}) or {}
    if not isinstance(family, dict):
        return _DEFAULT_ARYL_PREFIXES, _DEFAULT_ALKYL_PREFIXES
    aryl = family.get("aryl")
    alkyl = family.get("alkyl")
    aryl_prefixes = tuple(str(v) for v in aryl if isinstance(v, str) and v.strip()) if isinstance(aryl, list) else ()
    alkyl_prefixes = tuple(str(v) for v in alkyl if isinstance(v, str) and v.strip()) if isinstance(alkyl, list) else ()
    return (
        aryl_prefixes or _DEFAULT_ARYL_PREFIXES,
        alkyl_prefixes or _DEFAULT_ALKYL_PREFIXES,
    )


def _motif_ranking_rules() -> tuple[set[str], float, float]:
    logic = _load_molecule_logic()
    ranking = logic.get("motif_ranking", {}) or {}
    if not isinstance(ranking, dict):
        return (
            set(_DEFAULT_ORGANOMETAL_B_GROUPS),
            _DEFAULT_ORGANOMETAL_REACTIVITY_BONUS,
            _DEFAULT_REACTANT_MOTIF_BONUS,
        )
    raw_groups = ranking.get("organometal_group_b", [])
    groups = (
        {str(v) for v in raw_groups if isinstance(v, str) and v.strip()}
        if isinstance(raw_groups, list)
        else set()
    )
    organometal_bonus = ranking.get("organometal_bonus", _DEFAULT_ORGANOMETAL_REACTIVITY_BONUS)
    reactant_bonus = ranking.get("reactant_motif_bonus", _DEFAULT_REACTANT_MOTIF_BONUS)
    return (
        groups or set(_DEFAULT_ORGANOMETAL_B_GROUPS),
        float(organometal_bonus),
        float(reactant_bonus),
    )


def _is_inorganic_molecule(mol: Any) -> bool:
    has_carbon = any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if not has_carbon:
        return True
    for smarts in _inorganic_smarts():
        pattern = compile_smarts(smarts, validate=False)
        if pattern and mol.HasSubstructMatch(pattern):
            return True
    return False


def _build_payload(
    *,
    smiles: str,
    meta: Dict[str, Any],
    motifs: Optional[List[Dict[str, Any]]] = None,
    context_motifs: Optional[List[Dict[str, Any]]] = None,
    ranked_motifs: Optional[List[str]] = None,
    steric: Optional[Dict[str, Any]] = None,
    electronics: Optional[Dict[str, Any]] = None,
    nearby: Optional[List[Dict[str, Any]]] = None,
    aryl_analysis: Optional[Dict[str, Any]] = None,
    analyses: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    motif_list = motifs or []
    grouped_motifs = _group_motif_hits(motif_list)
    return {
        "schema_version": "v2",
        "smiles": smiles,
        "motifs": motif_list,
        "context_motifs": context_motifs or [],
        "ranked_motifs": ranked_motifs or [],
        "unique_motifs": [entry["compound_id"] for entry in grouped_motifs],
        "grouped_motifs": grouped_motifs,
        "steric": steric or {"aryl": [], "alkyl": []},
        "electronics": electronics or {"aryl": []},
        "nearby": nearby or [],
        "aryl_analysis": aryl_analysis or _empty_aryl_analysis(),
        "analyses": analyses or [],
        "meta": meta,
    }


def _group_motif_hits(motifs: Iterable[Dict[str, Any]]) -> List[Dict[str, Any]]:
    grouped: Dict[str, Dict[str, Any]] = {}
    for motif in motifs or []:
        if not isinstance(motif, dict):
            continue
        compound_id = str(motif.get("compound_id") or motif.get("id") or "").strip()
        if not compound_id:
            continue
        entry = grouped.setdefault(
            compound_id,
            {
                "compound_id": compound_id,
                "count": 0,
                "sites": [],
                "rank_score": motif.get("rank_score", 0.0),
                "undocumented": bool(motif.get("undocumented", False)),
            },
        )
        entry["count"] += 1
        if motif.get("rank_score", 0.0) > entry.get("rank_score", 0.0):
            entry["rank_score"] = motif.get("rank_score", 0.0)
        site = {
            key: motif.get(key)
            for key in ("a_atom_idx", "b_atom_idx", "bond")
            if motif.get(key) is not None
        }
        if site:
            entry["sites"].append(site)
        alt_ids = _normalize_alt_ids(motif.get("alt_compound_ids"))
        if alt_ids:
            existing = set(entry.get("alt_compound_ids") or [])
            entry["alt_compound_ids"] = sorted(existing | set(alt_ids))

    return sorted(
        grouped.values(),
        key=lambda item: (float(item.get("rank_score") or 0.0), item.get("compound_id", "")),
        reverse=True,
    )


def _organometal_rescue_smiles(smiles: str) -> Optional[str]:
    """Rewrite bracketed metal-halide-carbon notation into R-M-X connectivity."""
    stripped = (smiles or "").strip()
    match = _ORGANOMETAL_HALIDE_REWRITE_RE.match(stripped)
    if not match:
        return None
    metal, halide, remainder = match.groups()
    if not remainder:
        return None
    return f"{halide}[{metal}]{remainder}"


def _has_organometal_motif(motifs: Iterable[Dict[str, Any]]) -> bool:
    for motif in motifs or []:
        if not isinstance(motif, dict):
            continue
        group_b = str(motif.get("group_b") or "")
        compound_id = str(motif.get("compound_id") or motif.get("id") or "")
        if group_b in {"-Mg*", "-Zn*", "Mg", "Zn"}:
            return True
        if "-Mg*" in compound_id or "-Zn*" in compound_id:
            return True
    return False


def _normalize_alt_ids(alt_ids: Any) -> List[str]:
    if isinstance(alt_ids, set):
        return sorted(str(item) for item in alt_ids if item)
    if isinstance(alt_ids, list):
        return [str(item) for item in alt_ids if item]
    return []


def _normalize_target_groups(target_groups: Any) -> List[str]:
    if not target_groups:
        return []
    if isinstance(target_groups, str):
        return [target_groups]
    return [str(item) for item in target_groups if item]


def _matches_target_group(compound_id: str, target: str) -> bool:
    if compound_id == target:
        return True
    if compound_id.endswith("-" + target):
        return True
    return target.startswith("-") and compound_id.endswith(target)


def _filter_motifs_by_targets(
    motifs: List[Dict[str, Any]],
    target_groups: List[str],
) -> tuple[List[Dict[str, Any]], set[str]]:
    if not target_groups:
        return motifs, set()

    filtered: List[Dict[str, Any]] = []
    requested_ids: set[str] = set()

    for motif in motifs:
        compound_id = motif.get("compound_id", "")
        if not compound_id:
            continue
        if any(_matches_target_group(compound_id, target) for target in target_groups):
            filtered.append(motif)
            requested_ids.add(compound_id)

    if filtered:
        return filtered, requested_ids

    return motifs, set()


def _filter_background_motifs(
    motifs: List[Dict[str, Any]],
    *,
    include_h_motifs: bool,
    requested_ids: set[str],
) -> List[Dict[str, Any]]:
    if include_h_motifs:
        return motifs
    background_ids = _background_motif_ids()

    non_background = [m for m in motifs if m.get("compound_id") not in background_ids]
    if non_background:
        return [
            m
            for m in motifs
            if m.get("compound_id") not in background_ids
            or m.get("compound_id") in requested_ids
        ]
    return motifs


def _dedupe_background_motifs(motifs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    background_ids = _background_motif_ids()
    seen: set[str] = set()
    deduped: List[Dict[str, Any]] = []
    for motif in motifs:
        compound_id = motif.get("compound_id")
        if compound_id in background_ids:
            if compound_id in seen:
                continue
            seen.add(compound_id)
        deduped.append(motif)
    return deduped


def _motif_family(compound_id: str) -> Optional[str]:
    aryl_prefixes, alkyl_prefixes = _motif_family_prefixes()
    if compound_id.startswith(aryl_prefixes):
        return "aryl"
    if compound_id.startswith(alkyl_prefixes):
        return "alkyl"
    return None


@lru_cache(maxsize=1)
def _load_reaction_reactant_motifs() -> set[str]:
    try:
        from ..taxonomy import loader as taxonomy_loader
    except Exception:
        return set()

    reaction_payload = taxonomy_loader.load_reaction_types_raw()
    logic_data = taxonomy_loader.load_compound_logic()
    reaction_types = reaction_payload.get("reaction_types", []) or []
    if not reaction_types or not logic_data:
        return set()

    try:
        motif_sets = logic_data.get("motif_sets", {}) or {}
    except Exception:
        motif_sets = {}
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
    organometal_groups, organometal_bonus, reactant_motif_bonus = _motif_ranking_rules()
    reactivity = float(hit.get("reactivity_weight") or 0.0)
    group_b = hit.get("group_b")
    if group_b in organometal_groups:
        reactivity += organometal_bonus
    compound_id = hit.get("compound_id")
    if compound_id and compound_id in _load_reaction_reactant_motifs():
        reactivity += reactant_motif_bonus
    priority = int(hit.get("priority") or 0)
    complexity = int(hit.get("complexity") or 0)
    return (reactivity * 100.0) + (priority * 10.0) + complexity


def _motif_rank_breakdown(hit: Dict[str, Any]) -> Dict[str, Any]:
    organometal_groups, organometal_bonus, reactant_motif_bonus = _motif_ranking_rules()
    compound_id = hit.get("compound_id")
    group_b = hit.get("group_b")
    base_reactivity = float(hit.get("reactivity_weight") or 0.0)
    organometal_bonus_value = organometal_bonus if group_b in organometal_groups else 0.0
    reactant_bonus = 0.0
    if compound_id and compound_id in _load_reaction_reactant_motifs():
        reactant_bonus = reactant_motif_bonus
    reactivity_total = base_reactivity + organometal_bonus_value + reactant_bonus
    priority = int(hit.get("priority") or 0)
    complexity = int(hit.get("complexity") or 0)
    score = (reactivity_total * 100.0) + (priority * 10.0) + complexity
    return {
        "compound_id": compound_id,
        "group_b": group_b,
        "base_reactivity": base_reactivity,
        "organometal_bonus": organometal_bonus_value,
        "reactant_bonus": reactant_bonus,
        "reactivity_total": reactivity_total,
        "priority": priority,
        "complexity": complexity,
        "score": score,
    }


def print_motif_rank_breakdown(motifs: List[Dict[str, Any]]) -> None:
    """
    Print a scoring breakdown for each motif hit.
    """
    if not motifs:
        print("No motifs to score.")
        return
    for idx, hit in enumerate(motifs, start=1):
        breakdown = _motif_rank_breakdown(hit)
        compound_id = breakdown["compound_id"] or "Unknown"
        group_b = breakdown["group_b"] or "?"
        print(
            f"[{idx}] {compound_id} (group_b={group_b}) "
            f"score={breakdown['score']:.2f} "
            f"reactivity={breakdown['reactivity_total']:.2f} "
            f"(base={breakdown['base_reactivity']:.2f} "
            f"+ organometal={breakdown['organometal_bonus']:.2f} "
            f"+ reactant={breakdown['reactant_bonus']:.2f}) "
            f"priority={breakdown['priority']} complexity={breakdown['complexity']}"
        )


@lru_cache(maxsize=1)
def _load_heterocycle_scaffold_ids() -> set[str]:
    try:
        from ..taxonomy import loader as taxonomy_loader
    except Exception:
        return set()
    payload = taxonomy_loader.load_scaffold_motifs()
    if not payload:
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


def _strip_scaffold_registry_paths(registry_paths: Dict[str, str | Path]) -> Dict[str, str | Path]:
    return {k: v for k, v in registry_paths.items() if k != "scaffold_motifs"}


def _scaffold_only_registry_paths(
    registry_paths: Dict[str, str | Path],
) -> Optional[Dict[str, str | Path]]:
    scaffold_path = registry_paths.get("scaffold_motifs")
    groups_path = registry_paths.get("groups")
    if not scaffold_path or not groups_path:
        return None
    return {"groups": groups_path, "compounds": scaffold_path}


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
        return _build_payload(smiles=smiles, meta=meta)

    mol = parse_smiles(smiles)
    if mol is None:
        meta["error"] = "invalid_smiles"
        return _build_payload(smiles=smiles, meta=meta)

    is_inorganic = _is_inorganic_molecule(mol)

    from rdkit import Chem
    mol = Chem.AddHs(mol)

    registry_paths = registry_paths or _default_registry_paths()
    core_registry_paths = _strip_scaffold_registry_paths(registry_paths)
    registry = build_compound_registry(core_registry_paths)
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
    rescue_smiles = _organometal_rescue_smiles(smiles)
    if rescue_smiles and not _has_organometal_motif(motifs):
        rescue_mol = parse_smiles(rescue_smiles)
        if rescue_mol is not None:
            rescue_mol = Chem.AddHs(rescue_mol)
            rescue_motifs = list(
                detect_motifs(
                    rescue_mol,
                    compiled,
                    max_hits_per_compound=max_hits,
                    registry=registry,
                    discovery_mode=discovery_mode,
                    site_filter=site_filter,
                )
            )
            if _has_organometal_motif(rescue_motifs):
                mol = rescue_mol
                motifs = rescue_motifs
                meta["organometal_rescue"] = {
                    "input_smiles": smiles,
                    "normalized_smiles": rescue_smiles,
                    "reason": "rewrote leading [Mg/Zn]X-C notation to X-[Mg/Zn]-C connectivity",
                }
    context_motifs: List[Dict[str, Any]] = []
    if options.get("include_scaffold_motifs", True):
        scaffold_paths = _scaffold_only_registry_paths(registry_paths)
        if scaffold_paths:
            scaffold_registry = build_compound_registry(scaffold_paths)
            scaffold_compiled = scaffold_registry["compiled_compounds"]
            context_motifs = detect_motifs(
                mol,
                scaffold_compiled,
                max_hits_per_compound=max_hits,
                # Context pass is scaffold-only; disable registry-driven Any-* fallbacks
                # so generic substituent labels do not leak into reaction delta/CRK.
                registry=None,
                discovery_mode=discovery_mode,
                site_filter=site_filter,
            )
    all_context_motifs = motifs + context_motifs
    heterocycle_types = _collect_heterocycle_types(all_context_motifs)
    ring_metrics = _aryl_ring_metrics(mol)
    aryl_analysis = {
        "is_heterocycle": bool(heterocycle_types) or ring_metrics.get("heteroaromatic", False),
        "heterocycle_types": heterocycle_types,
        **ring_metrics,
    }

    target_groups = _normalize_target_groups(options.get("target_groups"))
    motifs, requested_ids = _filter_motifs_by_targets(motifs, target_groups)
    motifs = _filter_background_motifs(
        motifs,
        include_h_motifs=bool(options.get("include_h_motifs", False)),
        requested_ids=requested_ids,
    )
    motifs = _dedupe_background_motifs(motifs)
    if is_inorganic and not motifs:
        return _build_payload(
            smiles=smiles,
            meta=meta,
            motifs=[{"compound_id": "Inorganic"}],
            ranked_motifs=["Inorganic"],
        )

    ranked_motifs = []
    if motifs:
        for m in motifs:
            m["rank_score"] = _motif_rank_score(m)
        motifs.sort(key=lambda m: (m.get("rank_score", 0.0), m.get("compound_id", "")), reverse=True)
        ranked_motifs = [m.get("compound_id") for m in motifs if m.get("compound_id")]

    analyses = []
    for hit in motifs:
        compound_id = hit["compound_id"]
        family = _motif_family(compound_id)
        if family == "aryl":
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
            nearby = analyze_nearby_groups(mol, hit, all_context_motifs, groups, compound_map)
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
        elif family == "alkyl":
            steric = analyze_alkyl_steric(mol, hit, include_details=include_steric_details)
            nearby = analyze_nearby_groups(mol, hit, all_context_motifs, groups, compound_map)
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
        alt_ids = _normalize_alt_ids(analysis.get("alt_compound_ids", set()))
        steric_entry = {
            "compound_id": analysis.get("compound_id"),
            "alt_compound_ids": alt_ids,
            "center": analysis.get("center"),
            "result": analysis.get("steric"),
            "undocumented": analysis.get("undocumented", False),
        }
        if _motif_family(analysis.get("compound_id", "")) == "aryl":
            steric_payload["aryl"].append(steric_entry)
        else:
            steric_payload["alkyl"].append(steric_entry)
        if "electronic" in analysis:
            electronic_payload["aryl"].append(
                {
                    "compound_id": analysis.get("compound_id"),
                    "alt_compound_ids": alt_ids,
                    "center": analysis.get("center"),
                    "result": analysis.get("electronic"),
                    "undocumented": analysis.get("undocumented", False),
                }
            )
        if "nearby_groups" in analysis:
            nearby_payload.append(
                {
                    "compound_id": analysis.get("compound_id"),
                    "alt_compound_ids": alt_ids,
                    "center": analysis.get("center"),
                    "result": analysis.get("nearby_groups"),
                    "undocumented": analysis.get("undocumented", False),
                }
            )

    return _build_payload(
        smiles=smiles,
        meta=meta,
        motifs=motifs,
        context_motifs=context_motifs,
        ranked_motifs=ranked_motifs,
        steric=steric_payload,
        electronics=electronic_payload,
        nearby=nearby_payload,
        aryl_analysis=aryl_analysis,
        analyses=analyses,
    )


def analyze_smiles(
    smiles: str,
    registry_paths: Optional[Dict[str, str | Path]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Compatibility wrapper for motif/steric/electronic analysis."""
    return featurize_molecule(smiles, registry_paths=registry_paths, options=options)


__all__ = ["featurize_molecule", "analyze_smiles", "print_motif_rank_breakdown"]

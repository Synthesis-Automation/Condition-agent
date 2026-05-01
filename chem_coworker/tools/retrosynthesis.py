"""
Retrosynthesis tools for ChemCoworker.

Nine ToolPlugin entries for the retrosynthesis pipeline:

  inspect_target              — enhanced molecular analysis for retrosynthesis
  identify_retrons            — SMARTS retron pattern matching
  generate_disconnections     — core retrosynthetic transform engine
  find_retro_precedent        — search notes + literature for synthesis precedent
  search_hte_precedent        — DRFP k-NN search in the 231k-reaction HTE database
  search_by_product_similarity — Morgan FP product-space search (data-driven,
                                  no retron patterns; runs before generate_disconnections)
  apply_hte_templates         — apply 35+ atom-mapped HTE retrosynthetic SMARTS templates
                                  (SNAr, Chan-Lam, CuAAC, HWE, Knoevenagel, Wacker, etc.)
                                  via AllChem.RunReactants(); enriched with HTE conditions
  plan_route                  — multi-step greedy BFS route planner (AiZynthFinder-inspired);
                                  recursively disconnects complex precursors until all
                                  fragments are simple building blocks; InChI key cycle
                                  detection; returns complete route in a single tool call
  plan_route_candidates       — beam-search route candidate generator for strategy-aware
                                  comparison before the LLM presents a final route

All follow the _success() / _error() contract. Registered as ToolPlugins
so they appear in REASON_PROMPT tool descriptions and can be called by the executor.

Usage (by ChemCoworker, not directly):
  When task_type == "retrosynthesis", the agent selects these tools automatically.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from ._helpers import _error, _success, _to_jsonable
from ._base import ToolPlugin


# ---------------------------------------------------------------------------
# Tool A: inspect_target
# Enhanced molecular analysis specialized for retrosynthesis planning
# ---------------------------------------------------------------------------

def _inspect_target(smiles: str = "", target_smiles: str = "") -> Dict[str, Any]:
    """Analyze a target molecule's complexity and key features for retrosynthesis.

    Returns molecular complexity (BertzCT), ring system count, stereocenters,
    functional group density, and a synthetic accessibility estimate.
    Always run this first — it tells the LLM HOW HARD the synthesis will be
    and which molecular features define the strategic challenge.

    Args:
        smiles: Target molecule SMILES (not a reaction SMILES).
        target_smiles: Alias for smiles (accepted for compatibility).

    Returns:
        dict with complexity, ring_count, aromatic_rings, stereocenters,
        fg_density, heavy_atoms, strategic_notes.
    """
    try:
        smiles = smiles or target_smiles
        from rdkit import Chem, rdBase
        from rdkit.Chem import Descriptors, rdMolDescriptors, GraphDescriptors

        # Strip reaction component if accidentally passed
        mol_smiles = smiles.split(">>")[0].split(">")[0].strip() if ">" in smiles else smiles
        with rdBase.BlockLogs():
            mol = Chem.MolFromSmiles(mol_smiles)
        if mol is None:
            return _error(f"Cannot parse target SMILES: {mol_smiles!r}")

        # Core metrics
        try:
            complexity = round(GraphDescriptors.BertzCT(mol), 1)
        except Exception:
            complexity = 0.0

        heavy_atoms = mol.GetNumHeavyAtoms()
        ring_count = rdMolDescriptors.CalcNumRings(mol)
        aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        stereocenters = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
        rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        sp3_frac = round(rdMolDescriptors.CalcFractionCSP3(mol), 3)
        mw = round(Descriptors.MolWt(mol), 1)

        # Functional group landscape (reuse existing utility)
        fg_categories: Dict[str, List[str]] = {}
        fg_total = 0
        try:
            from chemtools.util.functional_groups import get_functional_groups, get_group_categories
            fgs = get_functional_groups(mol_smiles)
            cats = get_group_categories(mol_smiles)
            fg_categories = {c: g for c, g in cats.items() if g}
            fg_total = len(fgs)
        except Exception:
            pass

        # Complexity tier for LLM guidance
        if complexity < 100:
            complexity_tier = "simple"
            tier_note = "Straightforward target. 1-2 step route likely sufficient."
        elif complexity < 300:
            complexity_tier = "moderate"
            tier_note = "Moderate complexity. Expect 2-3 step linear route."
        elif complexity < 600:
            complexity_tier = "complex"
            tier_note = "Complex target. Convergent strategy recommended."
        else:
            complexity_tier = "highly_complex"
            tier_note = "Highly complex. Multi-step with protecting groups likely needed."

        # Strategic challenge notes
        strategic_notes = []
        if stereocenters > 0:
            strategic_notes.append(
                f"{stereocenters} stereocenter(s): stereocontrol needed (asymmetric synthesis or chiral resolution)"
            )
        if ring_count > 2:
            strategic_notes.append(
                f"{ring_count} rings: ring-forming reactions may be key steps"
            )
        if sp3_frac > 0.4:
            strategic_notes.append(
                "High sp3 fraction: C–C cross-coupling with sp3 centers may be challenging"
            )
        if fg_total > 5:
            strategic_notes.append(
                f"{fg_total} functional groups: protecting group strategy likely required"
            )
        if not strategic_notes:
            strategic_notes.append("No unusual strategic challenges identified.")

        return _success({
            "smiles": mol_smiles,
            "molecular_weight": mw,
            "heavy_atoms": heavy_atoms,
            "complexity_bertz": complexity,
            "complexity_tier": complexity_tier,
            "complexity_note": tier_note,
            "ring_count": ring_count,
            "aromatic_rings": aromatic_rings,
            "stereocenters": stereocenters,
            "rotatable_bonds": rot_bonds,
            "sp3_fraction": sp3_frac,
            "functional_group_count": fg_total,
            "functional_groups_by_category": fg_categories,
            "strategic_notes": strategic_notes,
        })

    except Exception as exc:
        return _error(f"Target inspection failed: {exc}")


inspect_target_tool = ToolPlugin(
    name="inspect_target",
    category="retrosynthesis",
    description=(
        "Analyze a target molecule's complexity, ring systems, stereocenters, and functional groups "
        "for retrosynthesis planning. Run this FIRST for any retrosynthesis query to understand "
        "the synthetic challenge before selecting a disconnection strategy."
    ),
    prerequisites=["normalize_reaction"],
    fn=_inspect_target,
)


# ---------------------------------------------------------------------------
# Tool B: identify_retrons
# Match SMARTS retron patterns → ranked list of applicable reactions
# ---------------------------------------------------------------------------

def _identify_retrons(smiles: str = "", target_smiles: str = "") -> Dict[str, Any]:
    """Match retron SMARTS patterns to identify applicable retrosynthetic moves.

    Scans the target molecule for structural motifs (retrons) that imply specific
    synthetic reactions. For example, a biaryl C–C bond implies Suzuki coupling;
    an aryl amine implies Buchwald-Hartwig. Returns matches ranked by difficulty
    (easiest first) so the LLM can prioritize accessible disconnections.

    Args:
        smiles: Target molecule SMILES.
        target_smiles: Alias for smiles (accepted for compatibility).

    Returns:
        dict with matched_retrons (name, reaction_name, difficulty, description,
        match_count, notes, precursor_hints), total_matched.
    """
    try:
        smiles = smiles or target_smiles
        mol_smiles = smiles.split(">>")[0].split(">")[0].strip() if ">" in smiles else smiles

        from chemtools.retro.disconnector import find_retrons
        matches = find_retrons(mol_smiles)

        if not matches:
            return _success({
                "smiles": mol_smiles,
                "total_matched": 0,
                "matched_retrons": [],
                "message": (
                    "No standard retrons detected. Target may require non-standard "
                    "disconnections or is a simple molecule without characteristic "
                    "retron motifs. Consider C–H functionalization or protection strategies."
                ),
            })

        retrons_data = []
        for m in matches:
            # Difficulty display: dots (●●○○○)
            filled = round(m.difficulty * 5)
            difficulty_display = "●" * filled + "○" * (5 - filled)

            # Resolve canonical taxonomy family ID for this retron via registry.
            # Output key remains `taxonomy_id` for API compatibility.
            taxonomy_id: Optional[str] = None
            try:
                from chemtools.retro.reaction_registry import get_taxonomy_id_for_retron
                taxonomy_id = get_taxonomy_id_for_retron(m.retron_name)
            except Exception:
                pass

            entry: Dict[str, Any] = {
                "name": m.retron_name,
                "reaction_name": m.reaction_name,
                "difficulty": m.difficulty,
                "difficulty_display": difficulty_display,
                "description": m.description,
                "notes": m.notes,
                "precursor_hints": m.precursor_hints,
                "match_count": m.match_count,
            }
            if taxonomy_id:
                entry["taxonomy_id"] = taxonomy_id
                entry["taxonomy_family_id"] = taxonomy_id
            retrons_data.append(entry)

        return _success({
            "smiles": mol_smiles,
            "total_matched": len(retrons_data),
            "matched_retrons": retrons_data,
        })

    except Exception as exc:
        return _error(f"Retron identification failed: {exc}")


identify_retrons_tool = ToolPlugin(
    name="identify_retrons",
    category="retrosynthesis",
    description=(
        "Match retron SMARTS patterns against the target molecule to identify applicable "
        "retrosynthetic disconnections. Returns ranked list of retrons (biaryl → Suzuki, "
        "aryl-amine → Buchwald-Hartwig, etc.) sorted by difficulty. "
        "Run after inspect_target; provides the basis for generate_disconnections."
    ),
    prerequisites=["normalize_reaction"],
    fn=_identify_retrons,
)


# ---------------------------------------------------------------------------
# Tool C: generate_disconnections
# Core retrosynthesis: apply transforms and produce precursor SMILES
# ---------------------------------------------------------------------------

def _generate_disconnections(smiles: str, top_k: int = 3, top_n: int = 0) -> Dict[str, Any]:
    """Generate retrosynthetic precursor pairs for the top disconnections.

    The core retrosynthesis engine. For each matched retron, applies the
    retrosynthetic transform to split the target into two precursor SMILES.
    Scores disconnections by: complexity reduction, fragment balance, and
    reaction difficulty. Returns ranked disconnections ready for conditions
    recommendation and LLM synthesis.

    Args:
        smiles: Target molecule SMILES.
        top_k: Maximum number of disconnections to return (default 3).
        top_n: Alias for top_k (accepted for compatibility).

    Returns:
        dict with disconnections (rank, reaction_name, precursor_1, precursor_2,
        complexity_delta, fragment_balance, overall_score, description, notes).
    """
    try:
        if top_n > 0:
            top_k = top_n
        mol_smiles = smiles.split(">>")[0].split(">")[0].strip() if ">" in smiles else smiles

        from chemtools.retro.disconnector import rank_disconnections
        results = rank_disconnections(mol_smiles, top_k=top_k)

        if not results:
            return _success({
                "smiles": mol_smiles,
                "total_disconnections": 0,
                "disconnections": [],
                "message": (
                    "No retrosynthetic disconnections could be generated automatically. "
                    "This may be a simple molecule (consider commercial availability), "
                    "or requires non-standard logic (C–H activation, rearrangements). "
                    "Use search_notes and search_literature for precedent."
                ),
            })

        disconnections = []
        for rank, r in enumerate(results, 1):
            # Difficulty display
            filled = round(r.difficulty * 5)
            difficulty_display = "●" * filled + "○" * (5 - filled)

            # Complexity delta label
            if r.complexity_delta > 30:
                simplification = "strong"
            elif r.complexity_delta > 0:
                simplification = "moderate"
            elif r.complexity_delta > -20:
                simplification = "minimal"
            else:
                simplification = "none (consider alternative)"

            disconnections.append({
                "rank": rank,
                "reaction_name": r.reaction_name,
                "retron_name": r.retron_name,
                "precursor_1": r.precursor_1,
                "precursor_2": r.precursor_2,
                "description": r.description,
                "difficulty": r.difficulty,
                "difficulty_display": difficulty_display,
                "complexity_delta": round(r.complexity_delta, 1),
                "complexity_simplification": simplification,
                "fragment_balance": round(r.fragment_balance, 2),
                "overall_score": round(r.overall_score, 3),
                "notes": r.notes,
                "precursor_hints": r.precursor_hints,
            })

        return _success({
            "smiles": mol_smiles,
            "total_disconnections": len(disconnections),
            "disconnections": disconnections,
        })

    except Exception as exc:
        return _error(f"Disconnection generation failed: {exc}")


generate_disconnections_tool = ToolPlugin(
    name="generate_disconnections",
    category="retrosynthesis",
    description=(
        "Core retrosynthesis tool: apply retrosynthetic transforms to generate precursor SMILES "
        "for the top disconnections. Returns ranked disconnections with precursor_1, precursor_2, "
        "complexity_delta, and overall_score. Run after identify_retrons. "
        "The output feeds directly into recommend_conditions for each forward reaction step."
    ),
    prerequisites=["normalize_reaction", "identify_retrons"],
    fn=_generate_disconnections,
)


# ---------------------------------------------------------------------------
# Tool D: find_retro_precedent
# Search notes + literature for existing synthesis routes to similar targets
# ---------------------------------------------------------------------------

def _find_retro_precedent(smiles: str = "", reaction_name: str = "", target_smiles: str = "", reaction_type: str = "") -> Dict[str, Any]:
    """Search the knowledge base for synthesis precedent for this target.

    Queries both the reaction notes (for specific reaction-type precedent)
    and literature sources (for experimental procedures). Also searches the
    routes/ note type if available. Returns relevant excerpts so the LLM
    can integrate real experimental data into its retrosynthetic analysis.

    Args:
        smiles: Target molecule SMILES (used to derive FG search terms).
        target_smiles: Alias for smiles (accepted for compatibility).
        reaction_name: Optional reaction taxonomy name (e.g., "suzuki_miyaura")
                       to focus the search. Pass empty string for broad search.
        reaction_type: Alias for reaction_name (accepted for compatibility).

    Returns:
        dict with notes_hits, literature_hits, route_hits, search_terms_used.
    """
    try:
        smiles = smiles or target_smiles
        reaction_name = reaction_name or reaction_type  # alias
        mol_smiles = smiles.split(">>")[0].split(">")[0].strip() if ">" in smiles else smiles

        # Build smart search query from FGs + reaction name
        query_parts = []
        if reaction_name:
            query_parts.append(reaction_name.replace("_", " "))

        # Extract key FG keywords from molecule
        try:
            from chemtools.util.functional_groups import get_functional_groups
            fgs = get_functional_groups(mol_smiles)
            # Map FG names to useful search terms
            fg_search_map = {
                "aryl_bromide": "aryl bromide coupling",
                "aryl_chloride": "aryl chloride coupling",
                "aryl_iodide": "aryl iodide coupling",
                "boronic_acid": "boronic acid Suzuki",
                "primary_amine": "amine synthesis",
                "secondary_amine": "amine synthesis",
                "aldehyde": "aldehyde synthesis reductive amination",
                "ketone": "ketone synthesis reduction",
                "carboxylic_acid": "acid coupling esterification",
                "ester": "ester synthesis hydrolysis",
                "alcohol": "alcohol oxidation reduction",
                "nitrile": "nitrile reduction hydrolysis",
            }
            for fg in fgs[:4]:  # Use first 4 FGs
                if fg in fg_search_map:
                    query_parts.append(fg_search_map[fg])
                else:
                    query_parts.append(fg.replace("_", " "))
        except Exception:
            pass

        if not query_parts:
            query_parts = ["synthesis", "route", "preparation"]

        search_query = " ".join(query_parts[:6])  # Keep query focused

        # Search notes
        notes_hits = []
        try:
            from chem_coworker.tools.notes import _search_notes
            notes_result = _search_notes(
                query=search_query,
                note_types=["reactions", "routes"],
                top_k=3,
            )
            if notes_result.get("success") and notes_result.get("found", 0) > 0:
                notes_hits = notes_result.get("results", [])
        except Exception:
            pass

        # Search literature
        lit_hits = []
        try:
            from chem_coworker.tools.literature import _search_literature
            lit_result = _search_literature(query=search_query, top_k=3)
            if lit_result.get("success") and lit_result.get("found", 0) > 0:
                lit_hits = lit_result.get("results", [])
        except Exception:
            pass

        # Targeted reaction note lookup
        reaction_note = {}
        if reaction_name:
            try:
                from chem_coworker.tools.notes import _read_notes
                note = _read_notes(id=reaction_name, note_type="reactions", max_chars=3000)
                if note.get("success") and note.get("has_notes"):
                    reaction_note = {
                        "id": note.get("id"),
                        "content_preview": note.get("content", "")[:1500],
                    }
            except Exception:
                pass

        total_found = len(notes_hits) + len(lit_hits) + (1 if reaction_note else 0)

        return _success({
            "smiles": mol_smiles,
            "search_query": search_query,
            "search_terms_used": query_parts[:6],
            "total_found": total_found,
            "reaction_note": reaction_note,
            "notes_hits": notes_hits,
            "literature_hits": lit_hits,
            "message": (
                f"Found {total_found} relevant sources. "
                + ("No precedent found — route may be novel." if total_found == 0 else "")
            ),
        })

    except Exception as exc:
        return _error(f"Precedent search failed: {exc}")


find_retro_precedent_tool = ToolPlugin(
    name="find_retro_precedent",
    category="retrosynthesis",
    description=(
        "Search the reaction notes, literature, and route knowledge base for synthesis precedent "
        "relevant to the target molecule. Uses functional groups and reaction name to build a "
        "smart search query. Args: smiles (target SMILES), reaction_name (e.g. 'suzuki_miyaura'; "
        "NOT reaction_type). Run in parallel with identify_retrons to gather experimental context "
        "before generating disconnections."
    ),
    prerequisites=[],
    fn=_find_retro_precedent,
)


# ---------------------------------------------------------------------------
# Tool E: search_hte_precedent
# DRFP k-NN search in the HTE reaction database (~231k reactions)
# ---------------------------------------------------------------------------

# Retron/reaction name → HTE family string
# Resolved via the unified reaction_registry (canonical taxonomy family ID → hte_families).
# The legacy _EXTRA_RETRON_MAP has been removed; all lookups go through the
# registry which is built from retron/template taxonomy links
# (v2 `taxonomy_family_id` preferred, legacy `taxonomy_id` fallback).


def _map_reaction_to_family(reaction_name: str) -> Optional[str]:
    """Map a retron/reaction name to a canonical HTE family string, or None.

    Resolution order:
    1. Treat *reaction_name* as a retron name → registry → first HTE family.
    2. Treat *reaction_name* as a template name → registry → first HTE family.
    3. Fall back to core_utils normaliser (unchanged behaviour).
    """
    if not reaction_name:
        return None
    # 1. Retron name lookup
    try:
        from chemtools.retro.reaction_registry import get_hte_families_for_retron
        families = get_hte_families_for_retron(reaction_name)
        if families:
            return families[0]
    except Exception:
        pass
    # 2. Template name lookup
    try:
        from chemtools.retro.reaction_registry import get_hte_families_for_template
        families = get_hte_families_for_template(reaction_name)
        if families:
            return families[0]
    except Exception:
        pass
    # 3. Core-utils normaliser fallback
    try:
        from chemtools.precedent.core_utils import _family_text
        mapped = _family_text(reaction_name)
        if mapped != reaction_name:
            return mapped
    except Exception:
        pass
    return None


def _fmt_reagent_list(reagents: List[Any]) -> Optional[str]:
    if not reagents:
        return None
    parts = []
    for r in reagents:
        if isinstance(r, dict):
            name = r.get("name", "")
            role = r.get("role", "")
            if name:
                parts.append(f"{name} ({role})" if role else name)
        elif isinstance(r, str) and r:
            parts.append(r)
    return ", ".join(parts) if parts else None


def _fmt_solvent_list(solvents: List[Any]) -> Optional[str]:
    if not solvents:
        return None
    parts = []
    for s in solvents:
        if isinstance(s, dict):
            name = s.get("name", "")
            if name:
                parts.append(name)
        elif isinstance(s, str) and s:
            parts.append(s)
    return "/".join(parts) if parts else None


def _extract_catalyst_name(obj: Any) -> Optional[str]:
    if obj is None:
        return None
    if isinstance(obj, dict):
        return obj.get("name")
    return str(obj)


def _fast_load_hte_family(family: Optional[str]) -> List[Dict[str, Any]]:
    """Read HTE CSV files directly for a given family, skipping featurization.

    This is a lightweight alternative to _load_selective() that avoids the
    expensive featurize_pair() / normalize_reaction() calls in _make_row_from_csv().
    Returns only the fields needed by search_hte_precedent: yield, conditions,
    reaction_smiles, reference.  Cached per family so the first call is the
    only slow one.
    """
    import csv as _csv
    import os as _os
    try:
        from chemtools.precedent.loader import HTE_DB_DIR, _file_family_from_name, _clean_text, _parse_float, _split_items
    except ImportError:
        return []

    rows: List[Dict[str, Any]] = []
    family_lower = (family or "").lower()

    try:
        from chemtools.precedent.loader import _dataset_family_map as _dfm
    except ImportError:
        _dfm = None  # type: ignore

    for subdir_name in ("literature", "protocols", "rules", "experiments"):
        subdir = _os.path.join(HTE_DB_DIR, subdir_name)
        if not _os.path.isdir(subdir):
            continue
        for fname in _os.listdir(subdir):
            if not fname.lower().endswith(".csv"):
                continue
            file_family = _file_family_from_name(fname)
            # Resolve file stem → canonical family name so e.g.
            # "suzuki_miyaura" (from filename) maps to "Suzuki" (canonical)
            mapped_family = (_dfm(file_family, fallback=file_family)
                             if _dfm else file_family)
            # Apply family filter: match on either the raw stem or the mapped name
            if family and (file_family.lower() != family_lower
                           and mapped_family.lower() != family_lower):
                continue
            path = _os.path.join(subdir, fname)
            for enc in ("utf-8", "latin-1"):
                try:
                    with open(path, encoding=enc, newline="") as fh:
                        reader = _csv.DictReader(fh)
                        for row_rec in reader:
                            cat = _clean_text(row_rec.get("catalyst"))
                            lig = _clean_text(row_rec.get("ligand"))
                            base = _clean_text(row_rec.get("base"))
                            acid = _clean_text(row_rec.get("acid"))
                            add = _clean_text(row_rec.get("additive"))
                            ca = _clean_text(row_rec.get("condensation_agent"))
                            solv = _clean_text(row_rec.get("solvent"))
                            reag = [
                                {"name": n, "role": r}
                                for n, r in [(base, "BASE"), (acid, "ACID"),
                                             (add, "ADDITIVE"), (ca, "COUPLING_REAGENT")]
                                if n
                            ]
                            rows.append({
                                "rxn_type": file_family,
                                "yield_value": _parse_float(
                                    row_rec.get("yield") or row_rec.get("yield_pct")
                                    or row_rec.get("yield_percent")
                                ),
                                "reaction_smiles": _clean_text(row_rec.get("reaction_smiles")),
                                "condition_core": "/".join(p for p in [cat, lig] if p),
                                "catalyst": {"name": cat} if cat else None,
                                "base_uid": base or None,
                                "solvent_uid": _split_items(solv)[0] if solv else None,
                                "reagents": reag,
                                "solvents": [{"name": s} for s in _split_items(solv) if s],
                                "reference": _clean_text(row_rec.get("reference")),
                                "source_file": fname,
                            })
                    break  # successful read
                except Exception:
                    continue

    return rows


# LRU cache keyed by family name so each family is loaded only once per process.
from functools import lru_cache as _lru_cache

@_lru_cache(maxsize=16)
def _fast_load_hte_family_cached(family_key: str) -> tuple:
    """Cached wrapper around _fast_load_hte_family; returns tuple for hashability."""
    rows = _fast_load_hte_family(family_key if family_key != "__ALL__" else None)
    return tuple(rows)   # tuples are hashable → safe for lru_cache key


def _search_hte_precedent(
    target_smiles: str,
    reaction_name: str = "",
    precursor_1: str = "",
    precursor_2: str = "",
    top_k: int = 5,
) -> Dict[str, Any]:
    """Search the HTE reaction database for synthesis precedents via DRFP k-NN.

    Composes a forward reaction SMILES (precursor_1.precursor_2 >> target) and
    uses DRFP fingerprint similarity to rank ~231k HTE literature reactions.
    Returns the nearest neighbors with their experimental conditions so the LLM
    can ground its route recommendations in real data.

    When precursors are not yet known, pass only target_smiles + reaction_name
    to get the top-yielding precedents for that reaction family (no DRFP ranking).

    Strategy:
      - Load the reaction family from the HTE database via selective loading.
      - Pre-rank by yield (descending) so that high-quality reactions are at the top.
      - If precursor SMILES are available: compute DRFP fingerprint for the proposed
        forward reaction (p1.p2 >> target) and re-rank the top-N yield candidates
        by Tanimoto similarity. N defaults to 200, keeping latency ~2 s per call.
      - If no precursors: return top-K by yield directly.

    Args:
        target_smiles: Target molecule SMILES (product of the forward step).
        reaction_name: Reaction family / retron name from identify_retrons or
                       generate_disconnections (e.g. "suzuki_miyaura",
                       "aryl_amine_buchwald", "amide_direct").
        precursor_1: First precursor SMILES from generate_disconnections.
        precursor_2: Second precursor SMILES from generate_disconnections.
        top_k: Number of precedents to return (default 5).

    Returns:
        dict with family, search_mode, support_in_family, precedents[].
        Each precedent has: reaction_smiles, similarity, yield, condition_core,
        catalyst, base, solvent, reagents, solvents, reference, source_file.
    """
    _DRFP_CANDIDATE_CAP = 300   # max rows to DRFP-score (keeps latency ~2 s)
    import time as _time

    try:
        t0 = _time.perf_counter()
        t_parse_0 = _time.perf_counter()
        mol_smiles = (
            target_smiles.split(">>")[0].split(">")[0].strip()
            if ">" in target_smiles
            else target_smiles.strip()
        )

        # Resolve family string
        family = _map_reaction_to_family(reaction_name)
        t_parse_1 = _time.perf_counter()

        # ── Load family rows (fast path: no featurization) ────────────────
        try:
            t_load_0 = _time.perf_counter()
            family_key = family if family else "__ALL__"
            rows = list(_fast_load_hte_family_cached(family_key))
            t_load_1 = _time.perf_counter()
        except Exception as load_err:
            return _error(f"HTE loader failed: {load_err}")

        if not rows:
            t_end = _time.perf_counter()
            timing_ms = {
                "input_parse_ms": round((t_parse_1 - t_parse_0) * 1000, 2),
                "load_family_ms": round((t_load_1 - t_load_0) * 1000, 2),
                "sort_ms": 0.0,
                "drfp_rerank_ms": 0.0,
                "format_ms": 0.0,
                "total_ms": round((t_end - t0) * 1000, 2),
            }
            return _success({
                "family": family,
                "search_mode": "none",
                "forward_smiles_used": "",
                "support_in_family": 0,
                "precedent_count": 0,
                "precedents": [],
                "hte_search_timing_ms": timing_ms,
                "message": (
                    f"No HTE rows loaded for family '{family}'. "
                    "Reaction type may not be in the database."
                ),
            })

        # ── Sort by yield descending (high-quality first) ─────────────────
        t_sort_0 = _time.perf_counter()
        rows_sorted = sorted(
            rows,
            key=lambda r: (r.get("yield_value") or 0.0),
            reverse=True,
        )
        t_sort_1 = _time.perf_counter()
        support = len(rows_sorted)

        # ── Build forward SMILES & decide search mode ──────────────────────
        p1 = (precursor_1 or "").strip()
        p2 = (precursor_2 or "").strip()
        fwd_smiles = ""
        use_drfp = False

        if mol_smiles and (p1 or p2):
            parts = [x for x in [p1, p2] if x]
            fwd_smiles = ".".join(parts) + ">>" + mol_smiles
            use_drfp = True

        # ── DRFP re-ranking (when precursors available) ───────────────────
        t_drfp_0 = _time.perf_counter()
        scored = []
        if use_drfp:
            try:
                from chemtools import reaction_similarity as rs
                if rs and rs.drfp_available():
                    q_fp = rs.encode_drfp_cached(fwd_smiles)
                    if q_fp is not None:
                        # Only DRFP-score the top-N by yield to keep latency low
                        candidates = rows_sorted[:_DRFP_CANDIDATE_CAP]
                        for r in candidates:
                            rsmi = r.get("reaction_smiles") or ""
                            if rsmi:
                                r_fp = rs.encode_drfp_cached(rsmi)
                                if r_fp is not None:
                                    sim = float(rs.tanimoto(q_fp, r_fp))
                                    # Blend DRFP sim with yield to prefer high-yield analogs
                                    y_norm = (r.get("yield_value") or 0.0) / 100.0
                                    blended = 0.85 * sim + 0.15 * y_norm
                                    scored.append((blended, sim, r))
                        if scored:
                            scored.sort(key=lambda x: -x[0])
                            rows_sorted = [r for _, _, r in scored]
                        # else: keep yield-sorted order as fallback
            except Exception:
                use_drfp = False   # silently fall back to yield ranking
        t_drfp_1 = _time.perf_counter()

        # ── Format top-K results ──────────────────────────────────────────
        t_fmt_0 = _time.perf_counter()
        formatted = []
        scored_lookup = (
            {id(r): sim for _, sim, r in scored}
            if use_drfp and scored
            else {}
        )

        for r in rows_sorted[:top_k]:
            entry: Dict[str, Any] = {
                "reaction_smiles": r.get("reaction_smiles") or "",
                "yield": r.get("yield_value"),
                "condition_core": r.get("condition_core"),
                "catalyst": _extract_catalyst_name(r.get("catalyst")),
                "base": r.get("base_uid"),
                "solvent": r.get("solvent_uid"),
                "reagents": _fmt_reagent_list(r.get("reagents") or []),
                "solvents": _fmt_solvent_list(r.get("solvents") or []),
                "reference": r.get("reference"),
                "rxn_type": r.get("rxn_type"),
                "source_file": r.get("source_file"),
            }
            sim_val = scored_lookup.get(id(r))
            if sim_val is not None:
                entry["drfp_similarity"] = round(sim_val, 4)
            formatted.append({k: v for k, v in entry.items() if v is not None})
        t_fmt_1 = _time.perf_counter()
        t_end = _time.perf_counter()
        timing_ms = {
            "input_parse_ms": round((t_parse_1 - t_parse_0) * 1000, 2),
            "load_family_ms": round((t_load_1 - t_load_0) * 1000, 2),
            "sort_ms": round((t_sort_1 - t_sort_0) * 1000, 2),
            "drfp_rerank_ms": round((t_drfp_1 - t_drfp_0) * 1000, 2),
            "format_ms": round((t_fmt_1 - t_fmt_0) * 1000, 2),
            "total_ms": round((t_end - t0) * 1000, 2),
        }

        search_mode = "drfp_yield_blend" if use_drfp else "family_yield"
        return _success({
            "family": family,
            "search_mode": search_mode,
            "forward_smiles_used": fwd_smiles or "(family-level only)",
            "support_in_family": support,
            "drfp_candidates_scored": min(support, _DRFP_CANDIDATE_CAP) if use_drfp else 0,
            "precedent_count": len(formatted),
            "precedents": formatted,
            "hte_search_timing_ms": timing_ms,
            "message": (
                f"Found {len(formatted)} HTE precedents from {support} reactions "
                f"in '{family}' family. "
                + (
                    f"DRFP-ranked (vs top-{min(support, _DRFP_CANDIDATE_CAP)} by yield)."
                    if use_drfp
                    else "Ranked by yield (pass precursor_1/precursor_2 for DRFP ranking)."
                )
            ),
        })

    except Exception as exc:
        return _error(f"HTE precedent search failed: {exc}")


search_hte_precedent_tool = ToolPlugin(
    name="search_hte_precedent",
    category="retrosynthesis",
    description=(
        "Search the ~231k-reaction HTE literature database for synthesis precedents using "
        "DRFP fingerprint k-NN similarity. Takes the target SMILES plus precursor SMILES "
        "from generate_disconnections, composes the forward reaction (p1.p2>>target), and "
        "returns the most similar real reactions with conditions (catalyst, base, solvent, "
        "yield, reference). Run in G3 alongside recommend_conditions. "
        "If precursors are not yet known, pass only target_smiles + reaction_name for "
        "family-level yield ranking."
    ),
    prerequisites=["generate_disconnections"],
    fn=_search_hte_precedent,
)


# ---------------------------------------------------------------------------
# Tool F: search_by_product_similarity
# Product-space Morgan fingerprint search — data-driven retrosynthesis
# Searches HTE products directly: "who made something like this target?"
# No retron patterns needed; runs before generate_disconnections.
# ---------------------------------------------------------------------------

# Module-level product FP index: family_key → {"rows": list, "fps": ndarray}
# Populated lazily on first call; persists for the process lifetime.
_PRODUCT_FP_INDEX: Dict[str, Any] = {}


def _get_morgan_fp(smiles: str) -> "Optional[Any]":
    """Compute Morgan fingerprint (radius=2, 2048 bits) as a numpy uint8 array."""
    try:
        from rdkit import Chem, rdBase
        from rdkit.Chem import rdFingerprintGenerator
        from rdkit.DataStructs import ConvertToNumpyArray
        import numpy as np
        with rdBase.BlockLogs():
            mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        fp = gen.GetFingerprint(mol)
        arr = np.zeros(2048, dtype=np.uint8)
        ConvertToNumpyArray(fp, arr)
        return arr
    except Exception:
        return None


def _build_product_fp_index(family_key: str) -> Dict[str, Any]:
    """Build (or retrieve) the product Morgan FP index for a reaction family.

    Parses the product SMILES from each HTE reaction SMILES, computes a Morgan
    fingerprint, and stacks the results into a (N, 2048) uint8 matrix for fast
    vectorised Tanimoto computation.  Results are cached in _PRODUCT_FP_INDEX
    so the expensive build only happens once per family per process.

    For '__ALL__', caps at the top-20k rows by yield to keep memory < ~5 MB
    and first-call latency < ~3 s.
    """
    if family_key in _PRODUCT_FP_INDEX:
        return _PRODUCT_FP_INDEX[family_key]

    try:
        import numpy as np

        rows = list(_fast_load_hte_family_cached(family_key))

        # Cap the all-families search so FP matrix stays manageable
        _ALL_CAP = 20_000
        if family_key == "__ALL__" and len(rows) > _ALL_CAP:
            rows = sorted(rows, key=lambda r: r.get("yield_value") or 0.0, reverse=True)
            rows = rows[:_ALL_CAP]

        valid_rows: List[Dict[str, Any]] = []
        fp_list = []

        for r in rows:
            rxn_smi = r.get("reaction_smiles") or ""
            if ">>" not in rxn_smi:
                continue
            # Take only the first listed product (before any "." in the product block)
            product_smi = rxn_smi.split(">>")[-1].strip().split(".")[0]
            if not product_smi:
                continue
            fp = _get_morgan_fp(product_smi)
            if fp is not None:
                valid_rows.append(r)
                fp_list.append(fp)

        fp_matrix = np.stack(fp_list, axis=0) if fp_list else np.zeros((0, 2048), dtype=np.uint8)
        index: Dict[str, Any] = {"rows": valid_rows, "fps": fp_matrix}
        _PRODUCT_FP_INDEX[family_key] = index
        return index
    except Exception:
        return {"rows": [], "fps": None}


def _search_by_product_similarity(
    target_smiles: str = "",
    smiles: str = "",
    product_smiles: str = "",
    reaction_name: str = "",
    top_k: int = 5,
) -> Dict[str, Any]:
    """Search HTE reactions by product structural similarity to the target molecule.

    Unlike search_hte_precedent (which needs precursor SMILES from
    generate_disconnections), this tool searches purely from the target:
    "who has made something structurally similar to this target, and how?"

    Uses Morgan fingerprint (radius=2, 2048 bits) Tanimoto similarity on the
    product side of HTE reactions.  The 'precursor_1'/'precursor_2' fields in
    the results are the actual reactants used — compare them against
    generate_disconnections output for cross-validation.

    Args:
        target_smiles: Target molecule SMILES (aliases: smiles, product_smiles).
        smiles: Alias for target_smiles.
        product_smiles: Alias for target_smiles.
        reaction_name: Optional family filter (e.g. "amide_coupling",
                       "suzuki_miyaura"). Leave empty to search all families.
        top_k: Number of results to return (default 5).

    Returns:
        dict with family, search_mode, support_searched, precedents[].
        Each precedent has: product_similarity, reaction_smiles,
        precursor_1, precursor_2, yield, condition_core, catalyst, base,
        solvent, reagents, solvents, reference, rxn_type.
    """
    try:
        import numpy as np

        mol_smiles = (target_smiles or smiles or product_smiles).strip()
        if ">" in mol_smiles:
            mol_smiles = mol_smiles.split(">>")[0].split(">")[0].strip()

        if not mol_smiles:
            return _error("target_smiles is required")

        query_fp = _get_morgan_fp(mol_smiles)
        if query_fp is None:
            return _error(f"Could not parse SMILES: {mol_smiles!r}")

        family = _map_reaction_to_family(reaction_name)
        family_key = family if family else "__ALL__"


        index = _build_product_fp_index(family_key)
        rows: List[Dict[str, Any]] = index["rows"]
        fp_matrix = index.get("fps")

        if not rows or fp_matrix is None or fp_matrix.shape[0] == 0:
            return _success({
                "family": family or "all",
                "search_mode": "product_morgan",
                "support_searched": 0,
                "precedent_count": 0,
                "precedents": [],
                "message": (
                    f"No indexed products found for family '{family or 'all'}'. "
                    "Try a different reaction_name or leave it empty."
                ),
            })

        # Vectorised Tanimoto: (N, 2048) & (2048,) → (N,)
        intersections = np.sum(fp_matrix & query_fp, axis=1).astype(np.float32)
        unions = np.sum(fp_matrix | query_fp, axis=1).astype(np.float32)
        tanimotos = intersections / np.maximum(unions, 1.0)

        # Blend with yield (80/20) to prefer high-yield similar reactions
        yields_norm = np.array(
            [float(r.get("yield_value") or 0.0) / 100.0 for r in rows],
            dtype=np.float32,
        )
        blended = 0.80 * tanimotos + 0.20 * yields_norm

        top_indices = np.argsort(blended)[::-1][:top_k]

        formatted: List[Dict[str, Any]] = []
        for idx in top_indices:
            r = rows[int(idx)]
            sim = float(tanimotos[idx])
            rxn_smi = r.get("reaction_smiles") or ""

            # Parse reactants from reaction SMILES (reactants>>product)
            precursor_1 = ""
            precursor_2 = ""
            if ">>" in rxn_smi:
                reactant_block = rxn_smi.split(">>")[0]
                parts = [x.strip() for x in reactant_block.split(".") if x.strip()]
                if parts:
                    precursor_1 = parts[0]
                if len(parts) >= 2:
                    precursor_2 = ".".join(parts[1:])

            entry: Dict[str, Any] = {
                "product_similarity": round(sim, 4),
                "reaction_smiles": rxn_smi,
                "precursor_1": precursor_1,
                "precursor_2": precursor_2,
                "yield": r.get("yield_value"),
                "condition_core": r.get("condition_core"),
                "catalyst": _extract_catalyst_name(r.get("catalyst")),
                "base": r.get("base_uid"),
                "solvent": r.get("solvent_uid"),
                "reagents": _fmt_reagent_list(r.get("reagents") or []),
                "solvents": _fmt_solvent_list(r.get("solvents") or []),
                "reference": r.get("reference"),
                "rxn_type": r.get("rxn_type"),
            }
            formatted.append({k: v for k, v in entry.items() if v is not None})

        return _success({
            "family": family or "all",
            "search_mode": "product_morgan",
            "support_searched": len(rows),
            "precedent_count": len(formatted),
            "precedents": formatted,
            "message": (
                f"Found {len(formatted)} HTE precedents from {len(rows):,} reactions "
                f"in '{family or 'all families'}' via Morgan product-similarity search. "
                "precursor_1/precursor_2 are the actual reactants — use them as "
                "data-driven disconnection candidates."
            ),
        })

    except Exception as exc:
        return _error(f"Product similarity search failed: {exc}")


search_by_product_similarity_tool = ToolPlugin(
    name="search_by_product_similarity",
    category="retrosynthesis",
    description=(
        "Data-driven retrosynthesis: search the ~231k-reaction HTE database for reactions "
        "where the PRODUCT is structurally similar to the target molecule (Morgan FP "
        "Tanimoto, radius=2). No retron patterns needed — runs BEFORE or IN PARALLEL with "
        "identify_retrons. Returns real precedents with actual reactants as candidate "
        "precursor pairs, plus conditions and yield. Especially powerful for complex, "
        "multi-functional molecules where SMARTS retrons may not fire. "
        "Pass reaction_name to restrict to a family; leave empty to search all 231k reactions. "
        "Parameter: target_smiles (also accepted as smiles or product_smiles)."
    ),
    prerequisites=[],  # Runs in G1 alongside identify_retrons
    fn=_search_by_product_similarity,
)


# ---------------------------------------------------------------------------
# Tool G: apply_hte_templates
# Apply atom-mapped retrosynthetic SMARTS templates from the HTE template library.
# Covers 35+ HTE reaction families not fully handled by the 46 hardcoded retrons.
# Uses AllChem.RunReactants() for atom-precise precursor SMILES generation.
# Enriches each hit with sample HTE conditions from the 231k-reaction database.
# ---------------------------------------------------------------------------

def _apply_one_template(
    target_mol: "Any",
    retro_smarts: str,
    n_precursors: int,
) -> "List[tuple]":
    """Apply a single retrosynthetic SMARTS template to the target molecule.

    Returns a list of precursor SMILES tuples (each with n_precursors entries).
    Invalid or duplicate results are filtered; at most 8 results returned.
    """
    try:
        from rdkit import Chem, rdBase
        from rdkit.Chem import AllChem

        with rdBase.BlockLogs():
            rxn = AllChem.ReactionFromSmarts(retro_smarts)
        if rxn is None:
            return []

        with rdBase.BlockLogs():
            results = rxn.RunReactants((target_mol,))

        valid_pairs: List[tuple] = []
        seen: set = set()

        for product_tuple in results:
            if len(product_tuple) != n_precursors:
                continue
            smiles_parts: List[str] = []
            valid = True
            for mol in product_tuple:
                if mol is None:
                    valid = False
                    break
                try:
                    # Use MolToSmiles directly — SanitizeMol can reject valid
                    # charged species (e.g. azides) produced by RunReactants.
                    # Round-trip through SMILES parser validates the structure.
                    with rdBase.BlockLogs():
                        smi = Chem.MolToSmiles(mol)
                    if not smi:
                        valid = False
                        break
                    # Validate round-trip
                    with rdBase.BlockLogs():
                        check_mol = Chem.MolFromSmiles(smi)
                    if check_mol is None:
                        valid = False
                        break
                    smiles_parts.append(smi)
                except Exception:
                    valid = False
                    break
            if valid and len(smiles_parts) == n_precursors:
                key = tuple(sorted(smiles_parts))
                if key not in seen:
                    seen.add(key)
                    valid_pairs.append(tuple(smiles_parts))
                    if len(valid_pairs) >= 8:
                        break

        return valid_pairs

    except Exception:
        return []


def _apply_hte_templates(
    target_smiles: str = "",
    smiles: str = "",
    reaction_name: str = "",
    top_k: int = 5,
) -> Dict[str, Any]:
    """Apply HTE-backed retrosynthetic SMARTS templates to generate precursor pairs.

    Unlike identify_retrons (which uses 46 hardcoded retrons) or
    generate_disconnections (which applies transforms from the retron library),
    this tool applies the *HTE template library* — 35+ atom-mapped retrosynthetic
    SMARTS tuned to the HTE reaction families that are NOT fully covered by the
    existing retron patterns: SNAr, Chan-Lam, CuAAC, HWE, Knoevenagel, Wacker,
    reductions, deoxyfluorination, Sandmeyer, Giese radical, and more.

    Each successful template application returns:
    - Validated precursor SMILES from AllChem.RunReactants()
    - HTE family name(s) cross-referenced to the 231k-reaction database
    - Sample HTE conditions from the matching family (best-yield row)
    - Difficulty score and chemistry notes from the template

    Run in G1 (no prerequisites) alongside identify_retrons and
    search_by_product_similarity for maximum coverage.

    Args:
        target_smiles: Target molecule SMILES (alias: smiles).
        smiles: Alias for target_smiles.
        reaction_name: Optional filter — only apply templates whose hte_families
                       or name matches this string (e.g. "SNAr_amination",
                       "CuAAC", "michael_addition"). Leave empty to try all 35+.
        top_k: Maximum number of template hits to return (default 5).

    Returns:
        dict with disconnections[]: template_name, hte_families, precursor_1,
        precursor_2, reaction_smiles, description, difficulty, notes,
        hte_conditions (sample from database).
    """
    try:
        from rdkit import Chem, rdBase  # noqa: F811

        mol_smiles = _normalize_target_smiles_for_route(target_smiles=target_smiles, smiles=smiles)
        if not mol_smiles:
            return _error("target_smiles is required")

        with rdBase.BlockLogs():
            target_mol = Chem.MolFromSmiles(mol_smiles)
        if target_mol is None:
            return _error(f"Cannot parse target SMILES: {mol_smiles!r}")

        from chemtools.retro.hte_templates import HTE_TEMPLATES

        # Optionally filter templates by reaction_name / hte_family
        templates_to_try = HTE_TEMPLATES
        if reaction_name:
            rn_lower = reaction_name.lower()
            templates_to_try = [
                t for t in HTE_TEMPLATES
                if t["name"].lower() == rn_lower
                or any(rn_lower in f.lower() for f in t.get("hte_families", []))
            ]

        # Collect ALL matching templates first, then sort and cap.
        # Early stopping by top_k would skip low-difficulty templates that come
        # later in the library order (e.g. CuAAC at position 31 with difficulty=0.10).
        hits: List[Dict[str, Any]] = []

        for tmpl in templates_to_try:
            n_prec = tmpl.get("n_precursors", 2)
            pairs = _apply_one_template(target_mol, tmpl["retro_smarts"], n_prec)
            if not pairs:
                continue

            # Take the first (most canonical) precursor pair
            pair = pairs[0]
            p1 = pair[0] if pair else ""
            p2 = pair[1] if n_prec >= 2 and len(pair) >= 2 else ""
            if p2:
                rxn_smi = f"{p1}.{p2}>>{mol_smiles}"
            else:
                rxn_smi = f"{p1}>>{mol_smiles}"

            # Fetch sample HTE conditions from the database for the first matching family
            hte_conditions: Optional[Dict[str, Any]] = None
            for family_name in tmpl.get("hte_families", []):
                try:
                    rows = list(_fast_load_hte_family_cached(family_name))
                    if rows:
                        best = max(rows, key=lambda r: r.get("yield_value") or 0.0)
                        cond: Dict[str, Any] = {
                            "family": family_name,
                            "catalyst": _extract_catalyst_name(best.get("catalyst")),
                            "base": best.get("base_uid"),
                            "solvent": best.get("solvent_uid"),
                            "reagents": _fmt_reagent_list(best.get("reagents") or []),
                            "yield": best.get("yield_value"),
                            "reference": best.get("reference"),
                            "source_file": best.get("source_file"),
                        }
                        hte_conditions = {k: v for k, v in cond.items() if v is not None}
                        break
                except Exception:
                    pass

            taxonomy_family_id = (
                tmpl.get("taxonomy_family_id")
                or tmpl.get("taxonomy_id", "")
            )
            hit: Dict[str, Any] = {
                "template_name": tmpl["name"],
                "taxonomy_id": taxonomy_family_id,
                "taxonomy_family_id": taxonomy_family_id,
                "hte_families": tmpl.get("hte_families", []),
                "precursor_1": p1,
                "precursor_2": p2,
                "reaction_smiles": rxn_smi,
                "description": tmpl.get("description", ""),
                "difficulty": tmpl.get("difficulty", 0.5),
                "notes": tmpl.get("notes", ""),
            }
            if hte_conditions:
                hit["hte_conditions"] = hte_conditions
            hits.append(hit)

        # Sort by difficulty (easiest first) and cap to top_k
        hits.sort(key=lambda h: h["difficulty"])
        hits = hits[:top_k]

        has_conditions = any(h.get("hte_conditions") for h in hits)
        return _success({
            "smiles": mol_smiles,
            "templates_tried": len(templates_to_try),
            "template_hits": len(hits),
            "disconnections": hits,
            "message": (
                f"Applied {len(templates_to_try)} HTE templates: "
                f"{len(hits)} matched the target structure. "
                + (
                    "HTE conditions from database included for each hit."
                    if has_conditions
                    else "No HTE family conditions found (families may not exist in database)."
                )
            ),
        })

    except Exception as exc:
        return _error(f"HTE template application failed: {exc}")


apply_hte_templates_tool = ToolPlugin(
    name="apply_hte_templates",
    category="retrosynthesis",
    description=(
        "Apply the HTE-backed retrosynthetic template library (35+ atom-mapped SMARTS) to "
        "generate validated precursor SMILES via AllChem.RunReactants(). Covers HTE families "
        "NOT in the standard 46 retrons: SNAr amination, Chan-Lam N-arylation, CuAAC triazole, "
        "Knoevenagel, HWE olefination, Wacker oxidation, NaBH4/LAH reductions, deoxyfluorination, "
        "Sandmeyer, Giese radical addition, sulfonamide, urea, carbamate, C–S coupling, and more. "
        "Each hit returns precursor_1, precursor_2, reaction_smiles, and sample HTE conditions "
        "from the database. Run in G1 alongside identify_retrons and search_by_product_similarity "
        "for maximum coverage. Pass reaction_name to filter to a specific HTE family."
    ),
    prerequisites=[],  # Runs in G1 alongside identify_retrons
    fn=_apply_hte_templates,
)


# ---------------------------------------------------------------------------
# Tool H: plan_route
# Multi-step greedy BFS retrosynthesis — AiZynthFinder-inspired.
# Recursively disconnects each complex precursor until all fragments are simple
# (BertzCT < threshold) or max_depth is reached.
# Cycle detection via InChI key blacklist prevents infinite loops.
# ---------------------------------------------------------------------------


def _bertz_complexity_fast(smiles: str) -> float:
    """Compute BertzCT complexity for a SMILES string. Returns 9999.0 on failure."""
    try:
        from rdkit import Chem, rdBase
        from rdkit.Chem import GraphDescriptors
        with rdBase.BlockLogs():
            mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 9999.0
        return float(GraphDescriptors.BertzCT(mol))
    except Exception:
        return 9999.0


def _inchi_key(smiles: str) -> Optional[str]:
    """Return the InChI key for a SMILES, or None on failure."""
    try:
        from rdkit import Chem, rdBase
        with rdBase.BlockLogs():
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            # MolToInchiKey is available directly on Chem in modern RDKit
            return Chem.MolToInchiKey(mol)
    except Exception:
        return None


def _normalize_target_smiles_for_route(target_smiles: str = "", smiles: str = "") -> str:
    """Normalize route-planning target input while preserving existing tool behavior."""
    mol_smiles = (target_smiles or smiles).strip()
    if ">" in mol_smiles:
        mol_smiles = mol_smiles.split(">>")[0].split(">")[0].strip()
    return mol_smiles


def _route_score_and_grade(
    *,
    steps: List[Dict[str, Any]],
    leaves: List[str],
    simple_leaves: List[str],
    total_difficulty: float,
) -> tuple[float, str]:
    """Compute the shared route-level score used by single and candidate planners."""
    if not steps:
        return 1.0, "A"

    n_steps = len(steps)
    avg_diff = total_difficulty / max(n_steps, 1)
    diff_score = max(0.0, 1.0 - avg_diff)
    max_step_diff = max(float(s.get("difficulty", 1.0) or 1.0) for s in steps)
    bottleneck_score = max(0.0, 1.0 - max_step_diff)
    completeness = len(simple_leaves) / max(len(leaves), 1)
    step_penalty = max(0.0, 1.0 - max(0, n_steps - 3) * 0.15)
    max_depth_actual = max(int(s.get("depth", 0) or 0) for s in steps) + 1
    convergence = min(1.0, n_steps / max(max_depth_actual, 1)) * 0.1

    route_score = round(
        0.30 * diff_score
        + 0.20 * bottleneck_score
        + 0.30 * completeness
        + 0.15 * step_penalty
        + 0.05 + convergence,
        3,
    )
    if route_score >= 0.80:
        route_grade = "A"
    elif route_score >= 0.65:
        route_grade = "B"
    elif route_score >= 0.50:
        route_grade = "C"
    elif route_score >= 0.35:
        route_grade = "D"
    else:
        route_grade = "F"
    return route_score, route_grade


def _route_candidate_signature(steps: List[Dict[str, Any]]) -> str:
    """Stable signature for de-duplicating equivalent route candidates."""
    parts: List[str] = []
    for step in steps:
        parts.append(
            "|".join(
                [
                    str(step.get("target", "")),
                    str(step.get("reaction_name", "")),
                    str(step.get("precursor_1", "")),
                    str(step.get("precursor_2", "")),
                ]
            )
        )
    return "||".join(parts)


def _route_strategy_flags(route: List[Dict[str, Any]]) -> Dict[str, int]:
    """Extract route-level strategic flags from generated disconnection metadata."""
    text = " ".join(
        str(step.get(key, ""))
        for step in route
        for key in ("reaction_name", "retron_name", "description")
    ).lower()
    metal_terms = (
        "suzuki", "buchwald", "heck", "negishi", "stille", "kumada",
        "sonogashira", "ullmann", "chan-lam", "pd", "palladium",
        "nickel", "copper",
    )
    protection_terms = (
        "protect", "deprotect", "boc", "cbz", "fmoc", "tbs", "tbdms",
        "benzyl protecting", "silyl",
    )
    ring_terms = ("cyclization", "annulation", "ring", "macrocycl", "lactam")
    return {
        "metal_steps": sum(1 for term in metal_terms if term in text),
        "protecting_group_steps": sum(1 for term in protection_terms if term in text),
        "ring_forming_steps": sum(1 for term in ring_terms if term in text),
    }


def _strategy_alignment_score(strategy_query: str, route: List[Dict[str, Any]], route_score: float) -> Dict[str, Any]:
    """Lightweight route-to-query scoring for pre-ranking candidates before LLM review."""
    query = str(strategy_query or "").lower()
    flags = _route_strategy_flags(route)
    score = 0.5
    notes: List[str] = []

    if any(term in query for term in ("feasible", "robust", "high yield", "high-yield", "scalable")):
        score += 0.25 * route_score
        notes.append("query asks for feasibility/robustness; deterministic route score weighted")
    if any(term in query for term in ("short", "concise", "few step", "few-step", "minimum step")):
        step_count = len(route)
        score += 0.18 if step_count <= 3 else -0.10
        notes.append("query asks for a concise route")
    if any(term in query for term in ("convergent", "convergence")):
        max_depth = max((int(s.get("depth", 0) or 0) for s in route), default=0) + 1
        convergence_ratio = len(route) / max(max_depth, 1)
        score += min(0.16, 0.08 * convergence_ratio)
        notes.append("query asks for convergent disconnections")
    if any(term in query for term in ("avoid protecting", "no protecting", "without protecting", "protecting group free")):
        penalty = min(0.25, 0.08 * flags["protecting_group_steps"])
        score -= penalty
        notes.append("query discourages protecting-group operations")
    if any(term in query for term in ("avoid metal", "metal-free", "no metal", "without metal")):
        penalty = min(0.25, 0.05 * flags["metal_steps"])
        score -= penalty
        notes.append("query discourages transition-metal chemistry")
    if any(term in query for term in ("ring early", "early ring", "form ring early", "construct ring early")):
        first_ring_depth = None
        for step in route:
            if _route_strategy_flags([step])["ring_forming_steps"]:
                first_ring_depth = int(step.get("depth", 99) or 99)
                break
        if first_ring_depth is not None:
            score += 0.14 if first_ring_depth <= 1 else -0.06
            notes.append("query asks about early ring construction")

    score = round(max(0.0, min(1.0, score)), 3)
    return {
        "strategy_alignment_score": score,
        "strategy_flags": flags,
        "strategy_notes": notes,
    }


def _plan_route(
    target_smiles: str = "",
    smiles: str = "",
    max_depth: int = 4,
    max_branches: int = 2,
    complexity_threshold: float = 100.0,
    max_difficulty: float = 0.8,
) -> Dict[str, Any]:
    """Plan a complete multi-step retrosynthetic route via greedy BFS.

    Recursively applies retrosynthetic disconnections to each complex precursor
    until all leaf nodes are simple building blocks (BertzCT < complexity_threshold)
    or max_depth is reached.  InChI key cycle detection prevents infinite loops.

    At each node the best-scored disconnection (from rank_disconnections) is chosen
    and both precursors are expanded.  The route is returned as an ordered list of
    single-step dicts so the LLM can summarise or present a forward synthesis plan.

    Inspired by AiZynthFinder's MCTS state expansion heuristic, adapted for a
    deterministic greedy agent context.

    Args:
        target_smiles: Target molecule SMILES (alias: smiles).
        smiles: Alias for target_smiles.
        max_depth: Maximum retrosynthetic depth (default 4).
        max_branches: Top-K disconnections considered at each node (default 2;
                      the best-scored one is always chosen — this controls how
                      many alternatives rank_disconnections generates).
        complexity_threshold: BertzCT score below which a fragment is treated as
                              a simple / purchasable building block and recursion
                              stops (default 100.0 ≈ simple aromatic or aliphatic
                              fragment with ≤ ~8 heavy atoms).
        max_difficulty: Maximum disconnection difficulty passed to
                        rank_disconnections (default 0.8).

    Returns:
        dict with route[], total_steps, cumulative_difficulty, all_leaves_simple,
        leaves[], simple_leaves[], complex_leaves[], route_summary (text).
    """
    try:
        from chemtools.retro.disconnector import rank_disconnections

        mol_smiles = (target_smiles or smiles).strip()
        if ">" in mol_smiles:
            mol_smiles = mol_smiles.split(">>")[0].split(">")[0].strip()
        if not mol_smiles:
            return _error("target_smiles is required")

        steps: List[Dict[str, Any]] = []
        leaves: List[str] = []
        simple_leaves: List[str] = []
        complex_leaves: List[str] = []

        def _expand(
            smi: str,
            depth: int,
            visited_inchi: set,
            cumulative_diff: float,
        ) -> float:
            """Recursively expand smi; returns updated cumulative difficulty."""
            complexity = _bertz_complexity_fast(smi)

            # Stop: simple enough or too deep
            if complexity < complexity_threshold or depth >= max_depth:
                leaves.append(smi)
                if complexity < complexity_threshold:
                    simple_leaves.append(smi)
                else:
                    complex_leaves.append(smi)
                return cumulative_diff

            # Cycle detection via InChI key
            ik = _inchi_key(smi)
            if ik and ik in visited_inchi:
                leaves.append(smi)
                complex_leaves.append(smi)
                return cumulative_diff
            new_visited = visited_inchi | ({ik} if ik else set())

            # Generate disconnections
            disconnections = rank_disconnections(
                smi, top_k=max_branches, max_difficulty=max_difficulty
            )
            if not disconnections:
                leaves.append(smi)
                complex_leaves.append(smi)
                return cumulative_diff

            # Greedy: pick the best-scored disconnection
            best = disconnections[0]
            step: Dict[str, Any] = {
                "step": len(steps) + 1,
                "depth": depth,
                "target": smi,
                "target_complexity": round(complexity, 1),
                "reaction_name": best.reaction_name,
                "retron_name": best.retron_name,
                "description": best.description,
                "difficulty": best.difficulty,
                "precursor_1": best.precursor_1,
                "precursor_2": best.precursor_2,
                "complexity_delta": round(best.complexity_delta, 1),
                "overall_score": round(best.overall_score, 3),
            }
            steps.append(step)
            cumulative_diff += best.difficulty

            # Recurse into both precursors
            for prec in [best.precursor_1, best.precursor_2]:
                if prec:
                    cumulative_diff = _expand(prec, depth + 1, new_visited, cumulative_diff)
            return cumulative_diff

        # Start with empty visited set — the root is added by _expand itself
        # (seeding the root would immediately trigger its own cycle check).
        total_difficulty = _expand(
            mol_smiles, depth=0, visited_inchi=set(), cumulative_diff=0.0
        )

        all_simple = len(complex_leaves) == 0

        # Build human-readable route summary
        if steps:
            lines = [f"Route to {mol_smiles} — {len(steps)} disconnection step(s):"]
            for s in steps:
                p2_part = f" + {s['precursor_2']}" if s.get("precursor_2") else ""
                lines.append(
                    f"  Step {s['step']} (depth {s['depth']}): {s['target']}  →"
                    f"  {s['precursor_1']}{p2_part}"
                    f"  via {s['reaction_name']} (difficulty={s['difficulty']:.2f})"
                )
            bb_list = ", ".join(leaves[:12]) + ("..." if len(leaves) > 12 else "")
            lines.append(f"Building blocks ({len(leaves)}): {bb_list}")
            lines.append(
                f"Cumulative difficulty: {total_difficulty:.2f}  |  "
                f"All leaves simple: {all_simple}"
            )
            route_summary = "\n".join(lines)
        else:
            route_summary = (
                f"{mol_smiles} is already below the complexity threshold "
                f"(BertzCT < {complexity_threshold}) — treat as a purchasable building block."
            )

        # ── Route-level composite score (0..1, higher = better) ───────
        route_score, route_grade = _route_score_and_grade(
            steps=steps,
            leaves=leaves,
            simple_leaves=simple_leaves,
            total_difficulty=total_difficulty,
        )

        return _success({
            "smiles": mol_smiles,
            "total_steps": len(steps),
            "cumulative_difficulty": round(total_difficulty, 3),
            "route_score": route_score,
            "route_grade": route_grade,
            "all_leaves_simple": all_simple,
            "complexity_threshold": complexity_threshold,
            "max_depth": max_depth,
            "leaves": leaves,
            "simple_leaves": simple_leaves,
            "complex_leaves": complex_leaves,
            "route": steps,
            "route_summary": route_summary,
        })

    except Exception as exc:
        return _error(f"Route planning failed: {exc}")


def _plan_route_candidates(
    target_smiles: str = "",
    smiles: str = "",
    strategy_query: str = "",
    max_depth: int = 4,
    max_branches: int = 4,
    beam_width: int = 8,
    top_k: int = 5,
    complexity_threshold: float = 100.0,
    max_difficulty: float = 0.9,
) -> Dict[str, Any]:
    """Generate multiple route candidates for strategy-aware LLM evaluation.

    Unlike plan_route, which greedily follows the top-ranked disconnection at
    each node, this tool keeps a beam of plausible partial routes and returns
    several complete candidate routes.  The LLM should compare these candidates
    against the user's synthesis strategy before presenting a final answer.

    Args:
        target_smiles: Target molecule SMILES (alias: smiles).
        smiles: Alias for target_smiles.
        strategy_query: User's natural-language route preference or full query.
        max_depth: Maximum retrosynthetic depth.
        max_branches: Number of disconnections retained at each expansion.
        beam_width: Number of partial routes retained during search.
        top_k: Number of final candidates returned.
        complexity_threshold: BertzCT below which a fragment is treated as simple.
        max_difficulty: Maximum disconnection difficulty passed to rank_disconnections.

    Returns:
        dict with candidates[], best_candidate, and search metadata.
    """
    try:
        from chemtools.retro.disconnector import rank_disconnections

        mol_smiles = _normalize_target_smiles_for_route(target_smiles=target_smiles, smiles=smiles)
        if not mol_smiles:
            return _error("target_smiles is required")

        max_depth = max(0, min(int(max_depth or 4), 8))
        max_branches = max(1, min(int(max_branches or 4), 8))
        beam_width = max(1, min(int(beam_width or 8), 24))
        top_k = max(1, min(int(top_k or 5), 12))
        complexity_threshold = float(complexity_threshold or 100.0)
        max_difficulty = max(0.0, min(float(max_difficulty or 0.9), 1.0))

        root_state = {
            "route": [],
            "pending": [{"smiles": mol_smiles, "depth": 0}],
            "leaves": [],
            "simple_leaves": [],
            "complex_leaves": [],
            "visited_inchi": set(),
            "cumulative_difficulty": 0.0,
        }
        states: List[Dict[str, Any]] = [root_state]
        expansions = 0
        max_expansions = max(1, beam_width * max(1, max_depth + 1) * 2)

        def _copy_state(state: Dict[str, Any]) -> Dict[str, Any]:
            return {
                "route": [dict(step) for step in state["route"]],
                "pending": [dict(item) for item in state["pending"]],
                "leaves": list(state["leaves"]),
                "simple_leaves": list(state["simple_leaves"]),
                "complex_leaves": list(state["complex_leaves"]),
                "visited_inchi": set(state["visited_inchi"]),
                "cumulative_difficulty": float(state["cumulative_difficulty"]),
            }

        def _interim_sort_key(state: Dict[str, Any]) -> tuple[float, float, int]:
            score, _grade = _route_score_and_grade(
                steps=state["route"],
                leaves=state["leaves"] + [p["smiles"] for p in state["pending"]],
                simple_leaves=state["simple_leaves"],
                total_difficulty=float(state["cumulative_difficulty"]),
            )
            return (score, -float(state["cumulative_difficulty"]), -len(state["pending"]))

        while states and expansions < max_expansions:
            if all(not state["pending"] for state in states):
                break

            next_states: List[Dict[str, Any]] = []
            for state in states:
                if not state["pending"]:
                    next_states.append(state)
                    continue

                work = _copy_state(state)
                current = work["pending"].pop(0)
                smi = str(current.get("smiles", "") or "").strip()
                depth = int(current.get("depth", 0) or 0)
                complexity = _bertz_complexity_fast(smi)

                if complexity < complexity_threshold or depth >= max_depth:
                    work["leaves"].append(smi)
                    if complexity < complexity_threshold:
                        work["simple_leaves"].append(smi)
                    else:
                        work["complex_leaves"].append(smi)
                    next_states.append(work)
                    continue

                ik = _inchi_key(smi)
                if ik and ik in work["visited_inchi"]:
                    work["leaves"].append(smi)
                    work["complex_leaves"].append(smi)
                    next_states.append(work)
                    continue
                if ik:
                    work["visited_inchi"].add(ik)

                disconnections = rank_disconnections(
                    smi,
                    top_k=max_branches,
                    max_difficulty=max_difficulty,
                )
                expansions += 1
                if not disconnections:
                    work["leaves"].append(smi)
                    work["complex_leaves"].append(smi)
                    next_states.append(work)
                    continue

                for disconnection in disconnections[:max_branches]:
                    branch = _copy_state(work)
                    step = {
                        "step": len(branch["route"]) + 1,
                        "depth": depth,
                        "target": smi,
                        "target_complexity": round(complexity, 1),
                        "reaction_name": disconnection.reaction_name,
                        "retron_name": disconnection.retron_name,
                        "description": disconnection.description,
                        "difficulty": float(disconnection.difficulty),
                        "precursor_1": disconnection.precursor_1,
                        "precursor_2": disconnection.precursor_2,
                        "complexity_delta": round(disconnection.complexity_delta, 1),
                        "overall_score": round(disconnection.overall_score, 3),
                    }
                    branch["route"].append(step)
                    branch["cumulative_difficulty"] += float(disconnection.difficulty)
                    for precursor in (disconnection.precursor_1, disconnection.precursor_2):
                        if precursor:
                            branch["pending"].append({"smiles": precursor, "depth": depth + 1})
                    next_states.append(branch)

            states = sorted(next_states, key=_interim_sort_key, reverse=True)[:beam_width]

        # Any unfinished pending fragments are complex leaves for transparent completeness scoring.
        finalized: List[Dict[str, Any]] = []
        for state in states:
            final_state = _copy_state(state)
            for pending in final_state["pending"]:
                smi = str(pending.get("smiles", "") or "").strip()
                if smi:
                    final_state["leaves"].append(smi)
                    final_state["complex_leaves"].append(smi)
            final_state["pending"] = []
            finalized.append(final_state)

        candidates: List[Dict[str, Any]] = []
        seen_signatures = set()
        for state in finalized:
            route = state["route"]
            signature = _route_candidate_signature(route) or mol_smiles
            if signature in seen_signatures:
                continue
            seen_signatures.add(signature)

            leaves = state["leaves"]
            simple_leaves = state["simple_leaves"]
            complex_leaves = state["complex_leaves"]
            cumulative_difficulty = float(state["cumulative_difficulty"])
            route_score, route_grade = _route_score_and_grade(
                steps=route,
                leaves=leaves,
                simple_leaves=simple_leaves,
                total_difficulty=cumulative_difficulty,
            )
            strategy_eval = _strategy_alignment_score(strategy_query, route, route_score)
            combined_score = round(
                0.70 * route_score + 0.30 * float(strategy_eval["strategy_alignment_score"]),
                3,
            )
            candidates.append({
                "candidate_id": f"route_{len(candidates) + 1}",
                "smiles": mol_smiles,
                "total_steps": len(route),
                "cumulative_difficulty": round(cumulative_difficulty, 3),
                "route_score": route_score,
                "route_grade": route_grade,
                "strategy_alignment_score": strategy_eval["strategy_alignment_score"],
                "combined_score": combined_score,
                "strategy_flags": strategy_eval["strategy_flags"],
                "strategy_notes": strategy_eval["strategy_notes"],
                "all_leaves_simple": len(complex_leaves) == 0,
                "leaves": leaves,
                "simple_leaves": simple_leaves,
                "complex_leaves": complex_leaves,
                "route": route,
            })

        candidates.sort(
            key=lambda item: (
                item.get("combined_score", 0.0),
                item.get("route_score", 0.0),
                -item.get("cumulative_difficulty", 999.0),
            ),
            reverse=True,
        )
        for idx, candidate in enumerate(candidates, 1):
            candidate["rank"] = idx

        best = candidates[0] if candidates else None
        result = _success({
            "smiles": mol_smiles,
            "strategy_query": strategy_query,
            "candidate_count": len(candidates[:top_k]),
            "candidates": candidates[:top_k],
            "best_candidate": best,
            "search": {
                "max_depth": max_depth,
                "max_branches": max_branches,
                "beam_width": beam_width,
                "top_k": top_k,
                "complexity_threshold": complexity_threshold,
                "max_difficulty": max_difficulty,
                "expansions": expansions,
                "finalized_states": len(finalized),
            },
        })
        if not strategy_query.strip():
            result["_warnings"] = [
                "No strategy_query was supplied; strategy_alignment_score is generic. "
                "Pass the user's route preferences for route-to-query reranking."
            ]
        return result

    except Exception as exc:
        return _error(f"Route candidate planning failed: {exc}")


plan_route_tool = ToolPlugin(
    name="plan_route",
    category="retrosynthesis",
    description=(
        "Plan a complete multi-step retrosynthetic route in a single tool call using "
        "greedy BFS (AiZynthFinder-inspired). Recursively disconnects each complex "
        "precursor (BertzCT >= complexity_threshold) until all leaf fragments are simple "
        "building blocks or max_depth is reached. Uses InChI key cycle detection to "
        "prevent loops. Returns route[] (ordered disconnection steps with reaction_name, "
        "precursors, difficulty), cumulative_difficulty, route_score (0-1 composite "
        "quality metric blending difficulty, completeness, step count, convergence), "
        "route_grade (A-F), all_leaves_simple, simple_leaves[], "
        "and a human-readable route_summary. "
        "Use instead of manually chaining identify_retrons + generate_disconnections "
        "across multiple turns when a full route is needed. "
        "Key parameters: max_depth (default 4), max_branches (default 2, controls "
        "disconnection alternatives at each node), complexity_threshold (default 100.0, "
        "BertzCT below which a fragment is considered purchasable)."
    ),
    prerequisites=[],
    fn=_plan_route,
)


plan_route_candidates_tool = ToolPlugin(
    name="plan_route_candidates",
    category="retrosynthesis",
    description=(
        "Generate multiple retrosynthetic route candidates with beam search for "
        "strategy-aware comparison. Returns candidates[] with route_score, "
        "strategy_alignment_score, combined_score, leaves, warnings, and route[] "
        "steps. Use for full-route planning when the user specifies preferences "
        "such as feasible, concise, convergent, metal-free, no protecting groups, "
        "or early ring construction; the LLM should compare candidates before "
        "choosing a final route."
    ),
    prerequisites=[],
    fn=_plan_route_candidates,
    provides=["route_candidates"],
)


def _project_plan_route_candidates(result: dict) -> Dict[str, Any]:
    if not isinstance(result, dict) or not result.get("success"):
        return {}
    return {
        "route_candidates": result.get("candidates", []),
        "best_route_candidate": result.get("best_candidate"),
    }


plan_route_candidates_tool.structured_projection = _project_plan_route_candidates


# ---------------------------------------------------------------------------
# Module export
# ---------------------------------------------------------------------------

RETROSYNTHESIS_TOOLS = [
    inspect_target_tool,
    identify_retrons_tool,
    generate_disconnections_tool,
    find_retro_precedent_tool,
    search_hte_precedent_tool,
    search_by_product_similarity_tool,
    apply_hte_templates_tool,
    plan_route_tool,
    plan_route_candidates_tool,
]

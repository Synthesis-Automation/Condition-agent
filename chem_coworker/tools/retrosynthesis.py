"""
Retrosynthesis tools for ChemCoworker.

Five ToolPlugin entries for the retrosynthesis pipeline:

  inspect_target         — enhanced molecular analysis for retrosynthesis
  identify_retrons       — SMARTS retron pattern matching
  generate_disconnections — core retrosynthetic transform engine
  find_retro_precedent   — search notes + literature for synthesis precedent
  search_hte_precedent   — DRFP k-NN search in the 231k-reaction HTE database

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

def _inspect_target(smiles: str) -> Dict[str, Any]:
    """Analyze a target molecule's complexity and key features for retrosynthesis.

    Returns molecular complexity (BertzCT), ring system count, stereocenters,
    functional group density, and a synthetic accessibility estimate.
    Always run this first — it tells the LLM HOW HARD the synthesis will be
    and which molecular features define the strategic challenge.

    Args:
        smiles: Target molecule SMILES (not a reaction SMILES).

    Returns:
        dict with complexity, ring_count, aromatic_rings, stereocenters,
        fg_density, heavy_atoms, strategic_notes.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors, GraphDescriptors

        # Strip reaction component if accidentally passed
        mol_smiles = smiles.split(">>")[0].split(">")[0].strip() if ">" in smiles else smiles
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

def _identify_retrons(smiles: str) -> Dict[str, Any]:
    """Match retron SMARTS patterns to identify applicable retrosynthetic moves.

    Scans the target molecule for structural motifs (retrons) that imply specific
    synthetic reactions. For example, a biaryl C–C bond implies Suzuki coupling;
    an aryl amine implies Buchwald-Hartwig. Returns matches ranked by difficulty
    (easiest first) so the LLM can prioritize accessible disconnections.

    Args:
        smiles: Target molecule SMILES.

    Returns:
        dict with matched_retrons (name, reaction_name, difficulty, description,
        match_count, notes, precursor_hints), total_matched.
    """
    try:
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

            retrons_data.append({
                "name": m.retron_name,
                "reaction_name": m.reaction_name,
                "difficulty": m.difficulty,
                "difficulty_display": difficulty_display,
                "description": m.description,
                "notes": m.notes,
                "precursor_hints": m.precursor_hints,
                "match_count": m.match_count,
            })

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

def _generate_disconnections(smiles: str, top_k: int = 3) -> Dict[str, Any]:
    """Generate retrosynthetic precursor pairs for the top disconnections.

    The core retrosynthesis engine. For each matched retron, applies the
    retrosynthetic transform to split the target into two precursor SMILES.
    Scores disconnections by: complexity reduction, fragment balance, and
    reaction difficulty. Returns ranked disconnections ready for conditions
    recommendation and LLM synthesis.

    Args:
        smiles: Target molecule SMILES.
        top_k: Maximum number of disconnections to return (default 3).

    Returns:
        dict with disconnections (rank, reaction_name, precursor_1, precursor_2,
        complexity_delta, fragment_balance, overall_score, description, notes).
    """
    try:
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

def _find_retro_precedent(smiles: str, reaction_name: str = "") -> Dict[str, Any]:
    """Search the knowledge base for synthesis precedent for this target.

    Queries both the reaction notes (for specific reaction-type precedent)
    and literature sources (for experimental procedures). Also searches the
    routes/ note type if available. Returns relevant excerpts so the LLM
    can integrate real experimental data into its retrosynthetic analysis.

    Args:
        smiles: Target molecule SMILES (used to derive FG search terms).
        reaction_name: Optional reaction taxonomy name (e.g., "suzuki_miyaura")
                       to focus the search. Pass empty string for broad search.

    Returns:
        dict with notes_hits, literature_hits, route_hits, search_terms_used.
    """
    try:
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
        "relevant to the target molecule. Uses functional groups and reaction type to build a "
        "smart search query. Run in parallel with identify_retrons to gather experimental context "
        "before generating disconnections."
    ),
    prerequisites=[],
    fn=_find_retro_precedent,
)


# ---------------------------------------------------------------------------
# Tool E: search_hte_precedent
# DRFP k-NN search in the HTE reaction database (~231k reactions)
# ---------------------------------------------------------------------------

# Retron/reaction name → HTE family string (supplements _family_text in core_utils)
_EXTRA_RETRON_MAP = {
    "biaryl_suzuki": "Suzuki",
    "aryl_alkyl_negishi": "Suzuki",          # closest family available
    "alkyl_alkyl_kumada": "Suzuki",
    "aryl_alkene_heck": "HeckMizoroki_coupling",
    "heck_mizoroki": "HeckMizoroki_coupling",
    "aryl_alkyne_sonogashira": "Sonogashira_coupling",
    "aryl_amine_buchwald": "C_N_Coupling",
    "aryl_amine_ullmann": "C_N_Coupling",
    "alpha_amino_reductive_amination": "Reductive_amination",
    "reductive_amination": "Reductive_amination",
    "nitrile_reduction_amine": "C_N_Coupling",
    "aryl_ether_ullmann_o": "C_O_Coupling",
    "williamson_ether": "C_O_Coupling",
    "mitsunobu_inversion": "C_O_Coupling",
    "ester_from_acid_alcohol": "Amide_formation",   # closest; esterification rare in HTE
    "amide_direct": "Amide_formation",
    "sulfonamide_formation": "Amide_formation",
    "sn2_alkyl_amine": "C_N_Coupling",
    "heterocycle_buchwald_n_arylation": "C_N_Coupling",
}


def _map_reaction_to_family(reaction_name: str) -> Optional[str]:
    """Map a retron/reaction name to a canonical HTE family string, or None."""
    if not reaction_name:
        return None
    rl = reaction_name.lower().strip()
    # Try extra retron map first (most specific)
    if rl in _EXTRA_RETRON_MAP:
        return _EXTRA_RETRON_MAP[rl]
    # Fall back to core_utils normalizer
    try:
        from chemtools.precedent.core_utils import _family_text
        mapped = _family_text(reaction_name)
        # _family_text returns the input unchanged if it doesn't recognize it
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

    try:
        mol_smiles = (
            target_smiles.split(">>")[0].split(">")[0].strip()
            if ">" in target_smiles
            else target_smiles.strip()
        )

        # Resolve family string
        family = _map_reaction_to_family(reaction_name)

        # ── Load family rows (fast path: no featurization) ────────────────
        try:
            family_key = family if family else "__ALL__"
            rows = list(_fast_load_hte_family_cached(family_key))
        except Exception as load_err:
            return _error(f"HTE loader failed: {load_err}")

        if not rows:
            return _success({
                "family": family,
                "search_mode": "none",
                "forward_smiles_used": "",
                "support_in_family": 0,
                "precedent_count": 0,
                "precedents": [],
                "message": (
                    f"No HTE rows loaded for family '{family}'. "
                    "Reaction type may not be in the database."
                ),
            })

        # ── Sort by yield descending (high-quality first) ─────────────────
        rows_sorted = sorted(
            rows,
            key=lambda r: (r.get("yield_value") or 0.0),
            reverse=True,
        )
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
        if use_drfp:
            try:
                from chemtools import reaction_similarity as rs
                if rs and rs.drfp_available():
                    q_fp = rs.encode_drfp_cached(fwd_smiles)
                    if q_fp is not None:
                        # Only DRFP-score the top-N by yield to keep latency low
                        candidates = rows_sorted[:_DRFP_CANDIDATE_CAP]
                        scored: List[tuple] = []
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

        # ── Format top-K results ──────────────────────────────────────────
        formatted = []
        scored_lookup = (
            {id(r): sim for _, sim, r in scored}
            if use_drfp and "scored" in dir() and scored
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

        search_mode = "drfp_yield_blend" if use_drfp else "family_yield"
        return _success({
            "family": family,
            "search_mode": search_mode,
            "forward_smiles_used": fwd_smiles or "(family-level only)",
            "support_in_family": support,
            "drfp_candidates_scored": min(support, _DRFP_CANDIDATE_CAP) if use_drfp else 0,
            "precedent_count": len(formatted),
            "precedents": formatted,
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
# Module export
# ---------------------------------------------------------------------------

RETROSYNTHESIS_TOOLS = [
    inspect_target_tool,
    identify_retrons_tool,
    generate_disconnections_tool,
    find_retro_precedent_tool,
    search_hte_precedent_tool,
]

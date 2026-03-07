"""
Forward synthesis tools for ChemCoworker.

Eight ToolPlugin entries for the forward synthesis pipeline:

  inspect_reactants           — analyze reactant pair: FGs, electro/nucleophilicity,
                                 compatibility flags, competing-reaction warnings
  identify_reactions          — find all HTE forward templates compatible with the
                                 reactant functional group pair (mirrors identify_retrons)
  generate_products           — apply forward SMARTS templates via
                                 AllChem.RunReactants(); returns ranked product SMILES
                                 with regioisomers (mirrors generate_disconnections)
  rank_products               — score product candidates by HTE yield proxy,
                                 chemoselectivity, difficulty; applies confidence labels
  find_forward_precedent      — RAG search in knowledge base for reactant-pair precedent
                                 (mirrors find_retro_precedent)
  search_reactant_precedent   — DRFP k-NN search in 231k HTE reactions using
                                 reactant SMILES as the query (mirrors search_hte_precedent)
  recommend_forward_conditions — get optimal catalyst/ligand/base/solvent for the
                                  predicted reaction type (bridges to conditions.py)
  plan_forward_route          — multi-step forward synthesis planner: applies reactions
                                 sequentially on a scaffold; returns ordered step list

All follow the _success() / _error() contract. Registered as ToolPlugins
so they appear in REASON_PROMPT tool descriptions.

Usage (by ChemCoworker, not directly):
  When task_type == "forward_synthesis", the agent selects these tools automatically.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from ._helpers import _error, _success, _to_jsonable, _validate_reaction_smiles
from ._base import ToolPlugin


# ---------------------------------------------------------------------------
# Tool 1: inspect_reactants
# Analyze reactant pair for FGs, compatibility, and strategic notes
# ---------------------------------------------------------------------------

def _inspect_reactants(
    smiles_a: str = "",
    smiles_b: str = "",
    reactant_a: str = "",
    reactant_b: str = "",
) -> Dict[str, Any]:
    """Analyze a reactant pair for forward synthesis planning.

    Returns functional group inventory, electrophile/nucleophile classification,
    compatibility assessment, and strategic notes for both reactants.
    Run this FIRST for any forward synthesis query — it tells the LLM what
    reactions are chemically plausible and flags potential competing reactions.

    Args:
        smiles_a: First reactant SMILES.
        smiles_b: Second reactant SMILES (leave empty for unimolecular reactions).
        reactant_a: Alias for smiles_a.
        reactant_b: Alias for smiles_b.

    Returns:
        dict with reactant_a_profile, reactant_b_profile (each with fgs,
        electrophilic_sites, nucleophilic_sites, leaving_group, mw, heavy_atoms),
        compatibility_flags, strategic_notes.
    """
    try:
        smiles_a = smiles_a or reactant_a
        smiles_b = smiles_b or reactant_b
        if not smiles_a:
            return _error("smiles_a is required")

        from chemtools.forward.reactor import ReactantAnalyzer
        analyzer = ReactantAnalyzer()

        profile_a = analyzer.profile(smiles_a)
        profile_b = analyzer.profile(smiles_b) if smiles_b else None

        def _fmt_profile(p: Any) -> Dict[str, Any]:
            return {
                "smiles": p.smiles,
                "canonical_smiles": p.canonical_smiles,
                "molecular_weight": p.molecular_weight,
                "heavy_atoms": p.heavy_atoms,
                "functional_groups": p.functional_groups,
                "functional_groups_by_role": p.fg_categories,
                "electrophilic_sites": p.electrophilic_sites,
                "nucleophilic_sites": p.nucleophilic_sites,
                "leaving_group": p.leaving_group,
            }

        # Compatibility analysis
        compatibility_flags: List[str] = []
        strategic_notes: List[str] = []

        if profile_b:
            # Check electrophile + nucleophile pairing
            has_elec_a = bool(profile_a.electrophilic_sites)
            has_nuc_b = bool(profile_b.nucleophilic_sites)
            has_nuc_a = bool(profile_a.nucleophilic_sites)
            has_elec_b = bool(profile_b.electrophilic_sites)

            if (has_elec_a and has_nuc_b) or (has_nuc_a and has_elec_b):
                compatibility_flags.append("COMPATIBLE: electrophile–nucleophile pair detected")
            else:
                compatibility_flags.append(
                    "CAUTION: both reactants appear to be the same role "
                    "(both electrophilic or both nucleophilic); "
                    "cross-coupling may require organometallic activation"
                )

            # Competing FG warnings
            all_fgs_a = set(profile_a.functional_groups)
            all_fgs_b = set(profile_b.functional_groups)
            competing = []
            if "primary_amine" in all_fgs_a and "secondary_amine" in all_fgs_a:
                competing.append("Both primary and secondary amine in reactant A — chemoselectivity challenge")
            if "aldehyde" in all_fgs_a and "ketone" in all_fgs_a:
                competing.append("Both aldehyde and ketone in reactant A — chemoselective reaction may be difficult")
            if len(profile_a.electrophilic_sites) > 2:
                competing.append(
                    f"Reactant A has {len(profile_a.electrophilic_sites)} electrophilic sites — "
                    "multiple reaction sites possible"
                )
            if competing:
                compatibility_flags.extend(competing)

            # Steric notes
            if profile_a.heavy_atoms > 40 or (profile_b and profile_b.heavy_atoms > 40):
                strategic_notes.append(
                    "Large reactant(s) detected — steric effects may reduce yield; "
                    "consider bulky ligands for cross-coupling"
                )
            # Amine safety
            if ("primary_amine" in all_fgs_b or "secondary_amine" in all_fgs_b) and (
                "aryl_bromide" in all_fgs_a or "aryl_chloride" in all_fgs_a
            ):
                strategic_notes.append(
                    "Aryl halide + amine: Buchwald-Hartwig (Pd) or Chan-Lam (Cu) are likely; "
                    "also consider SNAr if arene is electron-poor"
                )
        else:
            strategic_notes.append(
                "Unimolecular reaction mode: consider reductions, oxidations, "
                "rearrangements, or intramolecular cyclizations"
            )

        result: Dict[str, Any] = {
            "reactant_a": _fmt_profile(profile_a),
            "compatibility_flags": compatibility_flags,
            "strategic_notes": strategic_notes,
        }
        if profile_b:
            result["reactant_b"] = _fmt_profile(profile_b)

        return _success(result)

    except Exception as exc:
        return _error(f"Reactant inspection failed: {exc}")


inspect_reactants_tool = ToolPlugin(
    name="inspect_reactants",
    category="forward_synthesis",
    description=(
        "Analyze a reactant pair for forward synthesis planning. Returns functional group "
        "inventory, electrophile/nucleophile classification, and compatibility flags for both "
        "reactants. Run FIRST for any forward synthesis query to understand what reactions are "
        "chemically plausible. Args: smiles_a (required), smiles_b (optional for unimolecular)."
    ),
    prerequisites=[],
    fn=_inspect_reactants,
)


# ---------------------------------------------------------------------------
# Tool 2: identify_reactions
# Find all forward templates applicable to the reactant FG pair
# ---------------------------------------------------------------------------

def _identify_reactions(
    smiles_a: str = "",
    smiles_b: str = "",
    reactant_a: str = "",
    reactant_b: str = "",
) -> Dict[str, Any]:
    """Identify forward reactions compatible with the reactant functional group pair.

    Scans the HTE forward template library for templates whose forward SMARTS
    match the reactant functional groups. For example, aryl bromide + boronic acid
    implies Suzuki-Miyaura; aryl halide + amine implies Buchwald-Hartwig or SNAr.
    Returns matches sorted by difficulty (easiest first).

    Args:
        smiles_a: First reactant SMILES.
        smiles_b: Second reactant SMILES (empty for unimolecular reactions).
        reactant_a: Alias for smiles_a.
        reactant_b: Alias for smiles_b.

    Returns:
        dict with compatible_reactions (name, category, difficulty, description,
        hte_families, notes), total_matched.
    """
    try:
        smiles_a = smiles_a or reactant_a
        smiles_b = smiles_b or reactant_b
        if not smiles_a:
            return _error("smiles_a is required")

        from chemtools.forward.reactor import ReactantAnalyzer
        analyzer = ReactantAnalyzer()
        matches = analyzer.find_compatible_templates(smiles_a, smiles_b)

        if not matches:
            return _success({
                "smiles_a": smiles_a,
                "smiles_b": smiles_b,
                "total_matched": 0,
                "compatible_reactions": [],
                "message": (
                    "No standard forward reactions detected for this reactant pair. "
                    "The reactants may lack standard functional groups, or require "
                    "non-classical activation. Consider: C–H functionalization, "
                    "radical reactions, or photoredox chemistry."
                ),
            })

        reactions_data = []
        for m in matches:
            filled = round(m.difficulty * 5)
            difficulty_display = "●" * filled + "○" * (5 - filled)
            reactions_data.append({
                "name": m.template_name,
                "taxonomy_id": m.taxonomy_id,
                "taxonomy_family_id": m.taxonomy_id,
                "category": m.category,
                "difficulty": m.difficulty,
                "difficulty_display": difficulty_display,
                "description": m.description,
                "hte_families": m.hte_families,
                "n_reactants": m.n_reactants,
                "forward_smarts": m.forward_smarts,
                "notes": m.notes,
            })

        return _success({
            "smiles_a": smiles_a,
            "smiles_b": smiles_b,
            "total_matched": len(reactions_data),
            "compatible_reactions": reactions_data,
        })

    except Exception as exc:
        return _error(f"Reaction identification failed: {exc}")


identify_reactions_tool = ToolPlugin(
    name="identify_reactions",
    category="forward_synthesis",
    description=(
        "Identify forward reactions compatible with the reactant functional group pair. "
        "Scans 47 HTE forward templates and returns those whose SMARTS match the supplied "
        "reactants. Results are sorted by difficulty (easiest first). "
        "Run after inspect_reactants; provides the basis for generate_products. "
        "Args: smiles_a (required), smiles_b (optional for unimolecular reactions)."
    ),
    prerequisites=[],
    fn=_identify_reactions,
)


# ---------------------------------------------------------------------------
# Tool 3: generate_products
# Apply forward SMARTS templates and return ranked product SMILES
# ---------------------------------------------------------------------------

def _generate_products(
    smiles_a: str = "",
    smiles_b: str = "",
    reactant_a: str = "",
    reactant_b: str = "",
    reaction_name: str = "",
    top_k: int = 5,
) -> Dict[str, Any]:
    """Apply forward reaction templates to generate product SMILES candidates.

    Core forward synthesis engine. For each compatible template, applies the
    forward SMARTS using AllChem.RunReactants() to produce validated product
    SMILES. Scores products by HTE yield proxy, chemoselectivity, and difficulty.
    Returns ranked predictions including regioisomers.

    Args:
        smiles_a: First reactant SMILES.
        smiles_b: Second reactant SMILES (empty for unimolecular reactions).
        reactant_a: Alias for smiles_a.
        reactant_b: Alias for smiles_b.
        reaction_name: Optional filter — only apply templates matching this
                       name or hte_family (e.g. "suzuki_miyaura", "amide_coupling").
        top_k: Maximum number of product predictions to return (default 5).

    Returns:
        dict with products (rank, product_smiles, reaction_smiles, template_name,
        description, difficulty, hte_yield_proxy, overall_score, confidence_label,
        new_stereocenters, all_product_smiles, competing_templates, notes).
    """
    try:
        smiles_a = smiles_a or reactant_a
        smiles_b = smiles_b or reactant_b
        if not smiles_a:
            return _error("smiles_a is required")

        from chemtools.forward.reactor import ForwardReactor
        from chemtools.forward.scoring import score_products

        reactor = ForwardReactor()
        predictions = reactor.generate(
            smiles_a=smiles_a,
            smiles_b=smiles_b,
            reaction_name=reaction_name,
            top_k=top_k * 3,   # over-generate before scoring trim
        )

        if not predictions:
            return _success({
                "smiles_a": smiles_a,
                "smiles_b": smiles_b,
                "total_products": 0,
                "products": [],
                "message": (
                    "No products could be generated automatically. "
                    "Reactants may lack standard HTE-compatible functional groups, "
                    "or the reaction may require non-standard conditions (C–H activation, "
                    "photoredox, organocatalysis). Use find_forward_precedent for literature search."
                ),
            })

        # Score and trim
        ranked = score_products(predictions, smiles_a=smiles_a, smiles_b=smiles_b)
        ranked = ranked[:top_k]

        products_out = []
        invalid_rxn_filtered = 0
        for rank, pred in enumerate(ranked, 1):
            rxn_smiles_raw = str(pred.reaction_smiles or "")
            rxn_smiles, rxn_err = _validate_reaction_smiles(rxn_smiles_raw, require_product=True)
            if rxn_err:
                invalid_rxn_filtered += 1
                continue
            filled = round(pred.difficulty * 5)
            difficulty_display = "●" * filled + "○" * (5 - filled)

            products_out.append({
                "rank": len(products_out) + 1,
                "product_smiles": pred.product_smiles,
                "reaction_smiles": rxn_smiles,
                "template_name": pred.template_name,
                "taxonomy_id": pred.taxonomy_id,
                "taxonomy_family_id": pred.taxonomy_id,
                "description": pred.description,
                "difficulty": pred.difficulty,
                "difficulty_display": difficulty_display,
                "hte_yield_proxy": pred.hte_yield_proxy,
                "overall_score": pred.overall_score,
                "confidence_label": pred.confidence_label,
                "new_stereocenters": pred.new_stereocenters,
                "all_product_smiles": pred.all_product_smiles,
                "competing_templates": pred.competing_templates[:3],
                "notes": pred.notes,
                "hte_families": pred.hte_families,
            })

        payload: Dict[str, Any] = {
            "smiles_a": smiles_a,
            "smiles_b": smiles_b,
            "total_products": len(products_out),
            "products": products_out,
        }
        if invalid_rxn_filtered:
            payload["invalid_reaction_smiles_filtered"] = int(invalid_rxn_filtered)
            if not products_out:
                payload["message"] = (
                    "Generated product candidates were filtered because their reaction SMILES "
                    "failed validity checks."
                )
        return _success(payload)

    except Exception as exc:
        return _error(f"Product generation failed: {exc}")


generate_products_tool = ToolPlugin(
    name="generate_products",
    category="forward_synthesis",
    description=(
        "Core forward synthesis tool: apply HTE forward SMARTS templates to generate ranked "
        "product SMILES from the reactant pair. Returns products with hte_yield_proxy, "
        "overall_score, confidence_label, new_stereocenters, and all regioisomers. "
        "Run after identify_reactions. Output feeds into recommend_forward_conditions. "
        "Args: smiles_a (required), smiles_b (optional), reaction_name (optional filter), "
        "top_k (default 5)."
    ),
    prerequisites=[],
    fn=_generate_products,
)


# ---------------------------------------------------------------------------
# Tool 4: rank_products
# Explicit re-ranking step — useful when generate_products returned many candidates
# ---------------------------------------------------------------------------

def _rank_products(
    products: Optional[List[Dict[str, Any]]] = None,
    smiles_a: str = "",
    smiles_b: str = "",
    reactant_a: str = "",
    reactant_b: str = "",
    top_k: int = 3,
) -> Dict[str, Any]:
    """Re-rank and filter product candidates from generate_products.

    Takes the products list from generate_products and applies additional
    chemistry filters: flags products with new stereocenters, highlights
    the most selective reaction pathway, and trims to top_k.

    Useful when generate_products returns many competing candidates and the LLM
    needs explicit guidance on which product to prioritize.

    Args:
        products: List of product dicts from generate_products output.
        smiles_a: First reactant (used to recompute scores if products is empty).
        smiles_b: Second reactant.
        reactant_a: Alias for smiles_a.
        reactant_b: Alias for smiles_b.
        top_k: Maximum products to return (default 3).

    Returns:
        dict with ranked_products, selectivity_recommendation,
        stereochemistry_flags, competing_reaction_summary.
    """
    try:
        smiles_a = smiles_a or reactant_a
        smiles_b = smiles_b or reactant_b

        # If products is None (not passed), regenerate from smiles_a/smiles_b.
        # An explicitly passed empty list [] means genuinely no products — skip
        # regeneration and return success immediately.
        if products is None:
            regen = _generate_products(smiles_a=smiles_a, smiles_b=smiles_b, top_k=top_k * 2)
            if not regen.get("success"):
                return regen
            products = regen.get("products", [])

        if not products:
            return _success({
                "ranked_products": [],
                "selectivity_recommendation": "No products to rank.",
                "message": "No candidates available for re-ranking.",
            })

        # Sort by overall_score descending (may already be sorted)
        sorted_products = sorted(products, key=lambda p: p.get("overall_score", 0), reverse=True)
        top = sorted_products[:top_k]

        # Selectivity analysis
        if len(top) == 1:
            selectivity_note = (
                f"Single dominant product predicted: {top[0]['product_smiles']} "
                f"via {top[0]['template_name']} "
                f"(confidence: {top[0].get('confidence_label','?')})."
            )
        else:
            best = top[0]
            second = top[1]
            score_gap = round(best.get("overall_score", 0) - second.get("overall_score", 0), 1)
            if score_gap >= 15:
                selectivity_note = (
                    f"Strong preference for {best['template_name']} "
                    f"(score gap {score_gap:.1f} over {second['template_name']}). "
                    f"Chemoselectivity should be manageable."
                )
            else:
                selectivity_note = (
                    f"Multiple competing reactions within {score_gap:.1f} score points "
                    f"({best['template_name']} vs {second['template_name']}). "
                    "Chemoselectivity will be challenging; consider protecting groups or "
                    "reaction order optimization."
                )

        # Stereochemistry flags
        stereo_flags = [
            f"{p['product_smiles']} introduces {p['new_stereocenters']} new stereocenter(s)"
            for p in top
            if p.get("new_stereocenters", 0) > 0
        ]

        # Competing reaction summary
        all_competing = set()
        for p in top:
            all_competing.update(p.get("competing_templates", []))
        competing_summary = sorted(all_competing)[:5]

        return _success({
            "ranked_products": top,
            "selectivity_recommendation": selectivity_note,
            "stereochemistry_flags": stereo_flags,
            "competing_reaction_summary": competing_summary,
        })

    except Exception as exc:
        return _error(f"Product ranking failed: {exc}")


rank_products_tool = ToolPlugin(
    name="rank_products",
    category="forward_synthesis",
    description=(
        "Re-rank product candidates from generate_products and provide selectivity guidance. "
        "Highlights the dominant product pathway, flags new stereocenters, and identifies "
        "competing reactions. Run when generate_products returns multiple competing candidates "
        "and you need explicit chemoselectivity analysis. "
        "Args: products (list from generate_products), smiles_a, smiles_b, top_k (default 3)."
    ),
    prerequisites=["generate_products"],
    fn=_rank_products,
)


# ---------------------------------------------------------------------------
# Tool 5: find_forward_precedent
# RAG search in knowledge base for reactant-pair precedent
# ---------------------------------------------------------------------------

def _find_forward_precedent(
    smiles_a: str = "",
    smiles_b: str = "",
    reactant_a: str = "",
    reactant_b: str = "",
    reaction_name: str = "",
    reaction_type: str = "",
) -> Dict[str, Any]:
    """Search the knowledge base for experimental precedent for this reactant pair.

    Queries reaction notes and literature for precedents involving the
    functional groups and reactions matching the reactant pair.
    Returns relevant excerpts so the LLM can integrate real experimental data.

    Args:
        smiles_a: First reactant SMILES.
        smiles_b: Second reactant SMILES.
        reactant_a: Alias for smiles_a.
        reactant_b: Alias for smiles_b.
        reaction_name: Optional reaction taxonomy name to focus search
                       (e.g. "suzuki_miyaura", "amide_coupling").
        reaction_type: Alias for reaction_name.

    Returns:
        dict with notes_hits, literature_hits, reaction_note, search_query.
    """
    try:
        smiles_a = smiles_a or reactant_a
        smiles_b = smiles_b or reactant_b
        reaction_name = reaction_name or reaction_type
        if not smiles_a:
            return _error("smiles_a is required")

        # Build search query from FGs + reaction name
        query_parts = []
        if reaction_name:
            query_parts.append(reaction_name.replace("_", " "))

        # FG to search term mapping
        _FG_SEARCH_MAP = {
            "aryl_bromide": "aryl bromide coupling",
            "aryl_chloride": "aryl chloride coupling",
            "aryl_iodide": "aryl iodide coupling",
            "arylboronic_acid": "boronic acid Suzuki coupling",
            "primary_amine": "amine synthesis N-arylation",
            "secondary_amine": "amine alkylation reductive amination",
            "aldehyde": "aldehyde reductive amination condensation",
            "carboxylic_acid": "amide coupling esterification",
            "alcohol": "esterification etherification",
            "thiol": "thioether C-S coupling",
            "terminal_alkyne": "CuAAC click triazole",
            "alkyl_azide": "CuAAC click triazole",
            "michael_acceptor": "Michael addition conjugate addition",
        }

        # Collect FGs from both reactants
        for smi in filter(None, [smiles_a, smiles_b]):
            try:
                from chemtools.forward.reactor import ReactantAnalyzer
                profile = ReactantAnalyzer().profile(smi)
                for fg in profile.functional_groups[:3]:
                    if fg in _FG_SEARCH_MAP:
                        query_parts.append(_FG_SEARCH_MAP[fg])
                    else:
                        query_parts.append(fg.replace("_", " "))
            except Exception:
                pass

        if not query_parts:
            query_parts = ["synthesis", "forward", "reaction"]

        search_query = " ".join(query_parts[:6])

        # Search notes
        notes_hits = []
        try:
            from chem_coworker.tools.notes import _search_notes
            notes_result = _search_notes(
                query=search_query,
                note_types=["reactions", "protocols"],
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
            "smiles_a": smiles_a,
            "smiles_b": smiles_b,
            "search_query": search_query,
            "search_terms_used": query_parts[:6],
            "total_found": total_found,
            "reaction_note": reaction_note,
            "notes_hits": notes_hits,
            "literature_hits": lit_hits,
            "message": (
                f"Found {total_found} relevant sources. "
                + ("No precedent found — reactant combination may be novel." if total_found == 0 else "")
            ),
        })

    except Exception as exc:
        return _error(f"Forward precedent search failed: {exc}")


find_forward_precedent_tool = ToolPlugin(
    name="find_forward_precedent",
    category="forward_synthesis",
    description=(
        "Search the reaction notes and literature knowledge base for experimental precedent "
        "relevant to the reactant pair. Uses functional groups and reaction name to build a "
        "smart search query. Args: smiles_a, smiles_b, reaction_name (e.g. 'suzuki_miyaura'). "
        "Run in parallel with identify_reactions to gather context before generating products."
    ),
    prerequisites=[],
    fn=_find_forward_precedent,
)


# ---------------------------------------------------------------------------
# Tool 6: search_reactant_precedent
# DRFP k-NN search in HTE database using reactant SMILES as query
# ---------------------------------------------------------------------------

def _search_reactant_precedent(
    smiles_a: str = "",
    smiles_b: str = "",
    reactant_a: str = "",
    reactant_b: str = "",
    reaction_name: str = "",
    product_smiles: str = "",
    top_k: int = 5,
) -> Dict[str, Any]:
    """Search the HTE database for reactions with similar reactants via DRFP k-NN.

    Composes a forward reaction SMILES (reactant_a.reactant_b >> product) and
    uses DRFP fingerprint similarity to search ~231k HTE reactions. When no
    product SMILES is available yet, passes only the reaction family for
    yield-ranked family search.

    Args:
        smiles_a: First reactant SMILES.
        smiles_b: Second reactant SMILES.
        reactant_a: Alias for smiles_a.
        reactant_b: Alias for smiles_b.
        reaction_name: Reaction family / template name (e.g. "suzuki_miyaura").
        product_smiles: Optional product SMILES for DRFP-ranked search; if empty,
                        falls back to family-level yield ranking.
        top_k: Number of precedents to return (default 5).

    Returns:
        dict with family, search_mode, support_in_family, precedents[].
        Each precedent has: reaction_smiles, similarity, yield, condition_core,
        catalyst, base, solvent, reagents, reference.
    """
    try:
        smiles_a = smiles_a or reactant_a
        smiles_b = smiles_b or reactant_b
        if not smiles_a:
            return _error("smiles_a is required")

        # Delegate to the existing _search_hte_precedent tool which already
        # handles DRFP k-NN — we just compose the query differently:
        # in retro: target>>precursors; in forward: reactants>>product
        from chem_coworker.tools.retrosynthesis import _search_hte_precedent

        # Use product_smiles as the "target" for the underlying call
        target = product_smiles or smiles_a

        result = _search_hte_precedent(
            target_smiles=target,
            reaction_name=reaction_name,
            precursor_1=smiles_a,
            precursor_2=smiles_b,
            top_k=top_k,
        )

        if result.get("success"):
            result["note"] = (
                "Searched HTE database with reactant pair as query. "
                "Precedents shown are reactions with similar reactants/products."
            )
        return result

    except Exception as exc:
        return _error(f"Reactant precedent search failed: {exc}")


search_reactant_precedent_tool = ToolPlugin(
    name="search_reactant_precedent",
    category="forward_synthesis",
    description=(
        "Search the ~231k-reaction HTE database for precedents with similar reactants "
        "using DRFP fingerprint k-NN. Composes a forward reaction SMILES and finds the "
        "nearest experimental neighbors with conditions (catalyst, base, solvent, yield). "
        "Args: smiles_a (required), smiles_b (optional), reaction_name (family filter), "
        "product_smiles (optional; enables DRFP-ranked mode), top_k (default 5). "
        "Run after generate_products to ground conditions in real experimental data."
    ),
    prerequisites=[],
    fn=_search_reactant_precedent,
)


# ---------------------------------------------------------------------------
# Tool 7: recommend_forward_conditions
# Get optimal catalyst/ligand/base/solvent for the predicted reaction type
# ---------------------------------------------------------------------------

def _recommend_forward_conditions(
    smiles_a: str = "",
    smiles_b: str = "",
    reactant_a: str = "",
    reactant_b: str = "",
    product_smiles: str = "",
    reaction_type: str = "",
    top_k: int = 5,
) -> Dict[str, Any]:
    """Recommend optimal reaction conditions for the predicted forward reaction.

    Uses the HTE-backed condition recommender to suggest catalyst, ligand,
    base, and solvent for the forward reaction given the reactant SMILES and
    reaction type. Ranks conditions by Z-score (statistical performance across
    similar reactions in the 66k-experiment HTE database).

    Args:
        smiles_a: First reactant SMILES.
        smiles_b: Second reactant SMILES.
        reactant_a: Alias for smiles_a.
        reactant_b: Alias for smiles_b.
        product_smiles: Required product SMILES for valid forward reaction construction.
        reaction_type: Reaction taxonomy name (e.g. "suzuki_miyaura",
                       "C_N_Coupling", "Amide_formation"). If empty, the
                       recommender will auto-detect from the reactants.
        top_k: Number of condition sets to return (default 5).

    Returns:
        dict with recommended conditions (catalyst, ligand, base, solvent,
        avg_yield, z_score, num_experiments, confidence_score).
    """
    try:
        smiles_a = smiles_a or reactant_a
        smiles_b = smiles_b or reactant_b
        if not smiles_a:
            return _error("smiles_a is required")

        # Build a validated forward reaction SMILES and delegate to recommend_conditions
        from chem_coworker.tools.conditions import _recommend_conditions

        # Compose forward reaction SMILES: reactants >> product (or empty placeholder)
        can_a = smiles_a
        can_b = smiles_b
        if can_b:
            reactants_str = f"{can_a}.{can_b}"
        else:
            reactants_str = can_a

        if not str(product_smiles or "").strip():
            return _error(
                "product_smiles is required for recommend_forward_conditions; "
                "cannot build a valid reaction SMILES without a product."
            )
        rxn_smiles_raw = f"{reactants_str}>>{product_smiles}"
        rxn_smiles, rxn_err = _validate_reaction_smiles(rxn_smiles_raw, require_product=True)
        if rxn_err:
            return _error(f"Invalid reaction SMILES for condition recommendation: {rxn_err}")

        result = _recommend_conditions(
            reaction_smiles=rxn_smiles,
            top_k=top_k,
        )
        return result

    except Exception as exc:
        return _error(f"Forward condition recommendation failed: {exc}")


recommend_forward_conditions_tool = ToolPlugin(
    name="recommend_forward_conditions",
    category="forward_synthesis",
    description=(
        "Recommend optimal catalyst, ligand, base, and solvent for the predicted forward "
        "reaction. Uses the HTE-backed condition ranker (66k experiments, Z-score ranked). "
        "Args: smiles_a (required), smiles_b (optional), product_smiles (required), "
        "reaction_type (e.g. 'suzuki_miyaura', 'C_N_Coupling'; leave empty for "
        "auto-detection), top_k (default 5). "
        "Run in the final group after generate_products and search_reactant_precedent."
    ),
    prerequisites=["generate_products"],
    fn=_recommend_forward_conditions,
)


# ---------------------------------------------------------------------------
# Tool 8: plan_forward_route
# Multi-step forward synthesis planner for sequential reactions on a scaffold
# ---------------------------------------------------------------------------

def _plan_forward_route(
    scaffold_smiles: str = "",
    reaction_sequence: Optional[List[str]] = None,
    reagent_smiles_list: Optional[List[str]] = None,
    max_steps: int = 4,
) -> Dict[str, Any]:
    """Plan a multi-step forward synthesis route from a scaffold + reagent list.

    Applies reactions sequentially: each step transforms the current intermediate
    using a new reagent. Returns an ordered list of steps with products and
    conditions for each step.

    Unlike plan_route (which greedily disconnects a target), plan_forward_route
    builds from the scaffold outward — it is the forward analog of walking the
    retrosynthetic tree in the synthesis direction.

    Args:
        scaffold_smiles: Starting material / scaffold SMILES.
        reaction_sequence: Optional list of reaction names to apply in order
                           (e.g. ["suzuki_miyaura", "amide_coupling"]).
                           If empty, auto-selects based on reactive FGs.
        reagent_smiles_list: List of reagent/coupling-partner SMILES, one per step.
                             If provided, must match length of reaction_sequence.
        max_steps: Maximum number of forward steps to plan (default 4).

    Returns:
        dict with route_steps (step_number, reactant, reagent, product_smiles,
        reaction_name, conditions_summary, difficulty), total_steps,
        final_product, bm_complexity_change.
    """
    try:
        if not scaffold_smiles:
            return _error("scaffold_smiles is required")

        from chemtools.forward.reactor import ForwardReactor, ReactantAnalyzer
        from chemtools.forward.scoring import score_products
        from rdkit.Chem.GraphDescriptors import BertzCT
        from rdkit import Chem, rdBase

        reactor = ForwardReactor()
        analyzer = ReactantAnalyzer()

        steps = []
        current_smiles = scaffold_smiles
        reagents = list(reagent_smiles_list or [])
        rxn_names = list(reaction_sequence or [])

        with rdBase.BlockLogs():
            start_mol = Chem.MolFromSmiles(scaffold_smiles)
        start_complexity = BertzCT(start_mol) if start_mol else 0.0

        for step_idx in range(min(max_steps, max(len(reagents), len(rxn_names), 1))):
            reagent = reagents[step_idx] if step_idx < len(reagents) else ""
            rxn_name = rxn_names[step_idx] if step_idx < len(rxn_names) else ""

            # Generate products for this step
            preds = reactor.generate(
                smiles_a=current_smiles,
                smiles_b=reagent,
                reaction_name=rxn_name,
                top_k=3,
            )

            if not preds:
                steps.append({
                    "step_number": step_idx + 1,
                    "reactant": current_smiles,
                    "reagent": reagent,
                    "product_smiles": None,
                    "reaction_name": rxn_name or "unknown",
                    "error": "No products generated — step failed or FGs exhausted",
                })
                break

            ranked = score_products(preds)
            best = ranked[0]
            rxn_smiles, rxn_err = _validate_reaction_smiles(str(best.reaction_smiles or ""), require_product=True)
            if rxn_err:
                steps.append({
                    "step_number": step_idx + 1,
                    "reactant": current_smiles,
                    "reagent": reagent,
                    "product_smiles": best.product_smiles,
                    "reaction_name": best.template_name or (rxn_name or "unknown"),
                    "error": f"Generated invalid reaction SMILES for step: {rxn_err}",
                })
                break

            # Attempt to get top conditions hint from HTE
            conditions_hint: Optional[str] = None
            try:
                from chem_coworker.tools.retrosynthesis import _fast_load_hte_family_cached, _extract_catalyst_name
                for fam in best.hte_families:
                    rows = list(_fast_load_hte_family_cached(fam))
                    if rows:
                        top_row = max(rows, key=lambda r: r.get("yield_value") or 0.0)
                        cat = _extract_catalyst_name(top_row.get("catalyst"))
                        base = top_row.get("base_uid")
                        solv = top_row.get("solvent_uid")
                        parts = [x for x in [cat, base, solv] if x]
                        if parts:
                            conditions_hint = " / ".join(parts)
                        break
            except Exception:
                pass

            steps.append({
                "step_number": step_idx + 1,
                "reactant": current_smiles,
                "reagent": reagent,
                "product_smiles": best.product_smiles,
                "reaction_name": best.template_name,
                "taxonomy_id": best.taxonomy_id,
                "taxonomy_family_id": best.taxonomy_id,
                "reaction_smiles": rxn_smiles,
                "description": best.description,
                "difficulty": best.difficulty,
                "overall_score": best.overall_score,
                "confidence_label": best.confidence_label,
                "conditions_summary": conditions_hint,
                "notes": best.notes,
            })

            current_smiles = best.product_smiles

            # Stop if we've run out of reagents and there's no more to do
            if step_idx + 1 >= len(reagents) and not rxn_names:
                break

        # Compute complexity change
        final_mol = None
        with rdBase.BlockLogs():
            final_mol = Chem.MolFromSmiles(current_smiles)
        end_complexity = BertzCT(final_mol) if final_mol else 0.0

        return _success({
            "scaffold_smiles": scaffold_smiles,
            "final_product": current_smiles,
            "total_steps": len(steps),
            "complexity_start": round(start_complexity, 1),
            "complexity_end": round(end_complexity, 1),
            "complexity_increase": round(end_complexity - start_complexity, 1),
            "route_steps": steps,
            "message": (
                f"Forward route from scaffold → final product in {len(steps)} step(s). "
                f"Complexity increased from {start_complexity:.0f} to {end_complexity:.0f} (BertzCT)."
            ),
        })

    except Exception as exc:
        return _error(f"Forward route planning failed: {exc}")


plan_forward_route_tool = ToolPlugin(
    name="plan_forward_route",
    category="forward_synthesis",
    description=(
        "Plan a multi-step forward synthesis route: apply reactions sequentially on a scaffold "
        "to build up the target. Each step uses generate_products internally. Returns ordered "
        "steps with product SMILES, reaction names, and condition summaries. "
        "Args: scaffold_smiles (required), reaction_sequence (list of reaction names, optional), "
        "reagent_smiles_list (list of coupling-partner SMILES, one per step; optional), "
        "max_steps (default 4). "
        "Use for building up complex targets from simple scaffolds."
    ),
    prerequisites=[],
    fn=_plan_forward_route,
)


# ---------------------------------------------------------------------------
# FORWARD_SYNTHESIS_TOOLS — list exported for REGISTRY registration
# ---------------------------------------------------------------------------

FORWARD_SYNTHESIS_TOOLS = [
    inspect_reactants_tool,
    identify_reactions_tool,
    generate_products_tool,
    rank_products_tool,
    find_forward_precedent_tool,
    search_reactant_precedent_tool,
    recommend_forward_conditions_tool,
    plan_forward_route_tool,
]

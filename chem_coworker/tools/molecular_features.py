"""
Molecular feature analysis tools for ChemCoworker.

Exposes pre-existing chemtools capabilities as ToolPlugins:

  featurize_molecule   : Full electronic + steric + ring profile for all
                         reactive sites.  Covers aryl electronics (EWG/EDG
                         scoring), aryl ortho-steric bulk, alkyl alpha-carbon
                         classification, and ring metrics — in one call.

  assess_snar_feasibility : SNAr-specific feasibility check derived from
                            taxonomy-defined thresholds.  Use when an aryl
                            halide is the electrophile and you need to know
                            whether the ring is electronically activated
                            enough for nucleophilic aromatic substitution.
"""
from __future__ import annotations

from typing import Any, Dict, List

from ._helpers import _error, _success, _to_jsonable
from ._base import ToolPlugin


# ---------------------------------------------------------------------------
# Tool 1 — featurize_molecule
# ---------------------------------------------------------------------------

def _featurize_molecule(smiles: str) -> Dict[str, Any]:
    """Profile all reactive sites in a molecule: electronics, sterics, and ring metrics.

    Wraps chemtools.featurizers.molecule.featurize_molecule to expose:
      - ranked_motifs       : ordered list of detected site types
                              (e.g. ['Ar-Br', 'Ar-OMe', 'Alkyl-NH2'])
      - electronics         : per aryl-site EWG/EDG score (0=electron-rich,
                              10=electron-poor) + description
      - sterics             : per site bulk score (0=unhindered, 10=highly
                              hindered); aryl → ortho bulk; alkyl → alpha-C
                              classification (methyl/primary/secondary/tertiary)
      - aryl_analysis       : ring count, heteroaromatic flag, ring sizes

    Use to inform catalyst/ligand selection:
      • Low electronic score (electron-rich ring) → harder oxidative addition
        → use Pd(0)/Pd(II) with bulky electron-rich ligands (SPhos, RuPhos)
        or switch to Ni catalysis.
      • High steric score → use bulky monodentate ligand to match steric
        demand; avoid bidentate ligands that clash.
      • Alkyl secondary/tertiary → higher coupling challenge; may need Ni.

    Args:
        smiles: Molecule SMILES (single compound, not a reaction SMILES).
                For reaction SMILES pass the reactant side only.

    Returns:
        dict with ranked_motifs, electronics (list), sterics (list),
        aryl_analysis, and per_site_summary for quick inspection.
    """
    try:
        from chemtools.featurizers.molecule import featurize_molecule as _ff

        raw = _ff(smiles)

        meta = raw.get("meta", {})
        if meta.get("error"):
            return _error(
                f"Molecule featurization failed: {meta['error']}. "
                "Check that SMILES is valid and RDKit is available."
            )

        # --- Electronics (aryl sites) ---
        electronics: List[Dict[str, Any]] = []
        for e in raw.get("electronics", {}).get("aryl", []):
            result = e.get("result") or {}
            electronics.append({
                "site":         e.get("compound_id", ""),
                "score_0_10":   result.get("score_0_10"),
                "description":  result.get("description", ""),
            })

        # --- Sterics (aryl + alkyl) ---
        sterics: List[Dict[str, Any]] = []
        for s in raw.get("steric", {}).get("aryl", []):
            result = s.get("result") or {}
            sterics.append({
                "site":         s.get("compound_id", ""),
                "type":         "aryl",
                "score_0_10":   result.get("score_0_10"),
                "description":  result.get("description", ""),
            })
        for s in raw.get("steric", {}).get("alkyl", []):
            result = s.get("result") or {}
            entry: Dict[str, Any] = {
                "site":           s.get("compound_id", ""),
                "type":           "alkyl",
                "score_0_10":     result.get("score_0_10"),
                "description":    result.get("description", ""),
                "classification": result.get("classification", ""),
                "beta_hydrogens": result.get("beta_hydrogens"),
            }
            sterics.append(entry)

        # --- Aryl ring metrics ---
        aa = raw.get("aryl_analysis") or {}
        aryl_analysis = {
            "aromatic_ring_count": aa.get("aromatic_ring_count", 0),
            "heteroaromatic":      aa.get("heteroaromatic", False),
            "heterocycle_types":   aa.get("heterocycle_types", []),
            "ring_sizes":          aa.get("ring_sizes", []),
        }

        # --- Quick per-site summary for the LLM ---
        site_scores: Dict[str, Dict[str, Any]] = {}
        for e in electronics:
            sid = e["site"]
            site_scores.setdefault(sid, {})["electronic_score"] = e["score_0_10"]
            site_scores[sid]["electronic_desc"]  = e["description"]
        for s in sterics:
            sid = s["site"]
            site_scores.setdefault(sid, {})["steric_score"] = s["score_0_10"]
            site_scores[sid]["steric_desc"]   = s["description"]
            if s.get("classification"):
                site_scores[sid]["alkyl_class"] = s["classification"]

        per_site_summary = [
            {"site": sid, **scores}
            for sid, scores in site_scores.items()
        ]

        return _success({
            "smiles":           smiles,
            "ranked_motifs":    raw.get("ranked_motifs", []),
            "aryl_analysis":    aryl_analysis,
            "electronics":      electronics,
            "sterics":          sterics,
            "per_site_summary": per_site_summary,
        })

    except Exception as exc:
        return _error(f"Molecule featurization failed: {exc}")


featurize_molecule_tool = ToolPlugin(
    name="featurize_molecule",
    category="molecular_features",
    description=(
        "Get full electronic + steric profile for all reactive sites in a molecule: "
        "aryl EWG/EDG score (0=electron-rich → 10=electron-poor), ortho steric bulk "
        "(0=open → 10=highly hindered), alkyl alpha-carbon classification "
        "(methyl/primary/secondary/tertiary), and ring metrics. "
        "Use to decide catalyst family, ligand bulk, and whether Pd or Ni is needed. "
        "Pass a single molecule SMILES (reactant or product, not full reaction SMILES)."
    ),
    prerequisites=[],
    fn=_featurize_molecule,
)


# ---------------------------------------------------------------------------
# Tool 2 — assess_snar_feasibility
# ---------------------------------------------------------------------------

def _assess_snar_feasibility(smiles: str) -> Dict[str, Any]:
    """Check electronic activation of every aryl halide site for SNAr reactions.

    Uses taxonomy-defined thresholds (loaded from reaction_types.json) to score
    each Ar-X site.  A site is 'feasible' only if its scaffold electronic score
    exceeds the minimum_activation threshold (~6.0 on a 0-10 scale).

    Interpretation:
      score < 6.0  → electron-rich ring; SNAr will not proceed without forcing
                     conditions; consider Pd-catalysed C–N coupling instead.
      6.0–7.0      → moderate activation; SNAr may work with strong nucleophile
                     (e.g., alkoxide or secondary amine) at elevated temperature.
      > 7.0        → highly activated; SNAr is favoured; mild conditions suffice.

    Args:
        smiles: Molecule SMILES containing one or more aryl halide groups
                (Ar-F, Ar-Cl, Ar-Br, Ar-I, AromN-halide).

    Returns:
        dict with sites (list per aryl halide), any_feasible (bool),
        and a plain-language summary.
    """
    try:
        from chemtools.featurizers.molecule import featurize_molecule as _ff
        from chemtools.featurizers.analysis.feasibility import (
            analyze_molecule_snar_feasibility,
        )

        mol_payload = _ff(smiles)

        meta = mol_payload.get("meta", {})
        if meta.get("error"):
            return _error(
                f"SNAr feasibility analysis failed: {meta['error']}. "
                "Check that SMILES is valid and RDKit is available."
            )

        sites = analyze_molecule_snar_feasibility(mol_payload)

        if not sites:
            return _success({
                "smiles":       smiles,
                "sites":        [],
                "any_feasible": False,
                "summary": (
                    "No aryl halide sites detected (Ar-F/Cl/Br/I or AromN-halide). "
                    "SNAr is not applicable to this molecule; consider Pd-catalysed "
                    "C–N or C–O coupling instead."
                ),
            })

        feasible = [s for s in sites if s.get("feasible")]
        n = len(sites)
        nf = len(feasible)

        summary = (
            f"{nf}/{n} aryl halide site(s) electronically activated for SNAr. "
        )
        if feasible:
            top = max(feasible, key=lambda s: s.get("score", 0))
            summary += (
                f"Best site: {top['motif']} "
                f"(score={top['score']:.1f}, {top['confidence']} confidence)."
            )
        else:
            summary += (
                "All sites are below the activation threshold — "
                "Pd-catalysed coupling is likely a better choice."
            )

        return _success({
            "smiles":       smiles,
            "sites":        _to_jsonable(sites),
            "any_feasible": bool(feasible),
            "summary":      summary,
        })

    except Exception as exc:
        return _error(f"SNAr feasibility analysis failed: {exc}")


assess_snar_feasibility_tool = ToolPlugin(
    name="assess_snar_feasibility",
    category="molecular_features",
    description=(
        "Check whether each aryl halide site (Ar-F/Cl/Br/I) in a molecule is "
        "electronically activated for nucleophilic aromatic substitution (SNAr). "
        "Returns per-site feasibility (PASS/FAIL), confidence, score 0-10, and reason. "
        "Use when you have an aryl halide and want to decide between SNAr and "
        "Pd-catalysed C–N / C–O coupling. Taxonomy-calibrated thresholds."
    ),
    prerequisites=[],
    fn=_assess_snar_feasibility,
)


# ---------------------------------------------------------------------------
# Module export
# ---------------------------------------------------------------------------

MOLECULAR_FEATURE_TOOLS = [
    featurize_molecule_tool,
    assess_snar_feasibility_tool,
]

"""
Holdout roundtrip benchmark for extracted retrosynthetic templates.

Evaluates template quality by:
  1. Splitting HTE data into train (80%) / holdout (20%) by reaction
  2. Extracting templates from the train set
  3. For each holdout reaction, applying ALL templates to the product
  4. Checking if any template produces precursors matching the actual reactants

Metrics:
  - Coverage:   % of holdout products with at least 1 valid disconnection
  - Precision:  % of holdout products where a template reproduces the actual reactants
  - Top-k hit:   % where the correct precursors appear in the top-k ranked outcomes
  - Reactant recovery: Tanimoto similarity between predicted and actual precursors

Usage:
    python scripts/benchmark_retro_templates.py                     # quick (3 files, 50/file)
    python scripts/benchmark_retro_templates.py --all --limit 500   # medium
    python scripts/benchmark_retro_templates.py --all --limit 0     # full (slow)
    python scripts/benchmark_retro_templates.py --use-existing      # use pre-extracted templates
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, FrozenSet, List, Optional, Set, Tuple

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PROJECT_ROOT))

from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import AllChem

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

HTE_DB_DIR = _PROJECT_ROOT / "data" / "HTE_db" / "literature"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _canonical(smiles: str) -> Optional[str]:
    """Return canonical SMILES or None."""
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def _canonical_set(smiles_block: str) -> Optional[FrozenSet[str]]:
    """Parse dot-separated SMILES into a frozenset of canonical SMILES."""
    parts = [s.strip() for s in smiles_block.split(".") if s.strip()]
    canonical_parts = []
    for p in parts:
        c = _canonical(p)
        if c is None:
            return None
        canonical_parts.append(c)
    if not canonical_parts:
        return None
    return frozenset(canonical_parts)


def _fingerprint(smiles: str):
    """Morgan fingerprint for Tanimoto similarity."""
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def _tanimoto_smiles(smi1: str, smi2: str) -> float:
    """Tanimoto similarity between two SMILES strings."""
    fp1 = _fingerprint(smi1)
    fp2 = _fingerprint(smi2)
    if fp1 is None or fp2 is None:
        return 0.0
    return DataStructs.TanimotoSimilarity(fp1, fp2)


def _max_reactant_similarity(predicted_parts: List[str], actual_parts: List[str]) -> float:
    """Best-match average Tanimoto between predicted and actual reactant sets.

    For each actual reactant, find the best-matching predicted reactant.
    Return the mean of these best matches.
    """
    if not predicted_parts or not actual_parts:
        return 0.0

    pred_fps = []
    for p in predicted_parts:
        fp = _fingerprint(p)
        if fp is not None:
            pred_fps.append(fp)
    if not pred_fps:
        return 0.0

    actual_fps = []
    for a in actual_parts:
        fp = _fingerprint(a)
        if fp is not None:
            actual_fps.append(fp)
    if not actual_fps:
        return 0.0

    best_matches = []
    for afp in actual_fps:
        best = max(DataStructs.TanimotoSimilarity(afp, pfp) for pfp in pred_fps)
        best_matches.append(best)
    return sum(best_matches) / len(best_matches)


def _deterministic_split(
    reactions: List[Dict[str, Any]],
    test_fraction: float = 0.2,
    seed_str: str = "retro_bench_42",
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Deterministic train/test split using reaction SMILES hash.

    Uses MD5 hash of reaction SMILES to assign each reaction to train or test,
    ensuring reproducibility without needing a random seed.
    """
    train, test = [], []
    threshold = int(test_fraction * 1000)
    for rxn in reactions:
        h = hashlib.md5((seed_str + rxn["reaction_smiles"]).encode()).hexdigest()
        bucket = int(h[:3], 16) % 1000  # 0-999
        if bucket < threshold:
            test.append(rxn)
        else:
            train.append(rxn)
    return train, test


# ---------------------------------------------------------------------------
# Template application
# ---------------------------------------------------------------------------

def _clean_template_atom_maps(smarts: str) -> str:
    """Remove atom maps from leaving-group atoms (present only on reactant side).

    rdchiral requires that all mapped atoms in the reactant template have
    corresponding atoms in the product template. Leaving-group atoms (e.g., Br,
    B(OH)₂ in Suzuki) are mapped in our extraction but absent from the product
    pattern, causing KeyError crashes. This function strips those maps.
    """
    if ">>" not in smarts:
        return smarts
    prod_side, rct_side = smarts.split(">>", 1)

    # Extract atom map numbers from each side
    prod_maps = set(int(m) for m in re.findall(r":(\d+)", prod_side))
    rct_maps = set(int(m) for m in re.findall(r":(\d+)", rct_side))

    # Maps that only appear on the reactant side = leaving-group atoms
    lg_maps = rct_maps - prod_maps
    if not lg_maps:
        return smarts

    # Remove leaving-group atom maps from reactant side
    def _strip_map(match):
        mapnum = int(match.group(1))
        if mapnum in lg_maps:
            return ""  # remove the :N token
        return match.group(0)

    cleaned_rct = re.sub(r":(\d+)", _strip_map, rct_side)
    return f"{prod_side}>>{cleaned_rct}"


def _apply_template_to_product(
    template_smarts: str,
    product_smiles: str,
) -> List[List[str]]:
    """Apply a retro template to a product, return list of precursor sets.

    Each outcome is a list of precursor SMILES strings.
    Returns empty list if template cannot be applied.
    """
    try:
        from rdchiral.main import rdchiralReaction, rdchiralRun
        from rdchiral.initialization import rdchiralReactants
    except ImportError:
        return []

    # Clean leaving-group atom maps that crash rdchiral
    cleaned_smarts = _clean_template_atom_maps(template_smarts)

    try:
        rxn = rdchiralReaction(cleaned_smarts)
        # Clean product SMILES (remove atom maps if any)
        clean_smi = re.sub(r":\d+", "", product_smiles)
        with rdBase.BlockLogs():
            mol = Chem.MolFromSmiles(clean_smi)
        if mol is None:
            return []
        canonical_prod = Chem.MolToSmiles(mol)
        rcts = rdchiralReactants(canonical_prod)
        outcomes = rdchiralRun(rxn, rcts)
    except Exception:
        return []

    results = []
    for outcome in outcomes:
        parts = [s.strip() for s in outcome.split(".") if s.strip()]
        # Validate all parts are valid SMILES
        valid_parts = []
        for p in parts:
            c = _canonical(p)
            if c is not None:
                valid_parts.append(c)
        if valid_parts and len(valid_parts) == len(parts):
            results.append(valid_parts)
    return results


def _check_reactant_match(
    predicted_parts: List[str],
    actual_reactant_set: FrozenSet[str],
) -> bool:
    """Check if predicted precursors match the actual reactants.

    Uses canonical SMILES comparison. Allows subset/superset matching
    because templates may not capture all reagent components.
    """
    pred_set = frozenset(predicted_parts)
    if pred_set == actual_reactant_set:
        return True
    # Allow subset match: all actual reactants are found in predicted
    if actual_reactant_set.issubset(pred_set):
        return True
    # Allow superset match: all predicted reactants are found in actual
    if pred_set.issubset(actual_reactant_set):
        return True
    return False


def _inchikey_set(smiles_list) -> FrozenSet[str]:
    """Convert SMILES list to frozenset of InChIKey first blocks (connectivity)."""
    from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
    keys = set()
    for smi in smiles_list:
        with rdBase.BlockLogs():
            mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        try:
            inchi = MolToInchi(mol)
            if inchi:
                ikey = InchiToInchiKey(inchi)
                # Use first block only (connectivity layer, ignores stereo)
                keys.add(ikey.split("-")[0])
        except Exception:
            pass
    return frozenset(keys)


def _check_inchikey_match(predicted_parts: List[str], actual_parts: List[str]) -> bool:
    """Check reactant match using InChIKey connectivity layer (stereo-insensitive)."""
    pred_keys = _inchikey_set(predicted_parts)
    actual_keys = _inchikey_set(actual_parts)
    if not pred_keys or not actual_keys:
        return False
    # At least one actual key appears in predicted
    return bool(actual_keys & pred_keys)


def _count_component_hits(predicted_parts: List[str], actual_parts: List[str]) -> int:
    """Count how many actual reactant components appear in the predicted set."""
    pred_keys = _inchikey_set(predicted_parts)
    actual_keys = _inchikey_set(actual_parts)
    return len(actual_keys & pred_keys)


# ---------------------------------------------------------------------------
# Benchmark core
# ---------------------------------------------------------------------------

def _precompile_templates(
    templates: List[Dict[str, Any]],
) -> List[Tuple[str, Any]]:
    """Pre-compile product SMARTS patterns for fast substructure screening.

    Returns list of (template_smarts, compiled_product_pattern_or_None).
    """
    compiled = []
    for t in templates:
        smarts = t["reaction_smarts"]
        prod_pat = None
        try:
            # Retro template: product >> reactants
            # Product pattern is the LEFT side
            if ">>" in smarts:
                prod_smarts = smarts.split(">>")[0]
                with rdBase.BlockLogs():
                    prod_pat = Chem.MolFromSmarts(prod_smarts)
        except Exception:
            pass
        compiled.append((smarts, prod_pat))
    return compiled


def run_benchmark(
    templates: List[Dict[str, Any]],
    holdout: List[Dict[str, Any]],
    max_templates: int = 0,
    quality_threshold: float = 0.0,
) -> Dict[str, Any]:
    """Run the holdout roundtrip benchmark.

    Args:
        templates: List of template dicts with 'reaction_smarts', 'quality_score', etc.
        holdout: List of reaction dicts with 'reaction_smiles'.
        max_templates: If >0, only use this many top templates.
        quality_threshold: Only use templates with quality_score >= this value.

    Returns dict with all metrics.
    """
    # Filter and sort templates by quality
    active_templates = [
        t for t in templates
        if t.get("quality_score", 0) >= quality_threshold
    ]
    active_templates.sort(key=lambda t: t.get("quality_score", 0), reverse=True)
    if max_templates > 0:
        active_templates = active_templates[:max_templates]

    logger.info(f"Benchmark: {len(active_templates)} templates, {len(holdout)} holdout reactions")

    # Pre-compile product SMARTS for fast substructure screening
    logger.info("Pre-compiling template product patterns...")
    compiled_templates = _precompile_templates(active_templates)
    n_compilable = sum(1 for _, p in compiled_templates if p is not None)
    logger.info(f"  {n_compilable}/{len(compiled_templates)} product patterns compiled for screening")

    # Metrics accumulators
    n_holdout = len(holdout)
    n_coverage = 0           # products with at least 1 valid disconnection
    n_exact_match = 0        # products where template reproduces actual reactants
    n_inchikey_match = 0     # products where InChIKey connectivity matches
    n_partial_match = 0      # products with high similarity (Tanimoto > 0.7)
    n_parseable = 0          # holdout reactions with valid product + reactants
    similarity_scores = []   # best Tanimoto similarity per holdout reaction
    component_hit_scores = [] # best component hit fraction per holdout
    templates_per_hit = []   # how many templates fire per product
    source_type_results: Dict[str, Dict[str, int]] = defaultdict(
        lambda: {"total": 0, "coverage": 0, "exact": 0, "inchikey": 0}
    )

    t0 = time.perf_counter()
    for i, rxn in enumerate(holdout):
        if (i + 1) % max(1, min(100, n_holdout // 10)) == 0:
            elapsed = time.perf_counter() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            logger.info(f"  Progress: {i+1}/{n_holdout} ({rate:.1f} rxn/s)")

        rxn_smi = rxn["reaction_smiles"]
        if ">>" not in rxn_smi:
            continue

        parts = rxn_smi.split(">>")
        product_smi = parts[1].strip()
        reactant_block = parts[0].strip()

        # Canonical product
        canonical_product = _canonical(product_smi)
        if canonical_product is None:
            continue

        # Canonical actual reactants
        actual_reactants = _canonical_set(reactant_block)
        if actual_reactants is None:
            continue
        actual_parts = list(actual_reactants)

        n_parseable += 1
        rxn_type = rxn.get("rxn_type", "unknown")
        source_type_results[rxn_type]["total"] += 1

        # Pre-screen: parse product molecule once
        with rdBase.BlockLogs():
            prod_mol = Chem.MolFromSmiles(canonical_product)
        if prod_mol is None:
            similarity_scores.append(0.0)
            continue

        # Apply templates with substructure pre-screening
        all_outcomes: List[Tuple[List[str], int]] = []
        for t_idx, (smarts, prod_pat) in enumerate(compiled_templates):
            # Fast reject: if product pattern doesn't match, skip rdchiral
            if prod_pat is not None:
                try:
                    if not prod_mol.HasSubstructMatch(prod_pat):
                        continue
                except Exception:
                    pass
            outcomes = _apply_template_to_product(smarts, canonical_product)
            for outcome in outcomes:
                all_outcomes.append((outcome, t_idx))

        if not all_outcomes:
            similarity_scores.append(0.0)
            continue

        # Coverage: at least 1 valid disconnection
        n_coverage += 1
        templates_per_hit.append(len(all_outcomes))
        source_type_results[rxn_type]["coverage"] += 1

        # Check for exact reactant match
        found_exact = False
        found_inchikey = False
        best_sim = 0.0
        best_comp_frac = 0.0
        n_actual_parts = len(actual_parts)
        for outcome_parts, _ in all_outcomes:
            if _check_reactant_match(outcome_parts, actual_reactants):
                found_exact = True
                best_sim = 1.0
                best_comp_frac = 1.0
                break
            if not found_inchikey and _check_inchikey_match(outcome_parts, actual_parts):
                found_inchikey = True
            comp_hits = _count_component_hits(outcome_parts, actual_parts)
            comp_frac = comp_hits / max(n_actual_parts, 1)
            if comp_frac > best_comp_frac:
                best_comp_frac = comp_frac
            sim = _max_reactant_similarity(outcome_parts, actual_parts)
            if sim > best_sim:
                best_sim = sim

        if found_exact:
            n_exact_match += 1
            found_inchikey = True  # exact implies inchikey
            source_type_results[rxn_type]["exact"] += 1
        if found_inchikey:
            n_inchikey_match += 1
            source_type_results[rxn_type]["inchikey"] += 1
        if best_sim >= 0.7:
            n_partial_match += 1
        similarity_scores.append(best_sim)
        component_hit_scores.append(best_comp_frac)

    elapsed = time.perf_counter() - t0

    # Compute metrics
    n_base = max(n_parseable, 1)
    coverage = n_coverage / n_base
    exact_match_rate = n_exact_match / n_base
    inchikey_match_rate = n_inchikey_match / n_base
    partial_match_rate = n_partial_match / n_base
    mean_similarity = sum(similarity_scores) / len(similarity_scores) if similarity_scores else 0
    median_sim = sorted(similarity_scores)[len(similarity_scores)//2] if similarity_scores else 0
    mean_comp_frac = sum(component_hit_scores) / len(component_hit_scores) if component_hit_scores else 0

    # Per-type breakdown
    type_breakdown = {}
    for rtype, counts in sorted(source_type_results.items(), key=lambda x: -x[1]["total"]):
        t_total = max(counts["total"], 1)
        type_breakdown[rtype] = {
            "total": counts["total"],
            "coverage": counts["coverage"],
            "coverage_pct": round(100 * counts["coverage"] / t_total, 1),
            "exact_match": counts["exact"],
            "exact_match_pct": round(100 * counts["exact"] / t_total, 1),
            "inchikey_match": counts["inchikey"],
            "inchikey_match_pct": round(100 * counts["inchikey"] / t_total, 1),
        }

    return {
        "config": {
            "n_templates": len(active_templates),
            "n_holdout": n_holdout,
            "n_parseable": n_parseable,
            "quality_threshold": quality_threshold,
            "max_templates": max_templates,
        },
        "metrics": {
            "coverage": round(coverage, 4),
            "coverage_pct": round(100 * coverage, 1),
            "exact_match_rate": round(exact_match_rate, 4),
            "exact_match_pct": round(100 * exact_match_rate, 1),
            "inchikey_match_rate": round(inchikey_match_rate, 4),
            "inchikey_match_pct": round(100 * inchikey_match_rate, 1),
            "partial_match_rate_70": round(partial_match_rate, 4),
            "partial_match_pct_70": round(100 * partial_match_rate, 1),
            "mean_similarity": round(mean_similarity, 4),
            "median_similarity": round(median_sim, 4),
            "mean_component_hit_fraction": round(mean_comp_frac, 4),
            "mean_templates_per_hit": round(
                sum(templates_per_hit) / len(templates_per_hit), 1
            ) if templates_per_hit else 0,
        },
        "per_reaction_type": type_breakdown,
        "elapsed_seconds": round(elapsed, 1),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Holdout roundtrip benchmark for retro templates")
    parser.add_argument("--all", action="store_true", help="Use all CSV files")
    parser.add_argument("--file", type=str, default="", help="Filter CSV by name")
    parser.add_argument("--limit", type=int, default=50, help="Max rows per CSV (0=unlimited)")
    parser.add_argument("--use-existing", action="store_true",
                        help="Use pre-extracted templates from results/extracted_templates.json")
    parser.add_argument("--min-confidence", type=float, default=0.5)
    parser.add_argument("--min-count", type=int, default=3)
    parser.add_argument("--quality-threshold", type=float, default=0.0,
                        help="Only use templates with quality_score >= threshold")
    parser.add_argument("--max-templates", type=int, default=0,
                        help="Max templates to use (0=all)")
    parser.add_argument("--output", type=str, default="")
    args = parser.parse_args()

    # Import extraction functions
    from extract_retro_templates import (
        load_reactions_from_csv,
        extract_templates_from_reactions,
        _quality_score,
    )

    # Select CSV files
    csv_files = sorted(HTE_DB_DIR.glob("*.csv"))
    if args.file:
        pattern = args.file.lower()
        csv_files = [f for f in csv_files if pattern in f.name.lower()]
    elif not args.all:
        test_names = ["suzuki_miyaura", "C_N_Coupling", "Amide_formation", "Sulfonamide_formation"]
        csv_files = [
            f for f in csv_files
            if any(tn.lower() in f.name.lower() for tn in test_names)
        ]

    if not csv_files:
        logger.error("No CSV files found")
        sys.exit(1)

    # Load all reactions
    all_reactions: List[Dict[str, Any]] = []
    for csv_path in csv_files:
        rows = load_reactions_from_csv(csv_path, limit=args.limit)
        logger.info(f"  {csv_path.name}: {len(rows)} reactions")
        all_reactions.append(rows)

    flat_reactions = [r for batch in all_reactions for r in batch]
    logger.info(f"Total reactions loaded: {len(flat_reactions)}")

    # Deterministic 80/20 split
    train_rxns, test_rxns = _deterministic_split(flat_reactions, test_fraction=0.2)
    logger.info(f"Train: {len(train_rxns)}, Holdout: {len(test_rxns)}")

    # Get templates
    if args.use_existing:
        # Load pre-extracted templates
        templates_path = _PROJECT_ROOT / "results" / "extracted_templates.json"
        with open(templates_path) as f:
            data = json.load(f)
        templates = data["templates"]
        logger.info(f"Loaded {len(templates)} pre-extracted templates")
    else:
        # Extract templates from train set only
        logger.info("Extracting templates from train set...")
        t0 = time.perf_counter()
        template_stats, counters = extract_templates_from_reactions(
            train_rxns,
            min_confidence=args.min_confidence,
            validate_retro=True,
            min_count=args.min_count,
        )
        elapsed = time.perf_counter() - t0
        logger.info(f"Extracted {len(template_stats)} templates in {elapsed:.1f}s")

        # Convert to dicts with quality scores
        templates = []
        for entry in template_stats.values():
            d = entry.to_dict()
            d["quality_score"] = round(_quality_score(entry), 3)
            d["avg_confidence"] = round(
                sum(getattr(entry, '_confidences', [0.5])) /
                max(len(getattr(entry, '_confidences', [0.5])), 1), 3
            )
            templates.append(d)

    # Run benchmark
    print()
    print("=" * 70)
    print("HOLDOUT ROUNDTRIP BENCHMARK")
    print("=" * 70)

    results = run_benchmark(
        templates,
        test_rxns,
        max_templates=args.max_templates,
        quality_threshold=args.quality_threshold,
    )

    # Print results
    m = results["metrics"]
    print()
    print(f"Templates used:          {results['config']['n_templates']}")
    print(f"Holdout reactions:       {results['config']['n_parseable']}")
    print()
    print(f"Coverage:                {m['coverage_pct']:.1f}%  (products with >= 1 disconnection)")
    print(f"Exact match (SMILES):    {m['exact_match_pct']:.1f}%  (predicted == actual reactants)")
    print(f"InChIKey match:          {m['inchikey_match_pct']:.1f}%  (any reactant connectivity match)")
    print(f"Partial match (Tan>0.7): {m['partial_match_pct_70']:.1f}%")
    print(f"Mean similarity:         {m['mean_similarity']:.4f}")
    print(f"Median similarity:       {m['median_similarity']:.4f}")
    print(f"Mean component hit frac: {m['mean_component_hit_fraction']:.4f}")
    print(f"Templates/hit (mean):    {m['mean_templates_per_hit']}")
    print(f"Time:                    {results['elapsed_seconds']:.1f}s")

    # Per-type breakdown
    print()
    print("Per reaction type:")
    print(f"  {'Type':<40s} {'Total':>6s} {'Cov%':>6s} {'Exact%':>7s} {'IKey%':>6s}")
    print("  " + "-" * 63)
    for rtype, info in sorted(
        results["per_reaction_type"].items(),
        key=lambda x: -x[1]["total"]
    )[:25]:
        print(f"  {rtype:<40s} {info['total']:>6d} {info['coverage_pct']:>5.1f}% "
              f"{info['exact_match_pct']:>6.1f}% {info['inchikey_match_pct']:>5.1f}%")

    # Save results
    output_path = args.output or str(_PROJECT_ROOT / "results" / "retro_template_benchmark.json")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    logger.info(f"Results saved to {output_path}")


if __name__ == "__main__":
    main()

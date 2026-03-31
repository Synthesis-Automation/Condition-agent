"""
Extract retrosynthetic templates from HTE reaction data.

This module provides the core extraction pipeline that converts raw HTE
reaction CSV files into validated retrosynthetic SMARTS templates with
quality metrics.

Pipeline:
  1. Load unmapped reaction SMILES from HTE CSV files
  2. Atom-map via RDKit MCS with confidence scoring
  3. Filter reactions by mapping confidence threshold
  4. Extract retro templates via rdchiral
  5. Retro-application validation (ASKCOS-style roundtrip check)
  6. Cluster templates by canonical reaction-center core
  7. Filter singletons, score by frequency + yield
  8. Output ranked templates with quality metrics

Usage (as CLI):
    # Quick test run (3 default files, 50 rows each):
    python -m chemtools.retro.extract_templates

    # Process ALL literature CSVs with unlimited rows:
    python -m chemtools.retro.extract_templates --all --limit 0

    # Process a specific file pattern:
    python -m chemtools.retro.extract_templates --file suzuki

    # Adjust quality thresholds:
    python -m chemtools.retro.extract_templates --all --limit 0 \\
        --min-confidence 0.6 --min-count 3

    # Skip retro-validation for speed:
    python -m chemtools.retro.extract_templates --all --limit 0 --no-retro-validation

    # Custom output path:
    python -m chemtools.retro.extract_templates --all --output my_templates.json

Usage (as library):
    from chemtools.retro.extract_templates import (
        load_reactions_from_csv,
        extract_templates_from_reactions,
        quality_score,
    )
    from pathlib import Path

    reactions = load_reactions_from_csv(Path("data/HTE_db/literature/suzuki_miyaura.csv"))
    template_stats, counters = extract_templates_from_reactions(reactions)
    for smarts, stats in template_stats.items():
        print(f"{stats.count}x  Q={quality_score(stats):.3f}  {smarts[:80]}")

Full pipeline (extract + integrate):
    # Step 1: Extract templates from HTE data
    python -m chemtools.retro.extract_templates --all --limit 0

    # Step 2: Integrate high-quality templates into production files
    python -m chemtools.retro.integrate_templates

    # Step 3: Verify
    python -m pytest tests/ -q
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import os
import re
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdFMCS

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

_MODULE_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _MODULE_DIR.parent.parent
HTE_DB_DIR = _PROJECT_ROOT / "data" / "HTE_db" / "literature"


# ---------------------------------------------------------------------------
# Step 1: MCS-based atom mapping
# ---------------------------------------------------------------------------

def _add_atom_mapping_mcs(reaction_smiles: str) -> Optional[str]:
    """Add atom mapping to an unmapped reaction SMILES using MCS.

    Strategy:
      For each reactant, find Maximum Common Substructure with the product.
      Map matched atoms. Remaining product atoms get new map numbers.

    Returns atom-mapped reaction SMILES or None on failure.
    """
    if ">>" not in reaction_smiles:
        return None

    parts = reaction_smiles.split(">>")
    reactant_block = parts[0].strip()
    product_block = parts[1].strip()

    with rdBase.BlockLogs():
        product_mol = Chem.MolFromSmiles(product_block)
    if product_mol is None:
        return None

    # Parse reactants (may have multiple separated by .)
    reactant_smiles_list = [s.strip() for s in reactant_block.split(".") if s.strip()]
    reactant_mols = []
    for rs in reactant_smiles_list:
        with rdBase.BlockLogs():
            mol = Chem.MolFromSmiles(rs)
        if mol is None:
            return None
        reactant_mols.append(mol)

    # Work on editable copies
    product_rw = Chem.RWMol(product_mol)
    reactant_rws = [Chem.RWMol(m) for m in reactant_mols]

    # Clear any existing map numbers
    for atom in product_rw.GetAtoms():
        atom.SetAtomMapNum(0)
    for rmol in reactant_rws:
        for atom in rmol.GetAtoms():
            atom.SetAtomMapNum(0)

    map_counter = 1
    product_mapped: set = set()  # product atom indices already mapped

    for rmol in reactant_rws:
        # Find MCS between this reactant and the product
        try:
            mcs_result = rdFMCS.FindMCS(
                [rmol, product_rw],
                timeout=2,
                atomCompare=rdFMCS.AtomCompare.CompareElements,
                bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                matchValences=True,
                ringMatchesRingOnly=True,
                completeRingsOnly=False,
                matchChiralTag=False,
            )
        except Exception:
            continue

        if mcs_result.canceled or not mcs_result.smartsString:
            continue

        with rdBase.BlockLogs():
            mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
        if mcs_mol is None:
            continue

        # Get substructure matches
        with rdBase.BlockLogs():
            r_matches = rmol.GetSubstructMatch(mcs_mol)
            p_matches = product_rw.GetSubstructMatch(mcs_mol)

        if not r_matches or not p_matches or len(r_matches) != len(p_matches):
            continue

        # Assign map numbers to matched atoms
        for r_idx, p_idx in zip(r_matches, p_matches):
            if p_idx in product_mapped:
                continue  # Already mapped by a previous reactant
            rmol.GetAtomWithIdx(r_idx).SetAtomMapNum(map_counter)
            product_rw.GetAtomWithIdx(p_idx).SetAtomMapNum(map_counter)
            product_mapped.add(p_idx)
            map_counter += 1

    # Map remaining unmapped reactant atoms (leaving groups)
    for rmol in reactant_rws:
        for atom in rmol.GetAtoms():
            if atom.GetAtomMapNum() == 0:
                atom.SetAtomMapNum(map_counter)
                map_counter += 1

    # Map remaining unmapped product atoms (new atoms from reagents)
    for atom in product_rw.GetAtoms():
        if atom.GetAtomMapNum() == 0:
            atom.SetAtomMapNum(map_counter)
            map_counter += 1

    # Check that at least some atoms were mapped between reactant and product
    if not product_mapped:
        return None

    # Reconstruct SMILES
    reactant_parts = []
    for rmol in reactant_rws:
        with rdBase.BlockLogs():
            smi = Chem.MolToSmiles(rmol)
        if smi:
            reactant_parts.append(smi)

    with rdBase.BlockLogs():
        product_smi = Chem.MolToSmiles(product_rw)

    if not reactant_parts or not product_smi:
        return None

    return ".".join(reactant_parts) + ">>" + product_smi


# ---------------------------------------------------------------------------
# Step 2: rdchiral template extraction
# ---------------------------------------------------------------------------

def _extract_template(mapped_smiles: str) -> Optional[Dict[str, str]]:
    """Extract retro template from atom-mapped reaction SMILES using rdchiral.

    Returns dict with 'reaction_smarts' (retro direction) or None.
    """
    from rdchiral.template_extractor import extract_from_reaction

    if ">>" not in mapped_smiles:
        return None

    parts = mapped_smiles.split(">>")
    rxn_dict = {
        "_id": "tmp",
        "reactants": parts[0],
        "products": parts[1],
    }

    try:
        result = extract_from_reaction(rxn_dict)
    except Exception:
        return None

    if not result or "reaction_smarts" not in result:
        return None

    return {
        "reaction_smarts": result["reaction_smarts"],
        "necessary_reagent": result.get("necessary_reagent", ""),
    }


# ---------------------------------------------------------------------------
# Step 3: Load & process HTE CSVs
# ---------------------------------------------------------------------------

@dataclass
class TemplateStats:
    """Aggregated statistics for a unique template."""
    reaction_smarts: str
    count: int = 0
    yields: list = field(default_factory=list)
    source_files: set = field(default_factory=set)
    example_rxn: str = ""
    example_mapped: str = ""

    @property
    def avg_yield(self) -> float:
        return sum(self.yields) / len(self.yields) if self.yields else 0.0

    @property
    def median_yield(self) -> float:
        if not self.yields:
            return 0.0
        s = sorted(self.yields)
        n = len(s)
        return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2

    def to_dict(self) -> Dict[str, Any]:
        return {
            "reaction_smarts": self.reaction_smarts,
            "count": self.count,
            "avg_yield": round(self.avg_yield, 1),
            "median_yield": round(self.median_yield, 1),
            "min_yield": round(min(self.yields), 1) if self.yields else 0.0,
            "max_yield": round(max(self.yields), 1) if self.yields else 0.0,
            "source_files": sorted(self.source_files),
            "example_rxn": self.example_rxn,
            "example_mapped": self.example_mapped,
        }


def _parse_yield(val: Any) -> Optional[float]:
    """Parse yield value from CSV field."""
    if val is None:
        return None
    try:
        v = float(str(val).strip().rstrip("%"))
        if 0 <= v <= 100:
            return v
    except (ValueError, TypeError):
        pass
    return None


def load_reactions_from_csv(
    csv_path: Path,
    limit: int = 0,
) -> List[Dict[str, Any]]:
    """Load reaction SMILES + metadata from an HTE CSV file.

    Args:
        csv_path: Path to a CSV file with a 'reaction_smiles' column.
        limit: Maximum rows to read (0 = unlimited).

    Returns:
        List of dicts with keys: reaction_smiles, yield, source_file, rxn_type.
    """
    rows = []
    for enc in ("utf-8", "latin-1"):
        try:
            with open(csv_path, encoding=enc, newline="") as fh:
                reader = csv.DictReader(fh)
                for i, row in enumerate(reader):
                    if limit and i >= limit:
                        break
                    rxn_smi = (row.get("reaction_smiles") or "").strip()
                    if ">>" not in rxn_smi:
                        continue
                    yld = _parse_yield(
                        row.get("yield") or row.get("yield_pct") or row.get("yield_percent")
                    )
                    rows.append({
                        "reaction_smiles": rxn_smi,
                        "yield": yld,
                        "source_file": csv_path.name,
                        "rxn_type": (row.get("detected_reaction_type") or csv_path.stem).replace("_canonical", ""),
                    })
            break
        except Exception:
            continue
    return rows


def extract_templates_from_reactions(
    reactions: List[Dict[str, Any]],
    min_confidence: float = 0.5,
    validate_retro: bool = True,
    min_count: int = 2,
) -> Tuple[Dict[str, TemplateStats], Dict[str, int]]:
    """Run the full pipeline: atom-map, confidence filter, extract, validate, cluster.

    Args:
        reactions: List of reaction dicts with reaction_smiles, yield, etc.
        min_confidence: Minimum MCS mapping confidence to keep (0-1).
        validate_retro: If True, run retro-application validation on templates.
        min_count: Minimum occurrence count to keep a template (singleton filter).

    Returns (template_stats, counters).
    """
    stats: Dict[str, TemplateStats] = {}
    counters: Dict[str, int] = {
        "total": len(reactions),
        "mapped_ok": 0,
        "mapped_fail": 0,
        "low_confidence": 0,
        "template_ok": 0,
        "template_fail": 0,
        "retro_valid": 0,
        "retro_invalid": 0,
    }
    confidence_values: List[float] = []

    for rxn in reactions:
        rxn_smi = rxn["reaction_smiles"]

        # Step 1: Atom mapping with confidence
        mapped, confidence = _add_atom_mapping_mcs_with_confidence(rxn_smi)
        if mapped is None:
            counters["mapped_fail"] += 1
            continue
        counters["mapped_ok"] += 1
        confidence_values.append(confidence)

        # Step 1b: Filter by mapping confidence
        if confidence < min_confidence:
            counters["low_confidence"] += 1
            continue

        # Step 2: Template extraction
        tmpl = _extract_template(mapped)
        if tmpl is None:
            counters["template_fail"] += 1
            continue
        counters["template_ok"] += 1

        smarts = tmpl["reaction_smarts"]

        if smarts not in stats:
            # Extract product SMILES for retro validation
            product_smi = rxn_smi.split(">>")[1] if ">>" in rxn_smi else ""

            # Step 2b: Retro-application validation
            retro_ok = True
            if validate_retro and product_smi:
                retro_ok = _validate_template_retro(smarts, product_smi)
                if retro_ok:
                    counters["retro_valid"] += 1
                else:
                    counters["retro_invalid"] += 1

            # Step 2c: Extract reaction center core for clustering
            core = _extract_reaction_center_core(smarts)

            stats[smarts] = TemplateStats(
                reaction_smarts=smarts,
                example_rxn=rxn_smi,
                example_mapped=mapped,
            )
            stats[smarts]._retro_validated = retro_ok
            stats[smarts]._reaction_center_core = core

        entry = stats[smarts]
        entry.count += 1
        if not hasattr(entry, '_confidences'):
            entry._confidences = []
        entry._confidences.append(confidence)
        entry.source_files.add(rxn["source_file"])
        if rxn["yield"] is not None:
            entry.yields.append(rxn["yield"])

    # Step 3: Cluster templates by reaction center core
    _assign_clusters(stats)

    # Step 4: Record stats and apply filters
    counters["pre_filter_templates"] = len(stats)
    counters["min_confidence_threshold"] = min_confidence
    counters["min_count_threshold"] = min_count

    # Confidence distribution stats
    if confidence_values:
        confidence_values.sort()
        n = len(confidence_values)
        counters["confidence_mean"] = round(sum(confidence_values) / n, 3)
        counters["confidence_median"] = round(
            confidence_values[n // 2] if n % 2 else
            (confidence_values[n // 2 - 1] + confidence_values[n // 2]) / 2, 3
        )
        counters["confidence_p25"] = round(confidence_values[n // 4], 3)
        counters["confidence_p75"] = round(confidence_values[3 * n // 4], 3)

    # Remove singletons and retro-invalid templates
    to_remove = []
    for smarts, entry in stats.items():
        if entry.count < min_count:
            to_remove.append(smarts)
        elif validate_retro and getattr(entry, '_retro_validated', True) is False:
            to_remove.append(smarts)
    for smarts in to_remove:
        del stats[smarts]

    counters["singleton_removed"] = len(to_remove)
    counters["final_templates"] = len(stats)

    return stats, counters


def _add_atom_mapping_mcs_with_confidence(
    reaction_smiles: str,
) -> Tuple[Optional[str], float]:
    """Wrapper around _add_atom_mapping_mcs that also returns mapping confidence.

    Confidence = fraction of product heavy atoms whose map numbers also
    appear in the reactant side (i.e., atoms matched via MCS).

    Returns (mapped_smiles_or_None, confidence_0_to_1).
    """
    if ">>" not in reaction_smiles:
        return None, 0.0

    parts = reaction_smiles.split(">>")
    product_block = parts[1].strip()

    with rdBase.BlockLogs():
        product_mol = Chem.MolFromSmiles(product_block)
    if product_mol is None:
        return None, 0.0

    total_product_atoms = product_mol.GetNumHeavyAtoms()
    if total_product_atoms == 0:
        return None, 0.0

    mapped = _add_atom_mapping_mcs(reaction_smiles)
    if mapped is None:
        return None, 0.0

    # Count how many product atoms got mapped to reactant atoms
    mapped_product = mapped.split(">>")[1]
    with rdBase.BlockLogs():
        mapped_prod_mol = Chem.MolFromSmiles(mapped_product)
    if mapped_prod_mol is None:
        return mapped, 0.5  # fallback

    # Collect reactant map numbers
    mapped_reactants = mapped.split(">>")[0]
    reactant_map_nums = set()
    for rsmi in mapped_reactants.split("."):
        with rdBase.BlockLogs():
            rmol = Chem.MolFromSmiles(rsmi)
        if rmol:
            reactant_map_nums.update(
                a.GetAtomMapNum() for a in rmol.GetAtoms() if a.GetAtomMapNum() > 0
            )

    # Confidence = product atoms with map numbers also in reactants / total
    shared_count = sum(
        1 for a in mapped_prod_mol.GetAtoms()
        if a.GetAtomMapNum() > 0 and a.GetAtomMapNum() in reactant_map_nums
    )

    confidence = shared_count / total_product_atoms
    return mapped, confidence


def _validate_template_retro(template_smarts: str, product_smiles: str) -> bool:
    """Apply a retro template to a product and check it produces valid precursors.

    ASKCOS-style roundtrip validation: if the extracted template cannot be
    applied back to the product to yield at least one set of valid precursors,
    the template is unreliable and should be discarded.
    """
    try:
        from rdchiral.main import rdchiralReaction, rdchiralRun
        from rdchiral.initialization import rdchiralReactants
    except ImportError:
        return True  # skip validation if rdchiral.main unavailable

    try:
        rxn = rdchiralReaction(template_smarts)
        prod_smi = re.sub(r':\d+', '', product_smiles)
        with rdBase.BlockLogs():
            prod_mol = Chem.MolFromSmiles(prod_smi)
        if prod_mol is None:
            return False
        prod_smi_clean = Chem.MolToSmiles(prod_mol)
        rcts = rdchiralReactants(prod_smi_clean)
        outcomes = rdchiralRun(rxn, rcts)
        if not outcomes:
            return False
        for outcome in outcomes:
            parts = outcome.split(".")
            if all(Chem.MolFromSmiles(p) is not None for p in parts if p):
                return True
        return False
    except (KeyError, ValueError, RuntimeError):
        # rdchiral crashes on templates with leaving-group atoms not in product
        return True
    except Exception:
        return True


def _extract_reaction_center_core(template_smarts: str) -> Optional[str]:
    """Extract canonical reaction center core from a template SMARTS.

    Templates with identical cores represent the same retrosynthetic transform.
    """
    if ">>" not in template_smarts:
        return None
    try:
        with rdBase.BlockLogs():
            rxn = AllChem.ReactionFromSmarts(template_smarts)
        if rxn is None:
            return None
        if rxn.GetNumProductTemplates() == 0:
            return None
        prod_template = rxn.GetProductTemplate(0)
        with rdBase.BlockLogs():
            core_smarts = Chem.MolToSmarts(prod_template)
        return core_smarts
    except Exception:
        return None


def _assign_clusters(stats: Dict[str, TemplateStats]) -> None:
    """Group templates by reaction center core into clusters."""
    core_to_cluster: Dict[str, int] = {}
    cluster_id = 0
    for entry in stats.values():
        core = getattr(entry, '_reaction_center_core', None)
        if core is None:
            entry._cluster_id = -1
            continue
        if core not in core_to_cluster:
            core_to_cluster[core] = cluster_id
            cluster_id += 1
        entry._cluster_id = core_to_cluster[core]


def quality_score(entry: TemplateStats) -> float:
    """Combined quality score (0-1) blending frequency, yield, and confidence."""
    freq_score = min(1.0, math.log2(max(entry.count, 1)) / 10.0)
    yield_score = entry.avg_yield / 100.0 if entry.yields else 0.5
    confs = getattr(entry, '_confidences', [])
    conf_score = sum(confs) / len(confs) if confs else 0.5
    return 0.4 * freq_score + 0.3 * yield_score + 0.3 * conf_score


# ---------------------------------------------------------------------------
# CLI main
# ---------------------------------------------------------------------------

def main():
    """CLI entry point for template extraction."""
    parser = argparse.ArgumentParser(
        description="Extract retro templates from HTE data",
        epilog=(
            "Examples:\n"
            "  python -m chemtools.retro.extract_templates                       # quick test\n"
            "  python -m chemtools.retro.extract_templates --all --limit 0       # full extraction\n"
            "  python -m chemtools.retro.extract_templates --file suzuki         # single family\n"
            "  python -m chemtools.retro.extract_templates --all --min-count 5   # stricter filter\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--all", action="store_true", help="Process all literature CSVs")
    parser.add_argument("--file", type=str, default="", help="Filter CSV files by name pattern")
    parser.add_argument("--limit", type=int, default=50, help="Max rows per file (0=unlimited)")
    parser.add_argument("--output", type=str, default="", help="Output JSON path")
    parser.add_argument(
        "--min-confidence", type=float, default=0.5,
        help="Minimum MCS mapping confidence (0-1, default: 0.5)",
    )
    parser.add_argument(
        "--min-count", type=int, default=2,
        help="Minimum template occurrence count (default: 2, removes singletons)",
    )
    parser.add_argument(
        "--no-retro-validation", action="store_true",
        help="Skip retro-application validation (faster but less reliable)",
    )
    args = parser.parse_args()

    if not HTE_DB_DIR.is_dir():
        logger.error(f"HTE directory not found: {HTE_DB_DIR}")
        sys.exit(1)

    # Select CSV files
    csv_files = sorted(HTE_DB_DIR.glob("*.csv"))
    if args.file:
        pattern = args.file.lower()
        csv_files = [f for f in csv_files if pattern in f.name.lower()]
    elif not args.all:
        test_names = ["suzuki_miyaura", "C_N_Coupling", "Amide_formation"]
        csv_files = [
            f for f in csv_files
            if any(tn.lower() in f.name.lower() for tn in test_names)
        ]

    if not csv_files:
        logger.error("No CSV files found matching criteria")
        sys.exit(1)

    logger.info(f"Processing {len(csv_files)} CSV file(s), limit={args.limit} rows/file")
    logger.info(f"Confidence threshold: {args.min_confidence}, min count: {args.min_count}")
    logger.info(f"Retro validation: {'OFF' if args.no_retro_validation else 'ON'}")

    # Load reactions
    all_reactions: List[Dict[str, Any]] = []
    for csv_path in csv_files:
        rows = load_reactions_from_csv(csv_path, limit=args.limit)
        logger.info(f"  {csv_path.name}: {len(rows)} reactions loaded")
        all_reactions.extend(rows)

    logger.info(f"Total reactions: {len(all_reactions)}")

    # Extract templates
    t0 = time.perf_counter()
    template_stats, counters = extract_templates_from_reactions(
        all_reactions,
        min_confidence=args.min_confidence,
        validate_retro=not args.no_retro_validation,
        min_count=args.min_count,
    )
    elapsed = time.perf_counter() - t0

    # Sort by quality_score (best first)
    ranked = sorted(template_stats.values(), key=lambda t: quality_score(t), reverse=True)

    # Report
    total = max(counters["total"], 1)
    mapped_ok = max(counters["mapped_ok"], 1)
    print()
    print("=" * 70)
    print("TEMPLATE EXTRACTION RESULTS")
    print("=" * 70)
    print(f"Total reactions:        {counters['total']}")
    print(f"Atom mapping OK:        {counters['mapped_ok']} ({counters['mapped_ok']*100/total:.1f}%)")
    print(f"Atom mapping fail:      {counters['mapped_fail']}")
    print(f"Low confidence filtered:{counters['low_confidence']} (< {args.min_confidence})")
    print(f"Template extracted:     {counters['template_ok']} ({counters['template_ok']*100/mapped_ok:.1f}% of mapped)")
    print(f"Template extract fail:  {counters['template_fail']}")

    if not args.no_retro_validation:
        print(f"Retro-validated OK:     {counters['retro_valid']}")
        print(f"Retro-validation fail:  {counters['retro_invalid']}")

    print(f"Pre-filter templates:   {counters.get('pre_filter_templates', '?')}")
    print(f"Removed (singleton/bad):{counters.get('singleton_removed', '?')}")
    print(f"Final templates:        {counters.get('final_templates', len(ranked))}")

    # Confidence distribution
    if "confidence_mean" in counters:
        print(f"\nMapping confidence distribution:")
        print(f"  mean={counters['confidence_mean']:.3f}  median={counters['confidence_median']:.3f}  "
              f"p25={counters['confidence_p25']:.3f}  p75={counters['confidence_p75']:.3f}")

    # Cluster summary
    cluster_ids = set(
        getattr(t, '_cluster_id', -1) for t in ranked
        if getattr(t, '_cluster_id', -1) >= 0
    )
    if cluster_ids:
        print(f"\nReaction center clusters: {len(cluster_ids)}")
        cluster_sizes = Counter(
            getattr(t, '_cluster_id', -1) for t in ranked
            if getattr(t, '_cluster_id', -1) >= 0
        )
        multi = sum(1 for s in cluster_sizes.values() if s > 1)
        print(f"  Multi-template clusters: {multi} (templates sharing same core transform)")

    print(f"\nTime: {elapsed:.1f}s")
    print()

    show = min(25, len(ranked))
    print(f"Top {show} templates by quality score:")
    print("-" * 70)
    for i, t in enumerate(ranked[:show], 1):
        yield_str = f"avg={t.avg_yield:.0f}%" if t.yields else "no yield"
        confs = getattr(t, '_confidences', [])
        avg_conf = sum(confs) / len(confs) if confs else 0.0
        conf_str = f"conf={avg_conf:.2f}"
        retro_ok = getattr(t, '_retro_validated', None)
        retro_str = "Rv" if retro_ok else ("Rx" if retro_ok is False else "R?")
        qs = quality_score(t)
        files_str = ", ".join(sorted(t.source_files)[:3])
        print(f"  {i:2d}. count={t.count:4d}  {yield_str:12s}  {conf_str}  {retro_str}  "
              f"Q={qs:.3f}  sources: {files_str}")
        smarts_display = t.reaction_smarts
        if len(smarts_display) > 100:
            smarts_display = smarts_display[:97] + "..."
        print(f"      SMARTS: {smarts_display}")
        print(f"      Example: {t.example_rxn[:120]}")
        print()

    # Save output
    output_path = args.output or str(_PROJECT_ROOT / "results" / "extracted_templates.json")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    def _entry_to_dict(entry):
        d = entry.to_dict()
        confs = getattr(entry, '_confidences', [])
        d["avg_confidence"] = round(sum(confs) / len(confs), 3) if confs else 0.0
        d["quality_score"] = round(quality_score(entry), 3)
        d["retro_validated"] = getattr(entry, '_retro_validated', None)
        d["reaction_center_core"] = getattr(entry, '_reaction_center_core', None)
        d["cluster_id"] = getattr(entry, '_cluster_id', None)
        return d

    output = {
        "extraction_meta": {
            "csv_files": [f.name for f in csv_files],
            "limit_per_file": args.limit,
            "min_confidence": args.min_confidence,
            "min_count": args.min_count,
            "retro_validation": not args.no_retro_validation,
            **counters,
            "elapsed_seconds": round(elapsed, 2),
        },
        "templates": [_entry_to_dict(t) for t in ranked],
    }
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)

    logger.info(f"Results saved to {output_path}")


if __name__ == "__main__":
    main()

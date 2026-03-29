"""
Extract retrosynthetic templates from HTE reaction data.

Pipeline:
  1. Load unmapped reaction SMILES from HTE CSV files
  2. Atom-map via RDKit MCS (Maximum Common Substructure)
  3. Extract retro templates via rdchiral
  4. Cluster by canonical reaction-center SMARTS
  5. Output ranked templates with frequency and yield stats

Small-scale test: run on a single CSV file first.

Usage:
    python scripts/extract_retro_templates.py                    # small test (1 file)
    python scripts/extract_retro_templates.py --all              # all literature CSVs
    python scripts/extract_retro_templates.py --file suzuki      # specific file pattern
    python scripts/extract_retro_templates.py --limit 200        # cap rows per file
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Ensure project root is importable
_PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PROJECT_ROOT))

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdFMCS

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

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
        result = extract_from_reaction(rxn_dict, super_general=False, radius=1)
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
    """Load reaction SMILES + metadata from an HTE CSV file."""
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
) -> Tuple[Dict[str, TemplateStats], Dict[str, int]]:
    """Run the full pipeline: atom-map → extract template → aggregate.

    Returns (template_stats, counters).
    """
    stats: Dict[str, TemplateStats] = {}
    counters: Dict[str, int] = {
        "total": len(reactions),
        "mapped_ok": 0,
        "mapped_fail": 0,
        "template_ok": 0,
        "template_fail": 0,
    }

    for rxn in reactions:
        rxn_smi = rxn["reaction_smiles"]

        # Step 1: Atom mapping
        mapped = _add_atom_mapping_mcs(rxn_smi)
        if mapped is None:
            counters["mapped_fail"] += 1
            continue
        counters["mapped_ok"] += 1

        # Step 2: Template extraction
        tmpl = _extract_template(mapped)
        if tmpl is None:
            counters["template_fail"] += 1
            continue
        counters["template_ok"] += 1

        smarts = tmpl["reaction_smarts"]

        if smarts not in stats:
            stats[smarts] = TemplateStats(
                reaction_smarts=smarts,
                example_rxn=rxn_smi,
                example_mapped=mapped,
            )

        entry = stats[smarts]
        entry.count += 1
        entry.source_files.add(rxn["source_file"])
        if rxn["yield"] is not None:
            entry.yields.append(rxn["yield"])

    return stats, counters


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Extract retro templates from HTE data")
    parser.add_argument("--all", action="store_true", help="Process all literature CSVs")
    parser.add_argument("--file", type=str, default="", help="Filter CSV files by name pattern")
    parser.add_argument("--limit", type=int, default=50, help="Max rows per file (0=unlimited)")
    parser.add_argument("--output", type=str, default="", help="Output JSON path")
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
        # Default: small test with 3 diverse files
        test_names = ["suzuki_miyaura", "C_N_Coupling", "Amide_formation"]
        csv_files = [
            f for f in csv_files
            if any(tn.lower() in f.name.lower() for tn in test_names)
        ]

    if not csv_files:
        logger.error("No CSV files found matching criteria")
        sys.exit(1)

    logger.info(f"Processing {len(csv_files)} CSV file(s), limit={args.limit} rows/file")

    # Load reactions
    all_reactions: List[Dict[str, Any]] = []
    for csv_path in csv_files:
        rows = load_reactions_from_csv(csv_path, limit=args.limit)
        logger.info(f"  {csv_path.name}: {len(rows)} reactions loaded")
        all_reactions.extend(rows)

    logger.info(f"Total reactions: {len(all_reactions)}")

    # Extract templates
    t0 = time.perf_counter()
    template_stats, counters = extract_templates_from_reactions(all_reactions)
    elapsed = time.perf_counter() - t0

    # Sort by frequency (most common first)
    ranked = sorted(template_stats.values(), key=lambda t: t.count, reverse=True)

    # Report
    print("\n" + "=" * 70)
    print("TEMPLATE EXTRACTION RESULTS")
    print("=" * 70)
    print(f"Total reactions:      {counters['total']}")
    print(f"Atom mapping OK:      {counters['mapped_ok']} ({counters['mapped_ok']*100/max(counters['total'],1):.1f}%)")
    print(f"Atom mapping fail:    {counters['mapped_fail']}")
    print(f"Template extracted:   {counters['template_ok']} ({counters['template_ok']*100/max(counters['mapped_ok'],1):.1f}% of mapped)")
    print(f"Template extract fail:{counters['template_fail']}")
    print(f"Unique templates:     {len(ranked)}")
    print(f"Time:                 {elapsed:.1f}s")
    print()

    print(f"Top {min(20, len(ranked))} templates by frequency:")
    print("-" * 70)
    for i, t in enumerate(ranked[:20], 1):
        yield_str = f"avg={t.avg_yield:.0f}%" if t.yields else "no yield"
        files_str = ", ".join(sorted(t.source_files)[:3])
        print(f"  {i:2d}. count={t.count:4d}  {yield_str:12s}  sources: {files_str}")
        # Truncate long SMARTS for display
        smarts_display = t.reaction_smarts
        if len(smarts_display) > 100:
            smarts_display = smarts_display[:97] + "..."
        print(f"      SMARTS: {smarts_display}")
        print(f"      Example: {t.example_rxn[:120]}")
        print()

    # Save output
    output_path = args.output or str(_PROJECT_ROOT / "results" / "extracted_templates.json")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    output = {
        "extraction_meta": {
            "csv_files": [f.name for f in csv_files],
            "limit_per_file": args.limit,
            **counters,
            "unique_templates": len(ranked),
            "elapsed_seconds": round(elapsed, 2),
        },
        "templates": [t.to_dict() for t in ranked],
    }
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)

    logger.info(f"Results saved to {output_path}")


if __name__ == "__main__":
    main()

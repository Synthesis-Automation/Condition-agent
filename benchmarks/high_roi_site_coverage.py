"""Measure newly added reactive-site observations on a bounded mixed corpus."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from reactive_taxonomy import featurize_molecule  # noqa: E402


REPORT_SCHEMA_VERSION = "1.0"
DEFAULT_SOURCE = PROJECT_ROOT / "examples" / "all_dataset_300"
_CATEGORY_ORDER = (
    "explicit_c_n_o_anion",
    "polarized_c_n",
    "strained_heterocycle",
    "silyl_ether_link",
    "li_cu_al_transfer",
)


def _reactant_side(reaction_smiles: str) -> str:
    if ">>" in reaction_smiles:
        return reaction_smiles.split(">>", 1)[0].strip()
    parts = reaction_smiles.split(">")
    return parts[0].strip() if len(parts) == 3 else ""


def _site_category(site: Any) -> str | None:
    if (
        site.site_type == "nucleophile_anion"
        and site.details.get("center_token") in {"C", "N", "O"}
    ):
        return "explicit_c_n_o_anion"
    if site.canonical_signature == "PI|PolarizedC=N":
        return "polarized_c_n"
    if site.details.get("center_family") == "StrainedRing":
        return "strained_heterocycle"
    if site.canonical_signature == "LG|O|SiR3":
        return "silyl_ether_link"
    if (
        site.site_type == "transfer_group"
        and site.details.get("handle_token") in {"Li", "Cu", "Al"}
    ):
        return "li_cu_al_transfer"
    return None


def evaluate_high_roi_site_coverage(
    source: str | Path = DEFAULT_SOURCE,
    *,
    max_rows_per_file: int = 50,
) -> dict[str, Any]:
    """Return row-level incidence without using source reaction labels."""
    source_path = Path(source)
    files = (
        tuple(sorted(source_path.glob("*.csv")))
        if source_path.is_dir()
        else (source_path,)
    )
    started = time.perf_counter()
    row_counts: Counter[str] = Counter()
    site_counts: Counter[str] = Counter()
    datasets_by_category: dict[str, set[str]] = {
        category: set() for category in _CATEGORY_ORDER
    }
    examples: dict[str, list[dict[str, str]]] = {
        category: [] for category in _CATEGORY_ORDER
    }
    analyzed_rows = 0
    invalid_rows = 0
    rows_with_any_new_site = 0
    for path in files:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            for source_index, row in enumerate(csv.DictReader(handle)):
                if source_index >= max_rows_per_file:
                    break
                reaction_smiles = str(row.get("reaction_smiles") or "").strip()
                reactants = _reactant_side(reaction_smiles)
                if not reactants:
                    invalid_rows += 1
                    continue
                analysis = featurize_molecule(reactants)
                if not analysis.valid:
                    invalid_rows += 1
                    continue
                analyzed_rows += 1
                categories = Counter(
                    category
                    for site in analysis.sites
                    for category in (_site_category(site),)
                    if category is not None
                )
                if categories:
                    rows_with_any_new_site += 1
                for category, count in categories.items():
                    row_counts[category] += 1
                    site_counts[category] += count
                    datasets_by_category[category].add(path.stem)
                    if len(examples[category]) < 5:
                        examples[category].append(
                            {
                                "dataset": path.stem,
                                "reaction_id": str(
                                    row.get("reaction_id") or ""
                                ),
                                "reactant_smiles": reactants,
                            }
                        )
    return {
        "schema_version": REPORT_SCHEMA_VERSION,
        "artifact_type": "high_roi_reactive_site_coverage",
        "source": str(source_path),
        "source_label_used_for_detection": False,
        "max_rows_per_file": max_rows_per_file,
        "file_count": len(files),
        "analyzed_row_count": analyzed_rows,
        "invalid_or_missing_reaction_count": invalid_rows,
        "rows_with_any_new_site": rows_with_any_new_site,
        "row_counts_by_category": {
            category: row_counts[category]
            for category in _CATEGORY_ORDER
        },
        "site_counts_by_category": {
            category: site_counts[category]
            for category in _CATEGORY_ORDER
        },
        "dataset_counts_by_category": {
            category: len(datasets_by_category[category])
            for category in _CATEGORY_ORDER
        },
        "examples": examples,
        "elapsed_seconds": round(time.perf_counter() - started, 3),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--max-rows-per-file", type=int, default=50)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.max_rows_per_file < 1:
        raise SystemExit("--max-rows-per-file must be positive")
    report = evaluate_high_roi_site_coverage(
        args.source,
        max_rows_per_file=args.max_rows_per_file,
    )
    payload = json.dumps(report, indent=2, sort_keys=True)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload + "\n", encoding="utf-8")
    print(payload)


if __name__ == "__main__":
    main()

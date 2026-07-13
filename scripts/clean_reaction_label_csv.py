"""Clean source reaction-label CSV files using reactive-taxonomy labels."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

from reactive_taxonomy import featurize_molecule, validate_taxonomy


@dataclass(frozen=True)
class ReactiveLabelRule:
    """One source-label crosswalk verified by a representative molecule."""

    target: str
    representative_smiles: str
    canonical_signature: str


FG_LABEL_RULES: dict[str, ReactiveLabelRule] = {
    "ArBr": ReactiveLabelRule("Ar-Br", "Brc1ccccc1", "LG|Ar|Br"),
    "ArCl": ReactiveLabelRule("Ar-Cl", "Clc1ccccc1", "LG|Ar|Cl"),
    "ArI": ReactiveLabelRule("Ar-I", "Ic1ccccc1", "LG|Ar|I"),
    "ArF": ReactiveLabelRule("Ar-F", "Fc1ccccc1", "LG|Ar|F"),
    "ArNH2": ReactiveLabelRule("Ar-NH2", "Nc1ccccc1", "XH|N|H2|Ar"),
    "ArNHR": ReactiveLabelRule("Ar-NH-R", "CNc1ccccc1", "XH|N|H1|Ar,Alkyl"),
    "Ar2NH": ReactiveLabelRule("Ar1Ar2-NH", "c1ccccc1Nc2ccccc2", "XH|N|H1|Ar,Ar"),
    "ArB(OR)2": ReactiveLabelRule("Ar-B(OR)2", "COB(OC)c1ccccc1", "TM|Ar|B(OR)2"),
    "ArB(OH)2": ReactiveLabelRule("Ar-B(OH)2", "OB(O)c1ccccc1", "TM|Ar|B(OH)2"),
    "ArBF3K": ReactiveLabelRule("Ar-BF3K", "[K+].[B-](F)(F)(F)c1ccccc1", "TM|Ar|BF3K"),
    "ArOH": ReactiveLabelRule("Ar-OH", "Oc1ccccc1", "XH|O|H1|Ar"),
    "ArOSO2R": ReactiveLabelRule("Ar-OSO2R", "CCS(=O)(=O)Oc1ccccc1", "LG|Ar|OSO2R"),
    "arom. NH": ReactiveLabelRule("AromN-H", "c1cc[nH]c1", "XH|N_aromatic|H1|HeteroAr"),
    "alkeneB(OR)2": ReactiveLabelRule("Alkenyl-B(OR)2", "C=CB(OC)OC", "TM|Alkenyl|B(OR)2"),
    "alkene-Br": ReactiveLabelRule("Alkenyl-Br", "BrC=C", "LG|Alkenyl|Br"),
    "alkene-I": ReactiveLabelRule("Alkenyl-I", "IC=C", "LG|Alkenyl|I"),
    "RNH2": ReactiveLabelRule("R-NH2", "CN", "XH|N|H2|Alkyl"),
    "R2NH": ReactiveLabelRule("R1R2-NH", "CCNCC", "XH|N|H1|Alkyl,Alkyl"),
    "RCONH2": ReactiveLabelRule("R-C(O)-NH2", "CC(N)=O", "XH|N|H2|C(O)R"),
    "RCONHR": ReactiveLabelRule("R-C(O)-NHR", "CC(=O)NC", "XH|N|H1|C(O)R,Alkyl"),
    "RCO2R": ReactiveLabelRule("R-C(O)OR", "CC(=O)OC", "EC|Acyl|Alkyl|OR|ester"),
    "RSO2Cl": ReactiveLabelRule("R-S(O)2Cl", "CS(=O)(=O)Cl", "EC|Sulfonyl|Alkyl|Cl|activated"),
    "RSO2F": ReactiveLabelRule("R-S(O)2F", "CS(=O)(=O)F", "EC|Sulfonyl|Alkyl|F|activated"),
    "RSO2NHR": ReactiveLabelRule("R-S(O)2-NHR", "CS(=O)(=O)NC", "XH|N|H1|SO2R,Alkyl"),
    "RSH": ReactiveLabelRule("R-SH", "CCS", "XH|S|H1|Alkyl"),
    "Alkyl-Br": ReactiveLabelRule("R-Br", "CCCBr", "LG|Alkyl|Br"),
    "Alkyl-Cl": ReactiveLabelRule("R-Cl", "CCCCl", "LG|Alkyl|Cl"),
    "Alkyl-I": ReactiveLabelRule("R-I", "CCCI", "LG|Alkyl|I"),
    "Alkyl-OSO2R": ReactiveLabelRule("R-OSO2R", "CCOS(=O)(=O)CC", "LG|Alkyl|OSO2R"),
    "Alkyl-B(OH)2": ReactiveLabelRule("R-B(OH)2", "CCB(O)O", "TM|Alkyl|B(OH)2"),
    "Alkyl-B(OR)2": ReactiveLabelRule("R-B(OR)2", "CCB(OC)OC", "TM|Alkyl|B(OR)2"),
    "Alkyl-BF3K": ReactiveLabelRule("R-BF3K", "[K+].CC[B-](F)(F)F", "TM|Alkyl|BF3K"),
}

FG_LABEL_MAP: dict[str, str] = {
    source: rule.target for source, rule in FG_LABEL_RULES.items()
}

def validate_label_rules(label_rules: Mapping[str, ReactiveLabelRule]) -> None:
    """Verify crosswalk targets against actual reactive-taxonomy output."""
    taxonomy_errors = validate_taxonomy()
    if taxonomy_errors:
        raise ValueError(f"Reactive taxonomy is invalid: {taxonomy_errors}")
    failures: list[str] = []
    for source, rule in label_rules.items():
        analysis = featurize_molecule(rule.representative_smiles, label_style="hte_legacy")
        matching_sites = [
            site
            for site in analysis.sites
            if site.canonical_signature == rule.canonical_signature
        ]
        if not analysis.valid or len(matching_sites) != 1:
            failures.append(f"{source}: expected one {rule.canonical_signature}")
        elif matching_sites[0].chemist_label != rule.target:
            failures.append(
                f"{source}: expected {rule.target!r}, got {matching_sites[0].chemist_label!r}"
            )
    if failures:
        raise ValueError("Invalid reactive-label mappings: " + "; ".join(failures))


def clean_rows(
    rows: Iterable[dict[str, str]],
    *,
    label_map: Mapping[str, str] = FG_LABEL_MAP,
) -> tuple[list[dict[str, str]], Counter[str]]:
    """Filter requested rows and replace only explicitly mapped FG labels."""
    cleaned: list[dict[str, str]] = []
    stats: Counter[str] = Counter()

    for row in rows:
        stats["input_rows"] += 1
        fg_a = (row.get("FG A") or "").strip()
        fg_b = (row.get("FG B") or "").strip()

        both_blank = not fg_a and not fg_b
        identical = fg_a == fg_b
        protecting_group = "Protecting Group" in {fg_a, fg_b}

        if both_blank:
            stats["matched_both_blank"] += 1
        if identical:
            stats["matched_identical"] += 1
        if protecting_group:
            stats["matched_protecting_group"] += 1
        if both_blank or identical or protecting_group:
            stats["removed_union"] += 1
            continue

        for column in ("FG A", "FG B"):
            source = (row.get(column) or "").strip()
            replacement = label_map.get(source)
            if replacement:
                row[column] = replacement
                stats[f"replaced:{source}->{replacement}"] += 1

        cleaned.append(row)
        stats["output_rows"] += 1

    return cleaned, stats


def clean_csv(source: Path, destination: Path) -> Counter[str]:
    """Clean ``source`` and write a new CSV to ``destination``."""
    validate_label_rules(FG_LABEL_RULES)
    with source.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"CSV has no header: {source}")
        rows, stats = clean_rows(dict(row) for row in reader)
        fieldnames = list(reader.fieldnames)

    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()

    stats = clean_csv(args.source, args.destination)
    for key in sorted(stats):
        print(f"{key}: {stats[key]}")


if __name__ == "__main__":
    main()

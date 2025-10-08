#!/usr/bin/env python
"""Comprehensive reaction type detector run using sample reactions.

This utility loads the curated sample reaction SMILES from
``tests/sample_reactions.py`` and evaluates the chemtools routing detector
(`chemtools.router.detect_family_from_reaction`) across the full data set.

For CN coupling examples we additionally simulate two catalyst scenarios:
Pd-containing agents (expected to map to ``Buchwald_CN``) and Cu-containing
agents (expected to map to ``Ullmann_CN``). The script reports the detector
outputs for the baseline reaction as well as the catalyst-augmented variants
and highlights any mismatches against the expected families.

Usage::

    python scripts/reaction_type_detector_report.py

The report is printed to stdout as a human-readable summary.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import importlib.util
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from chemtools.router import detect_family_from_reaction
from chemtools.smiles import normalize_reaction


def _load_sample_reactions() -> Sequence[str]:
    sample_path = ROOT / "tests" / "sample_reactions.py"
    spec = importlib.util.spec_from_file_location("sample_reactions", sample_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Unable to load sample reactions from {sample_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[call-arg]
    try:
        reactions = getattr(module, "SAMPLE_REACTIONS")
    except AttributeError as exc:
        raise ImportError("SAMPLE_REACTIONS not defined in sample_reactions.py") from exc
    return reactions


SAMPLE_REACTIONS = tuple(_load_sample_reactions())


@dataclass
class DetectorRun:
    """Container for a single detector execution."""

    catalysts: Sequence[str]
    family: str
    confidence: float
    status: str
    rxn_insight_used: bool
    rxn_insight_family: Optional[str] = None
    rxn_insight_class: Optional[str] = None
    rxn_insight_name: Optional[str] = None

    def to_dict(self) -> dict:
        data = asdict(self)
        data["catalysts"] = list(self.catalysts)
        return data


@dataclass
class ReactionReport:
    """Aggregated detector outputs for one sample reaction."""

    index: int
    raw_entry: str
    smiles: str
    label: str
    runs: List[DetectorRun]

    def to_dict(self) -> dict:
        return {
            "index": self.index,
            "entry": self.raw_entry,
            "smiles": self.smiles,
            "label": self.label,
            "runs": [run.to_dict() for run in self.runs],
        }


def _split_entry(entry: str) -> Tuple[str, str]:
    """Separate SMILES string from trailing label if present."""
    entry = (entry or "").strip()
    if not entry:
        return "", ""
    # Labels in the dataset are appended as " ... (Label text)" at the end.
    # Use rfind on " (" to avoid interfering with parentheses inside SMILES.
    idx = entry.rfind(" (")
    if idx == -1 or not entry.endswith(")"):
        return entry, ""
    smiles = entry[:idx].strip()
    label = entry[idx + 2 : -1].strip()
    return smiles, label


def _inject_catalysts(norm_rsmi: str, catalysts: Sequence[str]) -> str:
    parts = norm_rsmi.split(">")
    if len(parts) != 3:
        raise ValueError(f"Unexpected reaction SMILES format: {norm_rsmi}")
    reactants, agents, products = parts
    cat_block = ".".join(catalysts)
    combined_agents = ".".join(filter(None, [agents, cat_block]))
    return f"{reactants}>{combined_agents}>{products}"


def _run_detector(rsmi: str) -> DetectorRun:
    result = detect_family_from_reaction(rsmi)
    family = str(result.get("family") or "Unknown")
    confidence = float(result.get("confidence") or 0.0)
    status = str(result.get("status") or "rule_only")
    catalysts = tuple(sorted(result.get("catalysts", {}).get("metals", [])))
    auto = result.get("auto") or {}
    rxn = result.get("rxn") or {}
    auto_available = bool(auto.get("available"))
    rxn_family = rxn.get("family")
    if not rxn_family and isinstance(auto, dict):
        rxn_family = auto.get("mapped_family")
    rxn_class = auto.get("rxn_class") if isinstance(auto, dict) else None
    rxn_name = auto.get("rxn_name") if isinstance(auto, dict) else None
    return DetectorRun(
        catalysts=catalysts,
        family=family,
        confidence=confidence,
        status=status,
        rxn_insight_used=auto_available,
        rxn_insight_family=str(rxn_family) if rxn_family else None,
        rxn_insight_class=str(rxn_class) if rxn_class else None,
        rxn_insight_name=str(rxn_name) if rxn_name else None,
    )


def _is_cn_label(label: str) -> bool:
    lower = (label or "").lower()
    return any(token in lower for token in ("c-n", "buchwald", "b-h", "ullmann"))


def generate_reports() -> Tuple[List[ReactionReport], List[str]]:
    reports: List[ReactionReport] = []
    mismatch_notes: List[str] = []
    for idx, entry in enumerate(SAMPLE_REACTIONS):
        smiles_raw, label = _split_entry(entry)
        if not smiles_raw or smiles_raw.lower().startswith("select "):
            continue
        normalized = normalize_reaction(smiles_raw)
        norm_rsmi = normalized.get("normalized") or smiles_raw

        runs: List[DetectorRun] = []

        # Baseline run (without forced catalysts)
        runs.append(_run_detector(norm_rsmi))

        if _is_cn_label(label):
            # Pd scenario (expect Buchwald)
            pd_run = _run_detector(_inject_catalysts(norm_rsmi, ("[Pd]",)))
            runs.append(pd_run)
            if pd_run.family != "Buchwald_CN":
                mismatch_notes.append(
                    f"[Pd] expected Buchwald_CN but got {pd_run.family} (index {idx}, label='{label}')"
                )

            # Cu scenario (expect Ullmann)
            cu_run = _run_detector(_inject_catalysts(norm_rsmi, ("[Cu]",)))
            runs.append(cu_run)
            if cu_run.family != "Ullmann_CN":
                mismatch_notes.append(
                    f"[Cu] expected Ullmann_CN but got {cu_run.family} (index {idx}, label='{label}')"
                )

        reports.append(
            ReactionReport(index=idx, raw_entry=entry, smiles=smiles_raw, label=label, runs=runs)
        )

    return reports, mismatch_notes


def _format_report(reports: Sequence[ReactionReport]) -> str:
    lines: List[str] = []
    total = len(reports)
    lines.append(f"Reaction Type Detector Comprehensive Report ({total} reactions)")
    lines.append("=" * 72)

    cn_total = sum(1 for report in reports if _is_cn_label(report.label))
    lines.append(f"C-N coupling cases: {cn_total}")
    lines.append("")

    for report in reports:
        lines.append(f"[{report.index:03d}] {report.label or 'Unlabeled reaction'}")
        lines.append(f"    SMILES: {report.smiles}")
        for run in report.runs:
            cats = ", ".join(run.catalysts) if run.catalysts else "none"
            rxn_flag = "yes" if run.rxn_insight_used else "no"
            lines.append(
                f"    - catalysts: {cats:12s} | family: {run.family:12s} | confidence: {run.confidence:5.2f} | status: {run.status:11s} | rxn_insight: {rxn_flag}"
            )
            if run.rxn_insight_used and (run.rxn_insight_family or run.rxn_insight_class or run.rxn_insight_name):
                extra = []
                if run.rxn_insight_family:
                    extra.append(f"family={run.rxn_insight_family}")
                if run.rxn_insight_class:
                    extra.append(f"class={run.rxn_insight_class}")
                if run.rxn_insight_name:
                    extra.append(f"name={run.rxn_insight_name}")
                lines.append(f"      rxn_insight details: {', '.join(extra)}")
        lines.append("")

    return "\n".join(lines)


def _rows_for_csv(reports: Sequence[ReactionReport]) -> List[dict]:
    rows: List[dict] = []
    for report in reports:
        for run in report.runs:
            rows.append(
                {
                    "index": report.index,
                    "label": report.label,
                    "smiles": report.smiles,
                    "raw_entry": report.raw_entry,
                    "catalysts": ";".join(run.catalysts),
                    "family": run.family,
                    "confidence": run.confidence,
                    "status": run.status,
                    "rxn_insight_used": run.rxn_insight_used,
                    "rxn_insight_family": run.rxn_insight_family or "",
                    "rxn_insight_class": run.rxn_insight_class or "",
                    "rxn_insight_name": run.rxn_insight_name or "",
                }
            )
    return rows


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--json",
        action="store_true",
        help="Also emit the full results as JSON to stdout after the human-readable report",
    )
    parser.add_argument(
        "--out-json",
        type=Path,
        help="Write the structured results to the given JSON file",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        help="Write a flattened table of results (one row per detector run) to CSV",
    )
    args = parser.parse_args(argv)

    reports, mismatches = generate_reports()
    print(_format_report(reports))
    if mismatches:
        print("\nMismatches detected:")
        for note in mismatches:
            print(f"  - {note}")

    json_payload: Optional[List[dict]] = None
    if args.json or args.out_json:
        json_payload = [report.to_dict() for report in reports]

    if args.json and json_payload is not None:
        print("\nJSON summary:")
        print(json.dumps(json_payload, indent=2))

    if args.out_json and json_payload is not None:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(json_payload, indent=2) + "\n", encoding="utf-8")
        print(f"\nWrote JSON results to {args.out_json}")

    if args.out_csv:
        csv_rows = _rows_for_csv(reports)
        fieldnames = [
            "index",
            "label",
            "smiles",
            "raw_entry",
            "catalysts",
            "family",
            "confidence",
            "status",
            "rxn_insight_used",
            "rxn_insight_family",
            "rxn_insight_class",
            "rxn_insight_name",
        ]
        args.out_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.out_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(csv_rows)
        print(f"\nWrote CSV results to {args.out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

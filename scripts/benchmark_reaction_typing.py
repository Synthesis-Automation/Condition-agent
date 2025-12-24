"""
Benchmark reaction typing vs taxonomy-v2 detection on sample reactions.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from chemtools import detect_reaction
from chemtools.reaction_typing import classify_reaction_smiles
from chemtools.reaction_typing.reaction_parser import strip_label
from examples.sample_reactions import SAMPLE_REACTIONS


LABEL_RE = re.compile(r"\(([^()]*)\)\s*$")

FEATURIZER_FAMILY_MAP: Dict[str, str] = {
    "suzuki_miyaura": "Suzuki",
    "sonogashira": "Sonogashira",
    "heck": "Heck",
    "c_n_cross_coupling": "C-N",
    "chan_lam_cn": "C-N",
    "snar_cn": "SNAr",
    "c_o_coupling": "C-O",
    "amide_coupling": "Amide",
    "acylation_acyl_halide_amide": "Amide",
    "sn2_substitution": "SN2",
}

TYPING_FAMILY_MAP: Dict[str, str] = {
    "Suzuki": "Suzuki",
    "Sonogashira": "Sonogashira",
    "Heck": "Heck",
    "C-N coupling": "C-N",
    "C-O coupling": "C-O",
    "C-S coupling": "C-S",
    "Amide formation": "Amide",
    "SN2 substitution": "SN2",
    "Negishi": "Negishi",
    "Stille": "Stille",
    "Kumada": "Kumada",
}

DEFAULT_INCLUDE = ["Suzuki", "Sonogashira", "Heck", "C-N", "C-O", "Amide", "SN2"]


@dataclass(frozen=True)
class Record:
    smiles: str
    label: str
    expected: str
    featurizer: str
    featurizer_raw: str
    featurizer_error: Optional[str]
    typing: str
    typing_raw: str
    typing_error: Optional[str]


def extract_label(line: str) -> Optional[str]:
    match = LABEL_RE.search(line or "")
    if not match:
        return None
    return match.group(1).strip()


def expected_family(label: str) -> Optional[str]:
    base = label.split(" - ", 1)[0].strip()
    lower = base.lower()

    if lower.startswith("suzuki"):
        return "Suzuki"
    if lower.startswith("sonogashira"):
        return "Sonogashira"
    if lower.startswith("heck"):
        return "Heck"
    if lower.startswith(("buchwald-hartwig", "b-h", "c-n", "ullmann c-n", "chan-lam")):
        return "C-N"
    if lower.startswith(("ullmann ether", "c-o", "mitsunobu")):
        return "C-O"
    if lower.startswith("c-s"):
        return "C-S"
    if lower.startswith("sn2"):
        return "SN2"
    if lower.startswith("snar"):
        return "SNAr"
    if lower.startswith("amide"):
        return "Amide"
    if lower.startswith("stille"):
        return "Stille"
    if lower.startswith("negishi"):
        return "Negishi"
    if lower.startswith("kumada"):
        return "Kumada"
    return None


def normalize_featurizer_family(family: str) -> str:
    return FEATURIZER_FAMILY_MAP.get(family, family or "Unknown")


def normalize_typing_family(predicted: str) -> str:
    return TYPING_FAMILY_MAP.get(predicted, predicted or "Unknown")


def iter_samples() -> Iterable[Tuple[str, str, str]]:
    for entry in SAMPLE_REACTIONS:
        if not entry or ">>" not in entry:
            continue
        label = extract_label(entry)
        if not label:
            continue
        expected = expected_family(label)
        if not expected:
            continue
        yield entry, label, expected


def run_detection(entry: str) -> Tuple[str, Optional[str], str]:
    cleaned = strip_label(entry)
    result = detect_reaction(cleaned, use_ml=False)
    family = result.get("family") or "Unknown"
    error = (result.get("details") or {}).get("error")
    return normalize_featurizer_family(family), error, family


def run_typing(entry: str) -> Tuple[str, Optional[str], str]:
    result = classify_reaction_smiles(entry)
    predicted = result.get("predicted") or "Unknown"
    error = None if result.get("ok") else result.get("error")
    return normalize_typing_family(predicted), error, predicted


def build_records(include: List[str], limit: Optional[int]) -> Tuple[List[Record], Counter]:
    include_set = {item.strip() for item in include if item.strip()}
    skipped = Counter()
    records: List[Record] = []

    for entry, label, expected in iter_samples():
        if expected not in include_set:
            skipped[expected] += 1
            continue
        featurizer, fe_error, fe_raw = run_detection(entry)
        typing, ty_error, ty_raw = run_typing(entry)
        records.append(
            Record(
                smiles=strip_label(entry),
                label=label,
                expected=expected,
                featurizer=featurizer,
                featurizer_raw=fe_raw,
                featurizer_error=fe_error,
                typing=typing,
                typing_raw=ty_raw,
                typing_error=ty_error,
            )
        )
        if limit and len(records) >= limit:
            break

    return records, skipped


def accuracy(records: List[Record], key: str) -> Tuple[int, int]:
    total = len(records)
    correct = sum(1 for rec in records if getattr(rec, key) == rec.expected)
    return correct, total


def report(records: List[Record], skipped: Counter, show_mismatches: int) -> None:
    if not records:
        print("No records matched the include set.")
        if skipped:
            print("Skipped families:", dict(skipped))
        return

    fe_correct, total = accuracy(records, "featurizer")
    ty_correct, _ = accuracy(records, "typing")

    print("Benchmark Summary")
    print(f"Total evaluated: {total}")
    print(f"Featurizer accuracy: {fe_correct}/{total} ({fe_correct / total:.1%})")
    print(f"Reaction typing accuracy: {ty_correct}/{total} ({ty_correct / total:.1%})")

    errors = Counter()
    for rec in records:
        if rec.featurizer_error:
            errors["featurizer_error"] += 1
        if rec.typing_error:
            errors["typing_error"] += 1
    if errors:
        print("Errors:", dict(errors))

    print("\nPer-family accuracy")
    by_family: Dict[str, List[Record]] = defaultdict(list)
    for rec in records:
        by_family[rec.expected].append(rec)
    for family in sorted(by_family):
        subset = by_family[family]
        fe_ok = sum(1 for rec in subset if rec.featurizer == rec.expected)
        ty_ok = sum(1 for rec in subset if rec.typing == rec.expected)
        total_family = len(subset)
        print(
            f"- {family}: featurizer {fe_ok}/{total_family} "
            f"({fe_ok / total_family:.1%}), typing {ty_ok}/{total_family} "
            f"({ty_ok / total_family:.1%})"
        )

    if skipped:
        print("\nSkipped families:", dict(skipped))

    if show_mismatches <= 0:
        return

    print("\nSample mismatches (featurizer)")
    fe_mismatches = [rec for rec in records if rec.featurizer != rec.expected]
    for rec in fe_mismatches[:show_mismatches]:
        print(
            f"- expected={rec.expected}, predicted={rec.featurizer} "
            f"(raw={rec.featurizer_raw}) :: {rec.label}"
        )

    print("\nSample mismatches (reaction_typing)")
    ty_mismatches = [rec for rec in records if rec.typing != rec.expected]
    for rec in ty_mismatches[:show_mismatches]:
        print(
            f"- expected={rec.expected}, predicted={rec.typing} "
            f"(raw={rec.typing_raw}) :: {rec.label}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark reaction typing vs taxonomy detection")
    parser.add_argument(
        "--include",
        default=",".join(DEFAULT_INCLUDE),
        help="Comma-separated expected families to include",
    )
    parser.add_argument("--limit", type=int, default=None, help="Limit evaluated reactions")
    parser.add_argument(
        "--show-mismatches",
        type=int,
        default=8,
        help="How many mismatches to display per method",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    include = [item.strip() for item in args.include.split(",") if item.strip()]
    records, skipped = build_records(include=include, limit=args.limit)
    report(records, skipped, show_mismatches=args.show_mismatches)


if __name__ == "__main__":
    main()

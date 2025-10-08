from __future__ import annotations

from argparse import ArgumentParser
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from chemtools import recommend
from chemtools.selector_payloads import build_amide_selector_payload
from chemtools.scdb_matcher.loader import load_db
from chemtools.scdb_matcher.matcher import match
from chemtools.scdb_matcher.types import SelectorRuleDB


def _parse_samples(path: Path) -> List[Tuple[str, str]]:
    reactions: List[Tuple[str, str]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not (line.startswith("\"") and ">>" in line):
            continue
        entry = line.strip('",')
        label = ""
        if " (Amide:" in entry:
            entry, label = entry.split(" (Amide:", 1)
            label = label.rstrip(")").strip()
        reactions.append((entry, label))
    return reactions


def _extract_reactants(reaction: str) -> Tuple[str, str]:
    reactant_side, _ = reaction.split(">>", 1)
    parts = reactant_side.split(".")
    acid = parts[0].strip()
    amine = parts[1].strip() if len(parts) > 1 else ""
    return acid, amine


def _run(db_path: Path, sample_path: Path, show_examples: int) -> Dict[str, Iterable[Dict[str, object]]]:
    db = load_db(db_path)
    if not isinstance(db, SelectorRuleDB):
        raise RuntimeError(f"Database at {db_path} is not a selector rule DB")

    reactions = _parse_samples(sample_path)
    counts: Counter[str] = Counter()
    grouped = defaultdict(list)

    for rxn, label in reactions:
        acid, amine = _extract_reactants(rxn)
        feats = recommend.feat_molecular.featurize(acid, amine)
        rule_features = recommend._compose_rule_features(
            "Amide_Formation",
            feats,
            role_pack=None,
            reactants=[acid, amine],
            detection={"rule_family": "Amide_Formation"},
        )
        payload = build_amide_selector_payload(rule_features)
        result = match(db, features=payload)
        counts[result.entry_id] += 1
        grouped[result.entry_id].append(
            {
                "reaction": rxn,
                "label": label,
                "payload": payload,
                "match_type": result.match_type,
            }
        )

    print(f"Processed {len(reactions)} reactions from {sample_path}")
    for entry_id, count in counts.most_common():
        print(f"- {entry_id}: {count}")
        for row in grouped[entry_id][:show_examples]:
            label = f" ({row['label']})" if row['label'] else ""
            print(
                f"    Â· {row['reaction']}{label} â†?class={row['payload']['acid']['class']} |"
                f" amine={row['payload']['amine']['class']} ({row['payload']['amine']['nucleophilicity']})"
            )
    if "DEFAULT" in grouped:
        print("\nDefaulted reactions:")
        for row in grouped["DEFAULT"]:
            print(f"- {row['reaction']} :: payload={row['payload']}")

    print()

    return grouped


def main() -> None:
    parser = ArgumentParser(description="Summarize amide selector matches for sample reactions.")
    parser.add_argument(
        "--db",
        type=Path,
        default=Path("data/conditionDB/amide_formation_db.json"),
        help="Path to the amide selector database JSON",
    )
    parser.add_argument(
        "--samples",
        type=Path,
        default=Path("tests/amide_formation.md"),
        help="Sample reactions to evaluate",
    )
    parser.add_argument(
        "--examples",
        type=int,
        default=3,
        help="Number of example reactions to print per selector",
    )
    args = parser.parse_args()

    _run(args.db, args.samples, max(0, args.examples))


if __name__ == "__main__":
    main()

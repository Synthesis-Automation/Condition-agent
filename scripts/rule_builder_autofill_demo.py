"""
CLI helper to draft rule databases from protocol text using the LLM-assisted
autofill workflow.

Example:
    python scripts/rule_builder_autofill_demo.py ^
        --family Suzuki_Miyaura ^
        --metadata-id suzuki_auto_v3 ^
        --metadata-name "Suzuki-Miyaura Draft v3" ^
        --metadata-version v3-draft ^
        --reference "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1" ^
        --reference-file data/rule_db_v2/Suzuki_db.json ^
        --protocol-file docs/AUTOMATION_FORMAT_SUMMARY.md ^
        --applies-if-hint sp2_halide_present ^
        --applies-if-hint sp2_boron_present ^
        --output build/suzuki_auto_v3.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable, List

from chem_assistant.chemtools_wrapper import (
    RuleBuilderAutoInput,
    run_rule_builder_autofill,
)


def _read_reference_file(path: Path) -> List[str]:
    if not path.exists():
        raise FileNotFoundError(f"Reference file not found: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    refs = data.get("reaction", {}).get("reference_reactions", [])
    if not isinstance(refs, list):
        raise ValueError(f"reference_reactions missing or malformed in {path}")
    return [str(r).strip() for r in refs if str(r).strip()]


def _collect_references(args: argparse.Namespace) -> List[str]:
    refs: List[str] = []
    if args.reference:
        refs.extend(args.reference)
    if args.reference_file:
        refs.extend(_read_reference_file(Path(args.reference_file)))
    deduped = []
    for ref in refs:
        if ref and ref not in deduped:
            deduped.append(ref)
    return deduped


def _load_protocol_text(args: argparse.Namespace) -> str:
    if args.protocol_text:
        return args.protocol_text
    if args.protocol_file:
        path = Path(args.protocol_file)
        if not path.exists():
            raise FileNotFoundError(f"Protocol file not found: {path}")
        return path.read_text(encoding="utf-8")
    raise SystemExit("Please supply --protocol-text or --protocol-file.")


def _print_issues(issues: Iterable[dict]) -> None:
    issues = list(issues)
    if not issues:
        print("No validation issues detected.")
        return
    print("Validation findings:")
    for issue in issues:
        field = issue.get("field", "unknown")
        severity = issue.get("severity", "warning").upper()
        message = issue.get("message", "")
        print(f"  [{severity}] {field}: {message}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="LLM-assisted rule-builder autofill helper."
    )
    parser.add_argument("--family", required=True, help="Reaction family name.")
    parser.add_argument("--metadata-id", required=True, help="Metadata id slug.")
    parser.add_argument("--metadata-name", required=True, help="Metadata name.")
    parser.add_argument("--metadata-version", required=True, help="Metadata version tag.")
    parser.add_argument("--created-date", help="Creation date (YYYY-MM-DD). Defaults to today.")
    parser.add_argument("--status", default="draft", help="Metadata status label.")
    parser.add_argument("--tag", action="append", dest="tags", help="Metadata tag (repeatable).")
    parser.add_argument(
        "--reference",
        action="append",
        help="Reference reaction SMILES (repeatable, format reactants>>products).",
    )
    parser.add_argument(
        "--reference-file",
        help="Path to existing rule DB JSON from which to copy reference_reactions.",
    )
    parser.add_argument(
        "--protocol-text",
        help="Protocol or HTE description text. Exclusive with --protocol-file.",
    )
    parser.add_argument(
        "--protocol-file",
        help="File containing protocol text (Markdown or plain).",
    )
    parser.add_argument("--notes", help="Override for reaction notes.")
    parser.add_argument(
        "--desired-focus",
        help="Optional prompt hint describing special focus (e.g., 'stress aryl chlorides').",
    )
    parser.add_argument(
        "--applies-if-hint",
        action="append",
        dest="applies_if_hints",
        help="Feature token to force into applies_if.all (repeatable).",
    )
    parser.add_argument(
        "--modifier-hint",
        action="append",
        dest="modifier_hints",
        help="Symptom trigger to prioritize in modifiers (repeatable).",
    )
    parser.add_argument(
        "--max-base-rules",
        type=int,
        default=4,
        help="Maximum number of base rules to request from the LLM (default 4).",
    )
    parser.add_argument(
        "--output",
        help="Optional path to write the resulting JSON database.",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Pretty-print the resulting JSON to stdout.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    references = _collect_references(args)
    if not references:
        raise SystemExit("Provide at least one reference reaction via --reference or --reference-file.")
    protocol_text = _load_protocol_text(args)

    params = RuleBuilderAutoInput(
        family=args.family,
        metadata_id=args.metadata_id,
        metadata_name=args.metadata_name,
        metadata_version=args.metadata_version,
        created_date=args.created_date,
        status=args.status,
        tags=args.tags,
        reference_reactions=references,
        protocol_text=protocol_text,
        notes=args.notes,
        desired_focus=args.desired_focus,
        applies_if_hints=args.applies_if_hints,
        modifier_hints=args.modifier_hints,
        max_base_rules=args.max_base_rules,
    )

    result = run_rule_builder_autofill(params)
    if not result.get("success"):
        print("Autofill failed:")
        print(json.dumps(result, indent=2))
        raise SystemExit(1)

    issues = result.get("issues", [])
    _print_issues(issues)

    rule_db = result.get("rule_database")
    if args.pretty:
        print(json.dumps(rule_db, indent=2))

    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(rule_db, indent=2), encoding="utf-8")
        print(f"Wrote rule database to {output_path}")
    else:
        print("No output path supplied; skipping file write.")


if __name__ == "__main__":
    main()

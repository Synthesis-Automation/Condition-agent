"""
CLI for testing compound and reaction featurizations.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
import argparse
from typing import Any, Dict, Iterable, List, Optional

# Ensure repo root is on sys.path for direct execution.
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.featurizers.unified import featurize_molecule, featurize_reaction


def _print_json(payload: Dict[str, Any]) -> None:
    print(json.dumps(payload, indent=2, sort_keys=True))


def _format_list(items: Iterable[str]) -> str:
    return "\n".join(f"  - {item}" for item in items)


def _get_molecule_payload(payload: Dict[str, Any]) -> Dict[str, Any]:
    if "molecule" in payload and isinstance(payload.get("molecule"), dict):
        return payload["molecule"]
    return payload


def _get_reaction_payload(payload: Dict[str, Any]) -> Dict[str, Any]:
    if "reaction" in payload and isinstance(payload.get("reaction"), dict):
        return payload["reaction"]
    return payload


def _print_molecule_summary(payload: Dict[str, Any]) -> None:
    molecule = _get_molecule_payload(payload)
    print("\nSummary (Molecule)")
    print("-" * 72)
    print(f"SMILES: {molecule.get('smiles')}")

    motifs = molecule.get("motifs") or []
    if motifs:
        lines = []
        for motif in motifs:
            compound_id = motif.get("compound_id", "unknown")
            a_idx = motif.get("a_atom_idx")
            b_idx = motif.get("b_atom_idx")
            lines.append(f"{compound_id} (a={a_idx}, b={b_idx})")
        print("Motifs:")
        print(_format_list(lines))

    steric = molecule.get("steric", {})
    aryl = steric.get("aryl", []) if isinstance(steric, dict) else []
    alkyl = steric.get("alkyl", []) if isinstance(steric, dict) else []
    if aryl:
        print("Aryl Steric:")
        entries = []
        for entry in aryl:
            result = entry.get("result", {})
            entries.append(
                f"{entry.get('compound_id', 'unknown')}: score={result.get('score_0_10')}, "
                f"ortho_subs={result.get('ortho_substituent_count')}"
            )
        print(_format_list(entries))
    if alkyl:
        print("Alkyl Steric:")
        entries = []
        for entry in alkyl:
            result = entry.get("result", {})
            entries.append(
                f"{entry.get('compound_id', 'unknown')}: score={result.get('score_0_10')}"
            )
        print(_format_list(entries))

    electronics = molecule.get("electronics", {})
    aryl_elec = electronics.get("aryl", []) if isinstance(electronics, dict) else []
    if aryl_elec:
        print("Aryl Electronics:")
        entries = []
        for entry in aryl_elec:
            result = entry.get("result", {})
            include_score = result.get("including_group_score_0_10")
            suffix = f", including_group={include_score}" if include_score is not None else ""
            entries.append(
                f"{entry.get('compound_id', 'unknown')}: score={result.get('score_0_10')}{suffix}"
            )
        print(_format_list(entries))

    nearby = molecule.get("nearby", [])
    if nearby:
        print("Nearby Groups (per motif):")
        entries = []
        for entry in nearby:
            compound_id = entry.get("compound_id", "unknown")
            groups = entry.get("result") or []
            if groups:
                # groups is a list of dicts with 'name' key
                group_names = [g.get("name", "unknown") if isinstance(g, dict) else str(g) for g in groups]
                entries.append(f"{compound_id}: {', '.join(group_names)}")
            else:
                entries.append(f"{compound_id}: none")
        print(_format_list(entries))

    snar = molecule.get("snar_feasibility")
    if snar:
        print("SNAr Feasibility:")
        for item in snar:
            status = "YES" if item.get("feasible") else "NO"
            print(f"  - {item.get('motif')}: {status} (confidence: {item.get('confidence')})")
            print(f"    Reason: {item.get('reason')}")


def _print_reaction_summary(payload: Dict[str, Any]) -> None:
    reaction = _get_reaction_payload(payload)
    print("\nSummary (Reaction)")
    print("-" * 72)
    print(f"Reaction SMILES: {reaction.get('reaction_smiles')}")

    reaction_type = reaction.get("reaction_type") or {}
    if isinstance(reaction_type, dict):
        print(
            "Reaction Type: "
            f"{reaction_type.get('reaction_type') or reaction_type.get('name')}"
            f" (confidence {reaction_type.get('confidence')})"
        )

    detection = reaction.get("detection", {})
    matches = detection.get("matches") if isinstance(detection, dict) else None
    if matches:
        print(f"Detection Matches: {len(matches)}")

    reactants = reaction.get("reactants") or []
    print(f"Reactants: {len(reactants)}")

    aggregates = reaction.get("aggregates") or {}
    if aggregates:
        print("Aggregates:")
        for key, value in aggregates.items():
            print(f"  - {key}: {value}")

    feasibility = reaction.get("feasibility")
    if feasibility:
        print("Feasibility Analysis:")
        status = "YES" if feasibility.get("feasible") else "NO"
        print(f"  - Possible: {status} (confidence: {feasibility.get('confidence')})")
        print(f"  - Reason: {feasibility.get('reason')}")


def _print_readable(payload: Dict[str, Any]) -> None:
    kind = payload.get("kind")
    if kind == "molecule":
        _print_molecule_summary(payload)
    elif kind == "reaction":
        _print_reaction_summary(payload)
    else:
        print("\nSummary")
        print("-" * 72)
        print(f"Kind: {kind or 'unknown'}")


def _prompt(text: str) -> str:
    try:
        return input(text)
    except EOFError:
        return "q"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Interactive CLI for ChemTools featurization.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Formats:\n"
            "  summary     - readable summary only (default)\n"
            "  full-json   - full JSON payload only\n"
            "  both        - summary + full JSON\n"
        ),
    )
    parser.add_argument(
        "--format",
        default="summary",
        choices=["summary", "full-json", "both"],
        help="Output format (default: summary).",
    )
    parser.add_argument(
        "--include-ar-h",
        action="store_true",
        help="Include Ar-H motifs even if other motifs are present.",
    )
    parser.add_argument(
        "--target-groups",
        help="Focus on specific motifs (e.g., 'Br', 'H', 'CN'). Comma-separated.",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    
    target_groups = None
    if args.target_groups:
        target_groups = [tg.strip() for tg in args.target_groups.split(",")]

    options = {
        "include_ar_h": args.include_ar_h,
        "target_groups": target_groups,
    }
    print("ChemTools Featurization CLI")
    print("Enter 'q' to quit.")

    while True:
        mode = _prompt("Mode (compound/reaction/auto) [auto]: ").strip().lower()
        if not mode:
            mode = "auto"
        if mode in {"q", "quit", "exit"}:
            return 0
        if mode not in {"compound", "reaction", "auto", "c", "r"}:
            print("Invalid mode. Choose compound, reaction, or auto.")
            continue

        smiles = _prompt("SMILES: ").strip()
        if not smiles:
            print("Empty input. Try again.")
            continue
        if smiles.lower() in {"q", "quit", "exit"}:
            return 0

        targets = _prompt("Target Groups (optional, e.g. Br,CN): ").strip()
        current_options = dict(options)
        if targets:
            current_options["target_groups"] = [tg.strip() for tg in targets.split(",")]

        try:
            if mode in {"reaction", "r"} or (mode == "auto" and ">" in smiles):
                payload = featurize_reaction(smiles, options=current_options)
            else:
                payload = featurize_molecule(smiles, options=current_options)
        except Exception as exc:
            print(f"Error: {exc}")
            continue

        if args.format in {"summary", "both"}:
            _print_readable(payload)
        if args.format in {"full-json", "both"}:
            print("\nFull Payload")
            print("-" * 72)
            _print_json(payload)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

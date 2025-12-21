#!/usr/bin/env python3
"""
Rule-based condition recommendation CLI (rule_db_v2).

This is a lightweight tester for the feature-driven rule engine in
`chemtools.rule`, backed by the schema v2 rule packs in `data/rule_db_v2/`.

Examples:
  python app/rule_based_recommendation_cli.py --rxn "Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1"
  python app/rule_based_recommendation_cli.py --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1" --family c_n_cross_coupling
  python app/rule_based_recommendation_cli.py --rxn "..." --db Suzuki_db.json --json
  python app/rule_based_recommendation_cli.py --interactive
  python app/rule_based_recommendation_cli.py --list-dbs
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Ensure project-root imports work when running as a script.
sys.path.insert(0, str(Path(__file__).parent.parent))

from chemtools import detect_reaction
from chemtools.formatters.rule_to_protocol import rule_conditions_to_reaction_setup
from chemtools.rule import RuleEngine
from chemtools.taxonomy.rule_db import resolve_rule_db_v2
from chemtools.util.rdkit_helpers import rdkit_available

_RULE_DIR = Path("data") / "rule_db_v2"


def _list_rule_dbs(rule_dir: Path = _RULE_DIR) -> List[str]:
    if not rule_dir.exists():
        return []
    return sorted([p.name for p in rule_dir.glob("*.json") if p.is_file()])


def _parse_reactants(reaction_smiles: str) -> List[str]:
    if ">>" in reaction_smiles:
        reactants_part = reaction_smiles.split(">>", 1)[0]
    elif ">" in reaction_smiles:
        reactants_part = reaction_smiles.split(">", 1)[0]
    else:
        reactants_part = reaction_smiles
    return [frag.strip() for frag in reactants_part.split(".") if frag.strip()]


def _resolve_db_path(
    *,
    db: Optional[str],
    family: Optional[str],
    reaction_smiles: Optional[str],
    rule_dir: Path = _RULE_DIR,
) -> Tuple[Path, Dict[str, Any]]:
    detection: Dict[str, Any] = {}
    if reaction_smiles:
        detection = detect_reaction(reaction_smiles, use_ml=False) or {}

    if db:
        candidate = Path(db)
        if not candidate.is_absolute():
            by_name = rule_dir / (db if db.lower().endswith(".json") else f"{db}.json")
            if by_name.exists():
                candidate = by_name
        if not candidate.exists():
            raise FileNotFoundError(
                f"Rule DB not found: '{db}'. Try --list-dbs or pass a full path."
            )
        return candidate.resolve(), detection

    effective_family = (family or detection.get("family") or "").strip()
    if effective_family:
        stem = resolve_rule_db_v2(effective_family)
        if stem:
            candidate = rule_dir / f"{stem}.json"
            if candidate.exists():
                return candidate.resolve(), detection

    raise FileNotFoundError(
        "Could not resolve rule DB. Provide --db, or provide a resolvable --family, "
        "or pass a reaction that auto-detects to a supported family."
    )


def _build_payload(
    *,
    reaction_smiles: str,
    db_path: Path,
    detection: Dict[str, Any],
    recommendation: Any,
    include_automation: bool,
    scale_mmol: float,
) -> Dict[str, Any]:
    payload: Dict[str, Any] = {
        "input": {"reaction_smiles": reaction_smiles},
        "detection": detection,
        "database": {"path": str(db_path), "name": db_path.name},
        "recommendation": recommendation.to_dict(),
    }

    if include_automation:
        substrates = [
            {
                "name": f"Reactant {idx + 1}",
                "smiles": smiles,
                "role": "starting_material",
                "mmol": scale_mmol,
                "equivalents": 1.0,
            }
            for idx, smiles in enumerate(_parse_reactants(reaction_smiles))
        ]
        payload["automation"] = rule_conditions_to_reaction_setup(
            recommendation.base_rule.conditions,
            user_substrates=substrates or None,
            scale_mmol=scale_mmol,
            reaction_family=detection.get("family"),
        )

    return payload


def _run_once(args: argparse.Namespace, reaction_smiles: str) -> bool:
    db_path, detection = _resolve_db_path(
        db=args.db,
        family=args.family,
        reaction_smiles=reaction_smiles,
    )

    engine = RuleEngine.from_file(db_path)
    rec = engine.recommend(
        reaction_smiles,
        symptoms=args.symptoms,
        combine_method=args.combine,
        include_reasoning=args.verbose,
    )

    if args.json or args.automation:
        payload = _build_payload(
            reaction_smiles=reaction_smiles,
            db_path=db_path,
            detection=detection,
            recommendation=rec,
            include_automation=bool(args.automation),
            scale_mmol=float(args.scale_mmol),
        )
        print(json.dumps(payload, indent=2, ensure_ascii=False))
        return True

    print(f"Detected family: {detection.get('family', 'Unknown')} ({detection.get('confidence', 0.0):.2f})")
    print(f"Rule DB: {db_path}")
    print(rec.format_summary())
    return True


def _interactive(args: argparse.Namespace) -> None:
    print("=" * 80)
    print("Rule-based Condition Recommendation (rule_db_v2) - Interactive Mode")
    print("=" * 80)
    print("Enter reaction SMILES (reactants>>products).")
    print("Commands:")
    print("  /db <name_or_path>          - set database (e.g., Suzuki_db.json)")
    print("  /family <family_label>      - set family label (used to resolve db)")
    print("  /symptoms <a,b,c>           - set symptoms list")
    print("  /combine union|all|first    - set feature combine mode")
    print("  /json on|off                - toggle JSON output")
    print("  /automation on|off          - toggle automation output (JSON)")
    print("  /scale <float>              - set scale (mmol) for automation")
    print("  /verbose on|off             - toggle reasoning/features")
    print("  /list-dbs                    - list available `data/rule_db_v2/*.json`")
    print("  /show                        - show current settings")
    print("  /help                        - show this help")
    print("  /quit                        - exit")
    print("=" * 80)

    current_db = args.db
    current_family = args.family
    current_symptoms = list(args.symptoms or [])
    current_combine = args.combine
    current_json = bool(args.json)
    current_automation = bool(args.automation)
    current_verbose = bool(args.verbose)
    current_scale = float(args.scale_mmol)

    while True:
        try:
            raw = input("rxn> ").strip()
        except (KeyboardInterrupt, EOFError):
            print()
            break

        if not raw:
            continue

        lowered = raw.lower()
        if lowered in {"/quit", "/exit"}:
            break
        if lowered == "/help":
            print("Type a reaction SMILES or use /show. Use /quit to exit.")
            continue
        if lowered == "/list-dbs":
            for name in _list_rule_dbs():
                print(f"  - {name}")
            continue
        if lowered == "/show":
            print("-" * 60)
            print(f"db: {current_db or '(auto)'}")
            print(f"family: {current_family or '(auto from reaction)'}")
            print(f"symptoms: {current_symptoms or '[]'}")
            print(f"combine: {current_combine}")
            print(f"json: {current_json}")
            print(f"automation: {current_automation}")
            print(f"scale_mmol: {current_scale}")
            print(f"verbose: {current_verbose}")
            print("-" * 60)
            continue

        if raw.startswith("/db "):
            current_db = raw.split(" ", 1)[1].strip() or None
            print(f"db set to: {current_db or '(auto)'}")
            continue
        if raw.startswith("/family "):
            current_family = raw.split(" ", 1)[1].strip() or None
            print(f"family set to: {current_family or '(auto)'}")
            continue
        if raw.startswith("/symptoms "):
            tail = raw.split(" ", 1)[1].strip()
            current_symptoms = [s.strip() for s in tail.split(",") if s.strip()]
            print(f"symptoms set to: {current_symptoms}")
            continue
        if raw.startswith("/combine "):
            tail = raw.split(" ", 1)[1].strip().lower()
            if tail not in {"union", "all", "first"}:
                print("combine must be one of: union, all, first")
                continue
            current_combine = tail
            print(f"combine set to: {current_combine}")
            continue
        if raw.startswith("/json "):
            tail = raw.split(" ", 1)[1].strip().lower()
            if tail not in {"on", "off"}:
                print("usage: /json on|off")
                continue
            current_json = tail == "on"
            print(f"json: {current_json}")
            continue
        if raw.startswith("/automation "):
            tail = raw.split(" ", 1)[1].strip().lower()
            if tail not in {"on", "off"}:
                print("usage: /automation on|off")
                continue
            current_automation = tail == "on"
            print(f"automation: {current_automation}")
            continue
        if raw.startswith("/scale "):
            tail = raw.split(" ", 1)[1].strip()
            try:
                current_scale = float(tail)
            except ValueError:
                print("usage: /scale <float>")
                continue
            print(f"scale_mmol: {current_scale}")
            continue
        if raw.startswith("/verbose "):
            tail = raw.split(" ", 1)[1].strip().lower()
            if tail not in {"on", "off"}:
                print("usage: /verbose on|off")
                continue
            current_verbose = tail == "on"
            print(f"verbose: {current_verbose}")
            continue
        if raw.startswith("/"):
            print("Unknown command. Use /help.")
            continue

        if ">>" not in raw and ">" not in raw:
            print("Reaction must contain '>>' (reactants>>products).")
            continue

        namespace = argparse.Namespace(
            db=current_db,
            family=current_family,
            symptoms=current_symptoms,
            combine=current_combine,
            json=current_json,
            automation=current_automation,
            scale_mmol=current_scale,
            verbose=current_verbose,
        )

        try:
            _run_once(namespace, raw)
        except Exception as exc:
            print(f"Error: {exc}")
            if current_verbose:
                import traceback

                traceback.print_exc()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="CLI tester for chemtools.rule + data/rule_db_v2 rule packs",
    )
    parser.add_argument("--rxn", help="Reaction SMILES (reactants>>products)")
    parser.add_argument("--family", help="Override detected family label")
    parser.add_argument("--db", help="Rule DB filename in data/rule_db_v2/ or a full path")
    parser.add_argument("--symptoms", nargs="*", default=[], help="Observed symptoms")
    parser.add_argument(
        "--combine",
        choices=["union", "all", "first"],
        default="union",
        help="How to combine multi-reactant features",
    )
    parser.add_argument("--verbose", action="store_true", help="Include reasoning/features")
    parser.add_argument("--json", action="store_true", help="Output JSON payload")
    parser.add_argument(
        "--automation",
        action="store_true",
        help="Include automation format (JSON payload)",
    )
    parser.add_argument(
        "--scale-mmol",
        type=float,
        default=1.0,
        help="Scale (mmol) for automation output",
    )
    parser.add_argument("--interactive", "-i", action="store_true", help="Run REPL mode")
    parser.add_argument("--list-dbs", action="store_true", help="List available rule DBs")

    args = parser.parse_args()

    if args.list_dbs:
        dbs = _list_rule_dbs()
        if not dbs:
            print(f"No rule DBs found under: {_RULE_DIR}")
            sys.exit(1)
        print("Available rule DBs:")
        for name in dbs:
            print(f"  - {name}")
        return

    if not rdkit_available():
        print("RDKit is required for rule-based feature extraction, but is not available.")
        sys.exit(1)

    if args.interactive:
        _interactive(args)
        return

    if not args.rxn:
        parser.error("Provide --rxn, or use --interactive, or use --list-dbs.")

    _run_once(args, args.rxn)


if __name__ == "__main__":
    main()

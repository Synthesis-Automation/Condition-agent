"""Integrated command-line tester for standalone reactive-taxonomy features."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Sequence

from . import featurize_molecule, featurize_reaction, validate_taxonomy


_MOLECULE_SITE_TYPES = (
    "leaving_group",
    "pronucleophile_XH",
    "transfer_group",
    "electrophilic_center",
    "aromatic_CH",
    "unsaturated_bond",
    "dipolar_group",
    "heteroatom_bond",
)

_SELF_TEST_MOLECULES = (
    "Brc1ccc(N)cc1C#N",
    "CC(=O)Cl",
    "CS(=O)(=O)Cl",
    "OB(O)c1ccccc1",
    "CNC",
)

_SELF_TEST_REACTIONS = (
    "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
    "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
    "Brc1ccccc1.Oc1ccccc1>>c1ccc(Oc2ccccc2)cc1",
    "Brc1ccccc1.Sc1ccccc1>>c1ccc(Sc2ccccc2)cc1",
)


def _json_dump(value: Any, *, compact: bool = False) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        indent=None if compact else 2,
    )


def _joined(values: Iterable[Any]) -> str:
    return ", ".join(str(value) for value in values) or "none"


def _molecule_summary(result: Any) -> str:
    lines = [
        f"valid: {result.valid}",
        f"input: {result.input_smiles}",
        f"canonical: {result.canonical_smiles or '-'}",
        f"components: {len(result.components)}",
        f"reactive sites: {len(result.sites)}",
    ]
    for site in result.sites:
        environment = (site.context_features or {}).get("environment") or {}
        steric = (environment.get("steric") or {}).get("class", "-")
        electronic = (environment.get("electronic") or {}).get("class", "-")
        lines.append(
            f"  {site.site_id}: {site.chemist_label} [{site.site_type}; "
            f"availability={site.availability}; steric={steric}; electronic={electronic}]"
        )
        if site.details.get("anchor_subtype") in {"benzylic", "allylic", "propargylic"}:
            lines.append(f"    attachment context: {site.details['anchor_subtype']}")
        attached = (environment.get("steric") or {}).get("attached_groups") or []
        for group in attached:
            group_class = group.get("attachment_carbon_class") or group.get("context") or "Other"
            branching = ", alpha-branched" if group.get("alpha_branched") else ""
            lines.append(f"    attached group: {group.get('context', 'Other')} ({group_class}{branching})")
    lines.append(f"functional groups: {len(result.functional_groups)}")
    for group in result.functional_groups:
        lines.append(
            f"  component {group.component_index}: {group.chemist_label} "
            f"[{group.group_id}; atoms={_joined(group.atom_indices)}]"
        )
    if result.warnings:
        lines.append(f"warnings: {_joined(result.warnings)}")
    if result.error:
        lines.append(f"error: {result.error}")
    return "\n".join(lines)


def _reaction_summary(result: Any) -> str:
    lines = [
        f"valid: {result.valid}",
        f"input: {result.input_reaction_smiles}",
        f"evidence: {result.evidence_quality}",
        f"transformation: {result.transformation_class or '-'}",
        f"named family: {result.named_family or '-'}",
        f"compatible families: {_joined(result.compatible_named_families)}",
        f"reaction label: {result.reaction_label or '-'}",
        f"label status: {result.reaction_label_status}",
        f"candidates: {len(result.candidates)}",
        f"mapped bond changes: {len(result.mapped_bond_changes)}",
    ]
    if result.selected_candidate:
        lines.append(f"selected grammar: {result.selected_candidate.grammar_id}")
        for change in result.selected_candidate.predicted_bond_changes:
            lines.append(
                f"  edit: {change.change_type} {change.atom_1_role}-{change.atom_2_role} "
                f"{change.old_order or '-'}->{change.new_order or '-'}"
            )
    if result.product_connection:
        lines.append(
            f"product connection: {result.product_connection.concise_label} "
            f"[{result.product_connection.connection_type}; {result.product_connection.evidence}]"
        )
    lines.append(f"spectator groups: {len(result.spectator_groups)}")
    for group in result.spectator_groups:
        distance = "-" if group.graph_distance is None else group.graph_distance
        lines.append(
            f"  component {group.component_index}: {group.chemist_label} "
            f"[{group.group_id}; distance={distance}]"
        )
    if result.family_environment:
        lines.append(f"family environment: {result.family_environment.family_id}")
        for partner in result.family_environment.partners:
            lines.append(
                f"  {partner.role}: {partner.chemist_label} "
                f"[handle={partner.handle_token or '-'}; context={partner.anchor_context or '-'}; "
                f"flags={_joined(partner.flags)}]"
            )
    if result.warnings:
        lines.append(f"warnings: {_joined(result.warnings)}")
    if result.error:
        lines.append(f"error: {result.error}")
    return "\n".join(lines)


def _molecule_concise_summary(result: Any) -> str:
    lines = [
        f"Molecule: {result.canonical_smiles or result.input_smiles}",
        f"Status: {'valid' if result.valid else 'invalid'}",
    ]
    if result.sites:
        lines.append("Reactive sites:")
        for site in result.sites:
            lines.append(
                f"  {site.chemist_label} — {site.site_type}, {site.availability}"
            )
            environment = (site.context_features or {}).get("environment") or {}
            attached = (environment.get("steric") or {}).get("attached_groups") or []
            for group in attached:
                group_class = group.get("attachment_carbon_class") or group.get("context") or "Other"
                lines.append(f"    attached {group.get('context', 'group')}: {group_class}")
    else:
        lines.append("Reactive sites: none")
    if result.functional_groups:
        groups = sorted({group.chemist_label for group in result.functional_groups})
        lines.append(f"Functional groups: {_joined(groups)}")
    else:
        lines.append("Functional groups: none")
    if result.warnings:
        lines.append(f"Warnings: {_joined(result.warnings)}")
    if result.error:
        lines.append(f"Error: {result.error}")
    return "\n".join(lines)


def _reaction_concise_summary(result: Any) -> str:
    lines = [
        f"Reaction: {result.reaction_label or '-'}",
        f"Status: {'valid' if result.valid else 'invalid'}",
        f"Evidence: {result.evidence_quality}",
        f"Transformation: {result.transformation_class or 'unknown'}",
        f"Named family: {result.named_family or 'unknown'}",
    ]
    if result.product_connection:
        lines.append(
            f"Product connection: {result.product_connection.concise_label} "
            f"({result.product_connection.connection_type})"
        )
    else:
        lines.append("Product connection: not verified")
    if result.family_environment:
        partners = [
            f"{partner.role}={partner.chemist_label}"
            for partner in result.family_environment.partners
        ]
        lines.append(f"Reactive partners: {_joined(partners)}")
    if result.spectator_groups:
        spectators = sorted({group.chemist_label for group in result.spectator_groups})
        lines.append(f"Spectator groups: {_joined(spectators)}")
    if result.warnings:
        lines.append(f"Warnings: {_joined(result.warnings)}")
    if result.error:
        lines.append(f"Error: {result.error}")
    return "\n".join(lines)


def _print_result(result: Any, output_format: str, *, concise: bool = False) -> None:
    if output_format == "json":
        print(_json_dump(result.to_dict()))
    elif hasattr(result, "input_reaction_smiles"):
        print(_reaction_concise_summary(result) if concise else _reaction_summary(result))
    else:
        print(_molecule_concise_summary(result) if concise else _molecule_summary(result))


def _command_validate(args: argparse.Namespace) -> int:
    errors = validate_taxonomy()
    payload = {"valid": not errors, "error_count": len(errors), "errors": errors}
    if args.format == "json":
        print(_json_dump(payload))
    elif errors:
        print(f"taxonomy invalid: {len(errors)} error(s)")
        for error in errors:
            print(f"  - {error}")
    else:
        print("taxonomy valid: all definitions and SMARTS passed validation")
    return 0 if not errors else 1


def _command_molecule(args: argparse.Namespace) -> int:
    result = featurize_molecule(
        args.smiles,
        site_types=args.site_type or None,
        include_context_features=not args.no_context,
        label_style=args.label_style,
    )
    _print_result(result, args.format, concise=args.concise)
    return 0 if result.valid else 1


def _command_reaction(args: argparse.Namespace) -> int:
    result = featurize_reaction(
        args.reaction_smiles,
        label_style=args.label_style,
        max_candidates=args.max_candidates,
    )
    _print_result(result, args.format, concise=args.concise)
    return 0 if result.valid else 1


def _detect_column(fieldnames: Sequence[str], mode: str) -> str | None:
    candidates = (
        ("reaction_smiles", "rxn_smiles", "rxn_smiles_clean", "reaction")
        if mode == "reaction"
        else ("smiles", "canonical_smiles", "molecule_smiles")
    )
    normalized = {name.strip().lower(): name for name in fieldnames}
    return next((normalized[name] for name in candidates if name in normalized), None)


def _batch_summary(records: Sequence[dict[str, Any]], mode: str) -> dict[str, Any]:
    valid_count = sum(bool(record["analysis"].get("valid")) for record in records)
    summary: dict[str, Any] = {
        "mode": mode,
        "total": len(records),
        "valid": valid_count,
        "invalid": len(records) - valid_count,
    }
    if mode == "molecule":
        summary["site_type_counts"] = dict(sorted(Counter(
            site["site_type"]
            for record in records
            for site in record["analysis"].get("sites", [])
        ).items()))
        summary["functional_group_counts"] = dict(sorted(Counter(
            group["group_id"]
            for record in records
            for group in record["analysis"].get("functional_groups", [])
        ).items()))
    else:
        summary["evidence_counts"] = dict(sorted(Counter(
            record["analysis"].get("evidence_quality", "unavailable")
            for record in records
        ).items()))
        summary["family_counts"] = dict(sorted(Counter(
            record["analysis"].get("named_family") or "unknown"
            for record in records
        ).items()))
        summary["label_status_counts"] = dict(sorted(Counter(
            record["analysis"].get("reaction_label_status", "unavailable")
            for record in records
        ).items()))
    return summary


def _molecule_csv_columns(
    source_fields: Sequence[str],
    smiles_column: str,
) -> list[str]:
    columns = ["source_row"]
    for field in source_fields:
        if field == "source_row":
            continue
        columns.append(field)
        if field == smiles_column:
            columns.extend(("reactive_site_labels", "functional_group_labels"))
    columns.extend(["valid", "canonical_smiles", "component_count", "reactive_site_count"])
    for site_type in _MOLECULE_SITE_TYPES:
        columns.extend((f"{site_type}_count", f"{site_type}_labels"))
    columns.extend((
        "canonical_signatures", "functional_group_count", "functional_group_ids",
        "warnings", "error",
    ))
    return columns


def _molecule_csv_row(record: dict[str, Any]) -> dict[str, Any]:
    analysis = record["analysis"]
    sites = list(analysis.get("sites") or [])
    groups = list(analysis.get("functional_groups") or [])
    row: dict[str, Any] = {
        "source_row": record["source_row"],
        **record["source"],
        "valid": bool(analysis.get("valid")),
        "canonical_smiles": analysis.get("canonical_smiles") or "",
        "component_count": len(analysis.get("components") or []),
        "reactive_site_count": len(sites),
        "reactive_site_labels": "; ".join(str(site.get("chemist_label") or "") for site in sites),
        "canonical_signatures": "; ".join(str(site.get("canonical_signature") or "") for site in sites),
        "functional_group_count": len(groups),
        "functional_group_ids": "; ".join(str(group.get("group_id") or "") for group in groups),
        "functional_group_labels": "; ".join(str(group.get("chemist_label") or "") for group in groups),
        "warnings": "; ".join(str(value) for value in analysis.get("warnings") or []),
        "error": analysis.get("error") or "",
    }
    for site_type in _MOLECULE_SITE_TYPES:
        selected = [site for site in sites if site.get("site_type") == site_type]
        row[f"{site_type}_count"] = len(selected)
        row[f"{site_type}_labels"] = "; ".join(str(site.get("chemist_label") or "") for site in selected)
    return row


def _reaction_csv_columns(
    source_fields: Sequence[str],
    reaction_smiles_column: str,
) -> list[str]:
    columns = ["source_row"]
    for field in source_fields:
        if field == "source_row":
            continue
        columns.append(field)
        if field == reaction_smiles_column:
            columns.extend(("reaction_label", "spectator_groups"))
    columns.extend((
        "valid", "evidence_quality", "transformation_class", "named_family",
        "reaction_label_status", "candidate_count", "warnings", "error",
    ))
    return columns


def _reaction_csv_row(record: dict[str, Any]) -> dict[str, Any]:
    analysis = record["analysis"]
    spectator_groups = analysis.get("spectator_groups") or []
    return {
        "source_row": record["source_row"],
        **record["source"],
        "valid": bool(analysis.get("valid")),
        "evidence_quality": analysis.get("evidence_quality") or "",
        "transformation_class": analysis.get("transformation_class") or "",
        "named_family": analysis.get("named_family") or "",
        "reaction_label": analysis.get("reaction_label") or "",
        "spectator_groups": "; ".join(
            str(group.get("group_id") or "") for group in spectator_groups
        ),
        "reaction_label_status": analysis.get("reaction_label_status") or "",
        "candidate_count": len(analysis.get("candidates") or []),
        "warnings": "; ".join(str(value) for value in analysis.get("warnings") or []),
        "error": analysis.get("error") or "",
    }


def _write_batch_csv(
    records: Sequence[dict[str, Any]],
    output_path: Path,
    mode: str,
    source_fields: Sequence[str],
    input_column: str,
    *,
    concise: bool = False,
) -> None:
    if concise:
        columns = ["reaction_smiles", "reaction_label", "spectator_groups"]
        row_builder = _reaction_csv_row
    else:
        columns = (
            _molecule_csv_columns(source_fields, input_column)
            if mode == "molecule"
            else _reaction_csv_columns(source_fields, input_column)
        )
        row_builder = _molecule_csv_row if mode == "molecule" else _reaction_csv_row
    with output_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for record in records:
            row = row_builder(record)
            if concise:
                row["reaction_smiles"] = record["source"].get(input_column) or ""
                row = {column: row[column] for column in columns}
            writer.writerow(row)


def _command_batch(args: argparse.Namespace) -> int:
    input_path = Path(args.input)
    with input_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = tuple(reader.fieldnames or ())
        column = args.column or _detect_column(fieldnames, args.mode)
        if not column or column not in fieldnames:
            print(
                f"error: could not find {args.mode} column; available columns: "
                f"{_joined(fieldnames)}",
                file=sys.stderr,
            )
            return 2
        rows = list(reader)

    records: list[dict[str, Any]] = []
    for source_row, row in enumerate(rows, start=2):
        value = str(row.get(column) or "").strip()
        analysis = (
            featurize_reaction(value, label_style=args.label_style, max_candidates=args.max_candidates)
            if args.mode == "reaction"
            else featurize_molecule(value, label_style=args.label_style)
        )
        records.append({
            "source_row": source_row,
            "source": row,
            "analysis": analysis.to_dict(),
        })

    summary = _batch_summary(records, args.mode)
    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_format = args.output_format or ("csv" if output_path.suffix.lower() == ".csv" else "jsonl")
        if args.concise and (args.mode != "reaction" or output_format != "csv"):
            print("error: --concise requires reaction mode and CSV output", file=sys.stderr)
            return 2
        if output_format == "csv":
            _write_batch_csv(
                records,
                output_path,
                args.mode,
                fieldnames,
                column,
                concise=args.concise,
            )
        else:
            with output_path.open("w", encoding="utf-8", newline="") as handle:
                for record in records:
                    handle.write(_json_dump(record, compact=True) + "\n")
        summary["output"] = str(output_path)
        summary["output_format"] = output_format
    if args.format == "json":
        print(_json_dump(summary))
    else:
        print(f"mode: {summary['mode']}")
        print(f"records: {summary['total']}")
        print(f"valid: {summary['valid']}")
        print(f"invalid: {summary['invalid']}")
        for key, value in summary.items():
            if key not in {"mode", "total", "valid", "invalid", "output", "output_format"}:
                print(f"{key}: {_json_dump(value, compact=True)}")
        if "output" in summary:
            print(f"{str(summary['output_format']).upper()} output: {summary['output']}")
    return 0 if summary["invalid"] == 0 else 1


def _command_self_test(args: argparse.Namespace) -> int:
    taxonomy_errors = validate_taxonomy()
    molecule_results = [featurize_molecule(value, label_style=args.label_style) for value in _SELF_TEST_MOLECULES]
    reaction_results = [featurize_reaction(value, label_style=args.label_style) for value in _SELF_TEST_REACTIONS]
    payload = {
        "taxonomy_valid": not taxonomy_errors,
        "taxonomy_errors": taxonomy_errors,
        "molecules": [result.to_dict() for result in molecule_results],
        "reactions": [result.to_dict() for result in reaction_results],
    }
    passed = not taxonomy_errors and all(result.valid for result in molecule_results + reaction_results)
    if args.format == "json":
        print(_json_dump(payload))
    else:
        print(f"taxonomy: {'PASS' if not taxonomy_errors else 'FAIL'}")
        print("\nmolecule feature checks")
        for result in molecule_results:
            print(f"  {'PASS' if result.valid else 'FAIL'} {result.input_smiles}: "
                  f"{len(result.sites)} sites, {len(result.functional_groups)} groups")
        print("\nreaction feature checks")
        for result in reaction_results:
            print(f"  {'PASS' if result.valid else 'FAIL'} {result.reaction_label or result.input_reaction_smiles}: "
                  f"{result.evidence_quality}, family={result.named_family or 'unknown'}")
        print(f"\noverall: {'PASS' if passed else 'FAIL'}")
    return 0 if passed else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m reactive_taxonomy.cli",
        description="Validate and exercise standalone compound/reaction featurization.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate_parser = subparsers.add_parser("validate", help="validate all taxonomy definitions")
    validate_parser.add_argument("--format", choices=("text", "json"), default="text")
    validate_parser.set_defaults(func=_command_validate)

    molecule_parser = subparsers.add_parser("molecule", help="featurize one molecule SMILES")
    molecule_parser.add_argument("smiles")
    molecule_parser.add_argument("--site-type", action="append", choices=_MOLECULE_SITE_TYPES)
    molecule_parser.add_argument("--no-context", action="store_true")
    molecule_parser.add_argument("--concise", action="store_true", help="show only key chemist-readable features")
    molecule_parser.add_argument("--label-style", choices=("unicode", "ascii", "hte_legacy"), default="unicode")
    molecule_parser.add_argument("--format", choices=("text", "json"), default="text")
    molecule_parser.set_defaults(func=_command_molecule)

    reaction_parser = subparsers.add_parser("reaction", help="featurize one reaction SMILES")
    reaction_parser.add_argument("reaction_smiles")
    reaction_parser.add_argument("--max-candidates", type=int, default=500)
    reaction_parser.add_argument("--concise", action="store_true", help="show only key chemist-readable features")
    reaction_parser.add_argument("--label-style", choices=("unicode", "ascii", "hte_legacy"), default="unicode")
    reaction_parser.add_argument("--format", choices=("text", "json"), default="text")
    reaction_parser.set_defaults(func=_command_reaction)

    batch_parser = subparsers.add_parser("batch", help="test every molecule or reaction in a CSV")
    batch_parser.add_argument("input")
    batch_parser.add_argument("--mode", choices=("molecule", "reaction"), required=True)
    batch_parser.add_argument("--column", help="input column; auto-detected when omitted")
    batch_parser.add_argument("--output", help="optional JSONL or CSV result path")
    batch_parser.add_argument("--output-format", choices=("jsonl", "csv"), help="infer from --output suffix when omitted")
    batch_parser.add_argument("--max-candidates", type=int, default=500)
    batch_parser.add_argument(
        "--concise",
        action="store_true",
        help=(
            "for reaction CSV output, write only reaction_smiles, "
            "reaction_label, and spectator_groups"
        ),
    )
    batch_parser.add_argument("--label-style", choices=("unicode", "ascii", "hte_legacy"), default="unicode")
    batch_parser.add_argument("--format", choices=("text", "json"), default="text")
    batch_parser.set_defaults(func=_command_batch)

    self_test_parser = subparsers.add_parser("self-test", help="exercise representative built-in features")
    self_test_parser.add_argument("--label-style", choices=("unicode", "ascii", "hte_legacy"), default="unicode")
    self_test_parser.add_argument("--format", choices=("text", "json"), default="text")
    self_test_parser.set_defaults(func=_command_self_test)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the integrated reactive-taxonomy tester."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())

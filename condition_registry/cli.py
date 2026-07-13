"""Command-line tester and audit tools for the condition registry."""

from __future__ import annotations

import argparse
import json
from typing import Any, Sequence

from .api import get_registry, resolve_substance
from .migration import run_migration_audit
from .models import ResolutionResult
from .validation import validate_registry


_SELF_TEST_CASES = (
    ("CAS resolution", {"cas": "14221-01-3"}, "resolved"),
    ("canonical-name resolution", {"name": "Hydrochloric Acid"}, "resolved"),
    ("alias resolution", {"name": "HCl"}, "resolved"),
    ("invalid CAS rejection", {"cas": "7732-18-4"}, "invalid_identifier"),
    ("unknown name handling", {"name": "not a registered substance"}, "unresolved"),
)


def _json_dump(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True)


def _resolution_summary(result: ResolutionResult) -> str:
    lines = [
        f"query: {result.query}",
        f"status: {result.status}",
        f"match: {result.match_kind or '-'}",
    ]
    if result.substance is not None:
        substance = result.substance
        lines.extend(
            (
                f"substance: {substance.canonical_name}",
                f"substance id: {substance.substance_id}",
                f"CAS: {substance.cas or '-'}",
                f"SMILES: {substance.smiles or '-'}",
            )
        )
        roles = ", ".join(
            f"{role.role_id}/{role.family_id or '-'}" for role in substance.roles
        )
        lines.append(f"roles: {roles or 'none'}")
        lines.append(f"aliases: {', '.join(substance.aliases) or 'none'}")
    if result.candidates:
        lines.append(f"candidates: {', '.join(result.candidates)}")
    if result.warnings:
        lines.append(f"warnings: {', '.join(result.warnings)}")
    return "\n".join(lines)


def _validation_summary(report: dict[str, Any]) -> str:
    lines = [
        f"registry substances: {get_registry().size}",
        f"rows checked: {report['total_rows']}",
        f"accepted rows: {report['accepted_rows']}",
        f"rows with issues: {report['issue_rows']}",
        f"duplicate CAS values: {report['duplicate_cas']}",
        f"duplicate normalized names: {report['duplicate_normalized_names']}",
    ]
    for issue, count in report["issue_counts"].items():
        lines.append(f"  {issue}: {count}")
    return "\n".join(lines)


def _run_self_test(*, output_format: str) -> int:
    results = []
    for label, query, expected_status in _SELF_TEST_CASES:
        result = resolve_substance(**query)
        passed = result.status == expected_status
        results.append(
            {
                "label": label,
                "query": result.query,
                "expected_status": expected_status,
                "actual_status": result.status,
                "passed": passed,
            }
        )
    passed_count = sum(bool(item["passed"]) for item in results)
    payload = {
        "passed": passed_count,
        "failed": len(results) - passed_count,
        "overall": "PASS" if passed_count == len(results) else "FAIL",
        "checks": results,
    }
    if output_format == "json":
        print(_json_dump(payload))
    else:
        for item in results:
            marker = "PASS" if item["passed"] else "FAIL"
            print(
                f"[{marker}] {item['label']}: {item['actual_status']} "
                f"(expected {item['expected_status']})"
            )
        print(f"overall: {payload['overall']} ({passed_count}/{len(results)})")
    return 0 if payload["overall"] == "PASS" else 1


def _add_format_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--format",
        choices=("text", "json"),
        default="text",
        help="Output format (default: text)",
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Test substance resolution and audit the condition registry"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    resolve_parser = subparsers.add_parser(
        "resolve", help="Resolve one exact CAS number, name, or alias"
    )
    query = resolve_parser.add_mutually_exclusive_group(required=True)
    query.add_argument("--cas", help="CAS registry number")
    query.add_argument("--name", help="Canonical name or curated alias")
    _add_format_argument(resolve_parser)

    validate_parser = subparsers.add_parser(
        "validate", help="Validate substance rows against the role/family taxonomy"
    )
    _add_format_argument(validate_parser)

    self_test_parser = subparsers.add_parser(
        "self-test", help="Run deterministic resolution smoke checks"
    )
    _add_format_argument(self_test_parser)

    audit_parser = subparsers.add_parser(
        "audit", help="Write accepted, quarantine, issue, and report files"
    )
    audit_parser.add_argument("output_dir", help="Directory for migration audit files")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the integrated condition-registry CLI tester."""
    args = _build_parser().parse_args(argv)

    if args.command == "resolve":
        result = resolve_substance(cas=args.cas, name=args.name)
        print(
            _json_dump(result.to_dict())
            if args.format == "json"
            else _resolution_summary(result)
        )
        return 0 if result.status == "resolved" else 1

    if args.command == "validate":
        report = validate_registry()
        print(_json_dump(report) if args.format == "json" else _validation_summary(report))
        return 0 if report["issue_rows"] == 0 else 1

    if args.command == "self-test":
        return _run_self_test(output_format=args.format)

    report = run_migration_audit(args.output_dir)
    print(_json_dump(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

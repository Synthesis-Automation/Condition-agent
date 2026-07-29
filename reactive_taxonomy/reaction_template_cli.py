"""CLI for authoring and querying the single-event reaction-template registry."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Sequence

from .reaction_templates import (
    DEFAULT_REACTION_TEMPLATE_REGISTRY_PATH,
    ReactionTemplateError,
    derive_reaction_template,
    load_reaction_template_registry,
    match_reaction_templates,
    upsert_reaction_template,
    validate_reaction_template_registry,
)


def _json_dump(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2)


def _registry_path(args: argparse.Namespace) -> Path:
    return Path(args.registry)


def _mapped_reaction(args: argparse.Namespace) -> str:
    if args.mapped_reaction:
        return str(args.mapped_reaction).strip()
    if args.input:
        return Path(args.input).read_text(encoding="utf-8").strip()
    raise ReactionTemplateError(
        "Provide --mapped-reaction or --input with one mapped reaction SMILES"
    )


def _command_import(args: argparse.Namespace) -> int:
    template = derive_reaction_template(
        _mapped_reaction(args),
        template_id=args.id,
        display_name=args.name,
        family_id=args.family,
        aliases=args.alias or (),
        transformation_class=args.transformation_class,
        status=args.status,
        provenance=args.provenance,
        notes=args.notes or "",
    )
    path = upsert_reaction_template(
        template,
        _registry_path(args),
        replace_existing=args.replace,
    )
    payload = {
        "saved": True,
        "registry": str(path.resolve()),
        "template": template.to_dict(),
    }
    if args.format == "json":
        print(_json_dump(payload))
    else:
        print(f"saved: {template.template_id}")
        print(f"registry: {path.resolve()}")
        print(f"status: {template.status}")
        print(f"event: {template.edit_archetype}; {len(template.edits)} edits")
        print(f"fingerprint: {template.edit_fingerprint}")
    return 0


def _command_list(args: argparse.Namespace) -> int:
    templates = load_reaction_template_registry(_registry_path(args))
    if not args.all:
        templates = tuple(
            template for template in templates if template.status != "retired"
        )
    if args.format == "json":
        print(_json_dump(
            {
                "registry": str(_registry_path(args).resolve()),
                "count": len(templates),
                "templates": [template.to_dict() for template in templates],
            }
        ))
    elif not templates:
        print("no reaction templates")
    else:
        for template in templates:
            family = template.family_id or "-"
            print(
                f"{template.template_id}\t{template.status}\t{family}\t"
                f"{template.display_name}"
            )
    return 0


def _command_show(args: argparse.Namespace) -> int:
    templates = load_reaction_template_registry(_registry_path(args))
    template = next(
        (
            item
            for item in templates
            if item.template_id == args.id
            or args.id.casefold()
            in {
                item.display_name.casefold(),
                *(alias.casefold() for alias in item.aliases),
            }
        ),
        None,
    )
    if template is None:
        raise ReactionTemplateError(f"Unknown reaction template: {args.id}")
    if args.format == "json":
        print(_json_dump(template.to_dict()))
    else:
        print(f"id: {template.template_id}")
        print(f"name: {template.display_name}")
        print(f"family: {template.family_id or '-'}")
        print(f"status: {template.status}")
        print(f"archetype: {template.edit_archetype}")
        print(f"transformation: {template.transformation_class or '-'}")
        print(f"edits: {len(template.edits)}")
        print(f"fingerprint: {template.edit_fingerprint}")
        print(f"definition hash: {template.definition_hash}")
        print(f"mapped reference: {template.mapped_reference_reaction}")
    return 0


def _command_validate(args: argparse.Namespace) -> int:
    errors = validate_reaction_template_registry(_registry_path(args))
    payload = {
        "valid": not errors,
        "registry": str(_registry_path(args).resolve()),
        "error_count": len(errors),
        "errors": errors,
    }
    if args.format == "json":
        print(_json_dump(payload))
    elif errors:
        print(f"registry invalid: {len(errors)} error(s)")
        for error in errors:
            print(f"  - {error}")
    else:
        print("registry valid")
    return 0 if not errors else 1


def _command_match(args: argparse.Namespace) -> int:
    result = match_reaction_templates(
        args.reaction_smiles,
        path=_registry_path(args),
        include_drafts=args.include_drafts,
    )
    if args.format == "json":
        print(_json_dump(result.to_dict()))
    else:
        print(f"valid: {result.valid}")
        print(f"evidence: {result.evidence}")
        print(f"signature: {result.signature_id or '-'}")
        print(f"edit fingerprint: {result.edit_fingerprint or '-'}")
        if result.matches:
            print("matches:")
            for match in result.matches:
                print(
                    f"  {match.template_id}: {match.display_name} "
                    f"[{match.status}; family={match.family_id or '-'}]"
                )
        else:
            print("matches: none")
        if result.warnings:
            print("warnings: " + ", ".join(result.warnings))
    return 0 if result.valid else 1


def build_parser() -> argparse.ArgumentParser:
    """Build the reaction-template authoring parser."""
    parser = argparse.ArgumentParser(
        prog="python -m reactive_taxonomy.reaction_template_cli",
        description=(
            "Import mapped single-event references into a versioned reaction-"
            "template registry and match query-derived edit signatures."
        ),
    )
    parser.add_argument(
        "--registry",
        default=str(DEFAULT_REACTION_TEMPLATE_REGISTRY_PATH),
        help="registry JSON path",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    import_parser = subparsers.add_parser(
        "import", help="derive and save a template from one mapped reference"
    )
    source = import_parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--mapped-reaction")
    source.add_argument("--input", help="UTF-8 file containing reaction SMILES")
    import_parser.add_argument("--id", required=True)
    import_parser.add_argument("--name", required=True)
    import_parser.add_argument("--family")
    import_parser.add_argument("--alias", action="append")
    import_parser.add_argument("--transformation-class")
    import_parser.add_argument(
        "--status", choices=("draft", "active", "retired"), default="draft"
    )
    import_parser.add_argument(
        "--provenance", default="manual_mapped_reference"
    )
    import_parser.add_argument("--notes")
    import_parser.add_argument(
        "--replace",
        action="store_true",
        help="replace an existing template with the same ID",
    )
    import_parser.add_argument(
        "--format", choices=("text", "json"), default="text"
    )
    import_parser.set_defaults(func=_command_import)

    list_parser = subparsers.add_parser("list", help="list registry templates")
    list_parser.add_argument(
        "--all", action="store_true", help="include retired templates"
    )
    list_parser.add_argument(
        "--format", choices=("text", "json"), default="text"
    )
    list_parser.set_defaults(func=_command_list)

    show_parser = subparsers.add_parser("show", help="show one template")
    show_parser.add_argument("id", help="template ID, display name, or alias")
    show_parser.add_argument(
        "--format", choices=("text", "json"), default="text"
    )
    show_parser.set_defaults(func=_command_show)

    validate_parser = subparsers.add_parser(
        "validate", help="validate the complete registry"
    )
    validate_parser.add_argument(
        "--format", choices=("text", "json"), default="text"
    )
    validate_parser.set_defaults(func=_command_validate)

    match_parser = subparsers.add_parser(
        "match", help="match a query-derived edit fingerprint"
    )
    match_parser.add_argument("reaction_smiles")
    match_parser.add_argument(
        "--include-drafts",
        action="store_true",
        help="include draft templates in matches",
    )
    match_parser.add_argument(
        "--format", choices=("text", "json"), default="text"
    )
    match_parser.set_defaults(func=_command_match)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the reaction-template registry CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except (OSError, ReactionTemplateError) as exc:
        if getattr(args, "format", "text") == "json":
            print(_json_dump({"valid": False, "error": str(exc)}))
        else:
            print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = ["build_parser", "main"]

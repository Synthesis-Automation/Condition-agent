#!/usr/bin/env python3

"""
Command-line reagent taxonomy generator.

Usage:
    python -m chemtools.reagent.taxonomy_cli --cas "14221-01-3" --name "Pd(PPh3)4"
    python -m chemtools.reagent.taxonomy_cli --list-families
    python -m chemtools.reagent.taxonomy_cli --cas "14221-01-3"  # Auto-resolve name
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

from .constants import ROLE_FILES
from .taxonomy_store import RoleHeuristics, TaxonomyStore
from .taxonomy_utils import (
    DEFAULT_RESOLVER_TIMEOUT,
    build_entry,
    dedupe_synonyms,
    normalize_cas,
    resolve_identity_from_cas,
    tokenize_all,
)

# Ensure UTF-8 encoding for output
for _stream in ("stdout", "stderr"):
    try:
        getattr(sys, _stream).reconfigure(encoding="utf-8")
    except Exception:
        pass


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Assign reagent role/family and update taxonomy.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Add reagent with auto-detection
  %(prog)s --cas "14221-01-3" --name "Pd(PPh3)4"
  
  # Auto-resolve name from CAS (requires internet)
  %(prog)s --cas "14221-01-3"
  
  # Specify role and family explicitly
  %(prog)s --cas "14221-01-3" --name "Pd(PPh3)4" \\
           --role metal_precursor --family pd_0_complexes
  
  # List all available families
  %(prog)s --list-families
  
  # Dry run (preview without saving)
  %(prog)s --cas "14221-01-3" --name "Pd(PPh3)4" --dry-run
        """
    )
    parser.add_argument("--cas", help="CAS number for the reagent.")
    parser.add_argument("--name", help="Primary name for the reagent.")
    parser.add_argument("--abbr", help="Short abbreviation.")
    parser.add_argument(
        "--synonym",
        action="append",
        default=[],
        help="Additional synonym (repeatable)."
    )
    parser.add_argument(
        "--role",
        choices=sorted(ROLE_FILES.keys()),
        help="Override detected role."
    )
    parser.add_argument("--family", help="Explicit family identifier to use.")
    parser.add_argument("--smiles", help="Optional SMILES string to store with the entry.")
    parser.add_argument(
        "--taxonomy-dir",
        default="data/compound_taxonomy",
        help="Path to taxonomy directory."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Infer role/family but do not write files."
    )
    parser.add_argument(
        "--list-families",
        action="store_true",
        help="Print known families and exit."
    )
    parser.add_argument(
        "--allow-default-family",
        action="store_true",
        help="Allow fallback default family when inference fails."
    )
    parser.add_argument(
        "--no-auto-resolve",
        action="store_true",
        help="Disable CAS-based identity lookup when name is missing."
    )
    parser.add_argument(
        "--resolver-timeout",
        type=float,
        default=DEFAULT_RESOLVER_TIMEOUT,
        help="Timeout (seconds) for resolver HTTP requests."
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Emit extra debug details."
    )
    return parser.parse_args()


def main() -> None:
    """Main CLI entry point."""
    args = parse_args()
    taxonomy_dir = Path(args.taxonomy_dir)
    base_dir: Optional[Path] = taxonomy_dir if taxonomy_dir.exists() else None
    if base_dir is None and args.verbose:
        print(f"# Using unified taxonomy registry (legacy directory '{taxonomy_dir}' not found).")
    store = TaxonomyStore(base_dir)
    heuristics = RoleHeuristics(store)

    # List families and exit
    if args.list_families:
        for role in sorted(ROLE_FILES.keys()):
            print(f"[{role}]")
            for _role, fid, label in store.list_families(role):
                label_str = f" - {label}" if label else ""
                print(f"  {fid}{label_str}")
        return

    if not args.dry_run and store.base_dir is None:
        raise SystemExit(
            "Cannot modify taxonomy because no legacy taxonomy directory is available. "
            "Provide --taxonomy-dir pointing to a writable directory."
        )

    # Validate CAS is provided
    if not args.cas:
        raise SystemExit("--cas is required unless --list-families is used.")

    # Normalize CAS
    cas = normalize_cas(args.cas)
    
    # Auto-resolve identity if name not provided
    resolved_identity: Optional[Dict[str, Any]] = None
    if not args.no_auto_resolve and not args.name:
        resolved_identity = resolve_identity_from_cas(cas, timeout=args.resolver_timeout)
        if args.verbose:
            if resolved_identity:
                print(f"# auto-resolved via {resolved_identity.get('source')}: {resolved_identity.get('name')}")
            else:
                print("# auto-resolve failed to supply a name")

    # Determine name
    name = args.name or (resolved_identity.get("name") if resolved_identity else None)
    if not name:
        raise SystemExit(
            "Unable to determine reagent name. Provide --name or use --no-auto-resolve to skip lookup."
        )

    # Build synonym list
    abbr = args.abbr or name
    resolved_synonyms = resolved_identity.get("synonyms", []) if resolved_identity else []
    synonyms = dedupe_synonyms([name, abbr, *args.synonym, *resolved_synonyms])
    input_tokens = tokenize_all([name, *synonyms])
    
    # Initialize variables for role/family detection
    role = args.role
    family_id = args.family
    used_default = False
    default_rejection_reason: Optional[str] = None
    family_reason: Optional[List[str]] = None
    role_reason: Optional[str] = None
    auto_resolve_source = resolved_identity.get("source") if resolved_identity else None
    resolved_smiles = resolved_identity.get("smiles") if resolved_identity else None

    # Validate family if provided
    if family_id:
        family_role = store.role_for_family(family_id)
        if not family_role:
            raise SystemExit(
                f"Unknown family '{family_id}'. Use --list-families to inspect available options."
            )
        if role and role != family_role:
            raise SystemExit(
                f"Provided role '{role}' conflicts with family '{family_id}' (expected role '{family_role}')."
            )
        role = family_role

    # Infer family if not provided
    inference = heuristics.infer_family(name, synonyms) if not family_id else None
    if inference:
        role, inferred_family, reason_tokens = inference
        family_id = family_id or inferred_family
        family_reason = reason_tokens

    # Infer role if not determined
    if not role:
        role_inference = heuristics.infer_role(name, synonyms)
        if role_inference:
            role, pattern = role_inference
            role_reason = pattern

    # Try default family if needed
    if not family_id:
        if role:
            default_family = heuristics.default_family_for_role(role)
            if default_family and args.allow_default_family:
                if store.family_token_overlap(role, default_family, input_tokens):
                    family_id = default_family
                    used_default = True
                else:
                    tokens_sample = ', '.join(sorted(input_tokens)[:6]) or 'none'
                    family_tokens = store.family_tokens.get((role, default_family), set())
                    family_sample = ', '.join(sorted(family_tokens)[:6]) or 'none'
                    default_rejection_reason = (
                        f"default family '{default_family}' rejected: no token overlap "
                        f"(input tokens: {tokens_sample}; family tokens sample: {family_sample})"
                    )
        if not family_id:
            message = (
                "Unable to determine family. Provide --family explicitly or use --allow-default-family."
            )
            if default_rejection_reason:
                message += f" Automatic fallback was skipped because {default_rejection_reason}."
            raise SystemExit(message)

    # Final role determination
    if not role:
        role = store.role_for_family(family_id)
    if not role:
        raise SystemExit("Unable to determine role; please supply --role.")

    # Check if entry already exists
    existing = store.find_by_cas(cas)
    if existing:
        existing_role, existing_family, data = existing
        result = {
            "cas": cas,
            "name": data.get("name"),
            "role": existing_role,
            "family_id": existing_family,
            "status": "exists",
        }
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return

    # Build entry
    family_data = store.family_data(role, family_id)
    numeric = store.numeric_baseline(role, family_id)
    entry = build_entry(
        role,
        family_data,
        cas,
        name,
        abbr,
        synonyms,
        args.smiles or resolved_smiles,
        numeric
    )

    # Prepare result
    result = {
        "cas": cas,
        "name": name,
        "role": role,
        "family_id": family_id,
        "taxonomy_file": str(store.file_for_role(role)),
        "dry_run": args.dry_run,
        "used_default_family": used_default,
    }
    if auto_resolve_source:
        result["auto_resolve_source"] = auto_resolve_source
    if resolved_smiles and not args.smiles:
        result["smiles_source"] = auto_resolve_source or "resolver"
    if family_reason:
        result["family_tokens"] = family_reason
    if role_reason:
        result["role_pattern"] = role_reason

    # Dry run or write
    if args.dry_run:
        result["status"] = "dry_run"
        result["entry_preview"] = entry
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return

    # Add entry and save
    store.add_entry(role, family_id, entry)
    path = store.save_role(role)
    result["status"] = "written"
    result["written_to"] = str(path)
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

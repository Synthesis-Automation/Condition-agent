#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Demonstrate attachment-point SMARTS templating for the Suzuki + Buchwald POC.

This script reads `smarts_templates.suzuki_buchwald.poc.json` and:
  1) Expands templates + fragments into SMARTS patterns
  2) Verifies the expansion reproduces `generated_atomic_features[].smarts`
  3) Optionally cross-checks that those patterns exist in the runtime atomic feature file
  4) Optionally compiles patterns with RDKit and shows attachment-point atom maps (e.g. ":1")

Run from repo root:
  python -m chemtools.taxonomy.new.demo_attachment_point_templating

Run from this folder:
  python demo_attachment_point_templating.py
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Optional


_PLACEHOLDER_RE = re.compile(r"{([a-zA-Z_][a-zA-Z0-9_]*)}")


def _ensure_repo_root_on_path() -> None:
    try:
        __import__("chemtools")
        return
    except ModuleNotFoundError:
        repo_root = Path(__file__).resolve().parents[3]
        sys.path.insert(0, str(repo_root))


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _infer_placeholders(pattern: str) -> list[str]:
    placeholders: list[str] = []
    seen: set[str] = set()
    for match in _PLACEHOLDER_RE.finditer(pattern):
        name = match.group(1)
        if name in seen:
            continue
        placeholders.append(name)
        seen.add(name)
    return placeholders


def _expand_pattern(pattern: str, placeholder_to_smarts: Mapping[str, str]) -> str:
    try:
        return pattern.format(**dict(placeholder_to_smarts))
    except KeyError as e:
        raise ValueError(f"Missing placeholder {e} for pattern: {pattern}") from e


@dataclass(frozen=True)
class GeneratedEntry:
    token: str
    template: str
    fragments: tuple[str, ...]
    expected_smarts: str


def _load_generated_entries(spec: Mapping[str, Any]) -> list[GeneratedEntry]:
    items = spec.get("generated_atomic_features", [])
    if not isinstance(items, list):
        raise ValueError("spec.generated_atomic_features must be a list")
    out: list[GeneratedEntry] = []
    for entry in items:
        if not isinstance(entry, dict):
            continue
        token = str(entry.get("token", "")).strip()
        template = str(entry.get("template", "")).strip()
        fragments = entry.get("fragments", [])
        smarts = str(entry.get("smarts", "")).strip()
        if not token or not template or not smarts:
            continue
        if not isinstance(fragments, list) or not all(isinstance(f, str) for f in fragments):
            raise ValueError(f"{token}: fragments must be a list[str]")
        out.append(GeneratedEntry(token=token, template=template, fragments=tuple(fragments), expected_smarts=smarts))
    return out


def _example_smiles(template: str, fragments: tuple[str, ...]) -> Optional[str]:
    if template == "core_halide":
        core, lg = fragments
        if core == "AromaticCore":
            return {"Cl": "Clc1ccccc1", "Br": "Brc1ccccc1", "I": "Ic1ccccc1"}.get(lg)
        if core == "VinylCore":
            return {"Cl": "C=CCl", "Br": "C=CBr", "I": "C=CI"}.get(lg)
        return None

    if template == "core_sulfonate":
        core, tail = fragments
        if core == "AromaticCore":
            return {
                "OTfTail": "FC(F)(F)S(=O)(=O)Oc1ccccc1",
                "OTsTail": "Cc1ccc(S(=O)(=O)Oc2ccccc2)cc1",
                "OMsTail": "CS(=O)(=O)Oc1ccccc1",
            }.get(tail)
        if core == "VinylCore":
            return {
                "OTfTail": "C=COS(=O)(=O)C(F)(F)F",
                "OTsTail": "C=COS(=O)(=O)c1ccc(C)cc1",
                "OMsTail": "C=COS(=O)(=O)C",
            }.get(tail)
        return None

    if template == "core_boron":
        core, boron = fragments
        if boron != "BoronicAcid":
            return None
        if core == "AromaticCore":
            return "OB(O)c1ccccc1"
        if core == "VinylCore":
            return "C=CB(O)O"
        return None

    return None


def main(argv: list[str] | None = None) -> int:
    _ensure_repo_root_on_path()

    from chemtools.util.rdkit_helpers import parse_smiles, rdkit_available
    from chemtools.util.smarts_cache import compile_smarts

    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", default="smarts_templates.suzuki_buchwald.poc.json", help="Template spec JSON file")
    ap.add_argument(
        "--atomic",
        default="calculable_features.atomic.suzuki_buchwald.poc.json",
        help="Runtime atomic feature JSON file (for cross-check)",
    )
    ap.add_argument("--limit", type=int, default=8, help="How many entries to print (0 = all)")
    ap.add_argument("--no-atomic-check", action="store_true", help="Skip cross-check against atomic feature file")
    ap.add_argument("--no-rdkit", action="store_true", help="Skip RDKit compilation/matching even if available")
    ap.add_argument("--no-examples", action="store_true", help="Skip example SMILES matching output")
    args = ap.parse_args(argv)

    folder = Path(__file__).resolve().parent
    spec_path = (folder / args.spec) if not Path(args.spec).is_absolute() else Path(args.spec)
    atomic_path = (folder / args.atomic) if not Path(args.atomic).is_absolute() else Path(args.atomic)

    spec = _read_json(spec_path)
    if not isinstance(spec, dict):
        raise ValueError(f"Expected object in {spec_path}, got {type(spec).__name__}")

    fragments = spec.get("fragments", {})
    templates = spec.get("templates", {})
    if not isinstance(fragments, dict) or not isinstance(templates, dict):
        raise ValueError("spec.fragments and spec.templates must be objects")

    fragment_smarts: dict[str, str] = {}
    for name, meta in fragments.items():
        if isinstance(meta, dict) and isinstance(meta.get("smarts"), str):
            fragment_smarts[str(name)] = meta["smarts"]

    template_patterns: dict[str, str] = {}
    for name, meta in templates.items():
        if isinstance(meta, dict) and isinstance(meta.get("pattern"), str):
            template_patterns[str(name)] = meta["pattern"]

    generated = _load_generated_entries(spec)
    if not generated:
        raise ValueError("No generated_atomic_features entries found in spec.")

    atomic_by_token: dict[str, dict[str, Any]] = {}
    if not args.no_atomic_check:
        atomic_payload = _read_json(atomic_path)
        if not isinstance(atomic_payload, list):
            raise ValueError(f"Expected list in {atomic_path}, got {type(atomic_payload).__name__}")
        for entry in atomic_payload:
            if isinstance(entry, dict) and entry.get("token"):
                atomic_by_token[str(entry["token"])] = entry

    can_rdkit = (not args.no_rdkit) and rdkit_available()
    if args.no_rdkit:
        can_rdkit = False

    mismatches: list[str] = []
    missing_in_atomic: list[str] = []
    compiled_failures: list[str] = []

    print(
        f"Loaded {len(fragment_smarts)} fragments, {len(template_patterns)} templates, "
        f"{len(generated)} generated entries from {spec_path.name}"
    )

    shown = 0
    for entry in generated:
        pattern = template_patterns.get(entry.template)
        if not pattern:
            mismatches.append(f"{entry.token}: unknown template '{entry.template}'")
            continue

        placeholders = _infer_placeholders(pattern)
        if len(placeholders) != len(entry.fragments):
            mismatches.append(
                f"{entry.token}: placeholder/fragment count mismatch "
                f"(placeholders={placeholders}, fragments={list(entry.fragments)})"
            )
            continue

        mapping: dict[str, str] = {}
        for name, frag_id in zip(placeholders, entry.fragments, strict=True):
            frag_smarts = fragment_smarts.get(frag_id)
            if frag_smarts is None:
                mismatches.append(f"{entry.token}: unknown fragment '{frag_id}'")
                frag_smarts = ""
            mapping[name] = frag_smarts

        expanded = _expand_pattern(pattern, mapping)
        if expanded != entry.expected_smarts:
            mismatches.append(
                f"{entry.token}: expansion mismatch (expanded='{expanded}', expected='{entry.expected_smarts}')"
            )

        if atomic_by_token:
            atomic_entry = atomic_by_token.get(entry.token)
            if not atomic_entry:
                missing_in_atomic.append(entry.token)
            else:
                detect = atomic_entry.get("detect") or {}
                smarts_any = detect.get("smarts_any") or []
                if expanded not in smarts_any:
                    missing_in_atomic.append(entry.token)

        attachment_maps: list[int] = []
        example_info = ""

        if can_rdkit:
            try:
                patt = compile_smarts(expanded, validate=True)
                if patt is not None:
                    attachment_maps = sorted(
                        {a.GetAtomMapNum() for a in patt.GetAtoms() if int(a.GetAtomMapNum() or 0) > 0}
                    )

                    if not args.no_examples:
                        example = _example_smiles(entry.template, entry.fragments)
                        if example:
                            mol = parse_smiles(example)
                            if mol is not None:
                                matches = mol.GetSubstructMatches(patt)
                                hit = bool(matches)
                                example_info = f" example={example} match={hit}"
                                if hit and attachment_maps:
                                    map_atoms = [
                                        (idx, a.GetAtomMapNum())
                                        for idx, a in enumerate(patt.GetAtoms())
                                        if int(a.GetAtomMapNum() or 0) > 0
                                    ]
                                    first = matches[0]
                                    attached = {m: int(first[i]) for i, m in map_atoms}
                                    example_info += f" attachment_atoms={attached}"
            except Exception as e:
                compiled_failures.append(f"{entry.token}: {e}")

        if args.limit != 0 and shown >= args.limit:
            continue

        shown += 1
        frag_str = ", ".join(entry.fragments)
        maps_str = f" maps={attachment_maps}" if attachment_maps else ""
        print(f"- {entry.token}: {entry.template} [{frag_str}] -> {expanded}{maps_str}{example_info}")

    if mismatches:
        print("\nExpansion mismatches:", file=sys.stderr)
        for m in mismatches:
            print(f"- {m}", file=sys.stderr)
    if missing_in_atomic:
        print("\nMissing or inconsistent in atomic feature file:", file=sys.stderr)
        for t in sorted(set(missing_in_atomic)):
            print(f"- {t}", file=sys.stderr)
    if compiled_failures:
        print("\nSMARTS compilation failures:", file=sys.stderr)
        for m in compiled_failures:
            print(f"- {m}", file=sys.stderr)

    if mismatches or missing_in_atomic or compiled_failures:
        return 1

    print("\nOK: templating expansion verified.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


#!/usr/bin/env python3
"""
Convert legacy chemtools taxonomy files into the unified schema.

Usage:
    python scripts/taxonomy/convert_legacy_taxonomies.py \
        --output-dir chemtools/taxonomy/data \
        --taxonomy-version 0.1.0
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

# Ensure repository root is importable (works when script launched directly)
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.reagent import constants as reagent_constants  # type: ignore  # noqa: E402
from chemtools.taxonomy.registry import TaxonomyRegistry  # type: ignore  # noqa: E402

LEGACY_REACTANT_TYPES = REPO_ROOT / "chemtools" / "reagent" / "data" / "reactant_types.json"
LEGACY_REACTION_TYPES = REPO_ROOT / "chemtools" / "reagent" / "data" / "reaction_types.json"
LEGACY_REAGENT_FAMILIES = REPO_ROOT / "chemtools" / "reagent" / "reagent_schema" / "families_registry.json"

ALIAS_ENTITY_REACTION_TYPE = "reaction_type"
ALIAS_ENTITY_REACTANT_TYPE = "reactant_type"
ALIAS_ENTITY_REAGENT_ROLE = "reagent_role"
ALIAS_ENTITY_REAGENT_FAMILY = "reagent_family"

_INVALID_CHAR_RE = re.compile(r"[^0-9a-zA-Z]+")


def slugify(value: str) -> str:
    """Convert an arbitrary string into a snake_case slug."""
    clean = _INVALID_CHAR_RE.sub("_", value.strip())
    clean = re.sub(r"_+", "_", clean).strip("_").lower()
    if not clean:
        clean = "unnamed"
    if clean[0].isdigit():
        clean = f"_{clean}"
    return clean


def ensure_unique(slug: str, used: set[str]) -> str:
    """Ensure slug uniqueness within ``used`` by appending numeric suffixes."""
    candidate = slug
    idx = 2
    while candidate in used:
        candidate = f"{slug}_{idx}"
        idx += 1
    used.add(candidate)
    return candidate


def read_json(path: Path) -> dict | list:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False, sort_keys=True)


def add_alias(
    aliases: Dict[str, dict],
    alias: Optional[str],
    *,
    entity_type: str,
    entity_id: str,
    source: Optional[str] = None,
    notes: Optional[str] = None,
) -> None:
    if not alias:
        return
    key = alias.strip()
    if not key:
        return
    lowered = key.lower()
    record = aliases.get(lowered)
    if record and record["entity_id"] != entity_id:
        # Conflicting alias mapping; keep the first mapping and annotate the clash.
        record.setdefault("notes", "")
        extras = f" clash:{entity_id}"
        if extras not in record["notes"]:
            record["notes"] = (record["notes"] + extras).strip()
        return
    aliases[lowered] = {
        "alias": key,
        "entity_type": entity_type,
        "entity_id": entity_id,
        "source": source,
        "notes": notes,
    }


@dataclass
class ReactantTypeBuilder:
    aliases: Dict[str, dict]
    used_ids: set = field(default_factory=set)
    reactant_entries: List[dict] = field(default_factory=list)
    token_map: Dict[str, str] = field(default_factory=dict)

    def load_legacy(self, payload: dict) -> None:
        for token, entry in payload.items():
            canonical_id = ensure_unique(slugify(token), self.used_ids)
            record = {
                "id": canonical_id,
                "name": entry.get("name") or token,
                "description": entry.get("description"),
                "category": entry.get("category"),
                "smarts": entry.get("smarts"),
                "members": [
                    {
                        "id": member.get("id"),
                        "name": member.get("name") or member.get("id"),
                        "smarts": member.get("smarts"),
                        "aliases": member.get("aliases", []),
                        "metadata": member.get("metadata", {}),
                    }
                    for member in entry.get("members", [])
                ],
                "aliases": entry.get("aliases", []),
                "metadata": entry.get("metadata", {}),
            }
            self.reactant_entries.append(record)
            self._register_token(token, canonical_id)
            add_alias(
                self.aliases,
                token,
                entity_type=ALIAS_ENTITY_REACTANT_TYPE,
                entity_id=canonical_id,
                source="legacy/reactant_types.json",
            )
            if entry.get("aliases"):
                for alias in entry["aliases"]:
                    add_alias(
                        self.aliases,
                        alias,
                        entity_type=ALIAS_ENTITY_REACTANT_TYPE,
                        entity_id=canonical_id,
                        source="legacy/reactant_types.json",
                        notes="top_level_alias",
                    )
            for member in entry.get("members", []):
                member_id = member.get("id")
                if member_id:
                    add_alias(
                        self.aliases,
                        member_id,
                        entity_type=ALIAS_ENTITY_REACTANT_TYPE,
                        entity_id=canonical_id,
                        source="legacy/reactant_types.json",
                        notes="member_id",
                    )
                    self._register_token(member_id, canonical_id)

    def ensure_token(self, token: str, *, source: str) -> str:
        if not token:
            return ensure_unique("unknown_reactant", self.used_ids)
        key = token.strip()
        if not key:
            return ensure_unique("unknown_reactant", self.used_ids)
        canonical = self.token_map.get(key)
        if canonical:
            return canonical
        canonical = self.token_map.get(key.lower())
        if canonical:
            self._register_token(key, canonical)
            return canonical

        canonical_base = slugify(key)
        canonical_id = ensure_unique(canonical_base, self.used_ids)
        record = {
            "id": canonical_id,
            "name": key,
            "description": None,
            "category": None,
            "smarts": None,
            "members": [],
            "aliases": [],
            "metadata": {"generated_from_token": key, "source": source},
        }
        self.reactant_entries.append(record)
        self._register_token(key, canonical_id)
        add_alias(
            self.aliases,
            key,
            entity_type=ALIAS_ENTITY_REACTANT_TYPE,
            entity_id=canonical_id,
            source=source,
            notes="generated_from_token",
        )
        return canonical_id

    def _register_token(self, token: str, canonical_id: str) -> None:
        self.token_map[token] = canonical_id
        self.token_map[token.lower()] = canonical_id

    def build(self) -> Tuple[List[dict], Dict[str, str]]:
        return self.reactant_entries, self.token_map


def convert_reagent_roles(additional_roles: Optional[Iterable[str]] = None) -> Tuple[List[dict], Dict[str, str]]:
    roles: List[dict] = []
    alias_map: Dict[str, str] = {}
    used: set[str] = set()
    base_roles = set(reagent_constants.ROLE_FILES.keys()) | set(reagent_constants.DEFAULT_FAMILY_BY_ROLE.keys()) | set(reagent_constants.ROLE_PRIORITY.keys())
    if additional_roles:
        base_roles |= set(additional_roles)

    for role in sorted(base_roles):
        canonical_id = ensure_unique(slugify(role), used)
        metadata = {}
        taxonomy_file = reagent_constants.ROLE_FILES.get(role)
        if taxonomy_file:
            metadata["taxonomy_file"] = taxonomy_file
        entry = {
            "id": canonical_id,
            "name": role.replace("_", " ").title(),
            "priority": int(reagent_constants.ROLE_PRIORITY.get(role, 100)),
            "default_family_id": reagent_constants.DEFAULT_FAMILY_BY_ROLE.get(role),
            "description": None,
            "metadata": metadata,
        }
        roles.append(entry)
        alias_map[role] = canonical_id
        alias_map[role.lower()] = canonical_id
        alias_map[slugify(role)] = canonical_id
        alias_map[canonical_id] = canonical_id
    return roles, alias_map


def convert_reagent_families(path: Path, alias_records: Dict[str, dict]) -> Tuple[List[dict], Dict[str, List[str]], set[str]]:
    payload = read_json(path)
    entries = payload.get("entries", []) if isinstance(payload, dict) else payload
    families: List[dict] = []
    role_to_families: Dict[str, List[str]] = {}
    roles_used: set[str] = set()
    for entry in entries:
        family_id = entry["family"]
        record = {
            "id": family_id,
            "role_id": entry["role"],
            "name": entry.get("definition") or family_id.replace("_", " "),
            "description": entry.get("definition"),
            "precedence": entry.get("precedence"),
            "include_smarts": entry.get("include_smarts", []),
            "exclude_smarts": entry.get("exclude_smarts", []),
            "required_props": entry.get("required_props", {}),
            "keywords": entry.get("keywords", []),
            "examples_pos": entry.get("examples_pos", []),
            "examples_neg": entry.get("examples_neg", []),
            "notes": entry.get("notes"),
        }
        families.append(record)
        role_to_families.setdefault(entry["role"], []).append(family_id)
        roles_used.add(entry["role"])
        add_alias(
            alias_records,
            family_id,
            entity_type=ALIAS_ENTITY_REAGENT_FAMILY,
            entity_id=family_id,
            source="legacy/families_registry.json",
        )
    return families, role_to_families, roles_used


def convert_reactions(
    payload: dict,
    reactant_builder: ReactantTypeBuilder,
    alias_records: Dict[str, dict],
) -> Tuple[List[dict], List[dict]]:
    categories: List[dict] = []
    reaction_types: List[dict] = []
    used_category_ids: set[str] = set()
    used_reaction_ids: set[str] = set()

    for category_key, category_data in payload.items():
        category_name = category_data.get("category") or category_key
        category_id = ensure_unique(slugify(category_name), used_category_ids)
        categories.append(
            {
                "id": category_id,
                "name": category_name,
                "description": category_data.get("description"),
            }
        )
        for reaction in category_data.get("reactions", []):
            base_slug = slugify(reaction.get("id") or reaction.get("name") or "reaction")
            reaction_id = ensure_unique(base_slug, used_reaction_ids)
            reactant_requirements = []
            for group in reaction.get("reactants", []):
                tokens = group if isinstance(group, list) else [group]
                canonical_id = None
                collected_tokens: List[str] = []
                for token in tokens:
                    canonical = reactant_builder.ensure_token(
                        token,
                        source=f"legacy/reaction_types.json:{reaction.get('id')}",
                    )
                    if canonical_id is None:
                        canonical_id = canonical
                    collected_tokens.append(token)
                reactant_requirements.append(
                    {
                        "reactant_type_id": canonical_id or "unknown_reactant",
                        "stoichiometry": None,
                        "notes": None,
                        "original_tokens": collected_tokens,
                    }
                )
            metadata = {}
            if reaction.get("typical_catalysts"):
                metadata["typical_catalysts"] = reaction["typical_catalysts"]
            if reaction.get("typical_conditions"):
                metadata["typical_conditions"] = reaction["typical_conditions"]
            if reaction.get("notes"):
                metadata["notes"] = reaction["notes"]
            reaction_types.append(
                {
                    "id": reaction_id,
                    "category_id": category_id,
                    "name": reaction.get("name") or reaction.get("id") or reaction_id,
                    "description": reaction.get("description"),
                    "aliases": reaction.get("aliases", []),
                    "reactants": reactant_requirements,
                    "required_roles": [],
                    "metadata": metadata,
                    "source_ids": [reaction["id"]] if reaction.get("id") else [],
                }
            )
            add_alias(
                alias_records,
                reaction.get("id"),
                entity_type=ALIAS_ENTITY_REACTION_TYPE,
                entity_id=reaction_id,
                source="legacy/reaction_types.json",
                notes="original_id",
            )
            add_alias(
                alias_records,
                reaction.get("name"),
                entity_type=ALIAS_ENTITY_REACTION_TYPE,
                entity_id=reaction_id,
                source="legacy/reaction_types.json",
                notes="display_name",
            )
            for alias in reaction.get("aliases", []):
                add_alias(
                    alias_records,
                    alias,
                    entity_type=ALIAS_ENTITY_REACTION_TYPE,
                    entity_id=reaction_id,
                    source="legacy/reaction_types.json",
                    notes="legacy_alias",
                )
    return categories, reaction_types


def build_manifest(schema_version: str, taxonomy_version: str) -> dict:
    now = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    return {
        "schema_version": schema_version,
        "taxonomy_version": taxonomy_version,
        "generated_at": now,
        "source": "scripts/taxonomy/convert_legacy_taxonomies.py",
        "notes": "Generated from legacy chemtools taxonomy assets.",
    }


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Convert legacy taxonomy assets into unified schema files.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory to write unified taxonomy JSON files (default: chemtools/taxonomy/data).",
    )
    parser.add_argument(
        "--taxonomy-version",
        default="0.1.0",
        help="Version string to embed in manifest (default: 0.1.0).",
    )
    parser.add_argument(
        "--schema-version",
        default="0.1.0",
        help="Schema version recorded in manifest (default: 0.1.0).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned updates without writing files.",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    output_dir = args.output_dir or (REPO_ROOT / "chemtools" / "taxonomy" / "data")

    alias_records: Dict[str, dict] = {}

    # Reactant types
    reactant_builder = ReactantTypeBuilder(aliases=alias_records)
    legacy_reactants = read_json(LEGACY_REACTANT_TYPES)
    if not isinstance(legacy_reactants, dict):
        raise RuntimeError("Expected reactant_types.json to contain an object mapping IDs to definitions.")
    reactant_builder.load_legacy(legacy_reactants)

    # Reagent families (collect role usage)
    reagent_families, role_to_families, roles_used = convert_reagent_families(LEGACY_REAGENT_FAMILIES, alias_records)

    # Reagent roles (including roles referenced by families)
    reagent_roles, role_alias_map = convert_reagent_roles(additional_roles=roles_used)
    for alias, canonical_id in role_alias_map.items():
        add_alias(
            alias_records,
            alias,
            entity_type=ALIAS_ENTITY_REAGENT_ROLE,
            entity_id=canonical_id,
            source="legacy/reagent_roles",
        )

    for family in reagent_families:
        role_alias = family["role_id"]
        canonical = (
            role_alias_map.get(role_alias)
            or role_alias_map.get(role_alias.lower())
            or role_alias_map.get(slugify(role_alias))
        )
        if canonical is None:
            raise RuntimeError(f"Unable to resolve reagent role '{role_alias}' referenced by family '{family['id']}'")
        family["role_id"] = canonical

    # Reaction categories/types
    legacy_reactions = read_json(LEGACY_REACTION_TYPES)
    if not isinstance(legacy_reactions, dict):
        raise RuntimeError("Expected reaction_types.json to contain category mapping.")
    reaction_categories, reaction_types = convert_reactions(
        legacy_reactions,
        reactant_builder=reactant_builder,
        alias_records=alias_records,
    )

    reactant_types, reactant_token_map = reactant_builder.build()

    manifest = build_manifest(args.schema_version, args.taxonomy_version)

    payloads = {
        "manifest.json": manifest,
        "reaction_categories.json": reaction_categories,
        "reaction_types.json": reaction_types,
        "reactant_types.json": reactant_types,
        "reagent_roles.json": reagent_roles,
        "reagent_families.json": reagent_families,
        "aliases.json": sorted(alias_records.values(), key=lambda item: item["alias"].lower()),
    }

    if args.dry_run:
        print("Dry-run conversion summary:")
        for name, payload in payloads.items():
            if isinstance(payload, list):
                print(f"  - {name}: {len(payload)} entries")
            elif isinstance(payload, dict):
                print(f"  - {name}: keys={list(payload.keys())[:5]}{'...' if len(payload) > 5 else ''}")
        return 0

    for filename, payload in payloads.items():
        write_json(output_dir / filename, payload)

    # Validate by reloading
    TaxonomyRegistry.from_path(output_dir)

    print("Unified taxonomy written to:", output_dir)
    print(f"  Categories: {len(reaction_categories)}")
    print(f"  Reaction types: {len(reaction_types)}")
    print(f"  Reactant types: {len(reactant_types)}")
    print(f"  Reagent roles: {len(reagent_roles)}")
    print(f"  Reagent families: {len(reagent_families)}")
    print(f"  Aliases: {len(payloads['aliases.json'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

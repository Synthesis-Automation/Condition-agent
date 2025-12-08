"""
Taxonomy registry implementation.
"""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from . import models

_MANIFEST_FILE = "manifest.json"
_REACTION_CATEGORIES_FILE = "reaction_categories.json"
_REACTION_TYPES_FILE = "reaction_types.json"
_REACTANT_TYPES_FILE = "reactant_types.json"
_REAGENT_ROLES_FILE = "reagent_roles.json"
_REAGENT_FAMILIES_FILE = "reagent_families.json"
_ALIASES_FILE = "aliases.json"


class TaxonomyRegistry:
    """
    In-memory representation of the unified taxonomy.
    """

    def __init__(
        self,
        *,
        root: Path,
        manifest: models.TaxonomyManifest,
        reaction_categories: Dict[str, models.ReactionCategory],
        reaction_types: Dict[str, models.ReactionType],
        reactant_types: Dict[str, models.ReactantType],
        reagent_roles: Dict[str, models.ReagentRole],
        reagent_families: Dict[str, models.ReagentFamily],
        aliases: Dict[str, models.AliasRecord],
    ) -> None:
        self.root = root
        self.manifest = manifest
        self.reaction_categories = reaction_categories
        self.reaction_types = reaction_types
        self.reactant_types = reactant_types
        self.reagent_roles = reagent_roles
        self.reagent_families = reagent_families
        self.aliases = aliases

    # ------------------------------------------------------------------ #
    # Construction
    # ------------------------------------------------------------------ #
    @classmethod
    def from_path(cls, root: Path) -> "TaxonomyRegistry":
        """
        Build a registry by reading taxonomy JSON files from ``root``.
        """
        root = root.resolve()
        if not root.exists():
            raise FileNotFoundError(f"Taxonomy directory not found: {root}")

        manifest = _load_manifest(root / _MANIFEST_FILE)
        reaction_categories = _load_reaction_categories(root / _REACTION_CATEGORIES_FILE)
        reactant_types = _load_reactant_types(root / _REACTANT_TYPES_FILE)
        reagent_roles = _load_reagent_roles(root / _REAGENT_ROLES_FILE)
        reagent_families = _load_reagent_families(root / _REAGENT_FAMILIES_FILE)
        reaction_types = _load_reaction_types(
            root / _REACTION_TYPES_FILE,
            reaction_categories=reaction_categories,
            reactant_types=reactant_types,
            reagent_roles=reagent_roles,
        )
        aliases = _load_aliases(root / _ALIASES_FILE)

        registry = cls(
            root=root,
            manifest=manifest,
            reaction_categories=reaction_categories,
            reaction_types=reaction_types,
            reactant_types=reactant_types,
            reagent_roles=reagent_roles,
            reagent_families=reagent_families,
            aliases=aliases,
        )
        registry._validate_integrity()
        return registry

    # ------------------------------------------------------------------ #
    # Lookups
    # ------------------------------------------------------------------ #
    def get_reaction_category(self, category_id: str) -> Optional[models.ReactionCategory]:
        return self.reaction_categories.get(category_id)

    def iter_reaction_categories(self) -> Iterable[models.ReactionCategory]:
        return self.reaction_categories.values()

    def get_reaction_type(self, reaction_id: str) -> Optional[models.ReactionType]:
        return self.reaction_types.get(reaction_id)

    def iter_reaction_types(self, *, category_id: Optional[str] = None) -> Iterable[models.ReactionType]:
        if category_id is None:
            return self.reaction_types.values()
        return (r for r in self.reaction_types.values() if r.category_id == category_id)

    def get_reactant_type(self, reactant_id: str) -> Optional[models.ReactantType]:
        return self.reactant_types.get(reactant_id)

    def iter_reactant_types(self, *, category: Optional[str] = None) -> Iterable[models.ReactantType]:
        if category is None:
            return self.reactant_types.values()
        return (r for r in self.reactant_types.values() if r.category == category)

    def get_reagent_role(self, role_id: str) -> Optional[models.ReagentRole]:
        return self.reagent_roles.get(role_id)

    def iter_reagent_roles(self) -> Iterable[models.ReagentRole]:
        return self.reagent_roles.values()

    def get_reagent_family(self, family_id: str) -> Optional[models.ReagentFamily]:
        return self.reagent_families.get(family_id)

    def iter_reagent_families(self, *, role_id: Optional[str] = None) -> Iterable[models.ReagentFamily]:
        if role_id is None:
            return self.reagent_families.values()
        return (f for f in self.reagent_families.values() if f.role_id == role_id)

    def resolve_alias(self, alias: str) -> Optional[models.AliasRecord]:
        return self.aliases.get(alias.lower())

    # ------------------------------------------------------------------ #
    # Export helpers
    # ------------------------------------------------------------------ #
    def to_dict(self) -> dict:
        """Return a nested dict of all taxonomy entities (mainly for testing/exports)."""
        return {
            "manifest": asdict(self.manifest),
            "reaction_categories": [asdict(rc) for rc in self.reaction_categories.values()],
            "reaction_types": [asdict(rt) for rt in self.reaction_types.values()],
            "reactant_types": [asdict(rt) for rt in self.reactant_types.values()],
            "reagent_roles": [asdict(rr) for rr in self.reagent_roles.values()],
            "reagent_families": [asdict(rf) for rf in self.reagent_families.values()],
            "aliases": [asdict(alias) for alias in self.aliases.values()],
        }

    # ------------------------------------------------------------------ #
    # Validation
    # ------------------------------------------------------------------ #
    def _validate_integrity(self) -> None:
        errors: List[str] = []

        for rt in self.reaction_types.values():
            if rt.category_id not in self.reaction_categories:
                errors.append(f"ReactionType '{rt.id}' references unknown category '{rt.category_id}'")
            for req in rt.reactants:
                if req.reactant_type_id not in self.reactant_types:
                    errors.append(
                        f"ReactionType '{rt.id}' references unknown reactant '{req.reactant_type_id}'"
                    )
            for role_req in rt.required_roles:
                if role_req.role_id not in self.reagent_roles:
                    errors.append(
                        f"ReactionType '{rt.id}' references unknown reagent role '{role_req.role_id}'"
                    )
                if (
                    role_req.default_family_id
                    and role_req.default_family_id not in self.reagent_families
                ):
                    errors.append(
                        f"ReactionType '{rt.id}' references unknown reagent family '{role_req.default_family_id}'"
                    )

        for family in self.reagent_families.values():
            if family.role_id not in self.reagent_roles:
                errors.append(
                    f"ReagentFamily '{family.id}' references unknown role '{family.role_id}'"
                )

        for alias in self.aliases.values():
            if alias.entity_type == "reaction_type" and alias.entity_id not in self.reaction_types:
                errors.append(f"Alias '{alias.alias}' targets unknown reaction type '{alias.entity_id}'")
            elif alias.entity_type == "reactant_type" and alias.entity_id not in self.reactant_types:
                errors.append(f"Alias '{alias.alias}' targets unknown reactant type '{alias.entity_id}'")
            elif alias.entity_type == "reagent_role" and alias.entity_id not in self.reagent_roles:
                errors.append(f"Alias '{alias.alias}' targets unknown reagent role '{alias.entity_id}'")
            elif alias.entity_type == "reagent_family" and alias.entity_id not in self.reagent_families:
                errors.append(f"Alias '{alias.alias}' targets unknown reagent family '{alias.entity_id}'")

        if errors:
            raise ValueError("Taxonomy integrity check failed:\n- " + "\n- ".join(errors))


# ---------------------------------------------------------------------- #
# JSON loaders
# ---------------------------------------------------------------------- #
def _load_json(path: Path) -> List[dict] | Dict[str, dict]:
    if not path.exists():
        raise FileNotFoundError(f"Taxonomy file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _load_manifest(path: Path) -> models.TaxonomyManifest:
    payload = _load_json(path)
    return models.TaxonomyManifest(
        schema_version=payload.get("schema_version", "0.0.0"),
        taxonomy_version=payload.get("taxonomy_version", "0.0.0"),
        generated_at=payload.get("generated_at", ""),
        source=payload.get("source"),
        notes=payload.get("notes"),
    )


def _load_reaction_categories(path: Path) -> Dict[str, models.ReactionCategory]:
    categories: Dict[str, models.ReactionCategory] = {}
    for entry in _load_json(path):
        obj = models.ReactionCategory(
            id=entry["id"],
            name=entry.get("name", entry["id"]),
            description=entry.get("description"),
        )
        categories[obj.id] = obj
    return categories


def _load_reactant_types(path: Path) -> Dict[str, models.ReactantType]:
    reactants: Dict[str, models.ReactantType] = {}
    for entry in _load_json(path):
        members = [
            models.ReactantTypeMember(
                id=item["id"],
                name=item.get("name", item["id"]),
                smarts=item.get("smarts"),
                aliases=item.get("aliases", []),
                metadata=item.get("metadata") or {},
            )
            for item in entry.get("members", [])
        ]
        obj = models.ReactantType(
            id=entry["id"],
            name=entry.get("name", entry["id"]),
            description=entry.get("description"),
            category=entry.get("category"),
            smarts=entry.get("smarts"),
            members=members,
            aliases=entry.get("aliases", []),
            metadata=entry.get("metadata") or {},
        )
        reactants[obj.id] = obj
    return reactants


def _load_reagent_roles(path: Path) -> Dict[str, models.ReagentRole]:
    roles: Dict[str, models.ReagentRole] = {}
    for entry in _load_json(path):
        obj = models.ReagentRole(
            id=entry["id"],
            name=entry.get("name", entry["id"]),
            priority=int(entry.get("priority", 100)),
            default_family_id=entry.get("default_family_id"),
            description=entry.get("description"),
            metadata=entry.get("metadata") or {},
        )
        roles[obj.id] = obj
    return roles


def _load_reagent_families(path: Path) -> Dict[str, models.ReagentFamily]:
    families: Dict[str, models.ReagentFamily] = {}
    for entry in _load_json(path):
        obj = models.ReagentFamily(
            id=entry["id"],
            role_id=entry["role_id"],
            name=entry.get("name", entry["id"]),
            description=entry.get("description"),
            precedence=entry.get("precedence"),
            include_smarts=entry.get("include_smarts", []),
            exclude_smarts=entry.get("exclude_smarts", []),
            required_props=entry.get("required_props") or {},
            keywords=entry.get("keywords", []),
            examples_pos=entry.get("examples_pos", []),
            examples_neg=entry.get("examples_neg", []),
            notes=entry.get("notes"),
        )
        families[obj.id] = obj
    return families


def _load_reaction_types(
    path: Path,
    *,
    reaction_categories: Dict[str, models.ReactionCategory],
    reactant_types: Dict[str, models.ReactantType],
    reagent_roles: Dict[str, models.ReagentRole],
) -> Dict[str, models.ReactionType]:
    reaction_types: Dict[str, models.ReactionType] = {}
    for entry in _load_json(path):
        # Support multiple formats for reactants:
        # - Simplified: ["ArX*", "RB*"]
        # - Old object: [{"type": "ArX*"}, {"reactant_type_id": "ArB*"}]
        reactant_reqs = []
        for item in entry.get("reactants", []):
            if isinstance(item, str):
                # Ultra-simplified format: just a string
                reactant_reqs.append(
                    models.ReactionReactantRequirement(
                        reactant_type_id=item,
                        stoichiometry=None,
                        notes=None,
                        original_tokens=[],
                    )
                )
            else:
                # Old format: dict with type/reactant_type_id
                type_id = item.get("type") or item.get("reactant_type_id", "")
                tokens = item.get("aliases") or item.get("original_tokens", [])
                reactant_reqs.append(
                    models.ReactionReactantRequirement(
                        reactant_type_id=type_id,
                        stoichiometry=item.get("stoichiometry"),
                        notes=item.get("notes"),
                        original_tokens=tokens,
                    )
                )
        
        role_reqs = [
            models.ReactionRoleRequirement(
                role_id=item["role_id"],
                required=bool(item.get("required", True)),
                default_family_id=item.get("default_family_id"),
                notes=item.get("notes"),
            )
            for item in entry.get("required_roles", [])
        ]
        
        # Load SMARTS - support multiple formats:
        # - Ultra-simplified: top-level "smarts" string
        # - Old format: "patterns" array with objects
        smarts = entry.get("smarts", "")
        patterns = []
        if entry.get("patterns"):
            for p in entry["patterns"]:
                pat_smarts = p.get("smarts") or p.get("smirks", "")
                pattern_id = p.get("id") or p.get("pattern_id", "")
                patterns.append(
                    models.ReactionPattern(
                        id=pattern_id,
                        smarts=pat_smarts,
                        source=p.get("source"),
                        scope=p.get("scope"),
                    )
                )
        elif smarts:
            # Create a single pattern from top-level smarts
            patterns.append(
                models.ReactionPattern(
                    id=f"{entry['id']}_pattern",
                    smarts=smarts,
                    source=None,
                    scope=None,
                )
            )
        
        # Load examples - support multiple formats:
        # - Ultra-simplified: "reference_reactions" list of full reaction SMILES (reactants>>product)
        # - Old format: "examples" or "example_reactants" with objects
        examples = []
        if entry.get("reference_reactions"):
            # Ultra-simplified format: list of full reaction SMILES
            for rxn in entry["reference_reactions"]:
                # Store full reaction SMILES as reactant1, product info extracted if needed
                examples.append(
                    models.ReactionExample(
                        reactant1=rxn,  # Full reaction SMILES
                        reactant2=None,
                        source=None,
                    )
                )
        else:
            # Old format: examples or example_reactants
            example_data = entry.get("examples") or entry.get("example_reactants", [])
            for e in example_data:
                r1 = e.get("reactant1") or e.get("educt1_smiles", "")
                r2 = e.get("reactant2") or e.get("educt2_smiles")
                examples.append(
                    models.ReactionExample(
                        reactant1=r1,
                        reactant2=r2,
                        source=e.get("source"),
                    )
                )
        
        # Support both old format (category_id) and new format (category)
        category_id = entry.get("category") or entry.get("category_id", "")
        
        # Catalysts: new format has top-level, old format has in metadata
        catalysts = entry.get("catalysts", [])
        if not catalysts and entry.get("metadata"):
            catalysts = entry["metadata"].get("typical_catalysts", [])
        
        # Conditions: new format has top-level, old format has in metadata
        conditions = entry.get("conditions")
        if not conditions and entry.get("metadata"):
            conditions = entry["metadata"].get("typical_conditions")
        
        obj = models.ReactionType(
            id=entry["id"],
            category_id=category_id,
            name=entry.get("name", entry["id"]),
            description=entry.get("description"),
            aliases=entry.get("aliases", []),
            reactants=reactant_reqs,
            required_roles=role_reqs,
            patterns=patterns,
            examples=examples,
            catalysts=catalysts,
            conditions=conditions,
            metadata=entry.get("metadata") or {},
            source_ids=entry.get("source_ids", []),
        )
        reaction_types[obj.id] = obj

    return reaction_types


def _load_aliases(path: Path) -> Dict[str, models.AliasRecord]:
    aliases: Dict[str, models.AliasRecord] = {}
    for entry in _load_json(path):
        alias_key = entry["alias"].lower()
        aliases[alias_key] = models.AliasRecord(
            alias=entry["alias"],
            entity_type=entry["entity_type"],
            entity_id=entry["entity_id"],
            source=entry.get("source"),
            notes=entry.get("notes"),
        )
    return aliases

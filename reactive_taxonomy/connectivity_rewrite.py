"""Bounded declarative connectivity-rewrite compiler and executor.

The instruction language is deliberately small. Definitions may select atoms
from already validated reactive sites, change localized bond/H/charge states,
enumerate endpoint permutations, and explicitly authorize product projection.
They cannot name Python callables or execute arbitrary definition content.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from types import MappingProxyType
from typing import Any, Dict, Mapping, Sequence, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .reaction_graph_editing import (
    bond_type,
    capture_join_stereochemistry,
    restore_join_stereochemistry,
    set_total_hydrogens,
)
from .reaction_models import (
    BondChange,
    PredictedStereoChange,
    RewriteOutcome,
    ReactionComponent,
    ReactionSiteReference,
)
from .reaction_site_interfaces import (
    NormalizedSiteInterfaces,
    SITE_INTERFACE_SCHEMA_VERSION,
    normalize_site_assignment,
)


CONNECTIVITY_REWRITE_SCHEMA_VERSION = "3.0"
CONNECTIVITY_REWRITE_INSTRUCTION_SET_VERSION = "1.2"

_PATH = Path(__file__).with_name("definitions") / "connectivity_rewrites.v3.json"
_LEGACY_SELECTOR = re.compile(r"^[a-z][a-z0-9_]*\.[a-z][a-z0-9_]*$")
_NORMALIZED_SELECTOR = re.compile(
    r"^[a-z][a-z0-9_]*\."
    r"(reactive_link|bond_capacity|connection_endpoint)\."
    r"[a-z][a-z0-9_]*$"
)
_BINDING = re.compile(r"^[a-z][a-z0-9_]*$")
_BOND_STATES = {"NONE", "SINGLE", "DOUBLE", "TRIPLE"}
_STATE_RANK = {"NONE": 0, "SINGLE": 1, "DOUBLE": 2, "TRIPLE": 3}
_ALLOWED_TEMPLATES = {
    "change_bond_order",
    "release_and_connect",
    "split_and_distribute",
    "depart_and_unsaturate",
}
_ALLOWED_INSTRUCTIONS = {
    "change_localized_bond_state",
    "change_schema_hydrogen_count",
    "change_observed_formal_charge",
    "declare_product_seed",
    "declare_projection_discardable_attachment",
    "enumerate_endpoint_permutation",
    "set_tetrahedral_outcome",
}
_INTERFACE_PREDICATE_FIELDS = {
    "reactive_link": {
        "available_units",
        "availability",
        "before_order",
        "source_kind",
        "symmetry_class",
    },
    "bond_capacity": {
        "availability",
        "bond_class",
        "current_order",
        "maximum_decrement",
        "maximum_increment",
    },
    "connection_endpoint": {
        "availability",
        "required_bond_capacity_decrement",
        "required_formal_charge_delta",
        "required_hydrogen_delta",
    },
}
_INTERFACE_SELECTOR_MEMBERS = {
    "reactive_link": {"endpoint_a", "endpoint_b", "carrier"},
    "bond_capacity": {"endpoint_a", "endpoint_b"},
    "connection_endpoint": {"atom"},
}


@dataclass(frozen=True)
class CompiledRewriteVariant:
    """One statically validated rewrite variant."""

    variant_id: str
    outcome_id: str
    predicates: Tuple[Mapping[str, Any], ...]
    instructions: Tuple[Mapping[str, Any], ...]


@dataclass(frozen=True)
class CompiledConnectivityRewrite:
    """One compiled, grammar-independent graph rewrite."""

    rewrite_id: str
    template: str
    variants: Tuple[CompiledRewriteVariant, ...]
    schema_version: str = CONNECTIVITY_REWRITE_SCHEMA_VERSION
    instruction_set_version: str = CONNECTIVITY_REWRITE_INSTRUCTION_SET_VERSION
    site_interface_schema_version: str = SITE_INTERFACE_SCHEMA_VERSION


def _freeze(value: Any) -> Any:
    """Recursively freeze definition content held by public compiled contracts."""
    if isinstance(value, Mapping):
        return MappingProxyType(
            {str(key): _freeze(item) for key, item in value.items()}
        )
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    return value


def _validate_selector(selector: object, *, allow_binding: bool = True) -> str:
    value = str(selector or "")
    if allow_binding and value.startswith("$"):
        if not _BINDING.fullmatch(value[1:]):
            raise ValueError(f"Invalid rewrite binding selector: {value}")
        return value
    if _NORMALIZED_SELECTOR.fullmatch(value):
        _, interface, member = value.split(".", 2)
        if member not in _INTERFACE_SELECTOR_MEMBERS[interface]:
            raise ValueError(f"Invalid normalized selector member: {value}")
        return value
    if not _LEGACY_SELECTOR.fullmatch(value):
        raise ValueError(f"Invalid rewrite atom selector: {value}")
    return value


def _compile_variant(raw: Mapping[str, Any]) -> CompiledRewriteVariant:
    variant_id = str(raw.get("id") or "")
    if not _BINDING.fullmatch(variant_id):
        raise ValueError(f"Invalid rewrite variant id: {variant_id}")
    outcome_id = str(raw.get("outcome_id") or "default")
    if not outcome_id:
        raise ValueError(f"Missing outcome id for rewrite variant: {variant_id}")

    raw_predicates = raw.get("when") or ()
    if not isinstance(raw_predicates, (list, tuple)):
        raise ValueError(f"Invalid predicate collection: {variant_id}")
    predicates = tuple(raw_predicates)
    for predicate in predicates:
        if not isinstance(predicate, dict):
            raise ValueError(f"Invalid predicate in rewrite variant: {variant_id}")
        if not _BINDING.fullmatch(str(predicate.get("role") or "")):
            raise ValueError(f"Invalid predicate role in rewrite variant: {variant_id}")
        detail = str(predicate.get("detail") or "")
        interface = str(predicate.get("interface") or "")
        field = str(predicate.get("field") or "")
        if detail:
            if interface or field or not _BINDING.fullmatch(detail):
                raise ValueError(
                    f"Invalid predicate detail in rewrite variant: {variant_id}"
                )
        elif (
            interface not in _INTERFACE_PREDICATE_FIELDS
            or field not in _INTERFACE_PREDICATE_FIELDS[interface]
        ):
            raise ValueError(
                f"Invalid interface predicate in rewrite variant: {variant_id}"
            )
        if predicate.get("operator") not in {"eq", "gte"}:
            raise ValueError(
                f"Unsupported predicate operator in rewrite variant: {variant_id}"
            )
        if "value" not in predicate:
            raise ValueError(f"Missing predicate value in rewrite variant: {variant_id}")

    raw_instructions = raw.get("instructions") or ()
    if not isinstance(raw_instructions, (list, tuple)):
        raise ValueError(f"Invalid instruction collection: {variant_id}")
    instructions = tuple(raw_instructions)
    if not instructions:
        raise ValueError(f"Rewrite variant has no instructions: {variant_id}")
    permutation_count = 0
    seed_count = 0
    changed_pairs: set[tuple[str, str]] = set()
    broken_pairs: set[tuple[str, str]] = set()
    projection_pairs: set[tuple[str, str]] = set()
    for instruction in instructions:
        if not isinstance(instruction, dict):
            raise ValueError(f"Invalid instruction in rewrite variant: {variant_id}")
        operation = str(instruction.get("op") or "")
        if operation not in _ALLOWED_INSTRUCTIONS:
            raise ValueError(
                f"Unsupported rewrite instruction {operation!r}: {variant_id}"
            )
        if operation == "change_localized_bond_state":
            endpoints = tuple(instruction.get("endpoints") or ())
            if len(endpoints) != 2:
                raise ValueError(f"Bond change requires two endpoints: {variant_id}")
            selectors = tuple(_validate_selector(value) for value in endpoints)
            before = str(instruction.get("before") or "").upper()
            after = str(instruction.get("after") or "").upper()
            if (
                before not in _BOND_STATES
                or after not in _BOND_STATES
                or before == after
            ):
                raise ValueError(f"Invalid localized bond states: {variant_id}")
            pair = tuple(sorted(selectors))
            if pair in changed_pairs:
                raise ValueError(f"Bond pair changed more than once: {variant_id}")
            changed_pairs.add(pair)
            if after == "NONE":
                broken_pairs.add(pair)
        elif operation == "change_schema_hydrogen_count":
            _validate_selector(instruction.get("selector"))
            delta = instruction.get("delta")
            if not isinstance(delta, int) or isinstance(delta, bool) or delta == 0:
                raise ValueError(f"Invalid schema hydrogen delta: {variant_id}")
        elif operation == "change_observed_formal_charge":
            _validate_selector(instruction.get("selector"))
            before = instruction.get("before")
            after = instruction.get("after")
            if (
                not isinstance(before, int)
                or isinstance(before, bool)
                or not isinstance(after, int)
                or isinstance(after, bool)
                or before == after
            ):
                raise ValueError(f"Invalid formal-charge transition: {variant_id}")
        elif operation == "declare_product_seed":
            _validate_selector(instruction.get("selector"))
            seed_count += 1
        elif operation == "declare_projection_discardable_attachment":
            retained = _validate_selector(instruction.get("retained"))
            discarded = _validate_selector(instruction.get("discarded"))
            if retained == discarded:
                raise ValueError(f"Invalid projection attachment: {variant_id}")
            projection_pairs.add(tuple(sorted((retained, discarded))))
        elif operation == "set_tetrahedral_outcome":
            _validate_selector(instruction.get("selector"))
            if instruction.get("outcome") not in {
                "invert_if_defined",
                "retain_if_defined",
            }:
                raise ValueError(
                    f"Invalid tetrahedral outcome: {variant_id}"
                )
        elif operation == "enumerate_endpoint_permutation":
            permutation_count += 1
            cases = tuple(instruction.get("cases") or ())
            if len(cases) < 2:
                raise ValueError(f"Endpoint permutation needs two cases: {variant_id}")
            case_ids: set[str] = set()
            binding_keys: set[str] | None = None
            for case in cases:
                if not isinstance(case, dict):
                    raise ValueError(f"Invalid permutation case: {variant_id}")
                outcome = str(case.get("outcome_id") or "")
                if not outcome or outcome in case_ids:
                    raise ValueError(f"Invalid permutation outcome id: {variant_id}")
                case_ids.add(outcome)
                bindings = case.get("bindings") or {}
                if not isinstance(bindings, dict) or not bindings:
                    raise ValueError(f"Missing permutation bindings: {variant_id}")
                keys = set(bindings)
                if binding_keys is None:
                    binding_keys = keys
                elif binding_keys != keys:
                    raise ValueError(
                        f"Inconsistent permutation bindings: {variant_id}"
                    )
                for name, selector in bindings.items():
                    if not _BINDING.fullmatch(str(name)):
                        raise ValueError(f"Invalid binding name: {variant_id}")
                    _validate_selector(selector, allow_binding=False)
    if permutation_count > 1:
        raise ValueError(f"Multiple permutation instructions: {variant_id}")
    if seed_count == 0:
        raise ValueError(f"Rewrite variant has no product seed: {variant_id}")
    if not projection_pairs <= broken_pairs:
        raise ValueError(
            f"Projection must follow a declared broken attachment: {variant_id}"
        )
    return CompiledRewriteVariant(
        variant_id=variant_id,
        outcome_id=outcome_id,
        predicates=tuple(_freeze(value) for value in predicates),
        instructions=tuple(_freeze(value) for value in instructions),
    )


def compile_connectivity_rewrite_definitions(
    payload: Mapping[str, Any],
) -> Tuple[CompiledConnectivityRewrite, ...]:
    """Statically validate and compile one definition payload."""
    if not isinstance(payload, Mapping):
        raise ValueError("Connectivity rewrite definition must be an object")
    if payload.get("schema_version") != CONNECTIVITY_REWRITE_SCHEMA_VERSION:
        raise ValueError("Unsupported connectivity rewrite schema version")
    if (
        payload.get("instruction_set_version")
        != CONNECTIVITY_REWRITE_INSTRUCTION_SET_VERSION
    ):
        raise ValueError("Unsupported connectivity rewrite instruction set")
    if (
        payload.get("site_interface_schema_version")
        != SITE_INTERFACE_SCHEMA_VERSION
    ):
        raise ValueError("Unsupported reactive-site interface schema")
    raw_rewrites = payload.get("rewrites")
    if not isinstance(raw_rewrites, (list, tuple)) or not raw_rewrites:
        raise ValueError("Connectivity rewrite definition has no rewrites")
    rewrites = []
    seen_rewrite_ids: set[str] = set()
    for raw in raw_rewrites:
        if not isinstance(raw, dict):
            raise ValueError("Connectivity rewrite records must be objects")
        legacy_keys = {"grammar_ids", "role_bindings"}.intersection(raw)
        if legacy_keys:
            raise ValueError(
                "Graph operators cannot contain grammar metadata: "
                + ",".join(sorted(legacy_keys))
            )
        rewrite_id = str(raw.get("id") or "")
        if not _BINDING.fullmatch(rewrite_id) or rewrite_id in seen_rewrite_ids:
            raise ValueError(f"Invalid or duplicate connectivity rewrite id: {rewrite_id}")
        seen_rewrite_ids.add(rewrite_id)
        template = str(raw.get("template") or "")
        if template not in _ALLOWED_TEMPLATES:
            raise ValueError(f"Unsupported connectivity rewrite template: {template}")
        raw_variants = raw.get("variants")
        if not isinstance(raw_variants, (list, tuple)) or not raw_variants:
            raise ValueError(f"Invalid connectivity rewrite variants: {rewrite_id}")
        if any(not isinstance(value, dict) for value in raw_variants):
            raise ValueError(f"Invalid connectivity rewrite variants: {rewrite_id}")
        variants = tuple(_compile_variant(value) for value in raw_variants)
        if not variants or len({item.variant_id for item in variants}) != len(variants):
            raise ValueError(f"Invalid connectivity rewrite variants: {rewrite_id}")
        rewrites.append(
            CompiledConnectivityRewrite(
                rewrite_id=rewrite_id,
                template=template,
                variants=variants,
            )
        )
    return tuple(sorted(rewrites, key=lambda item: item.rewrite_id))


@lru_cache(maxsize=1)
def load_connectivity_rewrites() -> Tuple[CompiledConnectivityRewrite, ...]:
    """Load and statically compile the versioned connectivity rewrites."""
    with _PATH.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    return compile_connectivity_rewrite_definitions(payload)


def connectivity_rewrite_by_id(
    rewrite_id: str,
) -> CompiledConnectivityRewrite | None:
    """Return one generic graph operator by its stable rewrite ID."""
    return next(
        (
            rewrite
            for rewrite in load_connectivity_rewrites()
            if rewrite.rewrite_id == rewrite_id
        ),
        None,
    )


def _variant_matches(
    variant: CompiledRewriteVariant,
    assignment: Mapping[str, ReactionSiteReference],
    normalized_assignment: Mapping[str, NormalizedSiteInterfaces],
) -> bool:
    for predicate in variant.predicates:
        role = str(predicate["role"])
        site = assignment.get(role)
        normalized = normalized_assignment.get(role)
        if site is None or normalized is None:
            return False
        if predicate.get("detail"):
            observed = site.details.get(str(predicate["detail"]))
        else:
            interface = str(predicate["interface"])
            collection = {
                "reactive_link": normalized.reactive_links,
                "bond_capacity": normalized.bond_capacities,
                "connection_endpoint": normalized.connection_endpoints,
            }[interface]
            if len(collection) != 1:
                return False
            observed = getattr(collection[0], str(predicate["field"]))
        expected = predicate["value"]
        if predicate["operator"] == "eq" and observed != expected:
            return False
        if predicate["operator"] == "gte":
            try:
                if float(observed) < float(expected):
                    return False
            except (TypeError, ValueError):
                return False
    return True


def _component_by_index(
    components: Sequence[ReactionComponent], component_index: int
) -> ReactionComponent:
    return next(
        component
        for component in components
        if component.component_index == component_index
    )


def _combined_assignment(
    assignment: Mapping[str, ReactionSiteReference],
    components: Sequence[ReactionComponent],
) -> tuple[Any, Dict[int, int]] | None:
    from rdkit import Chem

    component_indices = sorted(
        {site.component_index for site in assignment.values()}
    )
    molecules = []
    offsets: Dict[int, int] = {}
    atom_count = 0
    for component_index in component_indices:
        component = _component_by_index(components, component_index)
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            return None
        offsets[component_index] = atom_count
        atom_count += molecule.GetNumAtoms()
        molecules.append(molecule)
    if not molecules:
        return None
    combined = molecules[0]
    for molecule in molecules[1:]:
        combined = Chem.CombineMols(combined, molecule)
    return combined, offsets


def _resolve_selector(
    selector: object,
    *,
    bindings: Mapping[str, str],
    assignment: Mapping[str, ReactionSiteReference],
    normalized_assignment: Mapping[str, NormalizedSiteInterfaces],
    offsets: Mapping[int, int],
    label_roles: Mapping[str, str],
) -> tuple[int, str] | None:
    value = str(selector)
    if value.startswith("$"):
        value = bindings.get(value[1:], "")
    if _NORMALIZED_SELECTOR.fullmatch(value):
        role, interface, member = value.split(".", 2)
        site = assignment.get(role)
        normalized = normalized_assignment.get(role)
        if site is None or normalized is None:
            return None
        atom_index: int | None = None
        source_atom_role = ""
        if interface == "reactive_link" and len(normalized.reactive_links) == 1:
            link = normalized.reactive_links[0]
            if member in {"endpoint_a", "endpoint_b"}:
                endpoint = getattr(link, member)
                atom_index = endpoint.atom_index
                source_atom_role = endpoint.source_atom_role
            elif member == "carrier":
                endpoint = link.endpoint_b
                atom_index = endpoint.carrier_atom_index
                source_atom_role = endpoint.source_atom_role
        elif (
            interface == "bond_capacity"
            and len(normalized.bond_capacities) == 1
            and member in {"endpoint_a", "endpoint_b"}
        ):
            endpoint = getattr(normalized.bond_capacities[0], member)
            atom_index = endpoint.atom_index
            source_atom_role = endpoint.source_atom_role
        elif (
            interface == "connection_endpoint"
            and len(normalized.connection_endpoints) == 1
            and member == "atom"
        ):
            endpoint = normalized.connection_endpoints[0].endpoint
            atom_index = endpoint.atom_index
            source_atom_role = endpoint.source_atom_role
        if (
            atom_index is None
            or not source_atom_role
            or site.component_index not in offsets
        ):
            return None
        return (
            offsets[site.component_index] + atom_index,
            f"{label_roles.get(role, role)}.{source_atom_role}",
        )
    if not _LEGACY_SELECTOR.fullmatch(value):
        return None
    role, atom_role = value.split(".", 1)
    site = assignment.get(role)
    if site is None:
        return None
    indices = site.atom_roles.get(atom_role) or ()
    if len(indices) != 1 or site.component_index not in offsets:
        return None
    return (
        offsets[site.component_index] + int(indices[0]),
        f"{label_roles.get(role, role)}.{atom_role}",
    )


def _bond_state(molecule: Any, atom_1: int, atom_2: int) -> str | None:
    bond = molecule.GetBondBetweenAtoms(atom_1, atom_2)
    if bond is None:
        return "NONE"
    value = float(bond.GetBondTypeAsDouble())
    states = {1.0: "SINGLE", 2.0: "DOUBLE", 3.0: "TRIPLE"}
    return states.get(value)


def _set_bond_state(
    molecule: Any,
    atom_1: int,
    atom_2: int,
    *,
    before: str,
    after: str,
) -> bool:
    from rdkit import Chem

    bond = molecule.GetBondBetweenAtoms(atom_1, atom_2)
    if before == "NONE":
        if bond is not None or after == "NONE":
            return False
        molecule.AddBond(atom_1, atom_2, bond_type(after))
        return True
    if bond is None:
        return False
    if after == "NONE":
        molecule.RemoveBond(atom_1, atom_2)
        return True
    bond.SetBondType(bond_type(after))
    bond.SetStereo(Chem.BondStereo.STEREONONE)
    bond.SetBondDir(Chem.BondDir.NONE)
    if after in {"DOUBLE", "TRIPLE"}:
        molecule.GetAtomWithIdx(atom_1).SetChiralTag(
            Chem.ChiralType.CHI_UNSPECIFIED
        )
        molecule.GetAtomWithIdx(atom_2).SetChiralTag(
            Chem.ChiralType.CHI_UNSPECIFIED
        )
    return True


def _bond_change(
    instruction: Mapping[str, Any],
    resolved: Sequence[tuple[int, str]],
) -> BondChange:
    before = str(instruction["before"]).upper()
    after = str(instruction["after"]).upper()
    if before == "NONE":
        change_type = "formed"
        old_order, new_order = None, after
    elif after == "NONE":
        change_type = "broken"
        old_order, new_order = before, None
    else:
        change_type = "order_changed"
        old_order, new_order = before, after
    return BondChange(
        change_type,
        resolved[0][1],
        resolved[1][1],
        old_order,
        new_order,
        "connectivity_rewrite",
    )


def _permutation_cases(
    variant: CompiledRewriteVariant,
) -> Tuple[tuple[str, Mapping[str, str]], ...]:
    instruction = next(
        (
            value
            for value in variant.instructions
            if value["op"] == "enumerate_endpoint_permutation"
        ),
        None,
    )
    if instruction is None:
        return ((variant.outcome_id, {}),)
    return tuple(
        (
            str(case["outcome_id"]),
            {str(key): str(value) for key, value in case["bindings"].items()},
        )
        for case in instruction["cases"]
    )


def _atom_stereo_descriptor(atom: Any) -> str | None:
    """Return an assigned CIP code, or a defined tetrahedral tag fallback."""
    from rdkit import Chem

    if atom.HasProp("_CIPCode"):
        return str(atom.GetProp("_CIPCode"))
    tag = atom.GetChiralTag()
    if tag == Chem.ChiralType.CHI_UNSPECIFIED:
        return None
    return str(tag).replace("CHI_", "")


def _permutation_is_odd(
    logical_order: Sequence[int],
    actual_order: Sequence[int],
) -> bool | None:
    """Return parity mapping a logical neighbor order to RDKit's actual order."""
    if len(logical_order) != len(actual_order) or set(logical_order) != set(
        actual_order
    ):
        return None
    positions = [actual_order.index(value) for value in logical_order]
    inversions = sum(
        positions[left] > positions[right]
        for left in range(len(positions))
        for right in range(left + 1, len(positions))
    )
    return bool(inversions % 2)


def _execute_variant_case(
    variant: CompiledRewriteVariant,
    *,
    outcome_id: str,
    bindings: Mapping[str, str],
    assignment: Mapping[str, ReactionSiteReference],
    normalized_assignment: Mapping[str, NormalizedSiteInterfaces],
    components: Sequence[ReactionComponent],
    label_roles: Mapping[str, str],
) -> RewriteOutcome | None:
    from rdkit import Chem

    combined_result = _combined_assignment(assignment, components)
    if combined_result is None:
        return None
    source, offsets = combined_result
    bond_operations = []
    hydrogen_deltas: Dict[int, int] = {}
    charge_operations = []
    projection_pairs = []
    tetrahedral_directives = []
    seeds: set[int] = set()
    changes = []
    for instruction in variant.instructions:
        operation = str(instruction["op"])
        if operation == "enumerate_endpoint_permutation":
            continue
        if operation == "change_localized_bond_state":
            resolved = tuple(
                _resolve_selector(
                    selector,
                    bindings=bindings,
                    assignment=assignment,
                    normalized_assignment=normalized_assignment,
                    offsets=offsets,
                    label_roles=label_roles,
                )
                for selector in instruction["endpoints"]
            )
            if any(value is None for value in resolved):
                return None
            typed_resolved = tuple(value for value in resolved if value is not None)
            atom_1, atom_2 = typed_resolved[0][0], typed_resolved[1][0]
            before = str(instruction["before"]).upper()
            after = str(instruction["after"]).upper()
            if atom_1 == atom_2 or _bond_state(source, atom_1, atom_2) != before:
                return None
            bond_operations.append(
                (instruction, typed_resolved, atom_1, atom_2, before, after)
            )
            changes.append(_bond_change(instruction, typed_resolved))
        elif operation == "change_schema_hydrogen_count":
            resolved = _resolve_selector(
                instruction["selector"],
                bindings=bindings,
                assignment=assignment,
                normalized_assignment=normalized_assignment,
                offsets=offsets,
                label_roles=label_roles,
            )
            if resolved is None:
                return None
            atom_index, label = resolved
            delta = int(instruction["delta"])
            hydrogen_deltas[atom_index] = hydrogen_deltas.get(atom_index, 0) + delta
            changes.append(
                BondChange(
                    "hydrogen_change",
                    label,
                    None,
                    "SINGLE" if delta < 0 else None,
                    "SINGLE" if delta > 0 else None,
                    "connectivity_rewrite",
                )
            )
        elif operation == "change_observed_formal_charge":
            resolved = _resolve_selector(
                instruction["selector"],
                bindings=bindings,
                assignment=assignment,
                normalized_assignment=normalized_assignment,
                offsets=offsets,
                label_roles=label_roles,
            )
            if resolved is None:
                return None
            atom_index, _ = resolved
            before = int(instruction["before"])
            after = int(instruction["after"])
            if source.GetAtomWithIdx(atom_index).GetFormalCharge() != before:
                return None
            charge_operations.append((atom_index, after))
        elif operation == "declare_product_seed":
            resolved = _resolve_selector(
                instruction["selector"],
                bindings=bindings,
                assignment=assignment,
                normalized_assignment=normalized_assignment,
                offsets=offsets,
                label_roles=label_roles,
            )
            if resolved is None:
                return None
            seeds.add(resolved[0])
        elif operation == "declare_projection_discardable_attachment":
            retained = _resolve_selector(
                instruction["retained"],
                bindings=bindings,
                assignment=assignment,
                normalized_assignment=normalized_assignment,
                offsets=offsets,
                label_roles=label_roles,
            )
            discarded = _resolve_selector(
                instruction["discarded"],
                bindings=bindings,
                assignment=assignment,
                normalized_assignment=normalized_assignment,
                offsets=offsets,
                label_roles=label_roles,
            )
            if retained is None or discarded is None:
                return None
            projection_pairs.append((retained[0], discarded[0]))
        elif operation == "set_tetrahedral_outcome":
            resolved = _resolve_selector(
                instruction["selector"],
                bindings=bindings,
                assignment=assignment,
                normalized_assignment=normalized_assignment,
                offsets=offsets,
                label_roles=label_roles,
            )
            if resolved is None:
                return None
            tetrahedral_directives.append(
                (
                    resolved[0],
                    resolved[1],
                    str(instruction["outcome"]),
                )
            )

    def failed_outcome() -> RewriteOutcome:
        return RewriteOutcome(
            outcome_id=outcome_id,
            predicted_product_smiles=None,
            predicted_bond_changes=tuple(changes),
        )

    for atom_index, delta in hydrogen_deltas.items():
        if (
            int(source.GetAtomWithIdx(atom_index).GetTotalNumHs(includeNeighbors=True))
            + delta
            < 0
        ):
            return failed_outcome()

    rw = Chem.RWMol(source)
    for _, _, atom_1, atom_2, before, after in bond_operations:
        if _STATE_RANK[after] >= _STATE_RANK[before]:
            continue
        if not _set_bond_state(
            rw, atom_1, atom_2, before=before, after=after
        ):
            return failed_outcome()

    removable: set[int] = set()
    negative_product = rw.GetMol()
    fragments = tuple(
        set(int(index) for index in fragment)
        for fragment in Chem.GetMolFrags(
            negative_product, asMols=False, sanitizeFrags=False
        )
    )
    for retained, discarded in projection_pairs:
        fragment = next(
            (value for value in fragments if discarded in value),
            None,
        )
        if fragment is None or retained in fragment:
            return failed_outcome()
        removable.update(fragment)
    if removable.intersection(seeds):
        return failed_outcome()

    formed_pairs = []
    for _, _, atom_1, atom_2, before, after in bond_operations:
        if _STATE_RANK[after] <= _STATE_RANK[before]:
            continue
        if atom_1 in removable or atom_2 in removable:
            return failed_outcome()
        if not _set_bond_state(
            rw, atom_1, atom_2, before=before, after=after
        ):
            return failed_outcome()
        if before == "NONE":
            formed_pairs.append((atom_1, atom_2))

    Chem.AssignStereochemistry(source, cleanIt=True, force=True)
    tetrahedral_states = []
    for center_index, center_label, requested_outcome in tetrahedral_directives:
        center = source.GetAtomWithIdx(center_index)
        if center.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
            continue
        removed_neighbors = [
            neighbor.GetIdx()
            for neighbor in center.GetNeighbors()
            if neighbor.GetIdx() in removable
        ]
        added_neighbors = [
            atom_2 if atom_1 == center_index else atom_1
            for atom_1, atom_2 in formed_pairs
            if center_index in {atom_1, atom_2}
            and (atom_2 if atom_1 == center_index else atom_1) not in removable
        ]
        if len(removed_neighbors) != 1 or len(added_neighbors) != 1:
            return failed_outcome()
        removed_neighbor = removed_neighbors[0]
        added_neighbor = added_neighbors[0]
        source_order = tuple(
            int(neighbor.GetIdx()) for neighbor in center.GetNeighbors()
        )
        logical_order = tuple(
            added_neighbor if index == removed_neighbor else index
            for index in source_order
        )
        tetrahedral_states.append(
            (
                center_index,
                center_label,
                requested_outcome,
                center.GetChiralTag(),
                _atom_stereo_descriptor(center),
                logical_order,
            )
        )

    for atom_index, delta in hydrogen_deltas.items():
        if atom_index in removable or not set_total_hydrogens(
            source, rw, atom_index, delta
        ):
            return failed_outcome()
    for atom_index, charge in charge_operations:
        if atom_index in removable:
            return failed_outcome()
        rw.GetAtomWithIdx(atom_index).SetFormalCharge(charge)

    captured = ()
    if len(formed_pairs) == 1 and removable:
        captured = capture_join_stereochemistry(
            source,
            removed_indices=removable,
            join_indices=formed_pairs[0],
        )
    removed_indices = sorted(removable, reverse=True)
    for atom_index in removed_indices:
        rw.RemoveAtom(atom_index)

    removed_ascending = tuple(sorted(removable))

    def shifted(index: int) -> int:
        return index - sum(removed < index for removed in removed_ascending)

    product = rw.GetMol()
    try:
        product.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(product)
        if captured:
            restore_join_stereochemistry(
                product,
                captured=captured,
                shifted=shifted,
            )
        for (
            center_index,
            _,
            requested_outcome,
            source_tag,
            _,
            logical_order,
        ) in tetrahedral_states:
            target_center = product.GetAtomWithIdx(shifted(center_index))
            shifted_logical_order = tuple(shifted(index) for index in logical_order)
            actual_order = tuple(
                int(neighbor.GetIdx()) for neighbor in target_center.GetNeighbors()
            )
            odd = _permutation_is_odd(shifted_logical_order, actual_order)
            if odd is None:
                return failed_outcome()
            invert = requested_outcome == "invert_if_defined"
            flip_tag = odd != invert
            target_tag = source_tag
            if flip_tag:
                target_tag = (
                    Chem.ChiralType.CHI_TETRAHEDRAL_CCW
                    if source_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW
                    else Chem.ChiralType.CHI_TETRAHEDRAL_CW
                )
            target_center.SetChiralTag(target_tag)
        Chem.AssignStereochemistry(product, cleanIt=True, force=True)
    except Exception:
        return failed_outcome()

    shifted_seeds = {shifted(index) for index in seeds}
    retained_fragments = Chem.GetMolFrags(
        product, asMols=False, sanitizeFrags=False
    )
    if any(not set(fragment).intersection(shifted_seeds) for fragment in retained_fragments):
        return failed_outcome()
    smiles = Chem.MolToSmiles(product, canonical=True, isomericSmiles=True)
    stereo_changes = []
    for (
        center_index,
        center_label,
        requested_outcome,
        _,
        old_descriptor,
        _,
    ) in tetrahedral_states:
        new_descriptor = _atom_stereo_descriptor(
            product.GetAtomWithIdx(shifted(center_index))
        )
        stereo_changes.append(
            PredictedStereoChange(
                stereo_type="atom",
                atom_1_role=center_label,
                atom_2_role=None,
                old_descriptor=old_descriptor,
                new_descriptor=new_descriptor,
                change_type=(
                    "inverted"
                    if requested_outcome == "invert_if_defined"
                    else "retained"
                ),
                evidence=f"connectivity_rewrite:{requested_outcome}",
            )
        )
    return RewriteOutcome(
        outcome_id=outcome_id,
        predicted_product_smiles=smiles,
        predicted_bond_changes=tuple(changes),
        predicted_stereo_changes=tuple(stereo_changes),
    )


def apply_reaction_operator(
    operator: CompiledConnectivityRewrite | str,
    assignment: Mapping[str, ReactionSiteReference],
    components: Sequence[ReactionComponent],
    *,
    output_role_labels: Mapping[str, str] | None = None,
) -> Tuple[RewriteOutcome, ...]:
    """Execute one graph operator without loading grammar or family metadata."""
    rewrite = (
        connectivity_rewrite_by_id(operator)
        if isinstance(operator, str)
        else operator
    )
    if rewrite is None:
        return ()
    try:
        normalized_assignment = normalize_site_assignment(
            assignment,
            components,
        )
    except (KeyError, StopIteration, ValueError):
        return ()
    outcomes = []
    seen_products: set[str] = set()
    for variant in rewrite.variants:
        if not _variant_matches(variant, assignment, normalized_assignment):
            continue
        for outcome_id, bindings in _permutation_cases(variant):
            outcome = _execute_variant_case(
                variant,
                outcome_id=outcome_id,
                bindings=bindings,
                assignment=assignment,
                normalized_assignment=normalized_assignment,
                components=components,
                label_roles=dict(output_role_labels or {}),
            )
            if outcome is None:
                continue
            if outcome.predicted_product_smiles is not None:
                if outcome.predicted_product_smiles in seen_products:
                    continue
                seen_products.add(outcome.predicted_product_smiles)
            outcomes.append(outcome)
    return tuple(sorted(outcomes, key=lambda outcome: outcome.outcome_id))


__all__ = [
    "CONNECTIVITY_REWRITE_INSTRUCTION_SET_VERSION",
    "CONNECTIVITY_REWRITE_SCHEMA_VERSION",
    "CompiledConnectivityRewrite",
    "CompiledRewriteVariant",
    "apply_reaction_operator",
    "compile_connectivity_rewrite_definitions",
    "connectivity_rewrite_by_id",
    "load_connectivity_rewrites",
]

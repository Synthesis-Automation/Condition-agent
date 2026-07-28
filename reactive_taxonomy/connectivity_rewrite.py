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
    OperatorOutcome,
    ReactionComponent,
    ReactionSiteReference,
)


CONNECTIVITY_REWRITE_SCHEMA_VERSION = "1.0"
CONNECTIVITY_REWRITE_INSTRUCTION_SET_VERSION = "1.0"

_PATH = Path(__file__).with_name("definitions") / "connectivity_rewrites.v1.json"
_SELECTOR = re.compile(r"^[a-z][a-z0-9_]*\.[a-z][a-z0-9_]*$")
_BINDING = re.compile(r"^[a-z][a-z0-9_]*$")
_BOND_STATES = {"NONE", "SINGLE", "DOUBLE", "TRIPLE"}
_STATE_RANK = {"NONE": 0, "SINGLE": 1, "DOUBLE": 2, "TRIPLE": 3}
_ALLOWED_TEMPLATES = {
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
    """One compiled rewrite shared by one or more reaction grammars."""

    rewrite_id: str
    template: str
    grammar_ids: Tuple[str, ...]
    variants: Tuple[CompiledRewriteVariant, ...]
    schema_version: str = CONNECTIVITY_REWRITE_SCHEMA_VERSION
    instruction_set_version: str = CONNECTIVITY_REWRITE_INSTRUCTION_SET_VERSION


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
    if not _SELECTOR.fullmatch(value):
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
        if not _BINDING.fullmatch(str(predicate.get("detail") or "")):
            raise ValueError(
                f"Invalid predicate detail in rewrite variant: {variant_id}"
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
    raw_rewrites = payload.get("rewrites")
    if not isinstance(raw_rewrites, (list, tuple)) or not raw_rewrites:
        raise ValueError("Connectivity rewrite definition has no rewrites")
    rewrites = []
    seen_rewrite_ids: set[str] = set()
    seen_grammar_ids: set[str] = set()
    for raw in raw_rewrites:
        if not isinstance(raw, dict):
            raise ValueError("Connectivity rewrite records must be objects")
        rewrite_id = str(raw.get("id") or "")
        if not _BINDING.fullmatch(rewrite_id) or rewrite_id in seen_rewrite_ids:
            raise ValueError(f"Invalid or duplicate connectivity rewrite id: {rewrite_id}")
        seen_rewrite_ids.add(rewrite_id)
        template = str(raw.get("template") or "")
        if template not in _ALLOWED_TEMPLATES:
            raise ValueError(f"Unsupported connectivity rewrite template: {template}")
        grammar_ids = tuple(str(value) for value in raw.get("grammar_ids") or ())
        if (
            not grammar_ids
            or len(grammar_ids) != len(set(grammar_ids))
            or any(not _BINDING.fullmatch(value) for value in grammar_ids)
            or seen_grammar_ids.intersection(grammar_ids)
        ):
            raise ValueError(f"Invalid connectivity rewrite grammar ids: {rewrite_id}")
        seen_grammar_ids.update(grammar_ids)
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
                grammar_ids=grammar_ids,
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


def connectivity_rewrite_for_grammar(
    grammar_id: str,
) -> CompiledConnectivityRewrite | None:
    """Return the registered rewrite for a grammar, if Phase 2 migrated it."""
    return next(
        (
            rewrite
            for rewrite in load_connectivity_rewrites()
            if grammar_id in rewrite.grammar_ids
        ),
        None,
    )


def _variant_matches(
    variant: CompiledRewriteVariant,
    assignment: Mapping[str, ReactionSiteReference],
) -> bool:
    for predicate in variant.predicates:
        site = assignment.get(str(predicate["role"]))
        if site is None:
            return False
        observed = site.details.get(str(predicate["detail"]))
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
    offsets: Mapping[int, int],
) -> tuple[int, str] | None:
    value = str(selector)
    if value.startswith("$"):
        value = bindings.get(value[1:], "")
    if not _SELECTOR.fullmatch(value):
        return None
    role, atom_role = value.split(".", 1)
    site = assignment.get(role)
    if site is None:
        return None
    indices = site.atom_roles.get(atom_role) or ()
    if len(indices) != 1 or site.component_index not in offsets:
        return None
    return offsets[site.component_index] + int(indices[0]), value


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
        "grammar_operator",
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


def _execute_variant_case(
    variant: CompiledRewriteVariant,
    *,
    outcome_id: str,
    bindings: Mapping[str, str],
    assignment: Mapping[str, ReactionSiteReference],
    components: Sequence[ReactionComponent],
) -> OperatorOutcome | None:
    from rdkit import Chem

    combined_result = _combined_assignment(assignment, components)
    if combined_result is None:
        return None
    source, offsets = combined_result
    bond_operations = []
    hydrogen_deltas: Dict[int, int] = {}
    charge_operations = []
    projection_pairs = []
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
                    offsets=offsets,
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
                offsets=offsets,
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
                    "grammar_operator",
                )
            )
        elif operation == "change_observed_formal_charge":
            resolved = _resolve_selector(
                instruction["selector"],
                bindings=bindings,
                assignment=assignment,
                offsets=offsets,
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
                offsets=offsets,
            )
            if resolved is None:
                return None
            seeds.add(resolved[0])
        elif operation == "declare_projection_discardable_attachment":
            retained = _resolve_selector(
                instruction["retained"],
                bindings=bindings,
                assignment=assignment,
                offsets=offsets,
            )
            discarded = _resolve_selector(
                instruction["discarded"],
                bindings=bindings,
                assignment=assignment,
                offsets=offsets,
            )
            if retained is None or discarded is None:
                return None
            projection_pairs.append((retained[0], discarded[0]))

    for atom_index, delta in hydrogen_deltas.items():
        if (
            int(source.GetAtomWithIdx(atom_index).GetTotalNumHs(includeNeighbors=True))
            + delta
            < 0
        ):
            return None

    rw = Chem.RWMol(source)
    for _, _, atom_1, atom_2, before, after in bond_operations:
        if _STATE_RANK[after] >= _STATE_RANK[before]:
            continue
        if not _set_bond_state(
            rw, atom_1, atom_2, before=before, after=after
        ):
            return None

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
            return None
        removable.update(fragment)
    if removable.intersection(seeds):
        return None

    formed_pairs = []
    for _, _, atom_1, atom_2, before, after in bond_operations:
        if _STATE_RANK[after] <= _STATE_RANK[before]:
            continue
        if atom_1 in removable or atom_2 in removable:
            return None
        if not _set_bond_state(
            rw, atom_1, atom_2, before=before, after=after
        ):
            return None
        if before == "NONE":
            formed_pairs.append((atom_1, atom_2))

    for atom_index, delta in hydrogen_deltas.items():
        if atom_index in removable or not set_total_hydrogens(
            source, rw, atom_index, delta
        ):
            return None
    for atom_index, charge in charge_operations:
        if atom_index in removable:
            return None
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
        else:
            Chem.AssignStereochemistry(product, cleanIt=True, force=True)
    except Exception:
        return None

    shifted_seeds = {shifted(index) for index in seeds}
    retained_fragments = Chem.GetMolFrags(
        product, asMols=False, sanitizeFrags=False
    )
    if any(not set(fragment).intersection(shifted_seeds) for fragment in retained_fragments):
        return None
    smiles = Chem.MolToSmiles(product, canonical=True, isomericSmiles=True)
    return OperatorOutcome(
        outcome_id=outcome_id,
        predicted_product_smiles=smiles,
        predicted_bond_changes=tuple(changes),
    )


def apply_connectivity_rewrite(
    grammar: Mapping[str, Any] | str,
    assignment: Mapping[str, ReactionSiteReference],
    components: Sequence[ReactionComponent],
) -> Tuple[OperatorOutcome, ...]:
    """Execute a registered connectivity rewrite for one grammar assignment.

    Unsupported grammars and chemistry-invalid cases return no outcomes. This
    keeps absence distinct from a failed product represented as an outcome.
    """
    grammar_id = str(grammar if isinstance(grammar, str) else grammar.get("id") or "")
    rewrite = connectivity_rewrite_for_grammar(grammar_id)
    if rewrite is None:
        return ()
    outcomes = []
    seen_products: set[str] = set()
    for variant in rewrite.variants:
        if not _variant_matches(variant, assignment):
            continue
        for outcome_id, bindings in _permutation_cases(variant):
            outcome = _execute_variant_case(
                variant,
                outcome_id=outcome_id,
                bindings=bindings,
                assignment=assignment,
                components=components,
            )
            if (
                outcome is None
                or outcome.predicted_product_smiles is None
                or outcome.predicted_product_smiles in seen_products
            ):
                continue
            seen_products.add(outcome.predicted_product_smiles)
            outcomes.append(outcome)
    return tuple(sorted(outcomes, key=lambda outcome: outcome.outcome_id))


__all__ = [
    "CONNECTIVITY_REWRITE_INSTRUCTION_SET_VERSION",
    "CONNECTIVITY_REWRITE_SCHEMA_VERSION",
    "CompiledConnectivityRewrite",
    "CompiledRewriteVariant",
    "apply_connectivity_rewrite",
    "compile_connectivity_rewrite_definitions",
    "connectivity_rewrite_for_grammar",
    "load_connectivity_rewrites",
]

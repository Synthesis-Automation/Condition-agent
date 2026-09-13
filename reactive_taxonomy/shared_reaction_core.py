"""Qualified, versioned graph projections shared by precedent search channels.

The first generalization operator supports one intermolecular single-bond join
with graph-verified departing or hydrogen ports. Other transformations retain
their observed-local projection and explicitly abstain from generalization.
No hypothetical molecules or atom correspondence are generated here.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping

from rdkit import Chem, rdBase

from .chemistry.rdkit_utils import parse_smiles
from .reaction_parser import parse_reaction_smiles
from .reaction_models import (
    REACTION_CORE_PROJECTION_ALGORITHM_VERSION,
    REACTION_CORE_PROJECTION_SCHEMA_VERSION,
    REACTION_SIGNATURE_SCHEMA_VERSION,
)

SCHEMA_VERSION = "1.0"
ALGORITHM_VERSION = "shared_reaction_core.v1"
LEVELS = ("observed_local", "retained_local", "retained_typed")


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def _digest(namespace: str, value: Any) -> str:
    return namespace + ":" + hashlib.sha256(_json(value).encode()).hexdigest()


@lru_cache(maxsize=1)
def _definition_text() -> str:
    path = Path(__file__).with_name("definitions") / "shared_reaction_core.v1.json"
    return path.read_text(encoding="utf-8")


def load_shared_core_rules() -> dict[str, Any]:
    """Load validated definitions as a fresh value, without mutable globals."""
    rules = json.loads(_definition_text())
    if (
        rules.get("definition_id") != ALGORITHM_VERSION
        or rules.get("schema_version") != SCHEMA_VERSION
        or rules.get("levels") != dict(zip(LEVELS, (1, 1, 0)))
        or rules.get("generalization_operator") != "single_join_with_departing_ports"
        or rules.get("status") != "experimental_pending_independent_review"
    ):
        raise ValueError("invalid shared reaction core definition")
    for name in ("departing_root_elements", "hydrogen_port_elements"):
        values = rules.get(name, [])
        if not values or len(set(values)) != len(values):
            raise ValueError(f"invalid {name}")
        for symbol in values:
            if Chem.GetPeriodicTable().GetAtomicNumber(symbol) < 1:
                raise ValueError(f"invalid element: {symbol}")
    if rules.get("hydrogen_port_carbon_hybridizations") != ["SP"] or set(
        rules.get("protected_state_changes", [])
    ) != {
        "formal_charge",
        "radical",
        "isotope",
        "aromaticity",
        "hybridization",
        "atom_stereochemistry",
        "bond_stereochemistry",
    }:
        raise ValueError("invalid protected-state or hydrogen-port definition")
    for name in ("maximum_fragment_atoms", "maximum_core_atoms"):
        if type(rules.get(name)) is not int or not 1 <= rules[name] <= 256:
            raise ValueError(f"invalid {name}")
    return rules


@lru_cache(maxsize=1)
def shared_core_definition_hash() -> str:
    """Return the chemistry and algorithm identity required by artifacts."""
    return _digest(
        "SCD1",
        [
            ALGORITHM_VERSION,
            SCHEMA_VERSION,
            load_shared_core_rules(),
            REACTION_CORE_PROJECTION_ALGORITHM_VERSION,
            REACTION_CORE_PROJECTION_SCHEMA_VERSION,
            REACTION_SIGNATURE_SCHEMA_VERSION,
            rdBase.rdkitVersion,
        ],
    )


@dataclass(frozen=True)
class SharedCoreLevel:
    """One canonical projection with original evidence references."""

    level: str
    key: str
    graph_payload: str
    product_side_key: str
    generalized_edit_indices: tuple[int, ...] = ()


@dataclass(frozen=True)
class SharedReactionCore:
    """Derived observation views; original signature/core remain authoritative."""

    definition_hash: str
    observation_core_id: str
    reaction_identity: str
    input_identity: str
    product_identity: str
    levels: tuple[SharedCoreLevel, ...]
    realization_key: str
    realization_details: tuple[str, ...]
    evidence_status: str
    warnings: tuple[str, ...]
    unavailable_reasons: tuple[str, ...]
    schema_version: str = SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Serialize without dropping eligibility or provenance."""
        return asdict(self)

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> SharedReactionCore:
        """Validate artifact versions and restore immutable sequence fields."""
        if (
            value.get("schema_version") != SCHEMA_VERSION
            or value.get("definition_hash") != shared_core_definition_hash()
        ):
            raise ValueError("incompatible shared core projection; rebuild artifact")
        data = dict(value)
        data["levels"] = tuple(
            SharedCoreLevel(
                **{
                    **level,
                    "generalized_edit_indices": tuple(
                        level.get("generalized_edit_indices", ())
                    ),
                },
            )
            for level in value["levels"]
        )
        names = tuple(level.level for level in data["levels"])
        if names != LEVELS[: len(names)]:
            raise ValueError("invalid shared core level ordering")
        for name in ("realization_details", "warnings", "unavailable_reasons"):
            data[name] = tuple(data[name])
        return cls(**data)


@dataclass(frozen=True)
class SharedCoreComparison:
    """Channel-independent qualification against the original query."""

    eligible: bool
    level: str | None
    relation: str
    differences: tuple[str, ...]
    reasons: tuple[str, ...]
    query_core_id: str
    precedent_core_id: str
    definition_hash: str


def compare_reaction_cores(
    query: SharedReactionCore,
    precedent: SharedReactionCore,
    *,
    maximum_level: str = "retained_typed",
) -> SharedCoreComparison:
    """Qualify matching graphs; analogues never silently replace fixed inputs."""
    if maximum_level not in LEVELS:
        raise ValueError("unknown shared core level")
    common = dict(
        query_core_id=query.observation_core_id,
        precedent_core_id=precedent.observation_core_id,
        definition_hash=query.definition_hash,
    )
    if query.definition_hash != precedent.definition_hash:
        return SharedCoreComparison(
            False, None, "incompatible", (), ("DEFINITION_MISMATCH",), **common
        )
    other = {item.level: item for item in precedent.levels}
    for level in query.levels:
        if LEVELS.index(level.level) > LEVELS.index(maximum_level):
            break
        match = other.get(level.level)
        if (
            match is None
            or level.key != match.key
            or level.graph_payload != match.graph_payload
        ):
            continue
        differences = []
        if query.realization_key != precedent.realization_key:
            differences.append("Departing/source attachments differ")
        if query.input_identity != precedent.input_identity:
            differences.append("Supplied molecular inputs differ from the precedent")
        if level.level == "retained_typed":
            differences.append("Local environment detail was generalized")
        same = query.reaction_identity == precedent.reaction_identity
        return SharedCoreComparison(
            True,
            "whole_reaction" if same else level.level,
            "same_setup" if same else "analogue_evidence",
            tuple(differences),
            tuple(dict.fromkeys((*query.warnings, *precedent.warnings))),
            **common,
        )
    reason = (
        "DIFFERENT_TRANSFORMATION_SAME_PRODUCT"
        if query.product_identity
        and query.product_identity == precedent.product_identity
        else "NO_QUALIFIED_SHARED_CORE"
    )
    return SharedCoreComparison(False, None, "incompatible", (), (reason,), **common)


def _ref(state: Mapping[str, Any]) -> tuple[int, int]:
    return int(state["component_index"]), int(state["atom_index"])


def _atom_label(atom: Chem.Atom) -> list[Any]:
    return [
        atom.GetSymbol(),
        atom.GetFormalCharge(),
        atom.GetIsAromatic(),
        str(atom.GetHybridization()),
        atom.GetIsotope(),
        atom.GetNumRadicalElectrons(),
        atom.GetTotalNumHs(),
        atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else "",
    ]


def _canonical_graph(labels: list[Any], edges: list[tuple[int, int, Any]]) -> str:
    """Canonicalize a colored graph through an atom/edge incidence graph.

    Chemical labels become a sorted vocabulary of isotope colors on dummy
    atoms. Edge nodes preserve bond labels. The vocabulary is included so
    colors from different graphs cannot accidentally compare equal.
    """
    encoded = [_json(["atom", label]) for label in labels]
    edge_labels = [_json(["edge", edge[2]]) for edge in edges]
    vocabulary = sorted(set(encoded + edge_labels))
    colors = {label: i + 1 for i, label in enumerate(vocabulary)}
    graph = Chem.RWMol()
    for label in encoded + edge_labels:
        atom = Chem.Atom(0)
        atom.SetIsotope(colors[label])
        graph.AddAtom(atom)
    for i, (left, right, _) in enumerate(edges, len(labels)):
        graph.AddBond(left, i, Chem.BondType.SINGLE)
        graph.AddBond(i, right, Chem.BondType.SINGLE)
    return _json([vocabulary, Chem.MolToSmiles(graph, canonical=True)])


def _local_graph(
    molecule: Chem.Mol,
    centers: Mapping[int, Any],
    radius: int,
    limit: int,
) -> str:
    selected = set(centers)
    for _ in range(radius):
        selected |= {
            n.GetIdx()
            for i in selected
            for n in molecule.GetAtomWithIdx(i).GetNeighbors()
        }
    # Close complete conjugated multiple-bond units, but not aromatic scaffolds.
    while True:
        expanded = selected | {
            b.GetOtherAtomIdx(i)
            for i in selected
            for b in molecule.GetAtomWithIdx(i).GetBonds()
            if b.GetBondType() in {Chem.BondType.DOUBLE, Chem.BondType.TRIPLE}
        }
        if len(expanded) > limit:
            raise ValueError("CORE_SIZE_LIMIT")
        if expanded == selected:
            break
        selected = expanded
    order = sorted(selected)
    positions = {value: i for i, value in enumerate(order)}
    labels = [[_atom_label(molecule.GetAtomWithIdx(i)), centers.get(i)] for i in order]
    edges = [
        (
            positions[b.GetBeginAtomIdx()],
            positions[b.GetEndAtomIdx()],
            [str(b.GetBondType()), str(b.GetStereo())],
        )
        for b in molecule.GetBonds()
        if b.GetBeginAtomIdx() in selected and b.GetEndAtomIdx() in selected
    ]
    return _canonical_graph(labels, edges)


def _canonical_molecule(molecule: Chem.Mol) -> str:
    copy = Chem.Mol(molecule)
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(copy, canonical=True, isomericSmiles=True)


def build_shared_reaction_core(
    reaction_smiles: str,
    signature: Mapping[str, Any],
    core: Mapping[str, Any],
) -> SharedReactionCore:
    """Derive qualified views using stored edits and their atom correspondence.

    Unavailable generalized levels are a supported result, not an invitation to
    erase edits. Core-quality warnings survive even when graph checks pass.
    """
    rules = load_shared_core_rules()
    definition_hash = shared_core_definition_hash()
    base: dict[str, Any] = dict(
        definition_hash=definition_hash,
        observation_core_id=str(core.get("core_id", "")),
        reaction_identity="",
        input_identity="",
        product_identity="",
        levels=(),
        realization_key="",
        realization_details=(),
        evidence_status=str(core.get("evidence_status", "unavailable")),
        warnings=tuple(core.get("warnings", ())),
        unavailable_reasons=(),
    )

    def unavailable(reason: str) -> SharedReactionCore:
        # Unsupported generalization retains a valid detailed observation.
        # Contradictions invalidate every projection, including L0.
        invalidate = any(
            word in reason for word in ("CONTRADICT", "INCONSISTENT", "AMBIGUOUS")
        )
        return SharedReactionCore(
            **{
                **base,
                "unavailable_reasons": (reason,),
                **({"levels": ()} if invalidate else {}),
            }
        )

    parsed = parse_reaction_smiles(
        reaction_smiles, include_molecular_interpretation=False
    )
    if not parsed.valid or not signature or not core:
        return unavailable("MISSING_VALID_OBSERVATION")
    if (
        signature.get("schema_version") != REACTION_SIGNATURE_SCHEMA_VERSION
        or core.get("schema_version") != REACTION_CORE_PROJECTION_SCHEMA_VERSION
        or core.get("algorithm_version") != REACTION_CORE_PROJECTION_ALGORITHM_VERSION
    ):
        return unavailable("INCOMPATIBLE_OBSERVATION_SCHEMA")
    if (
        (core.get("quality") or {}).get("status") == "blocked"
        or core.get("evidence_status") not in {"verified", "inferred", "external"}
        or any(
            "CONFLICT" in str(w) or "AMBIGU" in str(w) for w in core.get("warnings", ())
        )
    ):
        return unavailable("UNQUALIFIED_OR_CONFLICTING_CORRESPONDENCE")
    reactants = {
        c.component_index: parse_smiles(c.input_smiles) for c in parsed.reactants
    }
    products = {
        c.component_index: parse_smiles(c.input_smiles) for c in parsed.products
    }
    if any(m is None for m in (*reactants.values(), *products.values())):
        return unavailable("INVALID_MOLECULE")
    inputs = sorted(_canonical_molecule(m) for m in reactants.values())
    outputs = sorted(_canonical_molecule(m) for m in products.values())
    base.update(
        reaction_identity=_digest("SCR1", [inputs, outputs]),
        input_identity=_digest("SCI1", inputs),
        product_identity=_digest("SCP1", outputs),
    )
    transitions = {}
    product_references = set()
    try:
        for transition in core.get("atom_transitions", ()):
            before, after = (
                transition.get("before_state"),
                transition.get("after_state"),
            )
            if after:
                if _ref(after) in product_references:
                    return unavailable("AMBIGUOUS_ATOM_CORRESPONDENCE")
                product_references.add(_ref(after))
            if before:
                key = _ref(before)
                if key in transitions:
                    return unavailable("AMBIGUOUS_ATOM_CORRESPONDENCE")
                transitions[key] = transition
            for state, molecules in ((before, reactants), (after, products)):
                if state is None:
                    continue
                component, index = _ref(state)
                atom = molecules[component].GetAtomWithIdx(index)
                if (
                    atom.GetSymbol() != state["element"]
                    or atom.GetFormalCharge() != state["formal_charge"]
                    or atom.GetIsAromatic() != state["aromatic"]
                    or str(atom.GetHybridization()) != state["hybridization"]
                    or atom.GetTotalNumHs() != state["total_hydrogens"]
                    or atom.GetIsotope() != state["isotope"]
                    or atom.GetNumRadicalElectrons() != state["radical_electrons"]
                    or (
                        atom.GetAtomMapNum()
                        and atom.GetAtomMapNum() != state["atom_map_number"]
                    )
                ):
                    return unavailable("ATOM_REFERENCE_CONTRADICTS_GRAPH")
        edits = tuple(signature.get("edits", ()))
        if not edits:
            return unavailable("NO_OBSERVED_EDITS")
        for edit in edits:
            refs = [edit.get("atom_1"), edit.get("atom_2")]
            for ref in refs:
                if not ref:
                    continue
                if ref.get("side") != "reactant" or _ref(ref) not in transitions:
                    return unavailable("EDIT_REFERENCE_CONTRADICTS_OBSERVATION")
                atom = reactants[_ref(ref)[0]].GetAtomWithIdx(_ref(ref)[1])
                if atom.GetSymbol() != ref["element"]:
                    return unavailable("EDIT_REFERENCE_CONTRADICTS_GRAPH")
            if all(refs):
                left, right = map(_ref, refs)
                old_bond = (
                    reactants[left[0]].GetBondBetweenAtoms(left[1], right[1])
                    if left[0] == right[0]
                    else None
                )
                old_order = str(old_bond.GetBondType()) if old_bond else None
                if old_order != edit.get("old_order"):
                    return unavailable("EDIT_OLD_BOND_CONTRADICTS_GRAPH")
                states = [transitions[p].get("after_state") for p in (left, right)]
                if all(states):
                    p, q = map(_ref, states)
                    new_bond = (
                        products[p[0]].GetBondBetweenAtoms(p[1], q[1])
                        if p[0] == q[0]
                        else None
                    )
                    new_order = str(new_bond.GetBondType()) if new_bond else None
                    if new_order != edit.get("new_order"):
                        return unavailable("EDIT_NEW_BOND_CONTRADICTS_GRAPH")
        graphs = []
        for component, mol in products.items():
            centers = {
                _ref(t["after_state"])[1]: "observed_center"
                for t in transitions.values()
                if t.get("after_state") and _ref(t["after_state"])[0] == component
            }
            if centers:
                graphs.append(
                    _local_graph(mol, centers, 1, rules["maximum_core_atoms"])
                )
        exact_key = core.get("exact_core_key")
        if not exact_key or not graphs:
            return unavailable("MISSING_OBSERVED_CORE")
        payload = _json([exact_key, sorted(graphs)])
        product_key = _digest("SCPS1", [definition_hash, sorted(graphs)])
        base["levels"] = (
            SharedCoreLevel(
                LEVELS[0],
                _digest("SCL0", [definition_hash, payload]),
                payload,
                product_key,
            ),
        )
        if signature.get("event_count") != 1 or core.get("event_count") != 1:
            return unavailable("GENERALIZATION_REQUIRES_SINGLE_EVENT")
        formed = [(i, e) for i, e in enumerate(edits) if e["edit_type"] == "formed"]
        if (
            len(formed) != 1
            or formed[0][1].get("new_order") != "SINGLE"
            or any(
                e["edit_type"] not in {"formed", "broken", "hydrogen_change"}
                for e in edits
            )
        ):
            return unavailable("UNSUPPORTED_EDIT_PATTERN")
        # Hybridization transitions are retained in the before/after center
        # labels (for example amine SP3 -> SP2), rather than discarded.
        if any(
            c["change_type"]
            in set(rules["protected_state_changes"]) - {"hybridization"}
            for c in core.get("state_changes", ())
        ):
            return unavailable("PROTECTED_STATE_CHANGE")
        if (signature.get("topology") or {}).get("reaction_scope") != "intermolecular":
            return unavailable("GENERALIZATION_REQUIRES_INTERMOLECULAR_JOIN")
        endpoints = [_ref(formed[0][1][name]) for name in ("atom_1", "atom_2")]
        if len(set(endpoints)) != 2 or endpoints[0][0] == endpoints[1][0]:
            return unavailable("CONTRADICTORY_JOIN_TOPOLOGY")
        after = [transitions[p].get("after_state") for p in endpoints]
        if any(s is None for s in after) or _ref(after[0])[0] != _ref(after[1])[0]:
            return unavailable("JOIN_CORRESPONDENCE_UNAVAILABLE")
        product = products[_ref(after[0])[0]]
        product_indices = [_ref(s)[1] for s in after]
        bond = product.GetBondBetweenAtoms(*product_indices)
        if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
            return unavailable("FORMED_BOND_CONTRADICTS_GRAPH")
        ports: dict[tuple[int, int], list[str]] = {p: [] for p in endpoints}
        generalized = []
        for i, edit in enumerate(edits):
            if edit["edit_type"] == "formed":
                continue
            if edit["edit_type"] == "hydrogen_change":
                refs = [edit.get("atom_1") or {}, edit.get("atom_2") or {}]
                heavy = [r for r in refs if r and r.get("element") != "H"]
                if len(heavy) != 1 or _ref(heavy[0]) not in ports:
                    return unavailable("UNRESOLVED_HYDROGEN_EDIT")
                point = _ref(heavy[0])
                atom = reactants[point[0]].GetAtomWithIdx(point[1])
                if not (
                    atom.GetSymbol() in rules["hydrogen_port_elements"]
                    or (
                        atom.GetSymbol() == "C"
                        and str(atom.GetHybridization())
                        in rules["hydrogen_port_carbon_hybridizations"]
                    )
                ):
                    return unavailable("PROTECTED_CARBON_HYDROGEN_EDIT")
                before_state = transitions[point]["before_state"]
                after_state = transitions[point]["after_state"]
                if (
                    before_state["total_hydrogens"] - after_state["total_hydrogens"]
                    != 1
                ):
                    return unavailable("UNSUPPORTED_HYDROGEN_DELTA")
                ports[point].append("[H]")
            else:
                a, b = _ref(edit["atom_1"]), _ref(edit["atom_2"])
                if edit.get("old_order") != "SINGLE" or (a in ports) == (b in ports):
                    return unavailable("PROTECTED_BROKEN_BOND")
                point, root = (a, b) if a in ports else (b, a)
                if point[0] != root[0] or transitions.get(root, {}).get("after_state"):
                    return unavailable("DEPARTING_PORT_NOT_ESTABLISHED")
                mol = reactants[point[0]]
                root_atom = mol.GetAtomWithIdx(root[1])
                old_bond = mol.GetBondBetweenAtoms(point[1], root[1])
                if (
                    old_bond is None
                    or old_bond.GetBondType() != Chem.BondType.SINGLE
                    or root_atom.GetSymbol() not in rules["departing_root_elements"]
                ):
                    return unavailable("UNQUALIFIED_DEPARTING_ROOT")
                fragment, pending = set(), [root[1]]
                while pending:
                    current = pending.pop()
                    if current in fragment:
                        continue
                    fragment.add(current)
                    if len(fragment) > rules["maximum_fragment_atoms"]:
                        return unavailable("DEPARTING_FRAGMENT_SIZE_LIMIT")
                    for neighbor in mol.GetAtomWithIdx(current).GetNeighbors():
                        neighbor_id = neighbor.GetIdx()
                        if {current, neighbor_id} == {point[1], root[1]}:
                            continue
                        pending.append(neighbor_id)
                if point[1] in fragment or any(
                    transitions.get((point[0], atom), {}).get("after_state")
                    for atom in fragment
                ):
                    return unavailable("DEPARTING_FRAGMENT_OVERLAPS_RETAINED_CORE")
                copy = Chem.Mol(mol)
                for atom in copy.GetAtoms():
                    atom.SetAtomMapNum(0)
                ports[point].append(
                    Chem.MolFragmentToSmiles(copy, sorted(fragment), canonical=True)
                )
            generalized.append(i)
        if any(len(value) != 1 for value in ports.values()):
            return unavailable("JOIN_REQUIRES_ONE_QUALIFIED_PORT_PER_ENDPOINT")
        # A hydrogen change cannot disappear merely because it was absent from signature edits.
        hydrogen_changes = [
            c for c in core.get("state_changes", ()) if c["change_type"] == "hydrogen"
        ]
        if len(hydrogen_changes) != sum(
            e["edit_type"] == "hydrogen_change" for e in edits
        ):
            return unavailable("UNRECONCILED_HYDROGEN_OBSERVATION")
        center_labels = {}
        realization_labels = {}
        for point, state in zip(endpoints, after):
            before_atom = reactants[point[0]].GetAtomWithIdx(point[1])
            # Aromatic C-H remains protected above. On supported centers, the
            # qualified H/source attachment is recorded separately from identity.
            label = _atom_label(before_atom)
            label[6] = "qualified_port"
            center_labels[_ref(state)[1]] = ["join", label]
            realization_labels[_ref(state)[1]] = ["join", label, ports[point]]
        realization = _local_graph(
            product, realization_labels, 0, rules["maximum_core_atoms"]
        )
        base.update(
            realization_key=_digest("SCRP1", realization),
            realization_details=tuple(sorted(value[0] for value in ports.values())),
        )
        levels = list(base["levels"])
        for name in LEVELS[1:]:
            graph = _local_graph(
                product,
                center_labels,
                rules["levels"][name],
                rules["maximum_core_atoms"],
            )
            product_graph = _local_graph(
                product,
                dict.fromkeys(product_indices, "join"),
                rules["levels"][name],
                rules["maximum_core_atoms"],
            )
            levels.append(
                SharedCoreLevel(
                    name,
                    _digest(
                        "SCL1" if name == LEVELS[1] else "SCL2",
                        [definition_hash, graph],
                    ),
                    graph,
                    _digest("SCPS1", [definition_hash, product_graph]),
                    tuple(generalized),
                )
            )
        return SharedReactionCore(**{**base, "levels": tuple(levels)})
    except (KeyError, IndexError, TypeError, ValueError, RuntimeError):
        return unavailable("INCONSISTENT_OR_UNSUPPORTED_GRAPH_EVIDENCE")

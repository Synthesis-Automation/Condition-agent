"""Direction-neutral contracts and deterministic reaction-operator application.

This module owns only graph-level operator representation and execution.  It
does not retrieve, rank, or interpret named reaction families.  Operators are
stored with explicit forward and reverse SMARTS because reversing an arbitrary
serialized reaction template is not, by itself, evidence that it is executable.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from dataclasses import asdict, dataclass
from functools import lru_cache
from typing import Any, Mapping, Optional, Sequence, Tuple

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from .chemistry.smarts_cache import compile_smarts


REACTION_OPERATOR_SCHEMA_VERSION = "1.0"
REACTION_OPERATOR_APPLICATION_VERSION = "reaction_operator_application.v2"


def _digest(namespace: str, *values: str) -> str:
    payload = "\0".join(values).encode("utf-8")
    return f"{namespace}:{hashlib.sha256(payload).hexdigest()}"


def _canonical_molecule_collection(smiles: str) -> Optional[str]:
    molecule = Chem.MolFromSmiles(str(smiles))
    if molecule is None:
        return None
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    try:
        Chem.SanitizeMol(molecule)
        fragments = Chem.GetMolFrags(molecule, asMols=True, sanitizeFrags=True)
    except Exception:
        return None
    return ".".join(
        sorted(
            Chem.MolToSmiles(
                fragment,
                canonical=True,
                isomericSmiles=True,
            )
            for fragment in fragments
        )
    )


def canonical_molecule_collection(smiles: str) -> Optional[str]:
    """Canonicalize a molecule or dot-separated molecule collection."""

    return _canonical_molecule_collection(smiles)


@dataclass(frozen=True)
class ReactionOperatorPrecedent:
    """Minimal source evidence retained by a direction-neutral operator."""

    reaction_id: str
    reference_id: str
    precursor_smiles: str
    product_smiles: str


@dataclass(frozen=True)
class BidirectionalReactionOperator:
    """One explicitly executable graph transformation in both directions."""

    operator_id: str
    realization_id: str
    template_id: str
    abstraction_level: str
    forward_smarts: str
    reverse_smarts: str
    precursor_smarts: str
    product_smarts: str
    edit_tokens: Tuple[str, ...]
    operator_signature: str
    stereo_policy: str
    observation_support: int
    independent_reference_support: int
    precedents: Tuple[ReactionOperatorPrecedent, ...]
    named_annotations: Tuple[str, ...] = ()
    schema_version: str = REACTION_OPERATOR_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not self.operator_id.startswith("OP1:"):
            raise ValueError("operator ID must use the OP1 namespace")
        if self.realization_id and not self.realization_id.startswith("REAL2:"):
            raise ValueError("realization ID must use the REAL2 namespace")
        if self.forward_smarts != f"{self.precursor_smarts}>>{self.product_smarts}":
            raise ValueError("forward SMARTS contradicts its explicit sides")
        if self.reverse_smarts != f"{self.product_smarts}>>{self.precursor_smarts}":
            raise ValueError("reverse SMARTS contradicts its explicit sides")
        if min(self.observation_support, self.independent_reference_support) < 0:
            raise ValueError("operator support cannot be negative")
        if self.schema_version != REACTION_OPERATOR_SCHEMA_VERSION:
            raise ValueError("unsupported reaction-operator schema")

    @property
    def forward_operator_id(self) -> str:
        """Return a direction-specific identity without changing OP1 identity."""

        return _digest(
            "FOP1",
            self.operator_id,
            self.realization_id,
            self.abstraction_level,
            self.forward_smarts,
            self.schema_version,
            REACTION_OPERATOR_APPLICATION_VERSION,
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible operator representation."""

        return asdict(self)

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "BidirectionalReactionOperator":
        """Load and validate a serialized direction-neutral operator."""

        return cls(
            operator_id=str(value["operator_id"]),
            realization_id=str(value.get("realization_id") or ""),
            template_id=str(value.get("template_id") or ""),
            abstraction_level=str(value.get("abstraction_level") or ""),
            forward_smarts=str(value["forward_smarts"]),
            reverse_smarts=str(value["reverse_smarts"]),
            precursor_smarts=str(value["precursor_smarts"]),
            product_smarts=str(value["product_smarts"]),
            edit_tokens=tuple(str(item) for item in value.get("edit_tokens") or ()),
            operator_signature=str(value.get("operator_signature") or ""),
            stereo_policy=str(value.get("stereo_policy") or "exact"),
            observation_support=int(value.get("observation_support") or 0),
            independent_reference_support=int(
                value.get("independent_reference_support") or 0
            ),
            precedents=tuple(
                ReactionOperatorPrecedent(
                    reaction_id=str(item.get("reaction_id") or ""),
                    reference_id=str(item.get("reference_id") or ""),
                    precursor_smiles=str(item.get("precursor_smiles") or ""),
                    product_smiles=str(item.get("product_smiles") or ""),
                )
                for item in value.get("precedents") or ()
            ),
            named_annotations=tuple(
                str(item) for item in value.get("named_annotations") or ()
            ),
            schema_version=str(
                value.get("schema_version") or REACTION_OPERATOR_SCHEMA_VERSION
            ),
        )


@dataclass(frozen=True)
class OperatorAtomCorrespondence:
    """Reactant provenance for one atom in a generated product graph."""

    product_component_index: int
    product_atom_index: int
    precursor_component_index: int
    precursor_atom_index: int
    atom_map_number: int
    operator_map_number: Optional[int] = None
    precursor_instance_index: int = 0


@dataclass(frozen=True)
class OperatorApplicationOutcome:
    """One sanitized graph outcome and the input components that produced it."""

    product_smiles: str
    mapped_product_smiles: str
    participating_component_indices: Tuple[int, ...]
    participating_precursor_smiles: str
    mapped_participating_precursor_smiles: str
    assignment: Tuple[int, ...]
    atom_correspondence: Tuple[OperatorAtomCorrespondence, ...]
    application_id: str
    reactant_stoichiometry: Tuple[Tuple[int, int], ...] = ()
    uses_virtual_copies: bool = False
    warnings: Tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible application result."""

        return asdict(self)


@lru_cache(maxsize=20_000)
def _compiled_reaction(reaction_smarts: str) -> Any:
    reaction = rdChemReactions.ReactionFromSmarts(reaction_smarts)
    if reaction is None:
        raise ValueError("reaction SMARTS could not be compiled")
    reaction.Initialize()
    return reaction


def _component_molecules(smiles: str) -> tuple[tuple[str, Any], ...]:
    canonical = _canonical_molecule_collection(smiles)
    if canonical is None:
        raise ValueError("starting materials must be valid SMILES")
    values = []
    for component in canonical.split("."):
        molecule = Chem.MolFromSmiles(component)
        if molecule is None:
            raise ValueError("starting-material component could not be parsed")
        values.append((component, molecule))
    return tuple(values)


def _candidate_assignments(
    reaction: Any,
    components: Sequence[tuple[str, Any]],
    *,
    max_assignments: int,
    allow_self_reaction: bool,
) -> tuple[tuple[int, ...], ...]:
    template_count = int(reaction.GetNumReactantTemplates())
    if template_count < 1 or (
        template_count > len(components) and not allow_self_reaction
    ):
        return ()
    matches: list[tuple[int, ...]] = []
    for template_index in range(template_count):
        pattern = reaction.GetReactantTemplate(template_index)
        matching = tuple(
            index
            for index, (_, molecule) in enumerate(components)
            if molecule.HasSubstructMatch(pattern)
        )
        if not matching:
            return ()
        matches.append(matching)
    assignments = []
    injective_signatures = set()
    for assignment in itertools.product(*matches):
        if len(set(assignment)) != len(assignment):
            continue
        normalized = tuple(int(index) for index in assignment)
        assignments.append(normalized)
        injective_signatures.add(
            tuple(components[index][0] for index in normalized)
        )
        if len(assignments) >= max_assignments:
            return tuple(assignments)
    if not allow_self_reaction:
        return tuple(assignments)
    for assignment in itertools.product(*matches):
        if len(set(assignment)) == len(assignment):
            continue
        normalized = tuple(int(index) for index in assignment)
        molecular_signature = tuple(
            components[index][0] for index in normalized
        )
        if molecular_signature in injective_signatures:
            continue
        assignments.append(normalized)
        if len(assignments) >= max_assignments:
            break
    return tuple(assignments)


def apply_reaction_smarts(
    reaction_smarts: str,
    starting_materials: str,
    *,
    max_assignments: int = 256,
    max_outcomes: int = 1024,
    allow_self_reaction: bool = False,
) -> Tuple[OperatorApplicationOutcome, ...]:
    """Apply reaction SMARTS to deterministic compatible component tuples.

    When ``allow_self_reaction`` is true, one supplied component may fill more
    than one reactant-template role.  Each repeated role is treated as a
    separate stoichiometric molecular instance, while correspondence retains
    the original user-supplied component index.
    """

    if min(max_assignments, max_outcomes) < 1:
        raise ValueError("application limits must be positive")
    reaction = _compiled_reaction(reaction_smarts)
    components = _component_molecules(starting_materials)
    assignments = _candidate_assignments(
        reaction,
        components,
        max_assignments=max_assignments,
        allow_self_reaction=allow_self_reaction,
    )
    outcomes: dict[tuple[str, tuple[int, ...]], OperatorApplicationOutcome] = {}
    for assignment in assignments:
        reactants = tuple(Chem.Mol(components[index][1]) for index in assignment)
        atom_maps = {}
        next_map = 1
        for reactant_position, molecule in enumerate(reactants):
            for atom in molecule.GetAtoms():
                atom_maps[(reactant_position, int(atom.GetIdx()))] = next_map
                next_map += 1
        try:
            generated = reaction.RunReactants(reactants)
        except Exception:
            continue
        for product_tuple in generated:
            product_parts = []
            mapped_product_parts = []
            correspondence = []
            valid = True
            for product_component_index, product in enumerate(product_tuple):
                try:
                    Chem.SanitizeMol(product)
                    for atom in product.GetAtoms():
                        if not (
                            atom.HasProp("react_idx") and atom.HasProp("react_atom_idx")
                        ):
                            continue
                        reactant_position = int(atom.GetProp("react_idx"))
                        precursor_component_index = assignment[reactant_position]
                        precursor_atom_index = int(atom.GetProp("react_atom_idx"))
                        map_number = atom_maps[
                            (reactant_position, precursor_atom_index)
                        ]
                        operator_map_number = (
                            int(atom.GetProp("old_mapno"))
                            if atom.HasProp("old_mapno")
                            else None
                        )
                        atom.SetAtomMapNum(map_number)
                        correspondence.append(
                            OperatorAtomCorrespondence(
                                product_component_index=product_component_index,
                                product_atom_index=int(atom.GetIdx()),
                                precursor_component_index=precursor_component_index,
                                precursor_atom_index=precursor_atom_index,
                                atom_map_number=map_number,
                                operator_map_number=operator_map_number,
                                precursor_instance_index=reactant_position,
                            )
                        )
                    mapped_product_parts.append(
                        Chem.MolToSmiles(
                            product,
                            canonical=True,
                            isomericSmiles=True,
                        )
                    )
                    unmapped = Chem.Mol(product)
                    for atom in unmapped.GetAtoms():
                        atom.SetAtomMapNum(0)
                    product_parts.append(
                        Chem.MolToSmiles(
                            unmapped,
                            canonical=True,
                            isomericSmiles=True,
                        )
                    )
                except Exception:
                    valid = False
                    break
            if not valid or not product_parts:
                continue
            product_smiles = _canonical_molecule_collection(".".join(product_parts))
            if product_smiles is None:
                continue
            mapped_product_smiles = ".".join(sorted(mapped_product_parts))
            indices = tuple(sorted(set(assignment)))
            precursor_smiles = ".".join(
                sorted(components[index][0] for index in assignment)
            )
            mapped_precursor_parts = []
            for reactant_position, component_index in enumerate(assignment):
                precursor = Chem.Mol(components[component_index][1])
                for atom in precursor.GetAtoms():
                    atom.SetAtomMapNum(
                        atom_maps[(reactant_position, int(atom.GetIdx()))]
                    )
                mapped_precursor_parts.append(
                    Chem.MolToSmiles(
                        precursor,
                        canonical=True,
                        isomericSmiles=True,
                    )
                )
            mapped_precursor_smiles = ".".join(sorted(mapped_precursor_parts))
            stoichiometry = tuple(
                sorted(
                    (component_index, assignment.count(component_index))
                    for component_index in set(assignment)
                )
            )
            uses_virtual_copies = any(count > 1 for _, count in stoichiometry)
            key = (product_smiles, indices)
            outcomes[key] = OperatorApplicationOutcome(
                product_smiles=product_smiles,
                mapped_product_smiles=mapped_product_smiles,
                participating_component_indices=indices,
                participating_precursor_smiles=precursor_smiles,
                mapped_participating_precursor_smiles=mapped_precursor_smiles,
                assignment=assignment,
                atom_correspondence=tuple(
                    sorted(
                        correspondence,
                        key=lambda item: (
                            item.product_component_index,
                            item.product_atom_index,
                        ),
                    )
                ),
                application_id=_digest(
                    "OPA1",
                    reaction_smarts,
                    precursor_smiles,
                    product_smiles,
                    json.dumps(assignment),
                    REACTION_OPERATOR_APPLICATION_VERSION,
                ),
                reactant_stoichiometry=stoichiometry,
                uses_virtual_copies=uses_virtual_copies,
                warnings=(
                    ("SELF_REACTION_USES_VIRTUAL_STOICHIOMETRIC_COPIES",)
                    if uses_virtual_copies
                    else ()
                ),
            )
            if len(outcomes) >= max_outcomes:
                break
        if len(outcomes) >= max_outcomes:
            break
    return tuple(
        sorted(
            outcomes.values(),
            key=lambda value: (
                value.product_smiles,
                value.participating_component_indices,
                value.assignment,
            ),
        )
    )


def apply_forward_operator(
    operator: BidirectionalReactionOperator,
    starting_materials: str,
    *,
    max_assignments: int = 256,
    max_outcomes: int = 1024,
    allow_self_reaction: bool = False,
) -> Tuple[OperatorApplicationOutcome, ...]:
    """Apply one admitted operator from precursors to possible products."""

    return apply_reaction_smarts(
        operator.forward_smarts,
        starting_materials,
        max_assignments=max_assignments,
        max_outcomes=max_outcomes,
        allow_self_reaction=allow_self_reaction,
    )


def reverse_recovers_precursors(
    operator: BidirectionalReactionOperator,
    outcome: OperatorApplicationOutcome,
) -> bool:
    """Return whether reverse application recovers the contributing precursors."""

    reverse = apply_reaction_smarts(
        operator.reverse_smarts,
        outcome.product_smiles,
        max_assignments=16,
        max_outcomes=256,
    )
    expected = _canonical_molecule_collection(outcome.participating_precursor_smiles)
    return bool(
        expected and any(candidate.product_smiles == expected for candidate in reverse)
    )


def precursor_pattern_tokens(precursor_smarts: str) -> Tuple[str, ...]:
    """Return conservative element tokens for precursor-side retrieval."""

    tokens = []
    for component in precursor_smarts.split("."):
        pattern = compile_smarts(component, validate=False)
        if pattern is None:
            continue
        elements = sorted(
            atom.GetSymbol() for atom in pattern.GetAtoms() if atom.GetAtomicNum() > 0
        )
        tokens.append("component:" + ",".join(elements))
    return tuple(sorted(tokens))


__all__ = [
    "BidirectionalReactionOperator",
    "OperatorAtomCorrespondence",
    "OperatorApplicationOutcome",
    "REACTION_OPERATOR_APPLICATION_VERSION",
    "REACTION_OPERATOR_SCHEMA_VERSION",
    "ReactionOperatorPrecedent",
    "apply_forward_operator",
    "apply_reaction_smarts",
    "canonical_molecule_collection",
    "precursor_pattern_tokens",
    "reverse_recovers_precursors",
]

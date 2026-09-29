"""Optional target-derived search fragments; no disconnections or corpus calls.

Candidates can overlap. Cutting a peripheral bond defines a search query, not a
precursor or a proposed reaction. Atom IDs refer to the returned canonical target.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

from rdkit import Chem

from .fragment_search import QUERY_COMPILER_VERSION, compile_fragment_query, indexed_product


SCHEMA_VERSION = "search_fragment_suggestions.v1"
KINDS = ("ring_system", "functional_region", "scaffold", "core_with_context", "whole_target")


def search_fragment_policy() -> dict[str, Any]:
    """Load validated, versioned generation limits and ordering weights."""
    path = Path(__file__).with_name("definitions") / "search_fragments.v1.json"
    policy = json.loads(path.read_text("utf-8"))
    if policy.get("schema_version") != "search_fragment_policy.v1" or not policy.get("definition_version"):
        raise ValueError("Invalid search fragment policy version")
    for field in ("max_target_characters", "max_target_atoms", "max_candidates",
                  "max_proposals", "min_region_atoms", "context_radius", "region_radius"):
        if type(policy.get(field)) is not int or policy[field] < 1:
            raise ValueError(f"Invalid search fragment policy: {field}")
    if policy.get("kind_order") != list(KINDS):
        raise ValueError("Invalid fragment kind order")
    weights = policy.get("feature_weights", {})
    if set(weights) != {"ring_junctions", "heteroatoms", "stereocenters", "multiple_bonds"}:
        raise ValueError("Invalid fragment feature weights")
    if any(type(v) is not int or v < 0 for v in weights.values()):
        raise ValueError("Fragment weights must be nonnegative integers")
    return policy


@dataclass(frozen=True)
class FragmentAtom:
    """Canonical-target atom ID with a link back to the parsed input atom."""

    atom_id: int
    input_atom_index: int
    element: str
    isotope: int
    charge: int
    aromatic: bool
    chiral_tag: str


@dataclass(frozen=True)
class FragmentBoundary:
    """Observed target bond crossing the query boundary; not a disconnection."""

    retained_atom_id: int
    omitted_atom_id: int
    bond_type: str


@dataclass(frozen=True)
class FragmentFeatures:
    """Structural descriptors used only as transparent ordering heuristics."""

    atom_count: int
    ring_junctions: int
    heteroatoms: int
    stereocenters: int
    multiple_bonds: int


@dataclass(frozen=True)
class SearchFragmentCandidate:
    """A search-ready query verified against its selected target atoms."""

    candidate_id: str
    kind: str
    query: str
    query_format: str
    topology: str
    target_atom_ids: tuple[int, ...]
    query_atom_target_ids: tuple[int, ...]
    omitted_atom_ids: tuple[int, ...]
    boundaries: tuple[FragmentBoundary, ...]
    features: FragmentFeatures
    structural_priority: int
    reasons: tuple[str, ...]
    cautions: tuple[str, ...]
    matches_target: bool = True


@dataclass(frozen=True)
class SearchFragmentSuggestions:
    """Bounded optional candidates with explicit coverage and stable identities."""

    schema_version: str
    definition_version: str
    query_compiler_version: str
    target_id: str
    target_smiles: str
    target_atoms: tuple[FragmentAtom, ...]
    candidates: tuple[SearchFragmentCandidate, ...]
    selection_mode: str
    generated_count: int
    rejected_count: int
    generation_truncated: bool
    output_truncated: bool
    limitations: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        """Serialize tuples as JSON-compatible arrays."""
        return json.loads(json.dumps(asdict(self)))


def _identity(prefix: str, value: Any) -> str:
    return prefix + hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def _ring_systems(mol: Any) -> tuple[frozenset[int], ...]:
    # Connected components of ring bonds include fused, bridged and spiro rings.
    remaining = {a.GetIdx() for a in mol.GetAtoms() if a.IsInRing()}
    systems = []
    while remaining:
        pending, selected = [min(remaining)], set()
        while pending:
            index = pending.pop()
            if index in selected:
                continue
            selected.add(index)
            pending.extend(b.GetOtherAtomIdx(index) for b in mol.GetAtomWithIdx(index).GetBonds()
                           if b.IsInRing())
        remaining -= selected
        systems.append(frozenset(selected))
    return tuple(systems)


def _closure(mol: Any, selected: set[int], rings: tuple[frozenset[int], ...]) -> set[int]:
    """Retain full rings, multiple bonds, charged valence and specified stereo."""
    selected = set(selected)
    while True:
        before = set(selected)
        for ring in rings:
            if selected & ring:
                selected.update(ring)
        for index in tuple(selected):
            atom = mol.GetAtomWithIdx(index)
            for bond in atom.GetBonds():
                if (bond.GetBondTypeAsDouble() != 1 or atom.GetFormalCharge()
                        or atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
                        or atom.GetAtomicNum() == 1):
                    selected.add(bond.GetOtherAtomIdx(index))
                if not bond.GetIsAromatic() and bond.GetBondTypeAsDouble() > 1:
                    for endpoint in (bond.GetBeginAtom(), bond.GetEndAtom()):
                        selected.update(n.GetIdx() for n in endpoint.GetNeighbors())
                if bond.GetStereo() != Chem.BondStereo.STEREONONE:
                    selected.update(bond.GetStereoAtoms())
        if selected == before:
            return selected


def _neighborhood(mol: Any, selected: set[int], radius: int) -> set[int]:
    selected = set(selected)
    for _ in range(radius):
        selected |= {n.GetIdx() for index in selected for n in mol.GetAtomWithIdx(index).GetNeighbors()}
    return selected


def _scaffold_atoms(mol: Any) -> set[int]:
    """Prune acyclic leaves, retaining ring systems and their connecting paths."""
    selected = set(range(mol.GetNumAtoms()))
    while True:
        leaves = {i for i in selected if not mol.GetAtomWithIdx(i).IsInRing()
                  and sum(n.GetIdx() in selected for n in mol.GetAtomWithIdx(i).GetNeighbors()) <= 1}
        if not leaves:
            return selected
        selected -= leaves


def _materialize(
    mol: Any, selected: set[int], kind: str, target_id: str, policy: dict[str, Any],
) -> SearchFragmentCandidate:
    smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=sorted(selected), canonical=True, isomericSmiles=True)
    query = compile_fragment_query(smiles)
    parameters = Chem.SubstructMatchParameters()
    parameters.useChirality = True
    parameters.setExtraFinalCheck(lambda _mol, match: set(match) == selected)
    match = mol.GetSubstructMatch(query.molecule, parameters)
    if not match:
        raise ValueError("Extracted query does not match the selected target atoms")
    for query_index, target_index in enumerate(match):
        source = mol.GetAtomWithIdx(target_index)
        extracted = query.molecule.GetAtomWithIdx(query_index)
        if (source.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
                and extracted.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED):
            raise ValueError("Selection loses specified stereochemistry; retain more context")
    for bond in query.molecule.GetBonds():
        source = mol.GetBondBetweenAtoms(match[bond.GetBeginAtomIdx()], match[bond.GetEndAtomIdx()])
        if source.GetStereo() != Chem.BondStereo.STEREONONE and bond.GetStereo() == Chem.BondStereo.STEREONONE:
            raise ValueError("Selection loses specified bond stereochemistry; retain more context")
    atoms = [mol.GetAtomWithIdx(i) for i in sorted(selected)]
    features = FragmentFeatures(
        len(selected), sum(sum(b.IsInRing() for b in a.GetBonds()) > 2 for a in atoms),
        sum(a.GetAtomicNum() not in (1, 6) for a in atoms),
        sum(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in atoms),
        sum(b.GetBeginAtomIdx() in selected and b.GetEndAtomIdx() in selected
            and not b.GetIsAromatic() and b.GetBondTypeAsDouble() > 1 for b in mol.GetBonds()),
    )
    priority = sum(getattr(features, k) * v for k, v in policy["feature_weights"].items())
    boundaries = []
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if (a in selected) != (b in selected):
            boundaries.append(FragmentBoundary(a if a in selected else b, b if a in selected else a,
                                               str(bond.GetBondType())))
    descriptions = {
        "ring_system": "Retains a complete ring system and required valence/stereo context.",
        "core_with_context": "Adds nearby substituent context to distinguish a ring core.",
        "scaffold": "Retains ring systems and the paths connecting them.",
        "functional_region": "Retains a local heteroatom, unsaturation or stereochemical region.",
        "whole_target": "Retains the whole target as a specific reference query.",
        "selected_atoms": "Extracts exactly the requested canonical-target atom IDs.",
    }
    cautions = ["SEARCH_QUERY_NOT_A_PRECURSOR"]
    if priority == 0:
        cautions.append("LOW_STRUCTURAL_SPECIFICITY_MAY_BE_BROAD")
    if boundaries:
        cautions.append("PERIPHERAL_SUBSTITUTION_IS_UNCONSTRAINED")
    return SearchFragmentCandidate(
        _identity("SFC1:", [target_id, query.query_id, sorted(selected), policy["definition_version"]]),
        kind, smiles, "smiles", "preserve_rings", tuple(sorted(selected)), tuple(match),
        tuple(sorted(set(range(mol.GetNumAtoms())) - selected)),
        tuple(sorted(boundaries, key=lambda b: (b.retained_atom_id, b.omitted_atom_id))),
        features, priority, (descriptions[kind],), tuple(cautions),
    )


def suggest_search_fragments(
    target_smiles: str, limit: int = 5, selected_atom_ids: list[int] | None = None,
) -> SearchFragmentSuggestions:
    """Suggest optional graph-derived queries or extract one explicit selection.

    Selection IDs are zero-based IDs in the returned canonical target, not the
    user's original SMILES order. Explicit selections must be connected and keep
    full ring systems and required valence/stereo context; they are never expanded
    silently. No corpus, atom mapper, forward check or route planner is called.
    """
    policy = search_fragment_policy()
    if type(limit) is not int or not 1 <= limit <= policy["max_candidates"]:
        raise ValueError("limit must be an integer from 1 to 5")
    if not isinstance(target_smiles, str) or not target_smiles.strip() or len(target_smiles) > policy["max_target_characters"]:
        raise ValueError("Provide one target SMILES of at most 5000 characters")
    canonical, input_order = indexed_product(target_smiles)
    mol = Chem.MolFromSmiles(canonical)
    if len(Chem.GetMolFrags(mol)) != 1 or mol.GetNumAtoms() > policy["max_target_atoms"]:
        raise ValueError("Provide one connected target with at most 200 atoms")
    if any(a.GetAtomicNum() == 0 or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        raise ValueError("Wildcard and radical targets are unsupported")
    target_id = _identity("SFT1:", canonical)
    references = tuple(FragmentAtom(a.GetIdx(), input_order[a.GetIdx()], a.GetSymbol(), a.GetIsotope(),
                                    a.GetFormalCharge(), a.GetIsAromatic(), str(a.GetChiralTag()))
                       for a in mol.GetAtoms())
    rings = _ring_systems(mol)
    proposals: list[tuple[str, set[int], str]] = []
    if selected_atom_ids is not None:
        if (not isinstance(selected_atom_ids, list) or not selected_atom_ids
                or any(type(i) is not int or not 0 <= i < mol.GetNumAtoms() for i in selected_atom_ids)
                or len(set(selected_atom_ids)) != len(selected_atom_ids)):
            raise ValueError("selected_atom_ids must be unique canonical-target atom IDs")
        selected = set(selected_atom_ids)
        if _closure(mol, selected, rings) != selected:
            raise ValueError("Selection must retain complete ring systems, multiple bonds and valence/stereo context")
        proposals.append(("selected_atoms", selected, "selected"))
    else:
        for number, ring in enumerate(rings):
            proposals.append(("ring_system", _closure(mol, set(ring), rings), f"ring-{number}"))
            context = _neighborhood(mol, set(ring), policy["context_radius"])
            proposals.append(("core_with_context", _closure(mol, context, rings), f"ring-{number}"))
        scaffold = _scaffold_atoms(mol)
        if scaffold:
            proposals.append(("scaffold", _closure(mol, scaffold, rings), "scaffold"))
        for atom in mol.GetAtoms():
            if (not atom.IsInRing() and (atom.GetAtomicNum() not in (1, 6)
                    or atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
                    or any(b.GetBondTypeAsDouble() > 1 for b in atom.GetBonds()))):
                region = _closure(mol, _neighborhood(mol, {atom.GetIdx()}, policy["region_radius"]), rings)
                if len(region) >= policy["min_region_atoms"]:
                    proposals.append(("functional_region", region, "functional"))
        proposals.append(("whole_target", set(range(mol.GetNumAtoms())), "whole"))
    truncated = len(proposals) > policy["max_proposals"]
    candidates: dict[str, tuple[SearchFragmentCandidate, str]] = {}
    rejected = 0
    for kind, atoms, family in proposals[:policy["max_proposals"]]:
        try:
            candidate = _materialize(mol, atoms, kind, target_id, policy)
        except ValueError:
            if selected_atom_ids is not None:
                raise
            rejected += 1
            continue
        existing = candidates.get(candidate.query)
        if existing is None or KINDS.index(kind) < KINDS.index(existing[0].kind):
            candidates[candidate.query] = candidate, family
    ordered = sorted(candidates.values(), key=lambda item: (
        item[0].structural_priority == 0,
        KINDS.index(item[0].kind) if item[0].kind in KINDS else -1,
        -item[0].structural_priority, item[0].features.atom_count, item[0].query,
    ))
    chosen, families = [], set()
    for candidate, family in ordered:
        if family not in families:
            chosen.append(candidate)
            families.add(family)
    chosen.extend(c for c, _ in ordered if c not in chosen)
    return SearchFragmentSuggestions(
        SCHEMA_VERSION, policy["definition_version"], QUERY_COMPILER_VERSION, target_id, canonical, references,
        tuple(chosen[:limit]), "selected_atoms" if selected_atom_ids is not None else "suggested",
        len(candidates), rejected, truncated, len(candidates) > limit,
        ("Candidates are optional, overlapping search regions, not a synthetic partition or route.",
         "Structural priority is a heuristic, not rarity, synthetic difficulty or experimental evidence.",
         "No corpus was searched. Use search_fragment_precedents only for candidates you choose.",
         "Specified stereo is retained within each selected region; omitted atoms and their features are not constrained.",
         "Atom IDs refer to the returned canonical target; input_atom_index links to the parsed input."),
    )

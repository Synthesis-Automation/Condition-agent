"""Bounded structural comparison for precedent transfer, never reaction mapping."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from itertools import product
from typing import Any

from rdkit import Chem, rdBase
from rdkit.Chem import rdFMCS

from .chemistry.smarts_cache import compile_smarts
from .fragment_search import compile_fragment_query, fragment_embeddings
from .molecule_inspection import InspectedMolecule, _prepare_molecule, molecule_inspection_policy


def _atom_identity(atom: Any) -> tuple[Any, ...]:
    return (atom.GetAtomicNum(), atom.GetIsotope(), atom.GetFormalCharge(),
            atom.GetIsAromatic(), sum(b.IsInRing() for b in atom.GetBonds()))


class _StrictAtoms(rdFMCS.MCSAtomCompare):
    """Require exact chemistry and ring junctions while allowing substituents."""

    def __call__(self, parameters: Any, left: Any, i: int, right: Any, j: int) -> bool:
        return (_atom_identity(left.GetAtomWithIdx(i)) == _atom_identity(right.GetAtomWithIdx(j))
                and self.CheckAtomRingMatch(parameters, left, i, right, j))


@dataclass(frozen=True)
class CoreBoundary:
    """Bond from one common-core position to an unmatched substituent."""

    core_position: int
    core_atom_id: int
    substituent_atom_id: int
    bond_type: str


@dataclass(frozen=True)
class CoreAlignment:
    """One possible structural alignment, with unmatched atoms and attachments."""

    atom_pairs: tuple[tuple[int, int], ...]
    left_only_atom_ids: tuple[int, ...]
    right_only_atom_ids: tuple[int, ...]
    left_boundaries: tuple[CoreBoundary, ...]
    right_boundaries: tuple[CoreBoundary, ...]


@dataclass(frozen=True)
class MoleculeComparison:
    """Bounded evidence with unresolved alignment and stereo scope kept explicit."""

    status: str
    left: InspectedMolecule
    right: InspectedMolecule
    same_constitution: bool
    stereo_relationship: str
    core_method: str
    core_smarts: str | None
    core_atom_count: int
    left_coverage: float
    right_coverage: float
    alignments: tuple[CoreAlignment, ...]
    alignment_count_observed: int
    alignment_ambiguous: bool
    embeddings_truncated: bool
    alignments_truncated: bool
    search_timed_out: bool
    warnings: tuple[str, ...]
    definition_version: str
    rdkit_version: str
    limitations: tuple[str, ...] = (
        "Core alignments are structural hypotheses, not reaction atom maps or evidence of feasibility.",
        "Automatic search returns one largest core; alternative equally large cores are not enumerated.",
        "core_smarts is a search pattern; only returned atom_pairs have passed strict atom and induced-bond checks.",
        "Stereo is compared for the whole graph only when constitution, charge and isotopes match; CIP changes are not mechanistic inversion.",
        "No tautomer, protonation, salt or aromaticity relaxation is performed. Atom IDs use each returned canonical SMILES.",
    )
    schema_version: str = "molecule_comparison.v1"

    def to_dict(self) -> dict[str, Any]:
        """Return serializable structural evidence and all coverage limits."""
        return asdict(self)


def _constitution(mol: Any) -> str:
    copy = Chem.Mol(mol)
    Chem.RemoveStereochemistry(copy)
    return Chem.MolToSmiles(copy, canonical=True, isomericSmiles=True)


def _boundaries(mol: Any, match: tuple[int, ...]) -> tuple[CoreBoundary, ...]:
    positions = {atom: position for position, atom in enumerate(match)}
    result = []
    for core in match:
        for bond in mol.GetAtomWithIdx(core).GetBonds():
            other = bond.GetOtherAtomIdx(core)
            if other not in positions:
                result.append(CoreBoundary(positions[core], core, other, str(bond.GetBondType())))
    return tuple(sorted(result, key=lambda b: (b.core_position, b.substituent_atom_id)))


def _automatic_core(left: Any, right: Any, seconds: int) -> tuple[Any, str | None, bool]:
    parameters = rdFMCS.MCSParameters()
    parameters.Timeout = seconds
    parameters.MaximizeBonds = False
    parameters.AtomTyper = _StrictAtoms()
    parameters.BondTyper = rdFMCS.BondCompare.CompareOrderExact
    parameters.AtomCompareParameters.MatchFormalCharge = True
    # The comparator checks isotopes; MatchIsotope's SMARTS writer replaces
    # element constraints with isotope wildcards, broadening re-embeddings.
    parameters.AtomCompareParameters.MatchIsotope = False
    parameters.AtomCompareParameters.RingMatchesRingOnly = True
    parameters.AtomCompareParameters.CompleteRingsOnly = True
    parameters.BondCompareParameters.RingMatchesRingOnly = True
    parameters.BondCompareParameters.CompleteRingsOnly = True
    parameters.BondCompareParameters.MatchFusedRings = True
    parameters.BondCompareParameters.MatchFusedRingsStrict = True
    result = rdFMCS.FindMCS([left, right], parameters)
    smarts = result.smartsString or None
    # Copy cached query objects before RDKit ring information is initialized.
    query = Chem.Mol(compile_smarts(smarts, validate=True)) if smarts else None
    if query is not None:
        Chem.GetSymmSSSR(query)
    return query, smarts, result.canceled


def _embeddings(mol: Any, query: Any, maximum: int) -> tuple[tuple[tuple[int, ...], ...], bool]:
    parameters = Chem.SubstructMatchParameters()
    parameters.uniquify = False
    parameters.maxMatches = maximum + 1
    matches = mol.GetSubstructMatches(query, parameters)
    return tuple(sorted(matches[:maximum])), len(matches) > maximum


def _compatible_alignment(left: Any, right: Any, a: tuple[int, ...], b: tuple[int, ...]) -> bool:
    # Re-embedding an MCS SMARTS must not loosen atom/bond constraints. Require
    # induced graph agreement, including bonds omitted by a subgraph search.
    if any(_atom_identity(left.GetAtomWithIdx(i)) != _atom_identity(right.GetAtomWithIdx(j))
           for i, j in zip(a, b)):
        return False
    mapping = dict(zip(a, b))
    def edges(mol: Any, atoms: set[int], remap: dict[int, int]) -> set[tuple[Any, ...]]:
        return {(tuple(sorted((remap[e.GetBeginAtomIdx()], remap[e.GetEndAtomIdx()]))),
                 str(e.GetBondType()), e.IsInRing()) for e in mol.GetBonds()
                if e.GetBeginAtomIdx() in atoms and e.GetEndAtomIdx() in atoms}
    return edges(left, set(a), mapping) == edges(right, set(b), {i: i for i in b})


def compare_molecules(
    left_smiles: str, right_smiles: str, core_smiles: str | None = None,
    timeout_seconds: int = 2,
) -> MoleculeComparison:
    """Compare strict common cores, differing atoms, attachments and stereo gaps.

    Optionally anchor to an explicit core SMILES; otherwise use bounded MCS.
    Molecules must be connected. Symmetric alternatives and incomplete searches
    are reported, never converted into invented reaction correspondence.
    """
    policy = molecule_inspection_policy()
    if type(timeout_seconds) is not int or not 1 <= timeout_seconds <= policy["max_timeout_seconds"]:
        raise ValueError(f"timeout_seconds must be 1..{policy['max_timeout_seconds']}")
    left, left_info = _prepare_molecule(left_smiles, policy)
    right, right_info = _prepare_molecule(right_smiles, policy)
    same = _constitution(left) == _constitution(right)
    if not same:
        stereo = "not_compared_different_graphs"
    elif left_info.canonical_smiles == right_info.canonical_smiles:
        stereo = "same_specification"
    elif any(s.specification != "Specified" for m in (left_info, right_info) for s in m.stereochemistry):
        stereo = "different_or_incomplete_specification"
    else:
        stereo = "different_specified_stereochemistry"
    maximum = policy["max_embeddings"]
    if core_smiles is not None:
        core = compile_fragment_query(core_smiles)
        query, smarts, timed_out = core.molecule, Chem.MolToSmarts(core.molecule), False
        a, cut_a = fragment_embeddings(core, left, maximum=maximum)
        b, cut_b = fragment_embeddings(core, right, maximum=maximum)
        method = "supplied_core"
    elif same:
        query = Chem.Mol(left)
        Chem.RemoveStereochemistry(query)
        smarts, timed_out = Chem.MolToSmarts(query), False
        a, cut_a = _embeddings(left, query, maximum)
        b, cut_b = _embeddings(right, query, maximum)
        method = "exact_constitution"
    else:
        query, smarts, timed_out = _automatic_core(left, right, timeout_seconds)
        a, cut_a = _embeddings(left, query, maximum) if query is not None else ((), False)
        b, cut_b = _embeddings(right, query, maximum) if query is not None else ((), False)
        method = "bounded_mcs"
    pairs = sorted({tuple(sorted(zip(i, j))) for i, j in product(a, b)
                    if _compatible_alignment(left, right, i, j)})
    alignments = []
    for pair in pairs[:policy["max_alignments"]]:
        first, second = tuple(i for i, _ in pair), tuple(j for _, j in pair)
        alignments.append(CoreAlignment(
            pair, tuple(sorted(set(range(left.GetNumAtoms())) - set(first))),
            tuple(sorted(set(range(right.GetNumAtoms())) - set(second))),
            _boundaries(left, first), _boundaries(right, second),
        ))
    count = query.GetNumAtoms() if pairs else 0
    status = "partial_timeout" if timed_out else "completed" if pairs else "no_verified_core"
    if not pairs and (cut_a or cut_b):
        status = "unresolved_embedding_limit"
    warnings = list(dict.fromkeys((*left_info.warnings, *right_info.warnings)))
    if timed_out:
        warnings.append("MCS_TIMEOUT_CORE_NOT_PROVEN_MAXIMAL")
    if cut_a or cut_b:
        warnings.append("EMBEDDING_LIMIT_ALIGNMENT_COUNTS_ARE_LOWER_BOUNDS")
    if len(pairs) > 1:
        warnings.append("AMBIGUOUS_STRUCTURAL_ALIGNMENT")
    if not same:
        warnings.append("PARTIAL_CORE_STEREOCHEMISTRY_NOT_COMPARED")
    return MoleculeComparison(
        status, left_info, right_info, same, stereo, method, smarts, count,
        round(count / left.GetNumAtoms(), 4), round(count / right.GetNumAtoms(), 4),
        tuple(alignments), len(pairs), len(pairs) > 1 or cut_a or cut_b,
        cut_a or cut_b, len(pairs) > len(alignments), timed_out, tuple(warnings),
        policy["definition_version"], rdBase.rdkitVersion,
    )

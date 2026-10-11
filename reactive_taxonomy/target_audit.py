"""Compact deterministic target audit for chemistry-tool consumers."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Optional, Tuple

from .api import analyze_molecule, observe_molecular_structure
from .chemistry.rdkit_utils import parse_smiles
from .fragment_search import indexed_product
from .descriptors.registry import descriptor_definition_versions


TARGET_AUDIT_SCHEMA_VERSION = "target_audit.v2"


@dataclass(frozen=True)
class TargetComponent:
    """Component-local atom/bond index -> whole canonical molecule index."""

    component_index: int
    canonical_atom_indices: Tuple[int, ...]
    canonical_bond_indices: Tuple[int, ...]


@dataclass(frozen=True)
class CarbonylNitrogenObservation:
    """Observed C(=O)-N adjacency, independent of motif interpretation."""

    carbon_index: int
    oxygen_index: int
    nitrogen_index: int
    cn_bond_in_ring: bool
    carbon_aromatic: bool
    nitrogen_aromatic: bool


@dataclass(frozen=True)
class TargetStereocenter:
    """One tetrahedral stereocenter perceived by RDKit."""

    atom_index: int
    element: str
    assignment: str
    assigned: bool


@dataclass(frozen=True)
class TargetDoubleBondStereo:
    """One explicitly annotated double-bond stereo relationship."""

    bond_index: int
    atom_1_index: int
    atom_2_index: int
    assignment: str


@dataclass(frozen=True)
class TargetMotif:
    """Compact definition-derived molecular motif evidence."""

    motif_id: str
    chemist_label: str
    component_index: int
    atom_indices: Tuple[int, ...]
    tags: Tuple[str, ...]


@dataclass(frozen=True)
class TargetReactiveSite:
    """Compact definition-derived reactive-site hypothesis."""

    hypothesis_id: str
    site_type: str
    topology: str
    component_index: int
    atom_indices: Tuple[int, ...]
    bond_indices: Tuple[int, ...]
    chemist_label: str
    availability: str
    warnings: Tuple[str, ...]


@dataclass(frozen=True)
class TargetAudit:
    """Validated target identity, stereo, motifs, and reactive sites."""

    input_smiles: str
    canonical_smiles: Optional[str]
    valid: bool
    component_count: int = 0
    atom_count: int = 0
    heavy_atom_count: int = 0
    formal_charge: int = 0
    stereocenters: Tuple[TargetStereocenter, ...] = ()
    double_bond_stereo: Tuple[TargetDoubleBondStereo, ...] = ()
    motifs: Tuple[TargetMotif, ...] = ()
    reactive_sites: Tuple[TargetReactiveSite, ...] = ()
    warnings: Tuple[str, ...] = ()
    error: Optional[str] = None
    molecular_formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    canonical_atom_to_input_atom: Tuple[int, ...] = ()
    input_atom_map_numbers: Tuple[int, ...] = ()
    components: Tuple[TargetComponent, ...] = ()
    rings: Tuple[Tuple[int, ...], ...] = ()
    carbonyl_nitrogen_observations: Tuple[CarbonylNitrogenObservation, ...] = ()
    definition_versions: Tuple[Tuple[str, str], ...] = ()
    rdkit_version: str = ""
    atom_index_scope: str = "zero_based_canonical_smiles"
    motif_site_index_scope: str = "component_local_in_canonical_smiles"
    coverage: str = "registered_patterns_only_not_exhaustive"
    limitations: Tuple[str, ...] = (
        "No motif match is not evidence of absent functionality or inertness.",
        "Motifs and sites are observations/hypotheses, not experimental feasibility or selectivity.",
        "Tautomers, protonation states and salt components are not merged or stripped.",
        "Rings are the RDKit perceived ring basis, not an exhaustive set of cycles.",
        "Use components to project motif/site indices; all other atom indices use canonical SMILES.",
    )
    schema_version: str = TARGET_AUDIT_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible audit without RDKit objects."""

        return asdict(self)


def audit_target(smiles: str) -> TargetAudit:
    """Audit one target without proposing or changing molecular structure."""

    structure = observe_molecular_structure(smiles)
    if not structure.valid or not structure.canonical_smiles:
        return TargetAudit(
            input_smiles=structure.input_smiles,
            canonical_smiles=structure.canonical_smiles,
            valid=False,
            warnings=tuple(structure.warnings),
            error=structure.error,
        )
    # Use the same map-free indexing contract as focused inspection/disconnection.
    canonical, input_order = indexed_product(smiles)
    analysis = analyze_molecule(canonical, include_context_features=False)
    molecule = parse_smiles(canonical)
    if molecule is None:
        return TargetAudit(
            input_smiles=structure.input_smiles,
            canonical_smiles=structure.canonical_smiles,
            valid=False,
            warnings=tuple(structure.warnings),
            error="CANONICAL_SMILES_REPARSE_FAILED",
        )

    from rdkit import Chem, rdBase  # type: ignore
    from rdkit.Chem import Descriptors, rdMolDescriptors  # type: ignore

    original = parse_smiles(smiles)
    component_atoms: list[tuple[int, ...]] = []
    fragments = Chem.GetMolFrags(molecule, asMols=True, fragsMolAtomMapping=component_atoms)
    components = tuple(
        TargetComponent(
            index, tuple(atoms), tuple(
                molecule.GetBondBetweenAtoms(
                    atoms[b.GetBeginAtomIdx()], atoms[b.GetEndAtomIdx()]
                ).GetIdx() for b in fragment.GetBonds()
            ),
        ) for index, (fragment, atoms) in enumerate(zip(fragments, component_atoms))
    )
    carbonyl_nitrogen = tuple(
        CarbonylNitrogenObservation(
            carbon.GetIdx(), oxygen.GetIdx(), nitrogen.GetIdx(),
            molecule.GetBondBetweenAtoms(carbon.GetIdx(), nitrogen.GetIdx()).IsInRing(),
            carbon.GetIsAromatic(), nitrogen.GetIsAromatic(),
        )
        for carbon in molecule.GetAtoms() if carbon.GetAtomicNum() == 6
        for oxygen in carbon.GetNeighbors()
        if oxygen.GetAtomicNum() == 8 and molecule.GetBondBetweenAtoms(
            carbon.GetIdx(), oxygen.GetIdx()
        ).GetBondType() == Chem.BondType.DOUBLE
        for nitrogen in carbon.GetNeighbors() if nitrogen.GetAtomicNum() == 7
    )

    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    stereocenters = tuple(
        TargetStereocenter(
            atom_index=int(atom_index),
            element=str(molecule.GetAtomWithIdx(atom_index).GetSymbol()),
            assignment=str(assignment),
            assigned=str(assignment) != "?",
        )
        for atom_index, assignment in Chem.FindMolChiralCenters(
            molecule,
            includeUnassigned=True,
            includeCIP=True,
        )
    )
    double_bond_stereo = tuple(
        TargetDoubleBondStereo(
            bond_index=int(bond.GetIdx()),
            atom_1_index=int(bond.GetBeginAtomIdx()),
            atom_2_index=int(bond.GetEndAtomIdx()),
            assignment=str(bond.GetStereo()),
        )
        for bond in molecule.GetBonds()
        if str(bond.GetStereo()) not in {"STEREONONE", "STEREOANY"}
    )
    motifs = tuple(
        TargetMotif(
            motif_id=item.motif_id,
            chemist_label=item.chemist_label,
            component_index=item.component_index,
            atom_indices=tuple(item.atom_indices),
            tags=tuple(item.tags),
        )
        for item in analysis.motifs
    )
    reactive_sites = tuple(
        TargetReactiveSite(
            hypothesis_id=item.hypothesis_id,
            site_type=item.site_type,
            topology=item.topology,
            component_index=item.component_index,
            atom_indices=tuple(item.atom_indices),
            bond_indices=tuple(item.bond_indices),
            chemist_label=item.chemist_label,
            availability=item.availability,
            warnings=tuple(item.warnings),
        )
        for item in analysis.reactive_site_hypotheses
    )
    return TargetAudit(
        input_smiles=structure.input_smiles,
        canonical_smiles=canonical,
        valid=True,
        component_count=len(structure.components),
        atom_count=int(molecule.GetNumAtoms()),
        heavy_atom_count=int(molecule.GetNumHeavyAtoms()),
        formal_charge=sum(
            int(atom.GetFormalCharge()) for atom in molecule.GetAtoms()
        ),
        stereocenters=stereocenters,
        double_bond_stereo=double_bond_stereo,
        motifs=motifs,
        reactive_sites=reactive_sites,
        warnings=tuple(structure.warnings),
        molecular_formula=rdMolDescriptors.CalcMolFormula(molecule),
        molecular_weight=round(Descriptors.MolWt(molecule), 4),
        canonical_atom_to_input_atom=input_order,
        input_atom_map_numbers=tuple(a.GetAtomMapNum() for a in original.GetAtoms()),
        components=components,
        rings=tuple(sorted(tuple(sorted(ring)) for ring in molecule.GetRingInfo().AtomRings())),
        carbonyl_nitrogen_observations=carbonyl_nitrogen,
        definition_versions=descriptor_definition_versions(),
        rdkit_version=rdBase.rdkitVersion,
    )


__all__ = [
    "TARGET_AUDIT_SCHEMA_VERSION",
    "TargetAudit",
    "TargetComponent",
    "CarbonylNitrogenObservation",
    "TargetDoubleBondStereo",
    "TargetMotif",
    "TargetReactiveSite",
    "TargetStereocenter",
    "audit_target",
]

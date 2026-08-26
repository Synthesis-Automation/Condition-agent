"""Compact deterministic target audit for chemistry-tool consumers."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Optional, Tuple

from .api import analyze_molecule
from .chemistry.rdkit_utils import parse_smiles


TARGET_AUDIT_SCHEMA_VERSION = "target_audit.v1"


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
    schema_version: str = TARGET_AUDIT_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible audit without RDKit objects."""

        return asdict(self)


def audit_target(smiles: str) -> TargetAudit:
    """Audit one target without proposing or changing molecular structure."""

    analysis = analyze_molecule(smiles, include_context_features=False)
    structure = analysis.structure
    if not structure.valid or not structure.canonical_smiles:
        return TargetAudit(
            input_smiles=structure.input_smiles,
            canonical_smiles=structure.canonical_smiles,
            valid=False,
            warnings=tuple(structure.warnings),
            error=structure.error,
        )
    molecule = parse_smiles(structure.canonical_smiles)
    if molecule is None:
        return TargetAudit(
            input_smiles=structure.input_smiles,
            canonical_smiles=structure.canonical_smiles,
            valid=False,
            warnings=tuple(structure.warnings),
            error="CANONICAL_SMILES_REPARSE_FAILED",
        )

    from rdkit import Chem  # type: ignore

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
        canonical_smiles=structure.canonical_smiles,
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
    )


__all__ = [
    "TARGET_AUDIT_SCHEMA_VERSION",
    "TargetAudit",
    "TargetDoubleBondStereo",
    "TargetMotif",
    "TargetReactiveSite",
    "TargetStereocenter",
    "audit_target",
]

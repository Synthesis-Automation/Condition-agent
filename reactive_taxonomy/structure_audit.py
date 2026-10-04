"""Deterministic graph checks for supplied structures, without source-name inference."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors

MAX_STRUCTURE_CHARACTERS = 4000

@dataclass(frozen=True)
class StructureAudit:
    """Graph validity and descriptors; parsing never verifies a literature assignment."""

    input_smiles: str
    valid: bool
    canonical_smiles: str | None
    formula: str | None
    formal_charge: int | None
    component_count: int | None
    ring_sizes: tuple[int, ...]
    unspecified_stereo_count: int
    warnings: tuple[str, ...]
    rdkit_version: str
    schema_version: str = "structure_audit.v1"

    def to_dict(self) -> dict[str, Any]:
        """Serialize observations without changing their scope."""
        return asdict(self)


def audit_structure(smiles: str) -> StructureAudit:
    """Check an explicit graph, retaining salts, stereo, charges and invalid input.

    Atom-map labels are removed only for canonical identity; no reaction atom
    correspondence is inferred. Formula/ring checks cannot verify compound names,
    tautomer assignment, material form, or experimental feasibility.
    """
    if not isinstance(smiles, str) or not smiles.strip() or len(smiles) > MAX_STRUCTURE_CHARACTERS:
        raise ValueError("Supply SMILES containing 1 to 4000 characters")
    with rdBase.BlockLogs():
        molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return StructureAudit(smiles, False, None, None, None, None, (), 0,
                              ("INVALID_SMILES",), rdBase.rdkitVersion)
    warnings = []
    if any(atom.GetAtomicNum() == 0 for atom in molecule.GetAtoms()):
        warnings.append("UNRESOLVED_WILDCARD_STRUCTURE")
    if any(atom.GetNumRadicalElectrons() for atom in molecule.GetAtoms()):
        warnings.append("RADICAL_STRUCTURE")
    if any(atom.GetAtomMapNum() for atom in molecule.GetAtoms()):
        warnings.append("ATOM_MAP_LABELS_NOT_REACTION_CORRESPONDENCE")
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(0)
    stereo_count = sum(str(item.specified) != "Specified" for item in Chem.FindPotentialStereo(molecule))
    if stereo_count:
        warnings.append("UNSPECIFIED_STEREOCHEMISTRY")
    return StructureAudit(
        smiles, not any(item in warnings for item in ("UNRESOLVED_WILDCARD_STRUCTURE", "RADICAL_STRUCTURE")),
        Chem.MolToSmiles(molecule, isomericSmiles=True),
        rdMolDescriptors.CalcMolFormula(molecule, separateIsotopes=True, abbreviateHIsotopes=False),
        Chem.GetFormalCharge(molecule), len(Chem.GetMolFrags(molecule)),
        tuple(sorted(len(ring) for ring in molecule.GetRingInfo().AtomRings())),
        stereo_count, tuple(warnings), rdBase.rdkitVersion,
    )

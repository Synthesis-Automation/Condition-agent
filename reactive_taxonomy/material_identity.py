"""Exact material identity and molecular weight without chemical-form relaxation."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors

from .fragment_search import indexed_product


@dataclass(frozen=True)
class MaterialIdentity:
    """A complete chemical form; atom maps are metadata, not identity."""

    input_smiles: str
    canonical_smiles: str
    molecular_weight: float
    component_count: int
    warnings: tuple[str, ...]
    schema_version: str = "material_identity.v1"
    identity_method: str = "canonical_isomeric_smiles_without_maps.v1"
    molecular_weight_method: str = "rdkit.Chem.Descriptors.MolWt"
    molecular_weight_unit: str = "g/mol"
    rdkit_version: str = rdBase.rdkitVersion

    def to_dict(self) -> dict[str, Any]:
        """Serialize the identity and descriptor provenance."""
        return asdict(self)


def identify_material(smiles: str) -> MaterialIdentity:
    """Preserve stereo, isotopes, charge and all components; reject query atoms.

    No salt stripping, neutralization or tautomer enumeration is performed.
    Closed-shell molecular forms are supported; radicals require explicit review.
    """
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("Provide a nonempty material SMILES")
    canonical, _ = indexed_product(smiles)
    molecule = Chem.MolFromSmiles(canonical)
    if any(atom.GetAtomicNum() == 0 or atom.GetNumRadicalElectrons()
           for atom in molecule.GetAtoms()):
        raise ValueError("Material lookup requires a concrete closed-shell structure")
    components = len(Chem.GetMolFrags(molecule))
    warnings = []
    if components > 1:
        warnings.append("All components are retained; molecular weight includes the entire chemical form.")
    if any(info.specified == Chem.StereoSpecified.Unspecified
           for info in Chem.FindPotentialStereo(molecule)):
        warnings.append("Unspecified stereochemistry is retained; no stereoisomer identity is assumed.")
    return MaterialIdentity(
        smiles, canonical, float(Descriptors.MolWt(molecule)), components, tuple(warnings),
    )

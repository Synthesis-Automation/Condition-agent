"""Focused graph inspection with canonical atom IDs and explicit stereo gaps.

This module reuses taxonomy detectors and descriptors. A detected site is a
hypothesis, not a prediction of reactivity, selectivity or experimental success.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Any

from rdkit import Chem, rdBase

from .api import detect_reactive_site_hypotheses
from .environments import build_site_environment
from .fragment_search import indexed_product
from .models import MolecularMotifMatch, ReactiveSiteEnvironment, ReactiveSiteHypothesis
from .molecular_motifs import detect_molecular_motifs


def molecule_inspection_policy() -> dict[str, Any]:
    """Load validated bounds and the supported strict comparison semantics."""
    path = Path(__file__).with_name("definitions") / "molecule_inspection.v1.json"
    policy = json.loads(path.read_text("utf-8"))
    if (policy.get("schema_version") != "molecule_inspection_policy.v1"
            or not isinstance(policy.get("definition_version"), str)
            or not policy["definition_version"]):
        raise ValueError("Invalid molecule inspection policy version")
    for field in ("max_characters", "max_atoms", "default_timeout_seconds",
                  "max_timeout_seconds", "max_embeddings", "max_alignments", "max_radius"):
        if type(policy.get(field)) is not int or policy[field] < 1:
            raise ValueError(f"Invalid molecule inspection bound: {field}")
    if policy["default_timeout_seconds"] > policy["max_timeout_seconds"]:
        raise ValueError("Default timeout exceeds inspection bound")
    if (policy.get("atom_match_fields") != ["element", "isotope", "formal_charge", "aromatic", "ring_degree"]
            or policy.get("bond_match") != "exact_order_and_ring_membership"
            or policy.get("complete_rings_only") is not True
            or policy.get("stereo_in_automatic_core") is not False):
        raise ValueError("Unsupported molecular comparison semantics")
    return policy


@dataclass(frozen=True)
class InspectionAtom:
    """Canonical SMILES atom ID, linked to the parsed input atom order."""

    atom_id: int
    input_atom_index: int
    element: str
    isotope: int
    formal_charge: int
    aromatic: bool
    hydrogen_count: int
    hybridization: str
    ring_sizes: tuple[int, ...]


@dataclass(frozen=True)
class StereoFeature:
    """Potential stereo from RDKit; no mechanistic retention/inversion claim."""

    stereo_type: str
    atom_ids: tuple[int, ...]
    bond_id: int | None
    specification: str
    descriptor: str
    cip_label: str | None


@dataclass(frozen=True)
class InspectedMolecule:
    """Unchanged chemical identity with explicit atom indexing and stereo scope."""

    input_smiles: str
    canonical_smiles: str
    atoms: tuple[InspectionAtom, ...]
    stereochemistry: tuple[StereoFeature, ...]
    warnings: tuple[str, ...]
    atom_id_scope: str = "zero_based_canonical_smiles"


def _prepare_molecule(smiles: str, policy: dict[str, Any]) -> tuple[Any, InspectedMolecule]:
    if not isinstance(smiles, str) or not smiles.strip() or len(smiles) > policy["max_characters"]:
        raise ValueError(f"Provide SMILES of 1..{policy['max_characters']} characters")
    canonical, order = indexed_product(smiles)
    mol = Chem.MolFromSmiles(canonical)
    if len(Chem.GetMolFrags(mol)) != 1 or mol.GetNumAtoms() > policy["max_atoms"]:
        raise ValueError(f"Select one connected molecule with at most {policy['max_atoms']} atoms; salts are not stripped")
    if any(a.GetAtomicNum() == 0 or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        raise ValueError("Wildcard and radical molecules are unsupported")
    atoms = tuple(InspectionAtom(
        a.GetIdx(), order[a.GetIdx()], a.GetSymbol(), a.GetIsotope(), a.GetFormalCharge(),
        a.GetIsAromatic(), a.GetTotalNumHs(), str(a.GetHybridization()),
        tuple(sorted(len(ring) for ring in mol.GetRingInfo().AtomRings() if a.GetIdx() in ring)),
    ) for a in mol.GetAtoms())
    stereo = []
    for info in Chem.FindPotentialStereo(mol):
        is_atom = str(info.type).startswith("Atom_")
        if is_atom:
            atom = mol.GetAtomWithIdx(info.centeredOn)
            ids, bond_id = (info.centeredOn,), None
            cip = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None
        else:
            bond = mol.GetBondWithIdx(info.centeredOn)
            ids, bond_id = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()), info.centeredOn
            cip = None
        stereo.append(StereoFeature(str(info.type), ids, bond_id, str(info.specified), str(info.descriptor), cip))
    warnings = []
    if any(s.specification != "Specified" for s in stereo):
        warnings.append("UNSPECIFIED_STEREOCHEMISTRY")
    original = Chem.MolFromSmiles(smiles)
    if any(a.GetAtomMapNum() for a in original.GetAtoms()):
        warnings.append("ATOM_MAP_LABELS_IGNORED_NOT_REACTION_CORRESPONDENCE")
    return mol, InspectedMolecule(smiles, canonical, atoms, tuple(stereo), tuple(warnings))


@dataclass(frozen=True)
class SiteInspection:
    """Focused observations plus alternative detected sites, without ranking."""

    molecule: InspectedMolecule
    selected_atom_ids: tuple[int, ...]
    neighborhood_atom_ids: tuple[int, ...]
    radius: int
    sites: tuple[ReactiveSiteHypothesis, ...]
    environments: tuple[ReactiveSiteEnvironment, ...]
    motifs: tuple[MolecularMotifMatch, ...]
    other_sites: tuple[ReactiveSiteHypothesis, ...]
    warnings: tuple[str, ...]
    definition_version: str
    rdkit_version: str
    limitations: tuple[str, ...] = (
        "Detected sites are hypotheses; unrecognized sites may also react.",
        "Steric/electronic descriptors are graph heuristics, not rates or condition-dependent selectivity.",
        "Atom and bond IDs refer to the returned canonical SMILES; input_atom_index links to parsed input.",
    )
    schema_version: str = "site_inspection.v1"

    def to_dict(self) -> dict[str, Any]:
        """Return serializable observations and descriptor provenance."""
        return asdict(self)


def inspect_reactive_sites(
    smiles: str, selected_atom_ids: list[int] | None = None, radius: int = 2,
) -> SiteInspection:
    """Inspect chosen canonical atoms and nearby sites, or all sites if omitted.

    Selections are never interpreted as reaction centers or bonds to break. Other
    sites remain visible without asserting they compete under any conditions.
    """
    policy = molecule_inspection_policy()
    if type(radius) is not int or not 0 <= radius <= policy["max_radius"]:
        raise ValueError(f"radius must be 0..{policy['max_radius']}")
    mol, molecule = _prepare_molecule(smiles, policy)
    if selected_atom_ids is not None and (
        not isinstance(selected_atom_ids, list) or not selected_atom_ids
        or any(type(i) is not int or not 0 <= i < mol.GetNumAtoms() for i in selected_atom_ids)
        or len(set(selected_atom_ids)) != len(selected_atom_ids)
    ):
        raise ValueError("selected_atom_ids must be unique zero-based canonical atom IDs")
    selected = set(selected_atom_ids if selected_atom_ids is not None else range(mol.GetNumAtoms()))
    neighborhood = set(selected)
    for _ in range(radius):
        neighborhood |= {n.GetIdx() for i in neighborhood for n in mol.GetAtomWithIdx(i).GetNeighbors()}
    hypotheses = detect_reactive_site_hypotheses(mol)
    motifs = tuple(detect_molecular_motifs(mol, 0))
    sites = tuple(s for s in hypotheses if neighborhood.intersection(s.atom_indices))
    environments = tuple(build_site_environment(mol, s, motifs) for s in sites)
    warnings = list(molecule.warnings)
    if not sites:
        warnings.append("NO_DETECTED_SITE_IN_REGION_NOT_PROOF_OF_INERTNESS")
    return SiteInspection(
        molecule, tuple(sorted(selected)), tuple(sorted(neighborhood)), radius, sites, environments,
        tuple(m for m in motifs if neighborhood.intersection(m.atom_indices)),
        tuple(s for s in hypotheses if s not in sites), tuple(warnings),
        policy["definition_version"], rdBase.rdkitVersion,
    )

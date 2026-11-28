"""Inspect the oxime geometry required by a mapped cyclic ketoxime-to-lactam edit."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import json
from pathlib import Path
from typing import Literal

from rdkit import Chem

from reactive_taxonomy.chemistry.reaction_smiles import split_reaction_fields
from reactive_taxonomy.chemistry.smarts_cache import compile_smarts


_DEFINITION_ID = "cyclic_oxime_migration.v1"
GeometryStatus = Literal["unspecified", "anti_aligned", "syn_aligned"]


@dataclass(frozen=True)
class MigrationGeometryReview:
    """Compare the drawn migrating bond with the oxime O across C=N."""

    status: GeometryStatus
    summary: str
    alternative_product_smiles: str
    audit_id: str = _DEFINITION_ID

    @property
    def requires_caution(self) -> bool:
        """Return whether the drawn oxime leaves the intended migration in question."""
        return self.status != "anti_aligned"


@lru_cache(maxsize=1)
def _patterns() -> tuple[Chem.Mol, Chem.Mol]:
    path = Path(__file__).with_name("definitions") / f"{_DEFINITION_ID}.json"
    rule = json.loads(path.read_text(encoding="utf-8"))
    if rule["definition_id"] != _DEFINITION_ID or rule["schema_version"] != "1.0":
        raise ValueError("Unsupported cyclic oxime migration definition")
    patterns = []
    for key, roles in (("reactant_smarts", [4, 1, 2, 3, 5]), ("product_smarts", [5, 1, 6, 2, 4])):
        pattern = compile_smarts(rule[key], validate=True)
        if [a.GetAtomMapNum() for a in pattern.GetAtoms()] != roles:
            raise ValueError(f"Invalid migration atom roles in {key}")
        patterns.append(pattern)
    return patterns[0], patterns[1]


def _alternative(product: Chem.Mol, site: tuple[int, ...]) -> str | None:
    retained, carbon, _, nitrogen, migrant = site
    # Replacing a bond at a specified tetrahedral center needs an explicit
    # stereo projection. This audit only interchanges achiral endpoints.
    endpoints = [product.GetAtomWithIdx(i) for i in (retained, migrant)]
    if any(atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for atom in endpoints):
        return None
    if any(bond.GetStereo() != Chem.BondStereo.STEREONONE for atom in endpoints for bond in atom.GetBonds()):
        return None
    if product.GetBondBetweenAtoms(carbon, migrant) or product.GetBondBetweenAtoms(nitrogen, retained):
        return None
    editable = Chem.RWMol(product)
    editable.RemoveBond(nitrogen, migrant)
    editable.RemoveBond(carbon, retained)
    editable.AddBond(carbon, migrant, Chem.BondType.SINGLE)
    editable.AddBond(nitrogen, retained, Chem.BondType.SINGLE)
    if Chem.SanitizeMol(editable, catchErrors=True) != Chem.SanitizeFlags.SANITIZE_NONE:
        return None
    # Equivalent migrations in a symmetric ring do not give regioisomers.
    if Chem.MolToSmiles(editable, isomericSmiles=False) == Chem.MolToSmiles(product, isomericSmiles=False):
        return None
    return Chem.MolToSmiles(editable, isomericSmiles=True)


def _geometry(bond: Chem.Bond, carbon: int, migrant: int, oxygen: int) -> GeometryStatus:
    stereo = bond.GetStereo()
    atoms = list(bond.GetStereoAtoms())
    if len(atoms) != 2 or stereo not in (
        Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOCIS, Chem.BondStereo.STEREOTRANS,
    ):
        return "unspecified"
    if bond.GetBeginAtomIdx() != carbon:
        atoms.reverse()
    if atoms[1] != oxygen:
        return "unspecified"
    opposite = stereo in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOTRANS)
    return "anti_aligned" if opposite == (atoms[0] == migrant) else "syn_aligned"


@lru_cache(maxsize=2048)
def assess_migration_geometry(reaction_smiles: str) -> MigrationGeometryReview | None:
    """Report a local geometry observation for one verified cyclic N-insertion.

    Atom maps identify the original C=N and which ring carbon moved to N.
    Oxygen correspondence is not required. The other migration must give
    a distinct sanitized constitutional isomer. The audit assigns no yield,
    reaction-condition preference, or experimental product distribution.
    """
    sides = split_reaction_fields(reaction_smiles)
    if len(sides) != 3 or not sides[0] or not sides[2] or "." in sides[2]:
        return None
    reactant = Chem.MolFromSmiles(sides[0])
    precursor_pattern, lactam_pattern = _patterns()
    if reactant is None or not (oximes := reactant.GetSubstructMatches(precursor_pattern)):
        return None
    product = Chem.MolFromSmiles(sides[2])
    if product is None or not (lactams := product.GetSubstructMatches(lactam_pattern)):
        return None
    maps = []
    for molecule in (reactant, product):
        numbers = [a.GetAtomMapNum() for a in molecule.GetAtoms()]
        assigned = [n for n in numbers if n]
        if len(assigned) != len(set(assigned)):
            return None
        maps.append(numbers)
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(0)
        Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    left, right = maps
    matches = []
    for first, carbon, nitrogen, oxygen, second in oximes:
        neighbors = {left[first]: first, left[second]: second}
        if not all((left[carbon], left[nitrogen], left[first], left[second])):
            continue
        for site in lactams:
            retained, acyl, _, amide, migrant = site
            if (left[carbon], left[nitrogen]) != (right[acyl], right[amide]):
                continue
            if {right[retained], right[migrant]} != neighbors.keys():
                continue
            alternative = _alternative(product, site)
            if alternative:
                bond = reactant.GetBondBetweenAtoms(carbon, nitrogen)
                status = _geometry(bond, carbon, neighbors[right[migrant]], oxygen)
                matches.append((status, alternative))
    if len(matches) != 1:
        return None
    status, alternative = matches[0]
    summaries = {
        "unspecified": "The cyclic oxime has unspecified C=N geometry. The other migration gives a distinct lactam regioisomer.",
        "anti_aligned": "The drawn oxime O is anti to the ring carbon that migrates in this proposal. Experimental regioselectivity remains open.",
        "syn_aligned": "The drawn oxime O is syn to the ring carbon that migrates in this proposal. Conventional anti migration would give the alternative lactam.",
    }
    return MigrationGeometryReview(status, summaries[status], alternative)

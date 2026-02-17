"""
Reaction SMILES Analysis - Core Deterministic Functions

Simple POC implementation for analyzing reaction SMILES using:
- RDKit for canonicalization and spectator detection
- rxnmapper for atom mapping
- Bond-change extraction for reaction center identification

Design: docs/reaction_smiles_analysis_agent_simple_v1.md
"""

from dataclasses import dataclass, field
from typing import List, Optional, Dict, Set, Tuple
from rdkit import Chem
from rdkit.Chem import rdChemReactions
import logging

# Import automatic interpretation
from chemtools.reaction_interpreter import interpret_reaction_pattern, format_interpretation_report

logger = logging.getLogger(__name__)

# Common spectators (simple ions and salts)
SPECTATORS = {
    "Cl", "Br", "I", "F",
    "[Cl-]", "[Br-]", "[I-]", "[F-]",
    "[Na+]", "[K+]", "[Li+]", "[Cs+]",
    "[OH-]", "[NH4+]",
    "O",  # Water as O in SMILES
}


@dataclass
class CleanReport:
    """Report from cleaning and canonicalizing reaction SMILES."""
    rxn_smiles_raw: str
    rxn_smiles_clean: str
    reactants_clean: List[str] = field(default_factory=list)
    products_clean: List[str] = field(default_factory=list)
    spectators: List[str] = field(default_factory=list)
    parse_warnings: List[str] = field(default_factory=list)
    standardization_actions: List[str] = field(default_factory=list)


@dataclass
class MappingReport:
    """Report from atom mapping."""
    mapped_rxn_smiles: Optional[str] = None
    mapping_qc: Dict = field(default_factory=lambda: {"ok": False, "notes": []})


@dataclass
class BondChange:
    """Single bond change event."""
    id: str
    change: str  # "broken", "formed", "order_change"
    a1: int  # atom map number 1
    a2: int  # atom map number 2
    bond_type: str  # "single", "double", "triple", "aromatic"

    def to_dict(self) -> Dict:
        return {
            "id": self.id,
            "change": self.change,
            "a1": self.a1,
            "a2": self.a2,
            "bond": self.bond_type
        }


@dataclass
class BondChangeReport:
    """Report from bond change extraction."""
    bond_changes: List[BondChange] = field(default_factory=list)
    reaction_center_atoms: List[int] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


def _canonicalize_smiles(smiles: str) -> Optional[str]:
    """Canonicalize a SMILES string using RDKit."""
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception as e:
        logger.debug(f"Failed to canonicalize '{smiles}': {e}")
        return None


def clean_reaction_smiles(rxn_smiles: str, drop_spectators: bool = True) -> CleanReport:
    """
    Clean and canonicalize reaction SMILES.

    Args:
        rxn_smiles: Raw reaction SMILES (reactants>>products)
        drop_spectators: If True, remove simple ions/salts from mapping

    Returns:
        CleanReport with canonicalized SMILES and detected spectators

    Raises:
        ValueError: If reaction SMILES is invalid
    """
    warnings = []
    actions = []
    spectators = []

    # Validate format
    if rxn_smiles.count(">>") != 1:
        raise ValueError("Invalid reaction SMILES: expected exactly one '>>'.")

    lhs, rhs = rxn_smiles.split(">>")
    lhs_parts = [p.strip() for p in lhs.split(".") if p.strip()]
    rhs_parts = [p.strip() for p in rhs.split(".") if p.strip()]

    if not lhs_parts or not rhs_parts:
        raise ValueError("Invalid reaction SMILES: empty reactants or products.")

    def process_side(parts: List[str], side_name: str) -> List[str]:
        """Process one side of reaction, canonicalizing and optionally dropping spectators."""
        cleaned = []
        for part in parts:
            canon = _canonicalize_smiles(part)
            if canon is None:
                warnings.append(f"{side_name}: failed to parse '{part}'")
                continue

            # Check if spectator
            if drop_spectators and canon in SPECTATORS:
                spectators.append(canon)
                actions.append(f"Dropped spectator: {canon}")
                continue

            cleaned.append(canon)
        return cleaned

    reactants_clean = process_side(lhs_parts, "reactants")
    products_clean = process_side(rhs_parts, "products")

    if not reactants_clean or not products_clean:
        warnings.append("All reactants or products were dropped/failed to parse")

    rxn_clean = ".".join(reactants_clean) + ">>" + ".".join(products_clean)

    return CleanReport(
        rxn_smiles_raw=rxn_smiles,
        rxn_smiles_clean=rxn_clean,
        reactants_clean=reactants_clean,
        products_clean=products_clean,
        spectators=sorted(set(spectators)),
        parse_warnings=warnings,
        standardization_actions=actions
    )


def map_reaction(rxn_smiles_clean: str) -> MappingReport:
    """
    Perform atom mapping using rxnmapper.

    Args:
        rxn_smiles_clean: Cleaned reaction SMILES

    Returns:
        MappingReport with mapped SMILES and QC info
    """
    try:
        from rxnmapper import RXNMapper

        mapper = RXNMapper()
        results = mapper.get_attention_guided_atom_maps([rxn_smiles_clean])

        if not results or not results[0]:
            return MappingReport(
                mapped_rxn_smiles=None,
                mapping_qc={"ok": False, "notes": ["rxnmapper returned empty result"]}
            )

        mapped_smiles = results[0]["mapped_rxn"]
        confidence = results[0].get("confidence", 0.0)

        # Simple QC: check if mapping actually added atom map numbers
        if ":" not in mapped_smiles:
            return MappingReport(
                mapped_rxn_smiles=mapped_smiles,
                mapping_qc={
                    "ok": False,
                    "notes": ["No atom map numbers found in output"],
                    "confidence": confidence
                }
            )

        # QC: check confidence threshold
        qc_ok = confidence >= 0.5
        qc_notes = []
        if confidence < 0.5:
            qc_notes.append(f"Low mapping confidence: {confidence:.2f}")

        return MappingReport(
            mapped_rxn_smiles=mapped_smiles,
            mapping_qc={
                "ok": qc_ok,
                "confidence": confidence,
                "notes": qc_notes
            }
        )

    except ImportError:
        logger.warning("rxnmapper not available - mapping skipped")
        return MappingReport(
            mapped_rxn_smiles=None,
            mapping_qc={"ok": False, "notes": ["rxnmapper not installed"]}
        )
    except Exception as e:
        logger.exception(f"Mapping failed: {e}")
        return MappingReport(
            mapped_rxn_smiles=None,
            mapping_qc={"ok": False, "notes": [f"Mapping error: {str(e)}"]}
        )


def extract_bond_changes(mapped_rxn_smiles: str) -> BondChangeReport:
    """
    Extract bond changes from mapped reaction SMILES.

    Args:
        mapped_rxn_smiles: Atom-mapped reaction SMILES

    Returns:
        BondChangeReport with bond changes and reaction center atoms
    """
    warnings = []

    try:
        # Parse mapped reaction
        rxn = rdChemReactions.ReactionFromSmiles(mapped_rxn_smiles)
        if rxn is None:
            warnings.append("Failed to parse mapped reaction")
            return BondChangeReport(warnings=warnings)

        # Get reactants and products
        reactant_mols = [mol for mol in rxn.GetReactants()]
        product_mols = [mol for mol in rxn.GetProducts()]

        if not reactant_mols or not product_mols:
            warnings.append("No reactants or products in mapped reaction")
            return BondChangeReport(warnings=warnings)

        # Build bond sets by atom map numbers
        reactant_bonds = _extract_bonds_by_mapnum(reactant_mols)
        product_bonds = _extract_bonds_by_mapnum(product_mols)

        # Find differences
        bond_changes = []
        change_id = 1

        # Broken bonds (in reactants but not products)
        for bond_key, bond_type in reactant_bonds.items():
            if bond_key not in product_bonds:
                a1, a2 = bond_key
                bond_changes.append(BondChange(
                    id=f"BC{change_id}",
                    change="broken",
                    a1=min(a1, a2),
                    a2=max(a1, a2),
                    bond_type=bond_type
                ))
                change_id += 1

        # Formed bonds (in products but not reactants)
        for bond_key, bond_type in product_bonds.items():
            if bond_key not in reactant_bonds:
                a1, a2 = bond_key
                bond_changes.append(BondChange(
                    id=f"BC{change_id}",
                    change="formed",
                    a1=min(a1, a2),
                    a2=max(a1, a2),
                    bond_type=bond_type
                ))
                change_id += 1

        # Bond order changes (same bond, different type)
        common_bonds = set(reactant_bonds.keys()) & set(product_bonds.keys())
        for bond_key in common_bonds:
            if reactant_bonds[bond_key] != product_bonds[bond_key]:
                a1, a2 = bond_key
                bond_changes.append(BondChange(
                    id=f"BC{change_id}",
                    change="order_change",
                    a1=min(a1, a2),
                    a2=max(a1, a2),
                    bond_type=f"{reactant_bonds[bond_key]}->{product_bonds[bond_key]}"
                ))
                change_id += 1

        # Extract reaction center atoms (all atoms involved in changes)
        reaction_center = set()
        for bc in bond_changes:
            reaction_center.add(bc.a1)
            reaction_center.add(bc.a2)

        return BondChangeReport(
            bond_changes=bond_changes,
            reaction_center_atoms=sorted(reaction_center),
            warnings=warnings
        )

    except Exception as e:
        logger.exception(f"Bond change extraction failed: {e}")
        warnings.append(f"Bond change extraction error: {str(e)}")
        return BondChangeReport(warnings=warnings)


def _extract_bonds_by_mapnum(mols: List) -> Dict[Tuple[int, int], str]:
    """
    Extract bonds indexed by atom map numbers from molecule list.

    Args:
        mols: List of RDKit Mol objects with atom mapping

    Returns:
        Dict mapping (mapnum1, mapnum2) -> bond_type_str
    """
    bonds = {}

    bond_type_map = {
        Chem.BondType.SINGLE: "single",
        Chem.BondType.DOUBLE: "double",
        Chem.BondType.TRIPLE: "triple",
        Chem.BondType.AROMATIC: "aromatic",
    }

    for mol in mols:
        if mol is None:
            continue

        for bond in mol.GetBonds():
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()

            map1 = atom1.GetAtomMapNum()
            map2 = atom2.GetAtomMapNum()

            # Skip unmapped atoms
            if map1 == 0 or map2 == 0:
                continue

            # Canonicalize bond key (smaller map num first)
            bond_key = (min(map1, map2), max(map1, map2))
            bond_type = bond_type_map.get(bond.GetBondType(), "unknown")

            bonds[bond_key] = bond_type

    return bonds


def analyze_reaction_smiles(
    rxn_smiles: str,
    drop_spectators: bool = True,
    skip_mapping: bool = False
) -> Dict:
    """
    Full deterministic analysis pipeline for reaction SMILES.

    Args:
        rxn_smiles: Raw reaction SMILES
        drop_spectators: Remove common ions/salts
        skip_mapping: Skip atom mapping and bond changes (faster)

    Returns:
        Dict with cleaned reaction, mapping, and bond changes
    """
    # Step 1: Clean
    clean = clean_reaction_smiles(rxn_smiles, drop_spectators=drop_spectators)

    result = {
        "input": {
            "rxn_smiles_raw": clean.rxn_smiles_raw,
            "rxn_smiles_clean": clean.rxn_smiles_clean,
            "reactants_clean": clean.reactants_clean,
            "products_clean": clean.products_clean,
            "spectators": clean.spectators,
            "parse_warnings": clean.parse_warnings,
            "standardization_actions": clean.standardization_actions,
        }
    }

    if skip_mapping:
        result["tool_facts"] = {
            "mapped_rxn_smiles": None,
            "mapping_qc": {"ok": False, "notes": ["Mapping skipped"]},
            "bond_changes": [],
            "reaction_center_atoms": []
        }
        return result

    # Step 2: Map
    mapping = map_reaction(clean.rxn_smiles_clean)

    # Step 3: Extract bond changes if mapping succeeded
    if mapping.mapping_qc.get("ok", False) and mapping.mapped_rxn_smiles:
        bond_report = extract_bond_changes(mapping.mapped_rxn_smiles)
    else:
        bond_report = BondChangeReport()

    result["tool_facts"] = {
        "mapped_rxn_smiles": mapping.mapped_rxn_smiles,
        "mapping_qc": mapping.mapping_qc,
        "bond_changes": [bc.to_dict() for bc in bond_report.bond_changes],
        "reaction_center_atoms": bond_report.reaction_center_atoms,
        "bond_change_warnings": bond_report.warnings
    }

    return result

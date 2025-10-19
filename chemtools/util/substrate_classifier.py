"""
Substrate Classification for Chemical Molecules.

Classifies chemical substrates by their structural features, functional groups,
and reactive positions. Designed to be reusable across multiple chemtools modules.

Usage:
    from chemtools.util.substrate_classifier import SubstrateClassifier, SubstrateInfo
    
    # Classify a substrate
    classifier = SubstrateClassifier()
    info = classifier.classify("CCCCCCI")
    
    print(info.substrate_class)        # 'primary_alkyl_iodide'
    print(info.carbon_types)           # {0: 'sp3', 1: 'sp3', ...}
    print(info.functional_groups)      # ['alkyl_iodide']
    print(info.special_positions.benzylic)  # []
    
Who uses this:
    - SMARTS generator (protocol scope definition)
    - Featurizers (ML feature extraction)
    - Recommendation engine (substrate validation)
    - Reaction type detector (understanding reactants)
    - Dataset analytics (classify substrate distributions)
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Any
from .rdkit_helpers import rdkit_available, parse_smiles
from .functional_groups import detect_all, get_functional_groups


# ============================================================================
# Data Classes
# ============================================================================

@dataclass
class SpecialPositions:
    """Tracks special reactive positions in molecule"""
    benzylic: List[int] = field(default_factory=list)      # [CH2] next to aromatic
    allylic: List[int] = field(default_factory=list)       # [CH2] next to C=C
    propargylic: List[int] = field(default_factory=list)   # [CH2] next to C≡C
    alpha_to_carbonyl: List[int] = field(default_factory=list)  # [CH2] next to C=O
    ortho_to_heteroatom: List[int] = field(default_factory=list)  # ortho to N,O,S on aromatic


@dataclass
class SubstrateInfo:
    """Complete classification of a chemical substrate"""
    # Primary classification
    substrate_class: str  # 'primary_alkyl_iodide', 'aryl_bromide', 'aniline', etc.
    substrate_family: str  # 'halide', 'amine', 'alcohol', 'carbonyl', etc.
    
    # Molecular properties
    smiles: str
    mol: Any = None  # RDKit mol object (if available)
    
    # Carbon hybridization map
    carbon_types: Dict[int, str] = field(default_factory=dict)  # atom_idx → 'sp3'|'sp2'|'sp'|'aromatic'
    
    # Functional groups
    functional_groups: List[str] = field(default_factory=list)  # detected functional groups
    functional_group_atoms: Dict[str, List[int]] = field(default_factory=dict)  # FG → atom indices
    
    # Special positions
    special_positions: SpecialPositions = field(default_factory=SpecialPositions)
    
    # Reactive centers
    reactive_centers: List[int] = field(default_factory=list)  # primary reactive atom indices
    reactive_center_types: Dict[int, str] = field(default_factory=dict)  # atom_idx → type
    
    # Additional context
    has_aromatic: bool = False
    has_heteroaromatic: bool = False
    ring_count: int = 0
    aromatic_ring_count: int = 0


# ============================================================================
# Substrate Classifier
# ============================================================================

class SubstrateClassifier:
    """Classify chemical substrates into meaningful categories"""
    
    def __init__(self):
        """Initialize classifier"""
        self._rdkit_available = rdkit_available()
    
    def classify(self, mol_or_smiles: Any) -> SubstrateInfo:
        """Complete substrate classification
        
        Args:
            mol_or_smiles: SMILES string or RDKit mol object
            
        Returns:
            SubstrateInfo with complete classification
        """
        # Parse input
        if isinstance(mol_or_smiles, str):
            smiles = mol_or_smiles
            mol = parse_smiles(smiles) if self._rdkit_available else None
        else:
            mol = mol_or_smiles
            if mol and self._rdkit_available:
                from rdkit import Chem
                smiles = Chem.MolToSmiles(mol)
            else:
                smiles = str(mol_or_smiles) if mol_or_smiles else ""
        
        # Initialize result
        info = SubstrateInfo(
            substrate_class="unknown",
            substrate_family="unknown",
            smiles=smiles,
            mol=mol
        )
        
        if not mol and not smiles:
            return info
        
        # Detect functional groups
        fg_dict = detect_all(smiles)
        info.functional_groups = get_functional_groups(smiles)
        
        # RDKit-based analysis (if available)
        if mol and self._rdkit_available:
            self._analyze_with_rdkit(mol, info, fg_dict)
        else:
            # Fallback text-based analysis
            self._analyze_with_text(smiles, info, fg_dict)
        
        # Classify substrate based on analysis
        self._determine_substrate_class(info, fg_dict)
        
        return info
    
    def _analyze_with_rdkit(self, mol, info: SubstrateInfo, fg_dict: Dict[str, bool]):
        """RDKit-based detailed analysis"""
        from rdkit import Chem
        
        # Carbon hybridization analysis
        info.carbon_types = self.get_carbon_types(mol)
        
        # Find special positions
        info.special_positions = self.find_special_positions(mol)
        
        # Analyze aromatic systems
        info.has_aromatic = any(atom.GetIsAromatic() for atom in mol.GetAtoms())
        info.has_heteroaromatic = any(
            atom.GetIsAromatic() and atom.GetAtomicNum() not in {6, 1}
            for atom in mol.GetAtoms()
        )
        
        # Ring analysis
        ring_info = mol.GetRingInfo()
        info.ring_count = ring_info.NumRings()
        info.aromatic_ring_count = sum(
            1 for ring in ring_info.AtomRings()
            if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
        )
        
        # Find reactive centers
        info.reactive_centers = self._find_reactive_centers(mol, fg_dict)
        info.reactive_center_types = self._classify_reactive_centers(mol, info.reactive_centers)
        
        # Map functional groups to atoms
        info.functional_group_atoms = self._map_functional_groups_to_atoms(mol, fg_dict)
    
    def _analyze_with_text(self, smiles: str, info: SubstrateInfo, fg_dict: Dict[str, bool]):
        """Text-based fallback analysis"""
        s = smiles.lower()
        
        # Simple aromatic detection
        info.has_aromatic = 'c1' in s or any(c in s for c in ['c[', 'cc(', 'cn'])
        info.has_heteroaromatic = info.has_aromatic and any(x in s for x in ['n', 'o', 's', 'p'])
        
        # Rough ring count (count '1', '2', '3' ring closure numbers)
        info.ring_count = sum(s.count(str(i)) for i in range(1, 10)) // 2
    
    def get_carbon_types(self, mol) -> Dict[int, str]:
        """Map atom index to carbon hybridization type
        
        Args:
            mol: RDKit mol object
            
        Returns:
            Dict mapping atom_idx → 'sp3'|'sp2'|'sp'|'aromatic'
        """
        if not self._rdkit_available:
            return {}
        
        from rdkit import Chem
        carbon_types = {}
        
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 6:  # Only carbons
                continue
            
            idx = atom.GetIdx()
            
            if atom.GetIsAromatic():
                carbon_types[idx] = 'aromatic'
            else:
                hybrid = atom.GetHybridization()
                if hybrid == Chem.HybridizationType.SP3:
                    carbon_types[idx] = 'sp3'
                elif hybrid == Chem.HybridizationType.SP2:
                    carbon_types[idx] = 'sp2'
                elif hybrid == Chem.HybridizationType.SP:
                    carbon_types[idx] = 'sp'
                else:
                    carbon_types[idx] = 'unknown'
        
        return carbon_types
    
    def find_special_positions(self, mol) -> SpecialPositions:
        """Identify benzylic, allylic, propargylic positions
        
        Args:
            mol: RDKit mol object
            
        Returns:
            SpecialPositions object with atom indices
        """
        positions = SpecialPositions()
        
        if not self._rdkit_available:
            return positions
        
        from rdkit import Chem
        
        # SMARTS patterns for special positions
        patterns = {
            'benzylic': Chem.MolFromSmarts("[CX4;H1,H2,H3]-[c,$(C=C-c)]"),  # CH next to aromatic
            'allylic': Chem.MolFromSmarts("[CX4;H1,H2,H3]-[CX3]=[CX3]"),    # CH next to C=C
            'propargylic': Chem.MolFromSmarts("[CX4;H1,H2,H3]-[CX2]#[CX2]"),  # CH next to C≡C
            'alpha_to_carbonyl': Chem.MolFromSmarts("[CX4;H1,H2,H3]-[CX3]=[OX1]"),  # CH next to C=O
        }
        
        for pos_type, pattern in patterns.items():
            if pattern:
                matches = mol.GetSubstructMatches(pattern)
                # First atom in match is the special position
                atom_indices = [match[0] for match in matches]
                setattr(positions, pos_type, sorted(set(atom_indices)))
        
        # Ortho to heteroatom on aromatic ring
        ortho_pattern = Chem.MolFromSmarts("c1ccccc1")  # Simple benzene
        if ortho_pattern:
            ring_matches = mol.GetSubstructMatches(ortho_pattern)
            ortho_atoms = []
            for ring in ring_matches:
                for i, atom_idx in enumerate(ring):
                    # Check neighbors for heteroatoms
                    atom = mol.GetAtomWithIdx(atom_idx)
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() not in {1, 6} and neighbor.GetIsAromatic():
                            ortho_atoms.append(atom_idx)
                            break
            positions.ortho_to_heteroatom = sorted(set(ortho_atoms))
        
        return positions
    
    def _find_reactive_centers(self, mol, fg_dict: Dict[str, bool]) -> List[int]:
        """Find primary reactive atom indices
        
        Args:
            mol: RDKit mol object
            fg_dict: Detected functional groups
            
        Returns:
            List of reactive atom indices
        """
        if not self._rdkit_available:
            return []
        
        from rdkit import Chem
        reactive = []
        
        # Priority order: heteroatoms > special carbons
        heteroatom_patterns = [
            ("[I,Br,Cl,F]", "halogen"),
            ("[NX3;H1,H2;!$(NC=O)]", "amine"),
            ("[OX2H]", "alcohol"),
            ("[SX2H]", "thiol"),
            ("[BX3]", "boron"),
            ("[CX3]=[OX1]", "carbonyl"),
        ]
        
        for smarts, _ in heteroatom_patterns:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern:
                matches = mol.GetSubstructMatches(pattern)
                for match in matches:
                    reactive.append(match[0])
        
        return sorted(set(reactive))
    
    def _classify_reactive_centers(self, mol, atom_indices: List[int]) -> Dict[int, str]:
        """Classify each reactive center by type
        
        Args:
            mol: RDKit mol object
            atom_indices: Reactive atom indices
            
        Returns:
            Dict mapping atom_idx → type string
        """
        if not self._rdkit_available:
            return {}
        
        types = {}
        for idx in atom_indices:
            atom = mol.GetAtomWithIdx(idx)
            symbol = atom.GetSymbol()
            
            if symbol in ['I', 'Br', 'Cl', 'F']:
                types[idx] = f"{symbol.lower()}_halide"
            elif symbol == 'N':
                types[idx] = "amine_or_amide"
            elif symbol == 'O':
                types[idx] = "alcohol_or_ether"
            elif symbol == 'S':
                types[idx] = "thiol_or_sulfide"
            elif symbol == 'B':
                types[idx] = "boron"
            elif symbol == 'C':
                # Check if carbonyl
                if any(n.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(idx, n.GetIdx()).GetBondTypeAsDouble() == 2.0
                       for n in atom.GetNeighbors()):
                    types[idx] = "carbonyl"
                else:
                    types[idx] = "carbon"
            else:
                types[idx] = "other"
        
        return types
    
    def _map_functional_groups_to_atoms(self, mol, fg_dict: Dict[str, bool]) -> Dict[str, List[int]]:
        """Map functional groups to their atom indices
        
        Args:
            mol: RDKit mol object
            fg_dict: Detected functional groups
            
        Returns:
            Dict mapping FG name → atom indices
        """
        if not self._rdkit_available:
            return {}
        
        from rdkit import Chem
        from .functional_groups import FUNCTIONAL_GROUP_SMARTS
        
        fg_atoms = {}
        
        for fg_name, is_present in fg_dict.items():
            if not is_present:
                continue
            
            smarts = FUNCTIONAL_GROUP_SMARTS.get(fg_name)
            if not smarts:
                continue
            
            pattern = Chem.MolFromSmarts(smarts)
            if not pattern:
                continue
            
            matches = mol.GetSubstructMatches(pattern)
            if matches:
                # Flatten all matches
                atoms = sorted(set(atom_idx for match in matches for atom_idx in match))
                fg_atoms[fg_name] = atoms
        
        return fg_atoms
    
    def _determine_substrate_class(self, info: SubstrateInfo, fg_dict: Dict[str, bool]):
        """Determine substrate class and family based on analysis
        
        Args:
            info: SubstrateInfo to update
            fg_dict: Detected functional groups
        """
        # Priority classification order
        
        # 1. Halides (most common in cross-coupling)
        if any(fg_dict.get(f"alkyl_{x}") for x in ['iodide', 'bromide', 'chloride', 'fluoride']):
            info.substrate_family = 'halide'
            info.substrate_class = self._classify_halide(info, fg_dict)
        
        elif any(fg_dict.get(f"aryl_{x}") for x in ['iodide', 'bromide', 'chloride', 'fluoride']):
            info.substrate_family = 'halide'
            info.substrate_class = self._classify_aryl_halide(info, fg_dict)
        
        # 2. Amines
        elif fg_dict.get('aniline'):
            info.substrate_family = 'amine'
            info.substrate_class = 'aniline'
        
        elif any(fg_dict.get(f) for f in ['amine_primary', 'amine_secondary', 'amine_tertiary']):
            info.substrate_family = 'amine'
            if fg_dict.get('amine_primary'):
                info.substrate_class = 'primary_amine'
            elif fg_dict.get('amine_secondary'):
                info.substrate_class = 'secondary_amine'
            else:
                info.substrate_class = 'tertiary_amine'
        
        # 3. Amides (separate from amines!)
        elif fg_dict.get('amide'):
            info.substrate_family = 'amide'
            if fg_dict.get('amide_primary'):
                info.substrate_class = 'primary_amide'
            elif fg_dict.get('amide_secondary'):
                info.substrate_class = 'secondary_amide'
            else:
                info.substrate_class = 'tertiary_amide'
        
        # 4. Alcohols and phenols
        elif fg_dict.get('phenol'):
            info.substrate_family = 'alcohol'
            info.substrate_class = 'phenol'
        
        elif fg_dict.get('alcohol'):
            info.substrate_family = 'alcohol'
            if info.special_positions.benzylic:
                info.substrate_class = 'benzylic_alcohol'
            elif info.special_positions.allylic:
                info.substrate_class = 'allylic_alcohol'
            else:
                info.substrate_class = 'aliphatic_alcohol'
        
        # 5. Carbonyls
        elif fg_dict.get('carboxylic_acid'):
            info.substrate_family = 'carbonyl'
            info.substrate_class = 'carboxylic_acid'
        
        elif fg_dict.get('ester'):
            info.substrate_family = 'carbonyl'
            info.substrate_class = 'ester'
        
        elif fg_dict.get('aldehyde'):
            info.substrate_family = 'carbonyl'
            info.substrate_class = 'aldehyde'
        
        elif fg_dict.get('ketone'):
            info.substrate_family = 'carbonyl'
            info.substrate_class = 'ketone'
        
        # 6. Boron compounds
        elif 'B(' in info.smiles or 'B1' in info.smiles:
            info.substrate_family = 'boron'
            if 'B(O)O' in info.smiles or 'B(OH)' in info.smiles:
                info.substrate_class = 'boronic_acid'
            elif 'B1OC(C)(C)C(C)(C)O1' in info.smiles or 'Bpin' in info.smiles.lower():
                info.substrate_class = 'boronic_ester_pinacol'
            else:
                info.substrate_class = 'boron_compound'
        
        # 7. Other heteroatoms
        elif fg_dict.get('thiol'):
            info.substrate_family = 'sulfur'
            info.substrate_class = 'thiol'
        
        elif fg_dict.get('sulfonic_acid'):
            info.substrate_family = 'sulfur'
            info.substrate_class = 'sulfonic_acid'
        
        elif fg_dict.get('triflate'):
            info.substrate_family = 'leaving_group'
            info.substrate_class = 'triflate'
        
        # 8. Simple hydrocarbons
        elif info.has_aromatic:
            info.substrate_family = 'hydrocarbon'
            if info.has_heteroaromatic:
                info.substrate_class = 'heteroaromatic'
            else:
                info.substrate_class = 'aromatic'
        
        elif fg_dict.get('alkene'):
            info.substrate_family = 'hydrocarbon'
            info.substrate_class = 'alkene'
        
        elif fg_dict.get('alkyne'):
            info.substrate_family = 'hydrocarbon'
            info.substrate_class = 'alkyne'
        
        else:
            info.substrate_family = 'hydrocarbon'
            info.substrate_class = 'alkane'
    
    def _classify_halide(self, info: SubstrateInfo, fg_dict: Dict[str, bool]) -> str:
        """Classify alkyl halide by position and substitution"""
        # Determine halogen
        halogen = None
        for x in ['iodide', 'bromide', 'chloride', 'fluoride']:
            if fg_dict.get(f"alkyl_{x}"):
                halogen = x
                break
        
        if not halogen:
            return 'alkyl_halide'
        
        # Check for special positions
        if info.special_positions.benzylic:
            return f'benzylic_{halogen}'
        
        if info.special_positions.allylic:
            return f'allylic_{halogen}'
        
        if info.special_positions.propargylic:
            return f'propargylic_{halogen}'
        
        # Determine substitution (primary, secondary, tertiary)
        # This requires analyzing the carbon bearing the halogen
        if info.mol and self._rdkit_available:
            from rdkit import Chem
            # Find halogen atom
            halogen_symbol = {'iodide': 'I', 'bromide': 'Br', 'chloride': 'Cl', 'fluoride': 'F'}[halogen]
            pattern = Chem.MolFromSmarts(f"[CX4]-[{halogen_symbol}]")
            if pattern:
                matches = info.mol.GetSubstructMatches(pattern)
                if matches:
                    carbon_idx = matches[0][0]
                    carbon = info.mol.GetAtomWithIdx(carbon_idx)
                    h_count = carbon.GetTotalNumHs()
                    
                    if h_count >= 2:
                        return f'primary_alkyl_{halogen}'
                    elif h_count == 1:
                        return f'secondary_alkyl_{halogen}'
                    else:
                        return f'tertiary_alkyl_{halogen}'
        
        return f'alkyl_{halogen}'
    
    def _classify_aryl_halide(self, info: SubstrateInfo, fg_dict: Dict[str, bool]) -> str:
        """Classify aryl halide"""
        # Determine halogen
        halogen = None
        for x in ['iodide', 'bromide', 'chloride', 'fluoride']:
            if fg_dict.get(f"aryl_{x}"):
                halogen = x
                break
        
        if not halogen:
            return 'aryl_halide'
        
        # Check if heteroaryl
        if info.has_heteroaromatic:
            return f'heteroaryl_{halogen}'
        
        return f'aryl_{halogen}'


# ============================================================================
# Convenience Functions
# ============================================================================

def classify_substrate(mol_or_smiles: Any) -> SubstrateInfo:
    """Convenience function to classify a substrate
    
    Args:
        mol_or_smiles: SMILES string or RDKit mol object
        
    Returns:
        SubstrateInfo with complete classification
    """
    classifier = SubstrateClassifier()
    return classifier.classify(mol_or_smiles)


def get_substrate_class(mol_or_smiles: Any) -> str:
    """Get just the substrate class string
    
    Args:
        mol_or_smiles: SMILES string or RDKit mol object
        
    Returns:
        Substrate class string (e.g., 'primary_alkyl_iodide')
    """
    info = classify_substrate(mol_or_smiles)
    return info.substrate_class


def get_substrate_family(mol_or_smiles: Any) -> str:
    """Get just the substrate family string
    
    Args:
        mol_or_smiles: SMILES string or RDKit mol object
        
    Returns:
        Substrate family string (e.g., 'halide', 'amine')
    """
    info = classify_substrate(mol_or_smiles)
    return info.substrate_family

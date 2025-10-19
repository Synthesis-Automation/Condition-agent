"""
Context-Aware SMARTS Pattern Builder.

Builds chemically meaningful SMARTS patterns based on substrate classification.
Uses SubstrateClassifier to understand chemical context and generate appropriate patterns.

Usage:
    from chemtools.util.smarts_builders import SmartsBuilder
    from chemtools.util.substrate_classifier import classify_substrate
    
    # Build SMARTS for a substrate
    builder = SmartsBuilder()
    info = classify_substrate("CCCCCCI")
    smarts = builder.build_for_substrate(info)
    # Result: "[CX4;H2,H3]-[I]"  (primary alkyl iodide)
    
    # Or use convenience function
    smarts = builder.build_from_smiles("c1ccccc1Br")
    # Result: "c-[Br]"  (aryl bromide)

Who uses this:
    - SMARTS generator CLI (protocol scope definition)
    - Protocol database queries (find compatible protocols)
    - Substrate validation (check if substrate matches scope)
    - Test frameworks (validate patterns)
"""

from __future__ import annotations
from typing import Optional, List, Dict, Any
from .rdkit_helpers import rdkit_available, parse_smiles
from .substrate_classifier import SubstrateClassifier, SubstrateInfo


class SmartsBuilder:
    """Build chemically meaningful SMARTS patterns from substrate classification"""
    
    def __init__(self):
        """Initialize SMARTS builder"""
        self.classifier = SubstrateClassifier()
        self._rdkit_available = rdkit_available()
    
    def build_from_smiles(self, smiles: str) -> str:
        """Build SMARTS pattern from SMILES string
        
        Args:
            smiles: SMILES string
            
        Returns:
            Context-aware SMARTS pattern
        """
        info = self.classifier.classify(smiles)
        return self.build_for_substrate(info)
    
    def build_for_substrate(self, info: SubstrateInfo) -> str:
        """Build SMARTS pattern for substrate based on classification
        
        Args:
            info: SubstrateInfo from SubstrateClassifier
            
        Returns:
            Context-aware SMARTS pattern
        """
        # Route to appropriate builder based on substrate family
        if info.substrate_family == 'halide':
            return self.build_halide_smarts(info)
        
        elif info.substrate_family == 'amine':
            return self.build_amine_smarts(info)
        
        elif info.substrate_family == 'amide':
            return self.build_amide_smarts(info)
        
        elif info.substrate_family == 'alcohol':
            return self.build_alcohol_smarts(info)
        
        elif info.substrate_family == 'carbonyl':
            return self.build_carbonyl_smarts(info)
        
        elif info.substrate_family == 'boron':
            return self.build_boron_smarts(info)
        
        elif info.substrate_family == 'sulfur':
            return self.build_sulfur_smarts(info)
        
        elif info.substrate_family == 'leaving_group':
            return self.build_leaving_group_smarts(info)
        
        else:
            # Fallback: generic pattern
            return self._build_generic_pattern(info)
    
    def build_halide_smarts(self, info: SubstrateInfo) -> str:
        """Build SMARTS for halide substrates
        
        Handles:
        - Primary/secondary/tertiary alkyl halides
        - Aryl halides
        - Heteroaryl halides
        - Benzylic halides
        - Allylic halides
        - Propargylic halides
        """
        # Extract halogen symbol
        halogen = self._extract_halogen(info)
        if not halogen:
            return "[X]"  # Generic halogen
        
        # Aryl halides (aromatic carbon bonded to halogen)
        if 'aryl' in info.substrate_class:
            if info.has_heteroaromatic:
                # Heteroaryl: any aromatic carbon with halogen
                return f"[a]-[{halogen}]"
            else:
                # Simple aryl: aromatic carbon (lowercase c)
                return f"c-[{halogen}]"
        
        # Special positions (benzylic, allylic, propargylic)
        if 'benzylic' in info.substrate_class:
            # Benzylic: CH2 directly bonded to aromatic
            return f"[CH2;$([CH2][c])]-[{halogen}]"
        
        if 'allylic' in info.substrate_class:
            # Allylic: CH2 next to C=C
            return f"[CH2;$([CH2]C=C)]-[{halogen}]"
        
        if 'propargylic' in info.substrate_class:
            # Propargylic: CH2 next to C≡C
            return f"[CH2;$([CH2]C#C)]-[{halogen}]"
        
        # Alkyl halides by substitution
        if 'primary' in info.substrate_class:
            # Primary: sp3 carbon with 2 or 3 hydrogens, bonded to halogen
            # This ensures we only match the carbon directly bonded to halogen
            return f"[CX4;H2,H3]-[{halogen}]"
        
        elif 'secondary' in info.substrate_class:
            # Secondary: sp3 carbon with 1 hydrogen, bonded to halogen
            return f"[CX4;H1]-[{halogen}]"
        
        elif 'tertiary' in info.substrate_class:
            # Tertiary: sp3 carbon with 0 hydrogens, bonded to halogen
            return f"[CX4;H0]-[{halogen}]"
        
        else:
            # Generic alkyl halide (sp3 carbon)
            return f"[CX4]-[{halogen}]"
    
    def build_amine_smarts(self, info: SubstrateInfo) -> str:
        """Build SMARTS for amine substrates
        
        Handles:
        - Aniline (aromatic amine)
        - Primary/secondary/tertiary aliphatic amines
        """
        # Aniline (aromatic amine)
        if info.substrate_class == 'aniline':
            # Aromatic carbon bonded to NH2
            return "c-[NX3;H2;!$(NC=O)]"
        
        # Aliphatic amines by degree
        if 'primary' in info.substrate_class:
            # Primary: NH2 bonded to carbon, not part of amide
            return "[CX4]-[NX3;H2;!$(NC=O)]"
        
        elif 'secondary' in info.substrate_class:
            # Secondary: NH bonded to carbon
            return "[NX3;H1;!$(NC=O)]"
        
        elif 'tertiary' in info.substrate_class:
            # Tertiary: N bonded to carbons, no H
            return "[NX3;H0;!$(NC=O);!$(N=*)]"
        
        else:
            # Generic amine
            return "[NX3;H2,H1;!$(NC=O)]"
    
    def build_amide_smarts(self, info: SubstrateInfo) -> str:
        """Build SMARTS for amide substrates
        
        Critical: Must distinguish from simple amines!
        """
        if 'primary' in info.substrate_class:
            # Primary amide: NH2-C(=O)
            return "[NX3;H2]-[CX3](=O)"
        
        elif 'secondary' in info.substrate_class:
            # Secondary amide: NH-C(=O)
            return "[NX3;H1]-[CX3](=O)"
        
        elif 'tertiary' in info.substrate_class:
            # Tertiary amide: N-C(=O)
            return "[NX3;H0]-[CX3](=O)"
        
        else:
            # Generic amide
            return "[NX3]-[CX3](=O)"
    
    def build_alcohol_smarts(self, info: SubstrateInfo) -> str:
        """Build SMARTS for alcohol substrates
        
        Handles:
        - Phenol (aromatic OH)
        - Benzylic alcohol
        - Allylic alcohol
        - Aliphatic alcohol
        """
        # Phenol (aromatic carbon with OH)
        if info.substrate_class == 'phenol':
            return "c-[OX2H]"
        
        # Benzylic alcohol (CH2OH next to aromatic)
        if 'benzylic' in info.substrate_class:
            return "[CH2;$([CH2][c])]-[OX2H]"
        
        # Allylic alcohol (CH2OH next to C=C)
        if 'allylic' in info.substrate_class:
            return "[CH2;$([CH2]C=C)]-[OX2H]"
        
        # Aliphatic alcohol (sp3 carbon with OH)
        return "[CX4]-[OX2H]"
    
    def build_carbonyl_smarts(self, info: SubstrateInfo) -> str:
        """Build SMARTS for carbonyl substrates
        
        Handles:
        - Carboxylic acid
        - Ester
        - Aldehyde
        - Ketone
        """
        if info.substrate_class == 'carboxylic_acid':
            # Carboxylic acid: C(=O)OH
            return "[CX3](=O)-[OX2H]"
        
        elif info.substrate_class == 'ester':
            # Ester: C(=O)O-C
            return "[CX3](=O)-[OX2]-[C]"
        
        elif info.substrate_class == 'aldehyde':
            # Aldehyde: C(=O)H
            return "[CX3;H1](=O)"
        
        elif info.substrate_class == 'ketone':
            # Ketone: C(=O)C
            return "[C]-[CX3](=O)-[C]"
        
        else:
            # Generic carbonyl
            return "[CX3]=O"
    
    def build_boron_smarts(self, info: SubstrateInfo) -> str:
        """Build SMARTS for boron substrates
        
        Handles:
        - Boronic acid B(OH)2
        - Boronic ester (pinacol)
        - BF3K
        """
        if info.substrate_class == 'boronic_acid':
            # Boronic acid: B(OH)2
            return "[B]([OH])([OH])"
        
        elif 'pinacol' in info.substrate_class:
            # Boronic ester pinacol: B-O-C ring system
            return "[B]1OC(C)(C)C(C)(C)O1"
        
        elif 'bf3k' in info.substrate_class.lower():
            # Potassium trifluoroborate
            return "[B-]([F])([F])[F]"
        
        else:
            # Generic boron
            return "[BX3]"
    
    def build_sulfur_smarts(self, info: SubstrateInfo) -> str:
        """Build SMARTS for sulfur-containing substrates"""
        if info.substrate_class == 'thiol':
            return "[SX2H]"
        
        elif info.substrate_class == 'sulfonic_acid':
            return "[SX4](=O)(=O)-[OX2H]"
        
        else:
            return "[S]"
    
    def build_leaving_group_smarts(self, info: SubstrateInfo) -> str:
        """Build SMARTS for leaving groups (triflate, tosylate, etc.)"""
        if info.substrate_class == 'triflate':
            # Triflate: O-S(=O)(=O)-CF3
            return "[OX2]-[SX4](=O)(=O)-[CX4](F)(F)F"
        
        else:
            return "[O]-[S](=O)(=O)"
    
    def _build_generic_pattern(self, info: SubstrateInfo) -> str:
        """Fallback generic pattern builder"""
        if info.has_aromatic:
            return "a"  # Any aromatic atom
        else:
            return "[#6]"  # Any carbon
    
    def _extract_halogen(self, info: SubstrateInfo) -> Optional[str]:
        """Extract halogen symbol from substrate class and functional groups"""
        # First check functional groups (most reliable)
        halogen_map = {
            'alkyl_iodide': 'I',
            'aryl_iodide': 'I',
            'alkyl_bromide': 'Br',
            'aryl_bromide': 'Br',
            'alkyl_chloride': 'Cl',
            'aryl_chloride': 'Cl',
            'alkyl_fluoride': 'F',
            'aryl_fluoride': 'F',
        }
        
        for fg in info.functional_groups:
            if fg in halogen_map:
                return halogen_map[fg]
        
        # Then check substrate class (check longer strings first to avoid 'i' matching 'bromide')
        for halogen in ['Br', 'Cl', 'I', 'F']:  # Br first since 'i' in 'bromide' would match I
            if halogen.lower() + 'ide' in info.substrate_class.lower():
                return halogen
        
        # Last resort: check for halogen symbol directly
        for halogen in ['Br', 'Cl', 'I', 'F']:
            if halogen.lower() == info.substrate_class.lower()[-1]:
                return halogen
        
        return None
    
    def build_guard_patterns(self, info: SubstrateInfo) -> List[str]:
        """Generate guard patterns (forbidden) based on substrate class
        
        Args:
            info: SubstrateInfo from classification
            
        Returns:
            List of SMARTS patterns that should be forbidden
        """
        guards = []
        
        # For primary alkyl halides, exclude secondary/tertiary/benzylic/allylic
        if 'primary_alkyl' in info.substrate_class:
            halogen = self._extract_halogen(info)
            if halogen:
                guards.extend([
                    f"[CX4;H1]-[{halogen}]  # Exclude secondary",
                    f"[CX4;H0]-[{halogen}]  # Exclude tertiary",
                    f"[CH2;$([CH2][c])]-[{halogen}]  # Exclude benzylic",
                    f"[CH2;$([CH2]C=C)]-[{halogen}]  # Exclude allylic",
                ])
        
        # For secondary alkyl halides, exclude primary/tertiary
        elif 'secondary_alkyl' in info.substrate_class:
            halogen = self._extract_halogen(info)
            if halogen:
                guards.extend([
                    f"[CX4;H2,H3]-[{halogen}]  # Exclude primary",
                    f"[CX4;H0]-[{halogen}]  # Exclude tertiary",
                ])
        
        # For aryl halides, exclude alkyl
        elif 'aryl' in info.substrate_class:
            halogen = self._extract_halogen(info)
            if halogen:
                guards.extend([
                    f"[CX4]-[{halogen}]  # Exclude aliphatic halides",
                ])
        
        # For aniline, exclude aliphatic amines
        elif info.substrate_class == 'aniline':
            guards.extend([
                "[CX4]-[NX3;H2;!$(NC=O)]  # Exclude aliphatic primary amines",
                "[CX4]-[NX3;H1;!$(NC=O)]  # Exclude aliphatic secondary amines",
            ])
        
        # For primary aliphatic amine, exclude aniline
        elif 'primary_amine' in info.substrate_class:
            guards.extend([
                "c-[NX3;H2]  # Exclude aniline",
            ])
        
        # For amide, exclude simple amines
        elif info.substrate_family == 'amide':
            guards.extend([
                "[NX3;H2;!$(NC=O)]  # Exclude non-amide primary amines",
                "[NX3;H1;!$(NC=O)]  # Exclude non-amide secondary amines",
            ])
        
        # For phenol, exclude aliphatic alcohols
        elif info.substrate_class == 'phenol':
            guards.extend([
                "[CX4]-[OX2H]  # Exclude aliphatic alcohols",
            ])
        
        return guards


class SmartsPatternMatcher:
    """Match molecules against SMARTS patterns with explanations"""
    
    def __init__(self):
        """Initialize pattern matcher"""
        self._rdkit_available = rdkit_available()
    
    def match(self, mol_or_smiles: Any, smarts: str) -> bool:
        """Check if molecule matches SMARTS pattern
        
        Args:
            mol_or_smiles: SMILES string or RDKit mol object
            smarts: SMARTS pattern
            
        Returns:
            True if matches, False otherwise
        """
        if not self._rdkit_available:
            # Fallback: simple text matching
            smiles = str(mol_or_smiles)
            return self._text_based_match(smiles, smarts)
        
        from rdkit import Chem
        
        # Parse molecule
        if isinstance(mol_or_smiles, str):
            mol = parse_smiles(mol_or_smiles)
        else:
            mol = mol_or_smiles
        
        if mol is None:
            return False
        
        # Parse SMARTS
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            return False
        
        # Match
        return mol.HasSubstructMatch(pattern)
    
    def explain_match(self, mol_or_smiles: Any, smarts: str) -> str:
        """Explain why molecule matches or doesn't match pattern
        
        Args:
            mol_or_smiles: SMILES string or RDKit mol object
            smarts: SMARTS pattern
            
        Returns:
            Human-readable explanation
        """
        matches = self.match(mol_or_smiles, smarts)
        
        smiles = str(mol_or_smiles)
        
        if matches:
            return f"✅ Molecule '{smiles}' matches pattern '{smarts}'"
        else:
            return f"❌ Molecule '{smiles}' does not match pattern '{smarts}'"
    
    def find_matching_atoms(self, mol_or_smiles: Any, smarts: str) -> List[int]:
        """Find atom indices that match SMARTS pattern
        
        Args:
            mol_or_smiles: SMILES string or RDKit mol object
            smarts: SMARTS pattern
            
        Returns:
            List of matching atom indices
        """
        if not self._rdkit_available:
            return []
        
        from rdkit import Chem
        
        # Parse molecule
        if isinstance(mol_or_smiles, str):
            mol = parse_smiles(mol_or_smiles)
        else:
            mol = mol_or_smiles
        
        if mol is None:
            return []
        
        # Parse SMARTS
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            return []
        
        # Find matches
        matches = mol.GetSubstructMatches(pattern)
        if not matches:
            return []
        
        # Flatten all matches
        return sorted(set(atom_idx for match in matches for atom_idx in match))
    
    def _text_based_match(self, smiles: str, smarts: str) -> bool:
        """Fallback text-based matching (very approximate)"""
        # Very simple fallback - just check if key atoms present
        s_lower = smiles.lower()
        smarts_lower = smarts.lower()
        
        # Check for common patterns
        if '[i]' in smarts_lower and 'i' not in s_lower:
            return False
        if '[br]' in smarts_lower and 'br' not in s_lower:
            return False
        if '[cl]' in smarts_lower and 'cl' not in s_lower:
            return False
        
        return True  # Default to match if can't determine


# ============================================================================
# Convenience Functions
# ============================================================================

def build_smarts(smiles: str) -> str:
    """Convenience function to build SMARTS from SMILES
    
    Args:
        smiles: SMILES string
        
    Returns:
        Context-aware SMARTS pattern
    """
    builder = SmartsBuilder()
    return builder.build_from_smiles(smiles)


def build_smarts_with_guards(smiles: str) -> Dict[str, Any]:
    """Build SMARTS with guard patterns
    
    Args:
        smiles: SMILES string
        
    Returns:
        Dict with 'core' SMARTS and 'guards_forbid' list
    """
    builder = SmartsBuilder()
    classifier = SubstrateClassifier()
    
    info = classifier.classify(smiles)
    core = builder.build_for_substrate(info)
    guards = builder.build_guard_patterns(info)
    
    return {
        'core': core,
        'guards_forbid': guards,
        'substrate_class': info.substrate_class,
        'substrate_family': info.substrate_family,
    }


def match_smarts(smiles: str, smarts: str) -> bool:
    """Check if SMILES matches SMARTS pattern
    
    Args:
        smiles: SMILES string
        smarts: SMARTS pattern
        
    Returns:
        True if matches, False otherwise
    """
    matcher = SmartsPatternMatcher()
    return matcher.match(smiles, smarts)

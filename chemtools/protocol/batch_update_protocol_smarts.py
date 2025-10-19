"""
Batch processor for protocol JSON files - auto-generate reaction_smarts_applicability.

This tool scans all protocol JSON files in data/protocol_db/, extracts reaction_smiles,
generates chemistry-aware SMARTS patterns, and updates the files with the new patterns.

Usage:
    python -m chemtools.protocol.batch_update_protocol_smarts
    python -m chemtools.protocol.batch_update_protocol_smarts --dry-run
    python -m chemtools.protocol.batch_update_protocol_smarts --protocol-dir path/to/protocols
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, Any, List, Optional
from dataclasses import dataclass

# Add parent directory to path if running as script (not as module)
if __name__ == "__main__" and __package__ is None:
    _script_dir = Path(__file__).resolve().parent
    _project_root = _script_dir.parent.parent
    if str(_project_root) not in sys.path:
        sys.path.insert(0, str(_project_root))

from chemtools.util.smarts_builders import build_smarts_with_guards
from chemtools.util.substrate_classifier import classify_substrate
from chemtools.util.reaction_center_detector import identify_changed_atoms_from_mapped_smiles

# RXNUtils for professional reaction template generation (industry standard)
# NOTE: reaction-utils requires C++ build tools and CGRtools, which may not be available
# Our fallback manual pattern generation works excellently and produces clean, general patterns
try:
    from rxn_utils.chem.reaction import ChemicalReaction
    RXNUTILS_AVAILABLE = True
except ImportError:
    try:
        # Try alternative import path
        from rxnutils.chem.reaction import ChemicalReaction
        RXNUTILS_AVAILABLE = True
    except ImportError:
        RXNUTILS_AVAILABLE = False
        # Don't print warning - our fallback is excellent
        # print("[INFO] rxn-utils not available - using proven manual pattern generation")

# RDKit for SMARTS validation (optional)
try:
    from rdkit import Chem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


@dataclass
class ProcessingResult:
    """Result of processing a single protocol file"""
    filename: str
    success: bool
    message: str
    old_pattern: Optional[Dict[str, Any]] = None
    new_pattern: Optional[Dict[str, Any]] = None


class ProtocolSmartsUpdater:
    """Batch updater for protocol JSON files"""
    
    def __init__(self, protocol_dir: Path, dry_run: bool = False, radius: int = 0):
        """Initialize updater.
        
        Args:
            protocol_dir: Directory containing protocol JSON files
            dry_run: If True, don't write changes, just report what would happen
            radius: Radius for reaction template generation (0=minimal/general, 1-3=more specific)
        """
        self.protocol_dir = protocol_dir
        self.dry_run = dry_run
        self.radius = radius
        self.results: List[ProcessingResult] = []
    
    def process_all_protocols(self) -> List[ProcessingResult]:
        """Process all protocol JSON files in the directory.
        
        Returns:
            List of ProcessingResult objects
        """
        if not self.protocol_dir.exists():
            print(f"Error: Directory not found: {self.protocol_dir}")
            return []
        
        json_files = list(self.protocol_dir.glob("*.json"))
        
        if not json_files:
            print(f"Warning: No JSON files found in {self.protocol_dir}")
            return []
        
        print(f"\nFound {len(json_files)} protocol files in {self.protocol_dir}")
        print("=" * 80)
        
        for json_file in sorted(json_files):
            result = self.process_protocol_file(json_file)
            self.results.append(result)
        
        return self.results
    
    def process_protocol_file(self, filepath: Path) -> ProcessingResult:
        """Process a single protocol JSON file.
        
        Args:
            filepath: Path to protocol JSON file
            
        Returns:
            ProcessingResult object
        """
        filename = filepath.name
        
        try:
            # Load protocol JSON
            with open(filepath, 'r', encoding='utf-8') as f:
                protocol = json.load(f)
            
            # Extract reaction SMILES
            reaction_smiles = protocol.get("reaction", {}).get("reaction_smiles")
            
            if not reaction_smiles:
                return ProcessingResult(
                    filename=filename,
                    success=False,
                    message="No reaction_smiles found in protocol"
                )
            
            print(f"\n[Processing] {filename}")
            print(f"   Reaction: {reaction_smiles}")
            
            # Analyze reactants to show substrate information
            if ">>" in reaction_smiles:
                parts = reaction_smiles.split(">>")
            elif ">" in reaction_smiles:
                parts = reaction_smiles.split(">")
            else:
                parts = [reaction_smiles, ""]
            
            reactants = parts[0].strip()
            reactant_smiles_list = [s.strip() for s in reactants.split('.') if s.strip()]
            
            # Show substrate classification for each reactant
            for i, smiles in enumerate(reactant_smiles_list, 1):
                try:
                    substrate_info = classify_substrate(smiles)
                    substrate_class = substrate_info.substrate_class
                    family = substrate_info.substrate_family
                    special_positions = substrate_info.special_positions
                    
                    # Get special positions if any
                    special_list = []
                    if hasattr(special_positions, 'benzylic') and special_positions.benzylic:
                        special_list.append('benzylic')
                    if hasattr(special_positions, 'allylic') and special_positions.allylic:
                        special_list.append('allylic')
                    if hasattr(special_positions, 'propargylic') and special_positions.propargylic:
                        special_list.append('propargylic')
                    
                    special_str = f" [{', '.join(special_list)}]" if special_list else ""
                    print(f"   Reactant {i}: {substrate_class} ({family}){special_str}")
                except Exception as e:
                    print(f"   Reactant {i}: {smiles[:40]}... - {str(e)[:50]}")
            
            # Check if pattern already exists
            old_pattern = protocol.get("reaction", {}).get("reaction_smarts_applicability")
            if old_pattern:
                print(f"   Info: Existing pattern found (will be overwritten)")
            
            # Generate new SMARTS pattern
            try:
                new_pattern = self.generate_smarts_pattern(reaction_smiles, protocol)
                
                # Validate the generated pattern
                validation_result = self.validate_smarts_pattern(new_pattern)
                if not validation_result['valid']:
                    return ProcessingResult(
                        filename=filename,
                        success=False,
                        message=f"Generated pattern is invalid: {validation_result['error']}",
                        old_pattern=old_pattern,
                        new_pattern=new_pattern
                    )
                
                # Print pattern details
                print(f"   [OK] Generated pattern:")
                print(f"      Core: {new_pattern['core']}")
                if new_pattern.get('guards_forbid'):
                    print(f"      Guards: {len(new_pattern['guards_forbid'])} exclusion patterns")
                    for guard in new_pattern['guards_forbid'][:3]:  # Show first 3
                        print(f"        - {guard}")
                    if len(new_pattern['guards_forbid']) > 3:
                        print(f"        ... and {len(new_pattern['guards_forbid']) - 3} more")
                        
            except Exception as e:
                return ProcessingResult(
                    filename=filename,
                    success=False,
                    message=f"Pattern generation failed: {e}",
                    old_pattern=old_pattern
                )
            
            # Update protocol with new pattern
            if "reaction" not in protocol:
                protocol["reaction"] = {}
            
            protocol["reaction"]["reaction_smarts_applicability"] = new_pattern
            
            # Write updated protocol (unless dry-run)
            if not self.dry_run:
                with open(filepath, 'w', encoding='utf-8') as f:
                    json.dump(protocol, f, indent=2, ensure_ascii=False)
                print(f"   [SUCCESS] Updated successfully")
            else:
                print(f"   [DRY RUN] Would update file")
            
            return ProcessingResult(
                filename=filename,
                success=True,
                message="Pattern generated and applied",
                old_pattern=old_pattern,
                new_pattern=new_pattern
            )
            
        except json.JSONDecodeError as e:
            return ProcessingResult(
                filename=filename,
                success=False,
                message=f"JSON parsing error: {e}"
            )
        except Exception as e:
            return ProcessingResult(
                filename=filename,
                success=False,
                message=f"Unexpected error: {e}"
            )
    
    def generate_smarts_pattern(self, reaction_smiles: str, protocol: Dict = None) -> Dict[str, Any]:
        """Generate SMARTS pattern from reaction SMILES focusing on reaction centers.
        
        Prefers atom-mapped SMILES if available (most reliable), falls back to heuristics.
        
        Args:
            reaction_smiles: Reaction SMILES string
            protocol: Full protocol dict (optional, for checking mapped SMILES)
            
        Returns:
            Dictionary with 'core' and 'guards_forbid' keys
        """
        # PRIORITY 1: Check for atom-mapped SMILES (most reliable)
        if protocol:
            mapped_smiles = protocol.get("reaction", {}).get("reaction_smiles_mapped")
            mapping_confidence = protocol.get("reaction", {}).get("mapping_confidence", 0.0)
            
            # Use atom-mapped SMILES whenever available - the mapping itself is valuable
            # even if confidence is low. Only skip if confidence is extremely low (<30%)
            if mapped_smiles and mapping_confidence > 0.3:
                if mapping_confidence > 0.7:
                    print(f"      Using atom-mapped SMILES (high confidence: {mapping_confidence:.1%})")
                else:
                    print(f"      Using atom-mapped SMILES (low confidence: {mapping_confidence:.1%}, but better than heuristics)")
                try:
                    return self.generate_smarts_from_mapped_reaction(mapped_smiles, radius=self.radius)
                except Exception as e:
                    print(f"      Warning: Mapped SMILES failed, falling back to heuristics: {e}")
            elif mapped_smiles and mapping_confidence <= 0.7:
                print(f"      Warning: Low mapping confidence ({mapping_confidence:.1%}), using heuristics instead")
        
        # PRIORITY 2: Fall back to heuristics-based detection
        return self.generate_smarts_from_unmapped_reaction(reaction_smiles)
    
    def generate_smarts_from_mapped_reaction(self, mapped_smiles: str, radius: int = 0) -> Dict[str, Any]:
        """
        Generate SMARTS pattern from atom-mapped reaction SMILES using industry-standard methods.
        
        Uses rxnutils ChemicalReaction.generate_reaction_template() when available,
        which provides radius-based control over pattern specificity.
        Falls back to manual pattern building if rxnutils is not available.
        
        Args:
            mapped_smiles: Atom-mapped reaction SMILES (e.g., "[c:1]-[I:2].[C:3]#[C:4]>>[c:1]-[C:3]#[C:4]")
            radius: Radius for reaction template (0=minimal, 1=one neighbor, 2=two neighbors, etc.)
                    Smaller radius = more general pattern
                    Larger radius = more specific pattern
        
        Returns:
            Dict with 'core' (SMARTS pattern) and 'guards_forbid' (list)
        """
        print(f"      Analyzing atom-mapped reaction (radius={radius})...")
        
        # PRIORITY 1: Use rxnutils professional template generation (industry standard)
        if RXNUTILS_AVAILABLE:
            try:
                # Create ChemicalReaction object from atom-mapped SMILES
                rxn = ChemicalReaction(mapped_smiles)
                
                # Generate reaction template with specified radius
                # expand_ring=False: don't expand pattern to include full rings
                # expand_hetero=False: don't expand to include all heteroatoms
                templates = rxn.generate_reaction_template(
                    radius=radius, 
                    expand_ring=False, 
                    expand_hetero=False
                )
                
                if templates and len(templates) > 0:
                    # Get the forward template (reactants>>products)
                    forward_template = templates[0].smarts
                    
                    print(f"      Generated template (rxnutils): {forward_template}")
                    
                    # Try to simplify the template if radius=0 (minimal pattern)
                    if radius == 0:
                        simplified = self._simplify_rxnutils_template(forward_template)
                        if simplified != forward_template:
                            print(f"      Simplified to: {simplified}")
                            forward_template = simplified
                    
                    return {
                        'core': forward_template,
                        'guards_forbid': [],
                        'method': 'rxnutils',
                        'radius': radius
                    }
                else:
                    print(f"      [WARNING] rxnutils returned no templates, using fallback")
            except Exception as e:
                print(f"      [WARNING] rxnutils failed: {e}, using fallback")
        
        # PRIORITY 2: Fallback to manual pattern building
        print(f"      Using fallback manual pattern generation...")
        
        # Use reaction center detector to find changed atoms
        change_info = identify_changed_atoms_from_mapped_smiles(mapped_smiles)
        
        if not change_info:
            print(f"      [WARNING] Could not identify reaction center from mapping, falling back to heuristics")
            # Remove atom mapping and use heuristics
            import re
            unmapped = re.sub(r':\d+', '', mapped_smiles)
            return self.generate_smarts_from_unmapped_reaction(unmapped)
        
        changed_atoms = change_info.get('changed_atoms', [])
        broken_bonds = change_info.get('broken_bonds', [])
        formed_bonds = change_info.get('formed_bonds', [])
        
        print(f"      Changed atoms (by atom map): {sorted(changed_atoms)}")
        print(f"      Broken bonds: {len(broken_bonds)}")
        print(f"      Formed bonds: {len(formed_bonds)}")
        
        # Split reaction into reactants and products
        parts = mapped_smiles.split('>>')
        if len(parts) != 2:
            print(f"      [ERROR] Invalid reaction SMILES format")
            return {'core': 'INVALID', 'guards_forbid': []}
        
        reactant_smiles = parts[0]
        product_smiles = parts[1]
        
        # Build reactant pattern focusing on changed atoms
        reactant_pattern = self._build_mapped_reactant_pattern(reactant_smiles, changed_atoms, broken_bonds)
        
        # Build product pattern focusing on formed bonds
        product_pattern = self._build_mapped_product_pattern(product_smiles, changed_atoms, formed_bonds)
        
        # Combine into full pattern
        smarts_pattern = f"{reactant_pattern}>>{product_pattern}"
        
        print(f"      Generated SMARTS: {smarts_pattern}")
        
        # Return with basic guards (can be enhanced later)
        return {
            'core': smarts_pattern,
            'guards_forbid': [],  # Atom mapping provides exact specificity, minimal guards needed
            'method': 'fallback'
        }
    
    def _simplify_rxnutils_template(self, template: str) -> str:
        """
        Simplify rxnutils-generated template to more general patterns.
        
        Rxnutils templates can be very specific (e.g., [CH;D2;+0:1]).
        This method simplifies common patterns for better generalization.
        
        Args:
            template: SMARTS template from rxnutils
        
        Returns:
            Simplified SMARTS template
        """
        import re
        
        simplified = template
        
        # Simplify aromatic carbon specifications
        # [c;D2;+0:1] or [cH;D2;+0:1] -> [c:1]
        simplified = re.sub(r'\[c[H]?;[^\]]+:(\d+)\]', r'[c:\1]', simplified)
        
        # Simplify aliphatic carbon specifications
        # [CH2;D2;+0:1] -> [C:1] or [CH2:1]
        simplified = re.sub(r'\[CH?[0-3]?;[^\]]+:(\d+)\]', r'[C:\1]', simplified)
        
        # Simplify nitrogen specifications
        # [NH2;D1;+0:1] -> [NH2:1] or [N:1]
        simplified = re.sub(r'\[N[H]?[0-3]?;[^\]]+:(\d+)\]', r'[N:\1]', simplified)
        
        # Simplify oxygen specifications
        # [O;D2;+0:1] -> [O:1]
        simplified = re.sub(r'\[O;[^\]]+:(\d+)\]', r'[O:\1]', simplified)
        
        # Simplify halogen specifications
        # [I;D1;+0:1] -> [I:1]
        for halogen in ['I', 'Br', 'Cl', 'F']:
            simplified = re.sub(rf'\[{halogen};[^\]]+:(\d+)\]', rf'[{halogen}:\1]', simplified)
        
        # Optionally: Remove atom mapping numbers for ultra-general patterns
        # (Keep them for now to maintain exact atom correspondence)
        
        return simplified
    
    def _build_mapped_reactant_pattern(self, reactant_smiles: str, changed_atoms: List[int], 
                                       broken_bonds: List[tuple]) -> str:
        """
        Build reactant side SMARTS pattern focusing on atoms involved in the reaction.
        
        Args:
            reactant_smiles: Reactant SMILES with atom mapping
            changed_atoms: List of atom map numbers that change
            broken_bonds: List of (atom_map1, atom_map2) tuples for broken bonds
        
        Returns:
            SMARTS pattern for reactant side (e.g., "c-[I].C#C")
        """
        import re
        
        # Parse multiple reactants (separated by dots)
        reactant_mols = reactant_smiles.split('.')
        patterns = []
        
        for mol_smiles in reactant_smiles.split('.'):
            mol_smiles = mol_smiles.strip()
            
            # Find atom map numbers in this molecule
            atom_maps = [int(m) for m in re.findall(r':(\d+)', mol_smiles)]
            
            # Check if this molecule contains changed atoms
            has_changed_atoms = any(am in changed_atoms for am in atom_maps)
            
            if not has_changed_atoms:
                # Skip spectator molecules
                print(f"         Skipping spectator molecule (no changed atoms)")
                continue
            
            print(f"         Processing reactant with {len(atom_maps)} atoms, changed: {[am for am in atom_maps if am in changed_atoms]}")
            
            # Extract pattern by identifying the reactive functional group
            mol_pattern = self._identify_reactive_group_in_reactant(mol_smiles, changed_atoms)
            patterns.append(mol_pattern)
        
        # Combine with dots
        result = '.'.join(patterns) if patterns else re.sub(r':\d+', '', reactant_smiles)
        print(f"         Reactant pattern: {result}")
        return result
    
    def _identify_reactive_group_in_reactant(self, smiles: str, changed_atoms: List[int]) -> str:
        """
        Identify the reactive functional group in a reactant based on changed atoms.
        
        Args:
            smiles: SMILES with atom mapping
            changed_atoms: List of atom map numbers that participate in reaction
        
        Returns:
            SMARTS pattern for the reactive group
        """
        import re
        
        # Remove atom mapping to get clean SMILES
        clean = re.sub(r':\d+', '', smiles)
        
        # Check for common reactive groups
        
        # Aryl/alkyl halides
        if 'I' in clean or 'Br' in clean or 'Cl' in clean or 'F' in clean:
            # Check if aromatic or aliphatic
            if 'c' in clean or '[c' in clean:
                # Aryl halide
                if 'I' in clean:
                    return 'c-[I]'
                elif 'Br' in clean:
                    return 'c-[Br]'
                elif 'Cl' in clean:
                    return 'c-[Cl]'
                else:
                    return 'c-[F]'
            else:
                # Alkyl halide
                if 'I' in clean:
                    return '[CX4]-[I]'
                elif 'Br' in clean:
                    return '[CX4]-[Br]'
                elif 'Cl' in clean:
                    return '[CX4]-[Cl]'
                else:
                    return '[CX4]-[F]'
        
        # Terminal or internal alkyne - just return C#C (most general)
        # Check for # (triple bond) - works for C#C, [C]#[C], etc.
        if '#' in clean and 'N' not in clean:  # Has triple bond but not nitrile
            return 'C#C'
        
        # Boronic acid/ester
        if 'B' in clean:
            if 'OC(C)(C)C(C)(C)O' in clean:
                return '[B]1OC(C)(C)C(C)(C)O1'  # Pinacol
            elif 'B(O)O' in clean:
                return '[B](O)O'  # Boronic acid
            else:
                return '[B]'  # Generic boron
        
        # Amine
        if '[NH2]' in clean or '[NH]' in clean or 'N' in clean:
            if 'c' in clean:
                return 'c-[NX3]'
            else:
                return '[NX3]'
        
        # Default: simplified pattern
        simplified = re.sub(r'\[CH3\]', 'C', clean)
        simplified = re.sub(r'\[CH2\]', 'C', simplified)
        simplified = re.sub(r'\[CH\]', 'C', simplified)
        simplified = re.sub(r'\[cH\]', 'c', simplified)
        return simplified
    
    def _build_mapped_product_pattern(self, product_smiles: str, changed_atoms: List[int],
                                      formed_bonds: List[tuple]) -> str:
        """
        Build product side SMARTS pattern focusing on atoms involved in formed bonds.
        
        Args:
            product_smiles: Product SMILES with atom mapping
            changed_atoms: List of atom map numbers that change
            formed_bonds: List of (atom_map1, atom_map2) tuples for formed bonds
        
        Returns:
            SMARTS pattern for product side (e.g., "c-C#C")
        """
        import re
        
        if not formed_bonds:
            # No formed bonds - just remove atom mapping
            return re.sub(r':\d+', '', product_smiles)
        
        # Analyze formed bonds to determine reaction type
        # Remove mapping to get clean SMILES
        clean_smiles = re.sub(r':\d+', '', product_smiles)
        
        # Detect common product patterns (order matters - most specific first!)
        
        # Sonogashira: aryl-alkyne bond (c-C#C) - check FIRST
        # Look for triple bond (works for C#C, [C]#[C], C#[C], etc.)
        if '#' in clean_smiles and 'c' in clean_smiles:
            # Make sure it's alkyne, not nitrile - check what follows the #
            triple_bond_idx = clean_smiles.index('#')
            next_char = clean_smiles[triple_bond_idx + 1] if triple_bond_idx + 1 < len(clean_smiles) else ''
            if next_char != 'N':  # Not #N (nitrile)
                print(f"         Product pattern: Sonogashira (aryl-alkyne)")
                return 'c-C#C'
        
        # Cyanation: C-CN bond
        if '#N' in clean_smiles:
            if 'c' in clean_smiles:
                print(f"         Product pattern: Cyanation (aryl-cyano)")
                return 'c-C#N'
            else:
                print(f"         Product pattern: Cyanation (alkyl-cyano)")
                return 'C-C#N'
        
        # Suzuki: biaryl bond (c-c) - check AFTER Sonogashira/cyanation
        if clean_smiles.count('c') >= 2 and '#' not in clean_smiles:
            print(f"         Product pattern: Suzuki (biaryl)")
            return 'c-c'
        
        # Buchwald-Hartwig: aryl-amine (c-N) - check AFTER triple bonds
        if 'c' in clean_smiles and ('[NH' in clean_smiles or 'N' in clean_smiles) and '#' not in clean_smiles:
            # Make sure it's not an amide carbonyl spectator
            if 'C(=O)N' not in clean_smiles and 'C(N)=O' not in clean_smiles:
                print(f"         Product pattern: Buchwald-Hartwig (aryl-amine)")
                return 'c-[NX3]'
        
        # Default: extract substructure (will simplify)
        print(f"         Product pattern: Default extraction")
        product_pattern = self._extract_substructure_around_bonds(product_smiles, formed_bonds)
        
        return product_pattern
    
    def _extract_substructure_around_bonds(self, smiles: str, bonds: List[tuple]) -> str:
        """
        Extract SMARTS substructure around specified bonds (by atom map numbers).
        
        Args:
            smiles: SMILES with atom mapping
            bonds: List of (atom_map1, atom_map2) tuples
        
        Returns:
            SMARTS pattern focusing on these bonds
        """
        import re
        
        # Get all atom map numbers involved in bonds
        bond_atoms = set()
        for atom1, atom2 in bonds:
            bond_atoms.add(atom1)
            bond_atoms.add(atom2)
        
        # For each bond, extract a minimal pattern
        # This is simplified - ideally would use RDKit to extract substructure
        
        # Remove atom mapping to get base SMILES
        pattern = re.sub(r':\d+', '', smiles)
        
        # Identify common patterns and simplify
        # Aryl halide: c-[I], c-[Br], c-[Cl], c-[F]
        if re.search(r'[ICBrClF]\[?c\]?', pattern):
            # Extract aryl-halide pattern
            if '[I]' in pattern:
                return 'c-[I]'
            elif '[Br]' in pattern:
                return 'c-[Br]'
            elif '[Cl]' in pattern:
                return 'c-[Cl]'
            elif '[F]' in pattern:
                return 'c-[F]'
        
        # Terminal alkyne: C#C (most general)
        if 'C#C' in pattern or '#C' in pattern:
            return 'C#C'
        
        # Alkyl halide patterns
        if '[I]' in pattern:
            # Check context
            if re.search(r'\[CH3\].*\[I\]', pattern):
                return '[CX4;H3]-[I]'  # Primary
            elif re.search(r'\[CH2\].*\[I\]', pattern):
                return '[CX4;H2]-[I]'  # Secondary
            elif re.search(r'\[CH\].*\[I\]', pattern):
                return '[CX4;H1]-[I]'  # Tertiary
            else:
                return '[CX4]-[I]'  # General alkyl iodide
        
        # Boronic esters
        if 'OC' in pattern and 'CO' in pattern and 'B' in pattern:
            if pattern.count('C(C)(C)') >= 2:  # Pinacol pattern
                return '[B]1OC(C)(C)C(C)(C)O1'
        
        # Amine patterns
        if '[NH2]' in pattern or '[NH]' in pattern:
            if 'c' in pattern:
                return 'c-[NX3]'  # Aryl amine
            else:
                return '[NX3]'  # Alkyl amine
        
        # Default: return simplified pattern (remove detailed atom specifications)
        pattern = re.sub(r'\[CH3\]', 'C', pattern)
        pattern = re.sub(r'\[CH2\]', 'C', pattern)
        pattern = re.sub(r'\[CH\]', 'C', pattern)
        pattern = re.sub(r'\[cH\]', 'c', pattern)
        
        return pattern
    
    def generate_smarts_from_unmapped_reaction(self, reaction_smiles: str) -> Dict[str, Any]:
        """Generate SMARTS pattern using heuristics (no atom mapping).
        
        Args:
            reaction_smiles: Unmapped reaction SMILES
            
        Returns:
            Dictionary with 'core' and 'guards_forbid' keys
        """
        # Parse reaction SMILES
        if ">>" in reaction_smiles:
            parts = reaction_smiles.split(">>")
        elif ">" in reaction_smiles:
            parts = reaction_smiles.split(">")
        else:
            raise ValueError(f"Invalid reaction SMILES format: {reaction_smiles}")
        
        reactants = parts[0].strip()
        products = parts[-1].strip()
        
        # Generate patterns for each reactant with focus on reaction centers
        reactant_patterns = []
        reactant_smiles_list = [s.strip() for s in reactants.split('.') if s.strip()]
        reactant_info_list = []
        
        for smiles in reactant_smiles_list:
            try:
                # Get substrate classification to identify reaction center
                substrate_info = classify_substrate(smiles)
                reactant_info_list.append(substrate_info)
                
                # Generate pattern focused on reaction center
                result = self.build_reaction_center_pattern(smiles, substrate_info)
                reactant_patterns.append(result)
            except Exception as e:
                print(f"      Warning: Could not generate pattern for {smiles}: {e}")
                continue
        
        # Generate pattern for product (focus on main product)
        product_smiles_list = [s.strip() for s in products.split('.') if s.strip()]
        product_pattern = None
        if product_smiles_list:
            try:
                # Classify product to focus on reaction center
                product_info = classify_substrate(product_smiles_list[0])
                # Pass both classified info and raw SMILES for better detection
                product_pattern = self.build_product_pattern(product_smiles_list[0], product_info, 
                                                             reactant_info_list, reactant_smiles_list)
            except Exception as e:
                print(f"      Warning: Could not generate product pattern: {e}")
        
        # Combine into reaction SMARTS focusing on reaction centers
        if reactant_patterns:
            # Combine reactant patterns with dots
            reactant_cores = [p['core'] for p in reactant_patterns if p.get('core')]
            core_pattern = '.'.join(reactant_cores) if len(reactant_cores) > 1 else reactant_cores[0]
            
            # Add product pattern if available
            if product_pattern:
                core_pattern = f"{core_pattern}>>{product_pattern}"
            
            # Collect all guard patterns
            all_guards = []
            for rp in reactant_patterns:
                all_guards.extend(rp.get('guards_forbid', []))
            
            # Remove duplicates while preserving order
            unique_guards = []
            seen = set()
            for guard in all_guards:
                guard_clean = guard.split('#')[0].strip()  # Remove comments for comparison
                if guard_clean not in seen:
                    seen.add(guard_clean)
                    unique_guards.append(guard)
            
            return {
                "core": core_pattern,
                "guards_forbid": unique_guards
            }
        else:
            raise ValueError("No patterns could be generated from reactants")
    
    def build_reaction_center_pattern(self, smiles: str, substrate_info) -> Dict[str, Any]:
        """Build pattern focused on the reaction center, not spectator groups.
        
        Args:
            smiles: SMILES string
            substrate_info: SubstrateInfo from classification
            
        Returns:
            Dictionary with 'core' and 'guards_forbid'
        """
        substrate_class = substrate_info.substrate_class
        family = substrate_info.substrate_family
        
        # For halides: focus on C-X bond
        if family == 'halide':
            if 'aryl' in substrate_class or 'heteroaryl' in substrate_class:
                # Aryl halide: focus on aromatic-halogen bond
                if 'iodide' in substrate_class:
                    return {
                        'core': 'c-[I]',
                        'guards_forbid': ['[CX4]-[I]  # Exclude aliphatic halides']
                    }
                elif 'bromide' in substrate_class:
                    return {
                        'core': 'c-[Br]',
                        'guards_forbid': ['[CX4]-[Br]  # Exclude aliphatic halides']
                    }
                elif 'chloride' in substrate_class:
                    return {
                        'core': 'c-[Cl]',
                        'guards_forbid': ['[CX4]-[Cl]  # Exclude aliphatic halides']
                    }
            elif 'alkyl' in substrate_class:
                # Alkyl halide: use existing detailed pattern
                result = build_smarts_with_guards(smiles)
                return result
        
        # For aryl sulfonates (mesylates, tosylates, triflates): focus on C-O-S bond
        # These are pseudo-halides used in cross-coupling
        if 'OS(=O)(=O)' in smiles or 'OS(=O)' in smiles:
            # Detect the type of sulfonate
            if 'OS(=O)(=O)C' in smiles and 'OS(=O)(=O)Cc' not in smiles:
                # Mesylate (methanesulfonate): -OSO2Me
                return {
                    'core': 'c-[OX2]-[SX4](=O)(=O)-[CH3]',
                    'guards_forbid': ['[CX4]-[OX2]-[SX4](=O)(=O)  # Exclude alkyl sulfonates']
                }
            elif 'OS(=O)(=O)C(F)(F)F' in smiles or 'OS(=O)(=O)CF' in smiles:
                # Triflate (trifluoromethanesulfonate): -OSO2CF3
                return {
                    'core': 'c-[OX2]-[SX4](=O)(=O)-[CX4](F)(F)F',
                    'guards_forbid': []
                }
            elif 'OS(=O)(=O)Cc' in smiles or 'OS(=O)(=O)c1ccc(C)cc1' in smiles:
                # Tosylate (toluenesulfonate): -OSO2-p-Tol
                return {
                    'core': 'c-[OX2]-[SX4](=O)(=O)',
                    'guards_forbid': []
                }
            else:
                # General aryl sulfonate
                return {
                    'core': 'c-[OX2]-[SX4](=O)(=O)',
                    'guards_forbid': ['[CX4]-[OX2]-[SX4]  # Exclude alkyl sulfonates']
                }
        
        # For alkynes: focus on terminal alkyne
        if 'alkyne' in substrate_class or any('alkyne' in fg for fg in substrate_info.functional_groups):
            return {
                'core': 'C#C',
                'guards_forbid': []
            }
        
        # For amines: focus on N, not neighboring groups
        if family == 'amine':
            if 'aniline' in substrate_class:
                return {
                    'core': 'c-[NX3;H2;!$(NC=O)]',
                    'guards_forbid': ['[CX4]-[NX3;H2;!$(NC=O)]  # Exclude aliphatic primary amines']
                }
            else:
                # Use existing detailed pattern for aliphatic amines
                result = build_smarts_with_guards(smiles)
                return result
        
        # For boronic acids/esters: focus on B
        if family == 'boron':
            if 'boronic_ester' in substrate_class:
                return {
                    'core': '[B]1OC(C)(C)C(C)(C)O1',  # Pinacol ester pattern
                    'guards_forbid': []
                }
            elif 'boronic_acid' in substrate_class:
                return {
                    'core': '[B](O)O',
                    'guards_forbid': []
                }
            else:
                # Generic boron pattern
                return {
                    'core': '[B]',
                    'guards_forbid': []
                }
        
        # For other cases: fall back to detailed pattern
        result = build_smarts_with_guards(smiles)
        return result
    
    def build_product_pattern(self, smiles: str, product_info, reactant_info_list, reactant_smiles_list=None) -> str:
        """Build product pattern focusing on newly formed bonds.
        
        Args:
            smiles: Product SMILES
            product_info: SubstrateInfo for product
            reactant_info_list: List of SubstrateInfo for reactants
            reactant_smiles_list: List of reactant SMILES (optional, for better detection)
            
        Returns:
            SMARTS pattern string for product
        """
        # Identify what reaction type based on reactants
        reactant_families = [info.substrate_family for info in reactant_info_list]
        reactant_classes = [info.substrate_class for info in reactant_info_list]
        
        # Sonogashira: aryl halide + alkyne -> aryl alkyne
        if 'halide' in reactant_families and any('alkyne' in str(r) or any('alkyne' in fg for fg in info.functional_groups) for r, info in zip(reactant_classes, reactant_info_list)):
            return 'c-C#C'  # Aryl-alkyne bond
        
        # Suzuki: aryl halide/sulfonate + boronic acid/ester -> biaryl
        # Check for halide OR sulfonate (pseudo-halide) + boron
        has_halide_or_sulfonate = ('halide' in reactant_families or 
                                    any('OS(=O)(=O)' in str(info.functional_groups) for info in reactant_info_list) or
                                    any('sulfonate' in str(info.substrate_class).lower() for info in reactant_info_list))
        
        if has_halide_or_sulfonate and 'boron' in reactant_families:
            # Check if product has biaryl
            if product_info.has_aromatic:
                return 'c-c'  # Biaryl bond
            else:
                return 'c-C'  # Aryl-alkyl bond
        
        # Buchwald-Hartwig: aryl halide + amine -> aryl amine
        if 'halide' in reactant_families and 'amine' in reactant_families:
            if 'aniline' in any(str(r) for r in reactant_classes):
                return 'c-[NX3;H1]'  # Aryl-NH-aryl (secondary amine)
            else:
                return 'c-[NX3]'  # Aryl-amine bond
        
        # Borylation: alkyl halide + boron -> alkyl-B
        if 'halide' in reactant_families and 'boron' in reactant_families:
            if any('alkyl' in str(r) for r in reactant_classes):
                if 'boronic_ester_pinacol' in str(product_info.substrate_class):
                    return '[CX4]-[B]1OC(C)(C)C(C)(C)O1'
                else:
                    return '[CX4]-[B]'
        
        # Generic: use first reactant's pattern with modifications
        # This is a fallback - may need refinement
        result = build_smarts_with_guards(smiles)
        return result['core']
    
    def validate_smarts_pattern(self, pattern: Dict[str, Any]) -> Dict[str, Any]:
        """Validate a SMARTS pattern.
        
        Args:
            pattern: Dictionary with 'core' and 'guards_forbid' keys
            
        Returns:
            Dictionary with 'valid' (bool) and optional 'error' (str) keys
        """
        if not HAS_RDKIT:
            # Can't validate without RDKit, assume valid
            return {'valid': True, 'warning': 'RDKit not available, validation skipped'}
        
        try:
            # Validate core pattern
            core = pattern.get('core', '')
            if not core:
                return {'valid': False, 'error': 'Core pattern is empty'}
            
            # Split reaction SMARTS if present
            if '>>' in core:
                reactant_part, product_part = core.split('>>')
                # Validate reactant part
                if reactant_part.strip():
                    mol = Chem.MolFromSmarts(reactant_part.strip())
                    if mol is None:
                        return {'valid': False, 'error': f'Invalid reactant SMARTS: {reactant_part}'}
                # Validate product part
                if product_part.strip():
                    mol = Chem.MolFromSmarts(product_part.strip())
                    if mol is None:
                        return {'valid': False, 'error': f'Invalid product SMARTS: {product_part}'}
            else:
                # Single SMARTS pattern
                mol = Chem.MolFromSmarts(core)
                if mol is None:
                    return {'valid': False, 'error': f'Invalid SMARTS: {core}'}
            
            # Validate guard patterns
            guards = pattern.get('guards_forbid', [])
            for i, guard in enumerate(guards):
                # Remove comments
                guard_clean = guard.split('#')[0].strip()
                if guard_clean:
                    mol = Chem.MolFromSmarts(guard_clean)
                    if mol is None:
                        return {'valid': False, 'error': f'Invalid guard pattern {i+1}: {guard_clean}'}
            
            return {'valid': True}
            
        except Exception as e:
            return {'valid': False, 'error': f'Validation error: {str(e)}'}
    
    def print_summary(self):
        """Print summary of processing results"""
        print("\n" + "=" * 80)
        print("PROCESSING SUMMARY")
        print("=" * 80)
        
        successful = [r for r in self.results if r.success]
        failed = [r for r in self.results if not r.success]
        
        print(f"\n[OK] Successful: {len(successful)}/{len(self.results)}")
        print(f"[FAIL] Failed: {len(failed)}/{len(self.results)}")
        
        if failed:
            print(f"\nFailed files:")
            for result in failed:
                print(f"   - {result.filename}: {result.message}")
        
        if successful:
            print(f"\nSuccessfully processed files:")
            for result in successful:
                status = "would be updated" if self.dry_run else "updated"
                had_existing = "(replaced existing)" if result.old_pattern else "(new pattern)"
                print(f"   - {result.filename} - {status} {had_existing}")
        
        if self.dry_run:
            print("\n[DRY RUN] No files were modified")
            print("   Run without --dry-run to apply changes")


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        prog="batch-update-protocol-smarts",
        description="Batch generate and update reaction_smarts_applicability in protocol JSON files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all protocols in default directory (data/protocol_db)
  python -m chemtools.protocol.batch_update_protocol_smarts
  
  # Dry run - see what would change without modifying files
  python -m chemtools.protocol.batch_update_protocol_smarts --dry-run
  
  # Process protocols in custom directory
  python -m chemtools.protocol.batch_update_protocol_smarts --protocol-dir path/to/protocols
  
  # Dry run with custom directory
  python -m chemtools.protocol.batch_update_protocol_smarts --protocol-dir path/to/protocols --dry-run

Output format:
  Each protocol JSON file will have reaction_smarts_applicability added/updated:
  
  "reaction": {
    "reaction_smiles": "...",
    "reaction_smarts_applicability": {
      "core": "[CX4;H2,H3]-[I]>>[B]1OC(C)(C)C(C)(C)O1",
      "guards_forbid": [
        "[CX4;H1]-[I]  # Exclude secondary",
        "[CX4;H0]-[I]  # Exclude tertiary",
        ...
      ]
    }
  }
"""
    )
    
    parser.add_argument(
        "--protocol-dir",
        type=Path,
        default=Path("data/protocol_db"),
        help="Directory containing protocol JSON files (default: data/protocol_db)"
    )
    
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Dry run mode - show what would change without modifying files"
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Verbose output with detailed pattern information"
    )
    
    parser.add_argument(
        "--radius",
        type=int,
        default=0,
        choices=[0, 1, 2, 3],
        help="Radius for reaction template generation (0=minimal/general, 1-3=more specific). "
             "Only applies when using rxnutils with atom-mapped SMILES. "
             "Smaller radius = more general pattern, larger = more specific. (default: 0)"
    )
    
    args = parser.parse_args()
    
    # Show header
    print("\n" + "=" * 80)
    print("Protocol SMARTS Batch Updater")
    print("=" * 80)
    
    if args.dry_run:
        print("\nRunning in DRY RUN mode - no files will be modified")
    
    if args.radius > 0:
        print(f"\nUsing radius={args.radius} for template generation (more specific patterns)")
    else:
        print(f"\nUsing radius=0 for template generation (minimal/general patterns)")
    
    # Create updater and process
    updater = ProtocolSmartsUpdater(
        protocol_dir=args.protocol_dir,
        dry_run=args.dry_run,
        radius=args.radius
    )
    
    try:
        results = updater.process_all_protocols()
        updater.print_summary()
        
        # Exit with error code if any failed
        failed_count = len([r for r in results if not r.success])
        return 0 if failed_count == 0 else 1
        
    except KeyboardInterrupt:
        print("\n\n⚠️  Interrupted by user")
        return 1
    except Exception as e:
        print(f"\n❌ Fatal error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())

"""
Reaction center identification using atom mapping and change detection.

Identifies which atoms/bonds change in a reaction to focus SMARTS patterns
on the reaction center, not spectator groups.
"""

from typing import Dict, List, Set, Tuple, Optional, Any
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


def identify_changed_atoms_from_mapped_smiles(mapped_smiles: str) -> Dict[str, Any]:
    """Identify reaction center from atom-mapped reaction SMILES.
    
    Args:
        mapped_smiles: Reaction SMILES with atom mapping (e.g., [C:1]#[C:2].[c:3][I:4]>>[C:1]#[C:2]-[c:3])
        
    Returns:
        Dictionary with:
        - changed_atoms: Set of atom map numbers that changed
        - broken_bonds: List of (atom1, atom2) bonds that broke
        - formed_bonds: List of (atom1, atom2) bonds that formed
        - spectator_atoms: Set of atom map numbers that didn't change
    """
    if not HAS_RDKIT:
        return {'error': 'RDKit not available'}
    
    try:
        # Split reaction
        reactants_smiles, products_smiles = mapped_smiles.split('>>')
        
        # Parse molecules
        reactant_mols = [Chem.MolFromSmiles(s) for s in reactants_smiles.split('.')]
        product_mols = [Chem.MolFromSmiles(s) for s in products_smiles.split('.')]
        
        # Build atom map → atom index mapping for reactants
        reactant_map_to_atom = {}
        for mol in reactant_mols:
            for atom in mol.GetAtoms():
                map_num = atom.GetAtomMapNum()
                if map_num > 0:
                    reactant_map_to_atom[map_num] = (atom, mol)
        
        # Build atom map → atom index mapping for products
        product_map_to_atom = {}
        for mol in product_mols:
            for atom in mol.GetAtoms():
                map_num = atom.GetAtomMapNum()
                if map_num > 0:
                    product_map_to_atom[map_num] = (atom, mol)
        
        # Find changed atoms by comparing bonds
        changed_atoms = set()
        broken_bonds = []
        formed_bonds = []
        seen_broken_bonds: Set[Tuple[Any, Any]] = set()
        seen_formed_bonds: Set[Tuple[Any, Any]] = set()

        def _append_unique_bond(
            sink: List[Tuple[Any, Any]],
            seen: Set[Tuple[Any, Any]],
            a: Any,
            b: Any,
        ) -> None:
            """Append a bond once using canonical ordering for mapped atom pairs."""
            if isinstance(a, int) and isinstance(b, int):
                key = (min(a, b), max(a, b))
                record = key
            else:
                key = (a, b)
                record = (a, b)
            if key in seen:
                return
            seen.add(key)
            sink.append(record)
        
        # Check all mapped atoms
        all_map_nums = set(reactant_map_to_atom.keys()) | set(product_map_to_atom.keys())
        
        # Also track leaving groups and joining groups (unmapped atoms)
        leaving_groups = []  # Bonds to unmapped atoms in reactants
        joining_groups = []  # Bonds to unmapped atoms in products
        
        for map_num in all_map_nums:
            if map_num not in reactant_map_to_atom:
                # Atom only in product (formed)
                changed_atoms.add(map_num)
            elif map_num not in product_map_to_atom:
                # Atom only in reactant (consumed)
                changed_atoms.add(map_num)
            else:
                # Atom in both - check if bonds changed
                r_atom, r_mol = reactant_map_to_atom[map_num]
                p_atom, p_mol = product_map_to_atom[map_num]
                
                # Get neighbors by map number (mapped atoms only)
                r_neighbors = set()
                for bond in r_atom.GetBonds():
                    other = bond.GetOtherAtom(r_atom)
                    other_map = other.GetAtomMapNum()
                    if other_map > 0:
                        r_neighbors.add((other_map, bond.GetBondType()))
                    else:
                        # Unmapped neighbor = leaving group
                        leaving_groups.append((
                            map_num,
                            other.GetSymbol(),
                            str(bond.GetBondType())
                        ))
                
                p_neighbors = set()
                for bond in p_atom.GetBonds():
                    other = bond.GetOtherAtom(p_atom)
                    other_map = other.GetAtomMapNum()
                    if other_map > 0:
                        p_neighbors.add((other_map, bond.GetBondType()))
                    else:
                        # Unmapped neighbor = joining group
                        joining_groups.append((
                            map_num,
                            other.GetSymbol(),
                            str(bond.GetBondType())
                        ))
                
                # If neighbors changed, this atom is in reaction center
                if r_neighbors != p_neighbors:
                    changed_atoms.add(map_num)
                    
                    # Identify specific bond changes
                    r_neighbor_maps = {n[0] for n in r_neighbors}
                    p_neighbor_maps = {n[0] for n in p_neighbors}
                    
                    # Broken bonds
                    for neighbor_map in r_neighbor_maps - p_neighbor_maps:
                        _append_unique_bond(
                            broken_bonds,
                            seen_broken_bonds,
                            map_num,
                            neighbor_map,
                        )
                    
                    # Formed bonds
                    for neighbor_map in p_neighbor_maps - r_neighbor_maps:
                        _append_unique_bond(
                            formed_bonds,
                            seen_formed_bonds,
                            map_num,
                            neighbor_map,
                        )
        
        # Add leaving/joining groups to broken/formed bonds
        # Format: (mapped_atom, "X" where X is element symbol)
        for leaving in leaving_groups:
            map_num, element, bond_type = leaving
            _append_unique_bond(
                broken_bonds,
                seen_broken_bonds,
                map_num,
                f"{element} (leaving group)",
            )
            changed_atoms.add(map_num)
        
        for joining in joining_groups:
            map_num, element, bond_type = joining
            _append_unique_bond(
                formed_bonds,
                seen_formed_bonds,
                map_num,
                f"{element} (joining group)",
            )
            changed_atoms.add(map_num)
        
        spectator_atoms = all_map_nums - changed_atoms
        
        return {
            'changed_atoms': changed_atoms,
            'broken_bonds': broken_bonds,
            'formed_bonds': formed_bonds,
            'spectator_atoms': spectator_atoms,
            'leaving_groups': leaving_groups,  # [(map_num, element, bond_type), ...]
            'joining_groups': joining_groups,  # [(map_num, element, bond_type), ...]
            'success': True
        }
        
    except Exception as e:
        return {'error': str(e), 'success': False}


def generate_smarts_from_mapped_reaction(mapped_smiles: str) -> Dict[str, Any]:
    """Generate focused SMARTS pattern from atom-mapped reaction SMILES.
    
    Args:
        mapped_smiles: Atom-mapped reaction SMILES
        
    Returns:
        Dictionary with 'core' and 'guards_forbid' keys
    """
    # Identify changed atoms
    change_info = identify_changed_atoms_from_mapped_smiles(mapped_smiles)
    
    if not change_info.get('success'):
        raise ValueError(f"Could not identify reaction center: {change_info.get('error')}")
    
    changed_atoms = change_info['changed_atoms']
    broken_bonds = change_info['broken_bonds']
    formed_bonds = change_info['formed_bonds']
    
    # Extract reactant and product SMILES
    reactants_smiles, products_smiles = mapped_smiles.split('>>')
    
    # Build SMARTS pattern focusing on changed atoms
    # This is a simplified version - could be expanded
    # For now, we'll extract the immediate environment of changed atoms
    
    # TODO: Implement smart SMARTS extraction based on changed atoms
    # For now, return the changed atoms info
    return {
        'changed_atoms': list(changed_atoms),
        'broken_bonds': broken_bonds,
        'formed_bonds': formed_bonds,
        'core': f"Reaction involves {len(changed_atoms)} atoms",
        'guards_forbid': []
    }


def compare_unmapped_reaction_to_find_changes(reaction_smiles: str) -> Dict[str, Any]:
    """Identify likely reaction center by comparing reactants to products.
    
    This function now uses RXNMapper for automatic atom mapping if available.
    Falls back to a helpful error message if RXNMapper is not installed.
    
    Args:
        reaction_smiles: Unmapped reaction SMILES
        
    Returns:
        Dictionary with identified changes (same format as identify_changed_atoms_from_mapped_smiles)
    """
    if not HAS_RDKIT:
        return {'error': 'RDKit not available', 'success': False}
    
    try:
        # Try to use RXNMapper for automatic atom mapping
        try:
            from .._atom_mapping import add_atom_mapping
            
            mapping_result = add_atom_mapping(reaction_smiles)
            
            if mapping_result['success']:
                # Successfully mapped, now analyze
                mapped_smiles = mapping_result['mapped_smiles']
                analysis = identify_changed_atoms_from_mapped_smiles(mapped_smiles)
                
                # Add mapping metadata
                if analysis.get('success'):
                    analysis['auto_mapped'] = True
                    analysis['mapping_confidence'] = mapping_result.get('confidence')
                    analysis['original_smiles'] = reaction_smiles
                    analysis['mapped_smiles'] = mapped_smiles
                
                return analysis
            else:
                # Mapping failed
                return {
                    'success': False,
                    'message': f'Automatic atom mapping failed: {mapping_result.get("error")}',
                    'recommendation': 'Install RXNMapper: pip install rxnmapper'
                }
        
        except ImportError:
            # RXNMapper module not available
            return {
                'success': False,
                'message': 'Automatic reaction center detection without atom mapping requires RXNMapper.',
                'recommendation': 'Install RXNMapper: pip install rxnmapper'
            }
        
    except Exception as e:
        return {'error': str(e), 'success': False}


if __name__ == '__main__':
    # Test with Sonogashira example
    print("Testing Sonogashira coupling with atom mapping")
    print("="*80)
    
    # Manually mapped Sonogashira (simplified)
    mapped = "[c:1]1[cH:2][cH:3][cH:4][cH:5][c:6]1[I:7].[C:8]#[C:9]>>[c:1]1[cH:2][cH:3][cH:4][cH:5][c:6]1[C:8]#[C:9]"
    
    result = identify_changed_atoms_from_mapped_smiles(mapped)
    
    if result.get('success'):
        print(f"Changed atoms: {result['changed_atoms']}")
        print(f"Broken bonds: {result['broken_bonds']}")
        print(f"Formed bonds: {result['formed_bonds']}")
        print(f"Spectator atoms: {result['spectator_atoms']}")
        print("\nConclusion: Focus SMARTS on atoms 6,7,8,9 (the C-I breaking and C-C forming)")
    else:
        print(f"Error: {result.get('error')}")

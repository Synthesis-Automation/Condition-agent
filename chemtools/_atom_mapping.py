"""
Atom mapping helpers for reaction analysis.

This module provides automatic atom mapping using RXNMapper, which enables:
- Bond breaking/formation analysis
- Reaction center identification
- Enhanced reaction type detection

RXNMapper is an optional dependency. If not installed, functions gracefully
degrade to returning unmapped SMILES.

Also provides MCS (Maximum Common Substructure) based bond detection as a
fallback method that works without RXNMapper.
"""

from typing import Any, Dict, List, Optional, Tuple, Set
import logging
import re

logger = logging.getLogger(__name__)

# Import the low-level analysis function
from .util.reaction_center_detector import identify_changed_atoms_from_mapped_smiles

# Lazy-loaded RXNMapper instance (expensive to initialize)
_MAPPER_INSTANCE: Optional[Any] = None
_MAPPER_AVAILABLE: Optional[bool] = None


def is_available() -> bool:
    """
    Check if RXNMapper is available.
    
    Returns:
        True if rxnmapper package is installed and can be imported
    """
    global _MAPPER_AVAILABLE
    
    if _MAPPER_AVAILABLE is not None:
        return _MAPPER_AVAILABLE
    
    try:
        import rxnmapper  # type: ignore
        _ = rxnmapper
        _MAPPER_AVAILABLE = True
        logger.debug("RXNMapper is available")
        return True
    except Exception as e:
        _MAPPER_AVAILABLE = False
        logger.debug(f"RXNMapper not available: {e}")
        return False


def _get_mapper() -> Optional[Any]:
    """
    Get or create the RXNMapper instance.
    
    Lazy initialization to avoid overhead if not needed.
    
    Returns:
        RXNMapper instance or None if unavailable
    """
    global _MAPPER_INSTANCE
    
    if not is_available():
        return None
    
    if _MAPPER_INSTANCE is not None:
        return _MAPPER_INSTANCE
    
    try:
        from rxnmapper import RXNMapper  # type: ignore
        logger.info("Initializing RXNMapper (this may take a few seconds)...")
        _MAPPER_INSTANCE = RXNMapper()
        logger.info("RXNMapper initialized successfully")
        return _MAPPER_INSTANCE
    except Exception as e:
        logger.warning(f"Failed to initialize RXNMapper: {e}")
        return None


def _has_atom_mapping(smiles: str) -> bool:
    """
    Check if SMILES string contains atom mapping numbers.
    
    Atom mapping is indicated by :N notation in bracket atoms (e.g., [C:1])
    
    Args:
        smiles: SMILES string to check
    
    Returns:
        True if atom mapping detected
    """
    # Look for pattern like [C:1] or [N:2]
    return bool(re.search(r'\[\w+:\d+\]', smiles))


def add_atom_mapping(reaction_smiles: str, force: bool = False) -> Dict[str, Any]:
    """
    Add atom mapping to unmapped reaction SMILES using RXNMapper.
    
    Args:
        reaction_smiles: Unmapped reaction SMILES (reactants>>products)
        force: If True, remap even if already has atom mapping
    
    Returns:
        Dictionary with:
        - mapped_smiles: Atom-mapped reaction SMILES
        - original_smiles: Original input SMILES
        - confidence: Mapping confidence (0-1) if available
        - method: "rxnmapper" or "passthrough"
        - success: True if mapping succeeded
        - error: Error message if failed
    """
    # Quick check if already mapped (has :N atom mapping)
    if not force and _has_atom_mapping(reaction_smiles):
        # Appears to already have atom mapping
        logger.debug("Reaction appears to already have atom mapping")
        return {
            'mapped_smiles': reaction_smiles,
            'original_smiles': reaction_smiles,
            'method': 'passthrough',
            'success': True,
            'already_mapped': True
        }
    
    mapper = _get_mapper()
    if mapper is None:
        logger.debug("RXNMapper not available, returning unmapped SMILES")
        return {
            'mapped_smiles': reaction_smiles,
            'original_smiles': reaction_smiles,
            'method': 'unavailable',
            'success': False,
            'error': 'RXNMapper not installed or failed to initialize'
        }
    
    try:
        # RXNMapper expects a list of reaction SMILES
        results = mapper.get_attention_guided_atom_maps([reaction_smiles])
        
        if not results or len(results) == 0:
            return {
                'mapped_smiles': reaction_smiles,
                'original_smiles': reaction_smiles,
                'method': 'rxnmapper',
                'success': False,
                'error': 'RXNMapper returned no results'
            }
        
        result = results[0]
        
        # Extract mapped SMILES and confidence
        # RXNMapper returns dict with 'mapped_rxn' and optionally 'confidence'
        if isinstance(result, dict):
            mapped = result.get('mapped_rxn') or result.get('mapped_smiles')
            confidence = result.get('confidence')
        else:
            # Fallback if result is just a string
            mapped = str(result)
            confidence = None
        
        if not mapped:
            return {
                'mapped_smiles': reaction_smiles,
                'original_smiles': reaction_smiles,
                'method': 'rxnmapper',
                'success': False,
                'error': 'RXNMapper returned empty mapping'
            }
        
        logger.debug(f"Successfully mapped reaction (confidence: {confidence})")
        
        return {
            'mapped_smiles': mapped,
            'original_smiles': reaction_smiles,
            'confidence': confidence,
            'method': 'rxnmapper',
            'success': True
        }
        
    except Exception as e:
        logger.warning(f"RXNMapper failed to map reaction: {e}")
        return {
            'mapped_smiles': reaction_smiles,
            'original_smiles': reaction_smiles,
            'method': 'rxnmapper',
            'success': False,
            'error': str(e)
        }


def batch_add_atom_mapping(reaction_smiles_list: List[str], force: bool = False) -> List[Dict[str, Any]]:
    """
    Add atom mapping to multiple reactions in batch.
    
    More efficient than calling add_atom_mapping() individually when processing
    many reactions.
    
    Args:
        reaction_smiles_list: List of unmapped reaction SMILES
        force: If True, remap even if already has atom mapping
    
    Returns:
        List of mapping result dictionaries (same format as add_atom_mapping)
    """
    if not reaction_smiles_list:
        return []
    
    mapper = _get_mapper()
    if mapper is None:
        # Return passthrough results
        return [
            {
                'mapped_smiles': smiles,
                'original_smiles': smiles,
                'method': 'unavailable',
                'success': False,
                'error': 'RXNMapper not installed or failed to initialize'
            }
            for smiles in reaction_smiles_list
        ]
    
    try:
        # Process in batch
        results = mapper.get_attention_guided_atom_maps(reaction_smiles_list)
        
        output = []
        for i, (original, result) in enumerate(zip(reaction_smiles_list, results)):
            # Extract mapped SMILES and confidence
            if isinstance(result, dict):
                mapped = result.get('mapped_rxn') or result.get('mapped_smiles')
                confidence = result.get('confidence')
            else:
                mapped = str(result)
                confidence = None
            
            if not mapped:
                output.append({
                    'mapped_smiles': original,
                    'original_smiles': original,
                    'method': 'rxnmapper',
                    'success': False,
                    'error': f'RXNMapper returned empty mapping for reaction {i}'
                })
            else:
                output.append({
                    'mapped_smiles': mapped,
                    'original_smiles': original,
                    'confidence': confidence,
                    'method': 'rxnmapper',
                    'success': True
                })
        
        return output
        
    except Exception as e:
        logger.warning(f"Batch RXNMapper failed: {e}")
        # Return passthrough results
        return [
            {
                'mapped_smiles': smiles,
                'original_smiles': smiles,
                'method': 'rxnmapper',
                'success': False,
                'error': str(e)
            }
            for smiles in reaction_smiles_list
        ]


def analyze_bond_changes(reaction_smiles: str, auto_map: bool = True) -> Dict[str, Any]:
    """
    Analyze which bonds break and form in a reaction.
    
    This is a high-level convenience function that combines atom mapping
    with reaction center detection.
    
    Args:
        reaction_smiles: Reaction SMILES (mapped or unmapped)
        auto_map: If True, automatically add atom mapping if needed
    
    Returns:
        Dictionary with:
        - broken_bonds: List of (atom1, atom2) bonds that broke
        - formed_bonds: List of (atom1, atom2) bonds that formed
        - changed_atoms: Set of atom map numbers that changed
        - spectator_atoms: Set of atom map numbers unchanged
        - mapped_smiles: Atom-mapped SMILES used for analysis
        - mapping_confidence: Confidence of atom mapping (if auto-mapped)
        - success: True if analysis succeeded
        - error: Error message if failed
    """
    from .util.reaction_center_detector import identify_changed_atoms_from_mapped_smiles
    
    mapped_smiles = reaction_smiles
    mapping_confidence = None
    
    # Auto-map if requested and needed
    if auto_map:
        mapping_result = add_atom_mapping(reaction_smiles)
        if mapping_result['success']:
            mapped_smiles = mapping_result['mapped_smiles']
            mapping_confidence = mapping_result.get('confidence')
        elif mapping_result.get('already_mapped'):
            # Already mapped, use as-is
            pass
        else:
            # Mapping failed, try analysis anyway
            logger.warning(f"Atom mapping failed: {mapping_result.get('error')}")
    
    # Analyze bond changes
    result = identify_changed_atoms_from_mapped_smiles(mapped_smiles)
    
    if result.get('success'):
        result['mapped_smiles'] = mapped_smiles
        result['mapping_confidence'] = mapping_confidence
        return result
    else:
        return {
            'success': False,
            'error': result.get('error', 'Unknown error'),
            'mapped_smiles': mapped_smiles,
            'mapping_confidence': mapping_confidence
        }


def analyze_with_mcs(reaction_smiles: str) -> Dict[str, Any]:
    """
    Analyze bond changes using Maximum Common Substructure (MCS) approach.
    
    This is a fallback method that works without atom mapping tools.
    Uses RDKit's FindMCS to identify unchanged atoms, then infers bond changes.
    
    Args:
        reaction_smiles: Unmapped reaction SMILES
    
    Returns:
        Dictionary with:
        - changed_atoms_estimated: Estimated number of changed atoms
        - likely_broken_bonds: Estimated number of broken bonds
        - likely_formed_bonds: Estimated number of formed bonds
        - mcs_size: Size of maximum common substructure
        - method: "mcs"
        - success: True if analysis succeeded
        - confidence: Lower than RXNMapper (0.6-0.8 range)
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdFMCS
    except ImportError:
        return {
            'success': False,
            'error': 'RDKit not available',
            'method': 'mcs'
        }
    
    try:
        # Parse reaction
        if '>>' not in reaction_smiles:
            return {
                'success': False,
                'error': 'Invalid reaction SMILES (missing >>)',
                'method': 'mcs'
            }
        
        reactants_smiles, products_smiles = reaction_smiles.split('>>')
        
        # Remove atom mapping if present (MCS works on structure, not mapping)
        def strip_atom_mapping(smiles: str) -> str:
            """Remove :N atom mapping from SMILES."""
            return re.sub(r':\d+', '', smiles)
        
        reactants_smiles = strip_atom_mapping(reactants_smiles)
        products_smiles = strip_atom_mapping(products_smiles)
        
        # Combine all reactants and products
        reactant_mols = [Chem.MolFromSmiles(s) for s in reactants_smiles.split('.') if s]
        product_mols = [Chem.MolFromSmiles(s) for s in products_smiles.split('.') if s]
        
        # Filter out None (failed parses) - often happens with malformed atom-mapped SMILES
        reactant_mols = [m for m in reactant_mols if m is not None]
        product_mols = [m for m in product_mols if m is not None]
        
        if not reactant_mols or not product_mols:
            return {
                'success': False,
                'error': 'Could not parse reactants or products',
                'method': 'mcs'
            }
        
        # Combine reactants into single molecule (for MCS comparison)
        # This is a simplification - MCS works best on similar molecules
        reactant_mol = reactant_mols[0] if len(reactant_mols) == 1 else Chem.CombineMols(*reactant_mols)
        product_mol = product_mols[0] if len(product_mols) == 1 else Chem.CombineMols(*product_mols)
        
        # Find Maximum Common Substructure
        mcs_result = rdFMCS.FindMCS(
            [reactant_mol, product_mol],
            timeout=10,  # 10 seconds max
            bondCompare=rdFMCS.BondCompare.CompareOrder,  # Consider bond order
            atomCompare=rdFMCS.AtomCompare.CompareElements,  # Match elements
        )
        
        if not mcs_result or mcs_result.numAtoms == 0:
            # No common substructure - major rearrangement
            return {
                'success': True,
                'method': 'mcs',
                'changed_atoms_estimated': reactant_mol.GetNumAtoms(),
                'likely_broken_bonds': reactant_mol.GetNumBonds(),
                'likely_formed_bonds': product_mol.GetNumBonds(),
                'mcs_size': 0,
                'confidence': 0.3,  # Low confidence for major changes
                'interpretation': 'Major rearrangement or completely different structures',
                'warning': 'MCS found no common structure - results highly uncertain'
            }
        
        # Calculate changes
        total_reactant_atoms = reactant_mol.GetNumAtoms()
        total_product_atoms = product_mol.GetNumAtoms()
        mcs_atoms = mcs_result.numAtoms
        mcs_bonds = mcs_result.numBonds
        
        # Estimate changed atoms (atoms not in MCS)
        changed_atoms_r = total_reactant_atoms - mcs_atoms
        changed_atoms_p = total_product_atoms - mcs_atoms
        changed_atoms_est = max(changed_atoms_r, changed_atoms_p)
        
        # Estimate bond changes
        total_reactant_bonds = reactant_mol.GetNumBonds()
        total_product_bonds = product_mol.GetNumBonds()
        
        # Rough estimate: bonds not in MCS are likely changed
        likely_broken = max(0, total_reactant_bonds - mcs_bonds)
        likely_formed = max(0, total_product_bonds - mcs_bonds)
        
        # Confidence based on MCS coverage
        mcs_coverage = mcs_atoms / max(total_reactant_atoms, total_product_atoms)
        confidence = 0.5 + (mcs_coverage * 0.3)  # 0.5-0.8 range
        
        return {
            'success': True,
            'method': 'mcs',
            'changed_atoms_estimated': changed_atoms_est,
            'likely_broken_bonds': likely_broken,
            'likely_formed_bonds': likely_formed,
            'mcs_size': mcs_atoms,
            'mcs_coverage': mcs_coverage,
            'confidence': confidence,
            'interpretation': _interpret_mcs_changes(changed_atoms_est, likely_broken, likely_formed),
            'warning': 'MCS-based estimates are approximate - exact bonds not identified'
        }
        
    except Exception as e:
        logger.warning(f"MCS analysis failed: {e}")
        return {
            'success': False,
            'error': str(e),
            'method': 'mcs'
        }


def _interpret_mcs_changes(changed_atoms: int, broken: int, formed: int) -> str:
    """Interpret MCS-based change estimates."""
    if changed_atoms == 0:
        return "No structural changes detected"
    elif changed_atoms <= 3:
        return "Minor structural changes (1-3 atoms)"
    elif changed_atoms <= 6:
        return "Moderate structural changes (4-6 atoms)"
    else:
        return "Major structural changes (>6 atoms)"


def analyze_bond_changes_hybrid(
    reaction_smiles: str,
    use_rxnmapper: bool = True,
    use_mcs: bool = True,
    use_manual: bool = True,
    auto_map: bool = True
) -> Dict[str, Any]:
    """
    Hybrid approach: Combine RXNMapper, MCS, and manual mapping for best results.
    
    This provides the best of all worlds:
    - Manual mapping (if available): Ground truth from atom-mapped SMILES
    - RXNMapper: ML-based atom mapping with high accuracy
    - MCS: Graph-based validation and fallback
    - Combined results increase confidence
    
    Priority order:
    1. Manual mapping (if reaction is already atom-mapped)
    2. RXNMapper (ML-based, high accuracy)
    3. MCS (graph-based, approximate)
    
    Args:
        reaction_smiles: Reaction SMILES (mapped or unmapped)
        use_rxnmapper: If True, try RXNMapper (default: True)
        use_mcs: If True, also run MCS analysis (default: True)
        use_manual: If True, check for existing atom mapping first (default: True)
        auto_map: If True, auto-map with RXNMapper if needed (default: True)
    
    Returns:
        Dictionary with:
        - manual_result: Results from manual mapping (if available)
        - rxnmapper_result: Results from RXNMapper (if used)
        - mcs_result: Results from MCS (if used)
        - combined_confidence: Combined confidence score
        - agreement: Dict showing which methods agree
        - recommended_result: Best result to use
        - method: "manual", "hybrid", "rxnmapper_only", or "mcs_only"
    """
    results = {
        'manual_result': None,
        'rxnmapper_result': None,
        'mcs_result': None,
        'combined_confidence': 0.0,
        'agreement': {},
        'method': 'hybrid',
        'success': False
    }
    
    # Check for manual mapping first (highest priority)
    if use_manual and _has_atom_mapping(reaction_smiles):
        try:
            manual_result = identify_changed_atoms_from_mapped_smiles(reaction_smiles)
            results['manual_result'] = manual_result
            
            if manual_result.get('success'):
                logger.info("Manual atom mapping detected - using as ground truth")
                # Manual mapping gets highest confidence (1.0)
                manual_result['mapping_confidence'] = 1.0
                manual_result['method'] = 'manual'
        except Exception as e:
            logger.warning(f"Manual mapping analysis failed: {e}")
    
    # Try RXNMapper first (more accurate)
    if use_rxnmapper and is_available():
        try:
            rxn_result = analyze_bond_changes(reaction_smiles, auto_map=auto_map)
            results['rxnmapper_result'] = rxn_result
            
            if rxn_result.get('success'):
                logger.info(f"RXNMapper analysis successful (conf: {rxn_result.get('mapping_confidence', 'N/A')})")
        except Exception as e:
            logger.warning(f"RXNMapper analysis failed: {e}")
    
    # Also try MCS (validation/fallback)
    if use_mcs:
        try:
            mcs_result = analyze_with_mcs(reaction_smiles)
            results['mcs_result'] = mcs_result
            
            if mcs_result.get('success'):
                logger.info(f"MCS analysis successful (conf: {mcs_result.get('confidence', 'N/A')})")
        except Exception as e:
            logger.warning(f"MCS analysis failed: {e}")
    
    # Combine results with priority: Manual > RXNMapper > MCS
    manual_ok = results['manual_result'] and results['manual_result'].get('success')
    rxn_ok = results['rxnmapper_result'] and results['rxnmapper_result'].get('success')
    mcs_ok = results['mcs_result'] and results['mcs_result'].get('success')
    
    # Extract bond counts for comparison
    def get_bond_counts(result):
        if not result or not result.get('success'):
            return None, None
        broken = len(result.get('broken_bonds', []))
        formed = len(result.get('formed_bonds', []))
        return broken, formed
    
    manual_bonds = get_bond_counts(results['manual_result'])
    rxn_bonds = get_bond_counts(results['rxnmapper_result'])
    mcs_bonds = (
        results['mcs_result'].get('likely_broken_bonds'),
        results['mcs_result'].get('likely_formed_bonds')
    ) if mcs_ok else (None, None)
    
    # Check agreements (within 2 bonds tolerance)
    def bonds_agree(bonds1, bonds2, tolerance=2):
        if bonds1[0] is None or bonds2[0] is None:
            return None
        broken_match = abs(bonds1[0] - bonds2[0]) <= tolerance
        formed_match = abs(bonds1[1] - bonds2[1]) <= tolerance
        return broken_match and formed_match
    
    results['agreement'] = {
        'manual_vs_rxnmapper': bonds_agree(manual_bonds, rxn_bonds) if (manual_ok and rxn_ok) else None,
        'manual_vs_mcs': bonds_agree(manual_bonds, mcs_bonds) if (manual_ok and mcs_ok) else None,
        'rxnmapper_vs_mcs': bonds_agree(rxn_bonds, mcs_bonds) if (rxn_ok and mcs_ok) else None,
    }
    
    # Determine best result and confidence
    if manual_ok:
        # Manual mapping is ground truth - highest priority
        results['success'] = True
        results['method'] = 'manual'
        results['recommended_result'] = results['manual_result']
        results['combined_confidence'] = 1.0  # Manual is 100% confidence
        
        # Validate against other methods
        validations = []
        if rxn_ok:
            if results['agreement']['manual_vs_rxnmapper']:
                validations.append('✓ RXNMapper agrees')
            else:
                validations.append('⚠ RXNMapper disagrees - check RXNMapper accuracy')
        if mcs_ok:
            if results['agreement']['manual_vs_mcs']:
                validations.append('✓ MCS agrees')
            else:
                validations.append('⚠ MCS disagrees - expected for approximate method')
        
        results['validation'] = 'Manual mapping (ground truth). ' + ', '.join(validations) if validations else 'Manual mapping (ground truth)'
        
    elif rxn_ok and mcs_ok:
        # Both ML methods succeeded - check agreement
        agreement = results['agreement']['rxnmapper_vs_mcs']
        
        results['success'] = True
        results['method'] = 'hybrid'
        
        # Combined confidence
        rxn_conf = results['rxnmapper_result'].get('mapping_confidence', 0.5)
        mcs_conf = results['mcs_result'].get('confidence', 0.5)
        
        if agreement:
            # Boost confidence if both agree
            results['combined_confidence'] = min(1.0, (rxn_conf + mcs_conf) / 1.5)
        else:
            # Lower confidence if disagree
            results['combined_confidence'] = max(rxn_conf, mcs_conf) * 0.8
        
        # Recommend RXNMapper result (more precise)
        results['recommended_result'] = results['rxnmapper_result']
        results['validation'] = 'Both methods agree' if agreement else 'Methods disagree - review recommended'
        
    elif rxn_ok:
        # Only RXNMapper worked
        results['success'] = True
        results['method'] = 'rxnmapper_only'
        results['recommended_result'] = results['rxnmapper_result']
        results['combined_confidence'] = results['rxnmapper_result'].get('mapping_confidence', 0.5)
        results['validation'] = 'MCS validation unavailable'
        
    elif mcs_ok:
        # Only MCS worked
        results['success'] = True
        results['method'] = 'mcs_only'
        results['recommended_result'] = results['mcs_result']
        results['combined_confidence'] = results['mcs_result'].get('confidence', 0.5)
        results['validation'] = 'RXNMapper unavailable - using MCS estimates'
        
    else:
        # All methods failed
        results['success'] = False
        results['error'] = 'All bond analysis methods failed'
        results['validation'] = 'No methods available'
    
    return results


# Expose convenience function at module level
__all__ = [
    'is_available',
    'add_atom_mapping',
    'batch_add_atom_mapping',
    'analyze_bond_changes',
    'analyze_with_mcs',
    'analyze_bond_changes_hybrid',
]

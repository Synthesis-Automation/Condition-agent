"""
Protocol Recommendation System

Uses DRFP similarity to find the most relevant protocol for a user's reaction.
Similar to the ML-based recommendation system but for experimental protocols.

This module:
1. Loads the protocol index (with DRFP fingerprints)
2. Computes DRFP similarity between query reaction and all protocols
3. Returns top-k most similar protocols in standard JSON format
4. Optionally filters by reaction family and tags
5. Supports SMARTS-based pre-filtering for structural matching

Usage:
    from chemtools.protocol.recommend import ProtocolRecommender
    
    recommender = ProtocolRecommender()
    
    results = recommender.recommend(
        reaction_smiles='CCBr.c1ccccc1B(O)O>>CCc1ccccc1',
        k=5,
        use_smarts_filter=True  # Enable SMARTS-based filtering
    )
    
    # Standard output format with meta, input, detection, recommended_conditions
    for rec in results['recommended_conditions']:
        print(f"Rank {rec['rank']}: {rec['protocol']['title']}")
"""

import logging
import time
from pathlib import Path
from typing import Dict, Any, List, Optional
import json

from .indexer import ProtocolIndexer, ProtocolRecord

try:
    from ..output_formatter import format_meta, format_input, format_detection
    HAS_OUTPUT_FORMATTER = True
except ImportError:
    HAS_OUTPUT_FORMATTER = False

logger = logging.getLogger(__name__)


def match_reaction_smarts(reaction_smiles: str, smarts_patterns: List[str]) -> bool:
    """
    Check if a reaction SMILES matches any of the provided SMARTS patterns.
    
    Uses RDKit's reaction SMARTS matching when available. Falls back to returning
    True if RDKit is not available (permissive fallback).
    
    Args:
        reaction_smiles: Input reaction SMILES string (e.g., "CCBr.c1ccccc1B(O)O>>CCc1ccccc1")
        smarts_patterns: List of reaction SMARTS patterns to match against
    
    Returns:
        True if the reaction matches any pattern, False otherwise
    """
    if not smarts_patterns:
        # No patterns specified, always match
        return True
    
    if not reaction_smiles or '>>' not in reaction_smiles:
        # Invalid reaction SMILES
        return False
    
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        # Parse the input reaction
        try:
            rxn_mol = AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)
            if rxn_mol is None:
                logger.debug(f"Could not parse reaction SMILES: {reaction_smiles}")
                return True  # Permissive fallback
        except Exception as e:
            logger.debug(f"Error parsing reaction SMILES: {e}")
            return True  # Permissive fallback
        
        # Try to match against each SMARTS pattern
        for pattern in smarts_patterns:
            if not pattern or '>>' not in pattern:
                continue
            
            try:
                # Parse the SMARTS pattern as a reaction
                pattern_rxn = AllChem.ReactionFromSmarts(pattern)
                if pattern_rxn is None:
                    logger.debug(f"Could not parse reaction SMARTS: {pattern}")
                    continue
                
                # Check if the input reaction matches the pattern
                # We need to check if the pattern is a substructure of the reaction
                # This is done by comparing reactants and products separately
                
                # Get reactants and products from both
                rxn_reactants = rxn_mol.GetReactants()
                rxn_products = rxn_mol.GetProducts()
                pattern_reactants = pattern_rxn.GetReactants()
                pattern_products = pattern_rxn.GetProducts()
                
                # Sanitize molecules to ensure proper aromaticity perception
                try:
                    for mol in rxn_reactants:
                        Chem.SanitizeMol(mol)
                    for mol in rxn_products:
                        Chem.SanitizeMol(mol)
                except Exception as e:
                    logger.debug(f"Error sanitizing reaction molecules: {e}")
                    continue
                
                # Check if all pattern reactants match some rxn reactants
                reactants_match = all(
                    any(rxn_r.HasSubstructMatch(pat_r) for rxn_r in rxn_reactants)
                    for pat_r in pattern_reactants
                )
                
                # Check if all pattern products match some rxn products
                products_match = all(
                    any(rxn_p.HasSubstructMatch(pat_p) for rxn_p in rxn_products)
                    for pat_p in pattern_products
                )
                
                if reactants_match and products_match:
                    logger.debug(f"Reaction matched SMARTS pattern: {pattern}")
                    return True
                    
            except Exception as e:
                logger.debug(f"Error matching SMARTS pattern {pattern}: {e}")
                continue
        
        # No patterns matched
        return False
        
    except ImportError:
        logger.debug("RDKit not available for SMARTS matching, using permissive fallback")
        return True  # Permissive fallback when RDKit is not available


class ProtocolRecommender:
    """
    Recommend protocols using DRFP similarity search
    
    This class provides ML-based protocol matching:
    - Loads precomputed DRFP fingerprints from index
    - Computes query reaction DRFP
    - Finds most similar protocols using cosine similarity
    - Supports filtering by family and tags
    """
    
    def __init__(
        self,
        index_path: Optional[Path] = None,
        protocol_dir: Optional[Path] = None
    ):
        """
        Initialize recommender
        
        Args:
            index_path: Path to protocol index file
            protocol_dir: Directory containing protocol JSON files
        """
        # Load index
        self.indexer = ProtocolIndexer.load(index_path)
        
        # Check if DRFP is available
        if not self.indexer.metadata.get('has_drfp', False):
            logger.warning(
                "Index does not have DRFP fingerprints. "
                "Rebuild index with: python -m chemtools.protocol.cli build"
            )
        
        # Set protocol directory
        if protocol_dir is None:
            protocol_dir = self.indexer.protocol_dir
        self.protocol_dir = Path(protocol_dir)
        
        logger.info(f"Loaded protocol recommender with {len(self.indexer.records)} protocols")
    
    def recommend(
        self,
        reaction_smiles: str,
        k: int = 5,
        reaction_family: Optional[str] = None,
        tags: Optional[List[str]] = None,
        min_similarity: float = 0.0,
        use_standard_format: bool = True,
        use_smarts_filter: bool = True
    ) -> Dict[str, Any]:
        """
        Find top-k most similar protocols for a reaction
        
        Args:
            reaction_smiles: Query reaction SMILES
            k: Number of recommendations to return
            reaction_family: Optional family filter (e.g., 'Suzuki')
            tags: Optional tag filter (match ANY tag)
            min_similarity: Minimum similarity threshold (0.0-1.0)
            use_standard_format: Return standard JSON format (default: True)
            use_smarts_filter: Use reaction SMARTS for structural pre-filtering (default: True)
        
        Returns:
            Dictionary with standard format:
                - meta: Model metadata and timing
                - input: Query information
                - detection: Reaction type detection
                - recommended_conditions: List of protocol recommendations
                
            Or legacy format if use_standard_format=False:
                - matches: List of protocol matches with similarity scores
                - query: Query information
                - metadata: Search metadata
        """
        start_time = time.time()
        
        # Compute DRFP for query reaction
        try:
            query_drfp = self._compute_drfp(reaction_smiles)
        except Exception as e:
            logger.error(f"Failed to compute DRFP for query: {e}")
            
            if use_standard_format and HAS_OUTPUT_FORMATTER:
                return {
                    'meta': format_meta(
                        model_type='Protocol-DRFP',
                        status='error',
                        processing_time_ms=(time.time() - start_time) * 1000
                    ),
                    'input': format_input(reaction_smiles=reaction_smiles),
                    'detection': format_detection(
                        detected_type='unknown',
                        confidence=None,
                        method='protocol-similarity'
                    ),
                    'recommended_conditions': [],
                    'extras': {'error': str(e)}
                }
            else:
                return {
                    'matches': [],
                    'query': {'reaction_smiles': reaction_smiles},
                    'metadata': {'error': str(e)}
                }
        
        # Get candidate protocols (filter by family/tags if specified)
        candidates = self._get_candidates(reaction_family, tags)
        
        # Apply SMARTS-based filtering if requested
        num_before_smarts = len(candidates)
        smarts_filter_warning = None
        
        if use_smarts_filter:
            candidates = self._filter_by_smarts(reaction_smiles, candidates)
            num_after_smarts = len(candidates)
            
            if num_after_smarts < num_before_smarts:
                removed = num_before_smarts - num_after_smarts
                logger.warning(
                    f"SMARTS filtering removed {removed} protocol(s): "
                    f"{num_before_smarts} → {num_after_smarts} candidates"
                )
                smarts_filter_warning = (
                    f"SMARTS filtering removed {removed} protocol(s) that did not match "
                    f"the reaction structure. {num_after_smarts} of {num_before_smarts} protocols remain. "
                    f"Use --no-smarts-filter to see all protocols ranked by similarity."
                )
            
            if num_after_smarts == 0:
                logger.warning(
                    f"SMARTS filtering eliminated ALL {num_before_smarts} candidate protocol(s). "
                    f"No protocols match the reaction structure pattern. "
                    f"Consider using --no-smarts-filter for DRFP-only similarity matching."
                )
                smarts_filter_warning = (
                    f"No protocols found matching the reaction SMARTS pattern. "
                    f"All {num_before_smarts} candidate(s) were filtered out. "
                    f"Try --no-smarts-filter to see protocols ranked by chemical similarity only."
                )
        
        if not candidates:
            processing_time_ms = (time.time() - start_time) * 1000
            
            if use_standard_format and HAS_OUTPUT_FORMATTER:
                extras = {
                    'num_candidates': 0,
                    'num_total_protocols': len(self.indexer.records),
                    'message': 'No candidates found matching filters'
                }
                if smarts_filter_warning:
                    extras['smarts_filter_warning'] = smarts_filter_warning
                
                return {
                    'meta': format_meta(
                        model_type='Protocol-DRFP',
                        status='success',
                        processing_time_ms=processing_time_ms
                    ),
                    'input': format_input(
                        reaction_smiles=reaction_smiles,
                        requested_reaction_type=reaction_family
                    ),
                    'detection': format_detection(
                        detected_type=reaction_family or 'unknown',
                        confidence=0.0,
                        method='protocol-similarity'
                    ),
                    'recommended_conditions': [],
                    'extras': extras
                }
            else:
                metadata = {
                    'num_candidates': 0,
                    'num_total_protocols': len(self.indexer.records),
                    'message': 'No candidates found matching filters'
                }
                if smarts_filter_warning:
                    metadata['smarts_filter_warning'] = smarts_filter_warning
                
                return {
                    'matches': [],
                    'query': {
                        'reaction_smiles': reaction_smiles,
                        'family': reaction_family,
                        'tags': tags
                    },
                    'metadata': metadata
                }
        
        # Compute similarities
        similarities = []
        
        # Debug: Check query_drfp type
        logger.debug(f"Query DRFP type: {type(query_drfp)}, len: {len(query_drfp) if hasattr(query_drfp, '__len__') else 'N/A'}")
        
        for record in candidates:
            logger.debug(f"Processing candidate: {record.filename}")
            
            # Get DRFP fingerprint from NPZ file (lazy loaded)
            drfp_fp = self.indexer.get_drfp_fingerprint(record.filename)
            
            if drfp_fp is None:
                logger.debug(f"  DRFP is None, skipping")
                continue
            
            logger.debug(f"  Protocol DRFP type: {type(drfp_fp)}, shape: {drfp_fp.shape if hasattr(drfp_fp, 'shape') else 'N/A'}")
            
            try:
                similarity = self._cosine_similarity(query_drfp, drfp_fp)
                logger.debug(f"  Similarity: {similarity:.4f}, threshold: {min_similarity}")
                
                if similarity >= min_similarity:
                    logger.debug(f"  ✓ Added to matches (similarity >= threshold)")
                    similarities.append({
                        'record': record,
                        'similarity': similarity
                    })
                else:
                    logger.debug(f"  ✗ Rejected (similarity < threshold)")
            except Exception as e:
                logger.debug(f"  ERROR computing similarity: {e}")
                import traceback
                logger.debug(traceback.format_exc())
        
        # Sort by similarity (descending)
        similarities.sort(key=lambda x: x['similarity'], reverse=True)
        
        # Take top-k
        top_matches = similarities[:k]
        
        logger.debug(f"Total similarities: {len(similarities)}, top_matches: {len(top_matches)}")
        if top_matches:
            for i, match in enumerate(top_matches):
                logger.debug(f"  Match {i+1}: {match['record'].filename}, similarity: {match['similarity']:.4f}")
        
        processing_time_ms = (time.time() - start_time) * 1000
        
        # Return standard format or legacy format
        if use_standard_format and HAS_OUTPUT_FORMATTER:
            return self._format_standard_output(
                reaction_smiles=reaction_smiles,
                reaction_family=reaction_family,
                tags=tags,
                top_matches=top_matches,
                num_candidates=len(candidates),
                num_before_smarts=num_before_smarts,
                smarts_filter_warning=smarts_filter_warning,
                processing_time_ms=processing_time_ms
            )
        else:
            return self._format_legacy_output(
                reaction_smiles=reaction_smiles,
                reaction_family=reaction_family,
                tags=tags,
                top_matches=top_matches,
                num_candidates=len(candidates),
                num_before_smarts=num_before_smarts,
                smarts_filter_warning=smarts_filter_warning
            )
    
    def _format_standard_output(
        self,
        reaction_smiles: str,
        reaction_family: Optional[str],
        tags: Optional[List[str]],
        top_matches: List[Dict[str, Any]],
        num_candidates: int,
        num_before_smarts: int,
        smarts_filter_warning: Optional[str],
        processing_time_ms: float
    ) -> Dict[str, Any]:
        """Format output in standard JSON format"""
        # Build recommended_conditions array
        recommended_conditions = []
        
        for rank, item in enumerate(top_matches, start=1):
            record = item['record']
            similarity = item['similarity']
            
            # Load full protocol to get chemicals list
            protocol_data = self.get_protocol_details(record.filename)
            chemicals_list = []
            conditions_data = {}
            
            if protocol_data:
                # Extract chemicals from reaction_setup
                reaction_setup = protocol_data.get('reaction_setup', [])
                if reaction_setup and len(reaction_setup) > 0:
                    chemicals_list = reaction_setup[0].get('chemicals', [])
                
                # Extract conditions using existing method
                conditions_data = self.extract_conditions(protocol_data)
            
            # Build protocol entry in standard format
            protocol_rec = {
                'rank': rank,
                'confidence': round(similarity, 4),
                'chemicals': chemicals_list,  # Now includes full chemicals list
                'conditions': conditions_data,
                'source': 'protocol_database',
                'similarity': round(similarity, 4),
                'protocol_metadata': {  # Metadata moved to nested object
                    'filename': record.filename,
                    'title': record.source_title,
                    'journal': record.source_journal,
                    'year': record.source_year,
                    'doi': record.source_doi,
                    'url': record.source_url,
                    'reaction_smiles': record.reaction_smiles,
                    'reaction_smarts': record.reaction_smarts,
                    'reaction_family': record.reaction_family,
                    'tags': record.tags,
                    'notes': record.notes
                }
            }
            
            recommended_conditions.append(protocol_rec)
        
        # Determine detected type from top match or use requested family
        detected_type = reaction_family or 'unknown'
        detection_confidence = None
        
        if top_matches:
            # Use top match's family as detected type
            detected_type = top_matches[0]['record'].reaction_family or detected_type
            detection_confidence = top_matches[0]['similarity']
        
        extras = {
            'num_candidates': num_candidates,
            'num_total_protocols': len(self.indexer.records),
            'num_matches': len(top_matches),
            'smarts_filtering_enabled': True
        }
        
        # Add SMARTS filter warning if applicable
        if smarts_filter_warning:
            extras['smarts_filter_warning'] = smarts_filter_warning
            extras['num_before_smarts_filter'] = num_before_smarts
            extras['num_after_smarts_filter'] = num_candidates
        
        return {
            'meta': format_meta(
                model_type='Protocol-DRFP',
                status='success',
                processing_time_ms=processing_time_ms,
                model_version='1.0.0'
            ),
            'input': format_input(
                reaction_smiles=reaction_smiles,
                requested_reaction_type=reaction_family,
                options={'k': len(top_matches), 'tags': tags} if tags else {'k': len(top_matches)}
            ),
            'detection': format_detection(
                detected_type=detected_type,
                confidence=detection_confidence,
                method='protocol-similarity'
            ),
            'recommended_conditions': recommended_conditions,
            'extras': extras
        }
    
    def _format_legacy_output(
        self,
        reaction_smiles: str,
        reaction_family: Optional[str],
        tags: Optional[List[str]],
        top_matches: List[Dict[str, Any]],
        num_candidates: int,
        num_before_smarts: int,
        smarts_filter_warning: Optional[str]
    ) -> Dict[str, Any]:
        """Format output in legacy format (for backward compatibility)"""
        # Format results
        matches = []
        for item in top_matches:
            record = item['record']
            match = {
                'filename': record.filename,
                'similarity': item['similarity'],
                'reaction_smiles': record.reaction_smiles,
                'reaction_smarts': record.reaction_smarts,
                'reaction_family': record.reaction_family,
                'tags': record.tags,
                'notes': record.notes,
                'source_title': record.source_title,
                'source_journal': record.source_journal,
                'source_year': record.source_year,
                'source_doi': record.source_doi,
                'source_url': record.source_url
            }
            matches.append(match)
        
        metadata = {
            'num_candidates': num_candidates,
            'num_matches': len(top_matches),
            'num_total_protocols': len(self.indexer.records)
        }
        
        # Add SMARTS filter warning if applicable
        if smarts_filter_warning:
            metadata['smarts_filter_warning'] = smarts_filter_warning
            metadata['num_before_smarts_filter'] = num_before_smarts
            metadata['num_after_smarts_filter'] = num_candidates
        
        return {
            'matches': matches,
            'query': {
                'reaction_smiles': reaction_smiles,
                'family': reaction_family,
                'tags': tags
            },
            'metadata': metadata
        }
    
    def _get_candidates(
        self,
        reaction_family: Optional[str],
        tags: Optional[List[str]]
    ) -> List[ProtocolRecord]:
        """Get candidate protocols based on filters"""
        candidates = list(self.indexer.records.values())
        
        # Filter by family
        if reaction_family:
            candidates = [
                r for r in candidates
                if r.reaction_family.lower() == reaction_family.lower()
            ]
        
        # Filter by tags (match ANY tag)
        if tags:
            tag_matches = self.indexer.search_tags(tags)
            tag_filenames = {r.filename for r in tag_matches}
            candidates = [
                r for r in candidates
                if r.filename in tag_filenames
            ]
        
        return candidates
    
    def _filter_by_smarts(
        self,
        reaction_smiles: str,
        candidates: List[ProtocolRecord]
    ) -> List[ProtocolRecord]:
        """
        Filter candidates by matching reaction SMARTS patterns
        
        Args:
            reaction_smiles: Query reaction SMILES
            candidates: List of protocol records to filter
        
        Returns:
            Filtered list of protocols that match the reaction structure
        """
        filtered = []
        
        for record in candidates:
            # If protocol has no SMARTS patterns, include it (permissive)
            if not record.reaction_smarts:
                filtered.append(record)
                continue
            
            # Check if the reaction matches any of the protocol's SMARTS patterns
            if match_reaction_smarts(reaction_smiles, record.reaction_smarts):
                filtered.append(record)
        
        return filtered
    
    def _compute_drfp(self, reaction_smiles: str) -> List[float]:
        """
        Compute DRFP fingerprint for a reaction
        
        Args:
            reaction_smiles: Reaction SMILES string
        
        Returns:
            DRFP fingerprint as list of floats
        """
        try:
            from drfp import DrfpEncoder
            import numpy as np
            
            # Use the correct DrfpEncoder initialization
            encoder = DrfpEncoder()
            fp = encoder.encode([reaction_smiles])[0]
            
            if isinstance(fp, np.ndarray):
                return fp.tolist()
            return list(fp)
            
        except ImportError:
            raise ImportError(
                "DRFP package required for protocol recommendation. "
                "Install with: pip install drfp"
            )
    
    def _cosine_similarity(self, vec1: List[float], vec2: List[float]) -> float:
        """
        Compute cosine similarity between two vectors
        
        Args:
            vec1: First vector
            vec2: Second vector
        
        Returns:
            Cosine similarity (0.0 to 1.0)
        """
        try:
            import numpy as np
            
            v1 = np.array(vec1)
            v2 = np.array(vec2)
            
            dot_product = np.dot(v1, v2)
            norm1 = np.linalg.norm(v1)
            norm2 = np.linalg.norm(v2)
            
            if norm1 == 0 or norm2 == 0:
                return 0.0
            
            similarity = dot_product / (norm1 * norm2)
            
            # Clamp to [0, 1] range
            return max(0.0, min(1.0, float(similarity)))
            
        except ImportError:
            raise ImportError("NumPy required for similarity computation")
    
    def get_protocol_details(self, filename: str) -> Optional[Dict[str, Any]]:
        """
        Load full protocol details from JSON file
        
        Handles both formats:
        - Legacy: filename is the actual JSON file (e.g., 'Suzuki_protocols.json')
        - New: filename includes index (e.g., 'Suzuki_protocols[0]')
        
        Args:
            filename: Protocol filename or identifier (e.g., 'Suzuki_protocols.json' or 'Suzuki_protocols[0]')
        
        Returns:
            Full protocol dictionary or None if not found
        """
        # Parse filename to extract actual file and index
        if '[' in filename and ']' in filename:
            # New format: filename[index]
            base_name = filename.split('[')[0]
            index_str = filename.split('[')[1].split(']')[0]
            try:
                protocol_index = int(index_str)
            except ValueError:
                logger.error(f"Invalid protocol index in filename: {filename}")
                return None
            actual_filename = f"{base_name}.json"
        else:
            # Legacy format: just the filename
            if not filename.endswith('.json'):
                actual_filename = f"{filename}.json"
            else:
                actual_filename = filename
            protocol_index = None
        
        protocol_path = self.protocol_dir / actual_filename
        
        if not protocol_path.exists():
            logger.error(f"Protocol file not found: {actual_filename}")
            return None
        
        try:
            with open(protocol_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            # Handle both array and single protocol formats
            if isinstance(data, list):
                if protocol_index is not None:
                    # Return specific protocol from array
                    if 0 <= protocol_index < len(data):
                        return data[protocol_index]
                    else:
                        logger.error(f"Protocol index {protocol_index} out of range for {actual_filename} (has {len(data)} protocols)")
                        return None
                else:
                    # Return first protocol if no index specified
                    return data[0] if data else None
            else:
                # Single protocol format (legacy)
                return data
                
        except Exception as e:
            logger.error(f"Error loading protocol {filename}: {e}")
            return None
    
    def extract_conditions(self, protocol_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract key reaction conditions from protocol
        
        Args:
            protocol_data: Full protocol dictionary
        
        Returns:
            Simplified conditions dictionary
        """
        conditions = {
            'catalyst': None,
            'ligand': None,
            'base': None,
            'solvent': None,
            'additives': [],
            'temperature_C': None,
            'time_h': None,
            'atmosphere': None
        }
        
        reaction_setup = protocol_data.get('reaction_setup', [])
        if not reaction_setup:
            return conditions
        
        # Get chemicals from first setup step
        chemicals = reaction_setup[0].get('chemicals', [])
        
        for chem in chemicals:
            role = chem.get('role', '')
            name = chem.get('name', '')
            abbrev = chem.get('abbreviation')
            display_name = abbrev if abbrev else name
            
            if role in ('metal_precursor', 'preformed_catalyst', 'co_catalyst'):
                if conditions['catalyst'] is None:
                    conditions['catalyst'] = display_name
                else:
                    conditions['additives'].append(display_name)
            elif role == 'ligand':
                conditions['ligand'] = display_name
            elif role == 'base':
                conditions['base'] = display_name
            elif role == 'solvent':
                conditions['solvent'] = display_name
            elif role == 'additive':
                conditions['additives'].append(display_name)
        
        # Get conditions from first setup step
        cond_list = reaction_setup[0].get('conditions', [])
        if cond_list:
            # Use the main reaction conditions (usually the last or longest step)
            main_cond = cond_list[-1] if len(cond_list) > 1 else cond_list[0]
            conditions['temperature_C'] = main_cond.get('temperature_C')
            conditions['time_h'] = main_cond.get('time_h')
            conditions['atmosphere'] = main_cond.get('atmosphere')
        
        return conditions
    
    def recommend_with_details(
        self,
        reaction_smiles: str,
        k: int = 5,
        use_standard_format: bool = True,
        use_smarts_filter: bool = True,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Recommend protocols and include extracted condition details
        
        NOTE: When use_standard_format=True (default), chemicals and conditions
        are already included in the output, so this method is equivalent to recommend().
        This method is maintained for backward compatibility.
        
        Args:
            reaction_smiles: Query reaction SMILES
            k: Number of recommendations
            use_standard_format: Return standard JSON format (default: True)
            use_smarts_filter: Use reaction SMARTS for structural pre-filtering (default: True)
            **kwargs: Additional arguments for recommend()
        
        Returns:
            Same as recommend() - with chemicals and conditions already included
        """
        # In standard format mode, details are already included
        return self.recommend(
            reaction_smiles,
            k=k,
            use_standard_format=use_standard_format,
            use_smarts_filter=use_smarts_filter,
            **kwargs
        )


def recommend_protocol(
    reaction_smiles: str,
    k: int = 5,
    reaction_family: Optional[str] = None,
    tags: Optional[List[str]] = None,
    index_path: Optional[Path] = None,
    use_smarts_filter: bool = True
) -> Dict[str, Any]:
    """
    Standalone function to recommend protocols
    
    Args:
        reaction_smiles: Query reaction SMILES
        k: Number of recommendations
        reaction_family: Optional family filter
        tags: Optional tag filter
        index_path: Optional index file path
        use_smarts_filter: Use reaction SMARTS for structural pre-filtering (default: True)
    
    Returns:
        Recommendation results dictionary
    
    Example:
        from chemtools.protocol.recommend import recommend_protocol
        
        results = recommend_protocol(
            'CCBr.c1ccccc1B(O)O>>CCc1ccccc1',
            k=3,
            tags=['Suzuki'],
            use_smarts_filter=True
        )
        
        for match in results['matches']:
            print(f"{match['similarity']:.3f}: {match['source_title']}")
    """
    recommender = ProtocolRecommender(index_path=index_path)
    return recommender.recommend(
        reaction_smiles=reaction_smiles,
        k=k,
        reaction_family=reaction_family,
        tags=tags,
        use_smarts_filter=use_smarts_filter
    )

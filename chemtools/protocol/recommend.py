"""
Protocol Recommendation System

Uses DRFP similarity to find the most relevant protocol for a user's reaction.
Similar to the ML-based recommendation system but for experimental protocols.

This module:
1. Loads the protocol index (with DRFP fingerprints)
2. Computes DRFP similarity between query reaction and all protocols
3. Returns top-k most similar protocols
4. Optionally filters by reaction family and tags

Usage:
    from chemtools.protocol.recommend import ProtocolRecommender
    
    recommender = ProtocolRecommender()
    
    results = recommender.recommend(
        reaction_smiles='CCBr.c1ccccc1B(O)O>>CCc1ccccc1',
        k=5
    )
    
    for match in results['matches']:
        print(f"{match['similarity']:.3f}: {match['source_title']}")
"""

import logging
from pathlib import Path
from typing import Dict, Any, List, Optional
import json

from .indexer import ProtocolIndexer, ProtocolRecord

logger = logging.getLogger(__name__)


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
        min_similarity: float = 0.0
    ) -> Dict[str, Any]:
        """
        Find top-k most similar protocols for a reaction
        
        Args:
            reaction_smiles: Query reaction SMILES
            k: Number of recommendations to return
            reaction_family: Optional family filter (e.g., 'Suzuki')
            tags: Optional tag filter (match ANY tag)
            min_similarity: Minimum similarity threshold (0.0-1.0)
        
        Returns:
            Dictionary with:
                - matches: List of protocol matches with similarity scores
                - query: Query information
                - metadata: Search metadata
        """
        # Compute DRFP for query reaction
        try:
            query_drfp = self._compute_drfp(reaction_smiles)
        except Exception as e:
            logger.error(f"Failed to compute DRFP for query: {e}")
            return {
                'matches': [],
                'query': {'reaction_smiles': reaction_smiles},
                'metadata': {'error': str(e)}
            }
        
        # Get candidate protocols (filter by family/tags if specified)
        candidates = self._get_candidates(reaction_family, tags)
        
        if not candidates:
            return {
                'matches': [],
                'query': {
                    'reaction_smiles': reaction_smiles,
                    'family': reaction_family,
                    'tags': tags
                },
                'metadata': {
                    'num_candidates': 0,
                    'num_total_protocols': len(self.indexer.records),
                    'message': 'No candidates found matching filters'
                }
            }
        
        # Compute similarities
        similarities = []
        for record in candidates:
            if record.drfp_fingerprint is None:
                logger.debug(f"Skipping {record.filename} (no DRFP)")
                continue
            
            try:
                similarity = self._cosine_similarity(query_drfp, record.drfp_fingerprint)
                
                if similarity >= min_similarity:
                    similarities.append({
                        'record': record,
                        'similarity': similarity
                    })
            except Exception as e:
                logger.debug(f"Error computing similarity for {record.filename}: {e}")
        
        # Sort by similarity (descending)
        similarities.sort(key=lambda x: x['similarity'], reverse=True)
        
        # Take top-k
        top_matches = similarities[:k]
        
        # Format results
        matches = []
        for item in top_matches:
            record = item['record']
            match = {
                'filename': record.filename,
                'similarity': item['similarity'],
                'reaction_smiles': record.reaction_smiles,
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
        
        return {
            'matches': matches,
            'query': {
                'reaction_smiles': reaction_smiles,
                'family': reaction_family,
                'tags': tags
            },
            'metadata': {
                'num_candidates': len(candidates),
                'num_matches': len(similarities),
                'num_total_protocols': len(self.indexer.records),
                'min_similarity': min_similarity
            }
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
        
        Args:
            filename: Protocol filename (e.g., 'Suzuki_Cu_C(sp3)-C(sp2).json')
        
        Returns:
            Full protocol dictionary or None if not found
        """
        protocol_path = self.protocol_dir / filename
        
        if not protocol_path.exists():
            logger.error(f"Protocol file not found: {filename}")
            return None
        
        try:
            with open(protocol_path, 'r', encoding='utf-8') as f:
                return json.load(f)
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
        **kwargs
    ) -> Dict[str, Any]:
        """
        Recommend protocols and include extracted condition details
        
        Args:
            reaction_smiles: Query reaction SMILES
            k: Number of recommendations
            **kwargs: Additional arguments for recommend()
        
        Returns:
            Same as recommend() but with 'conditions' added to each match
        """
        results = self.recommend(reaction_smiles, k=k, **kwargs)
        
        # Add detailed conditions to each match
        for match in results['matches']:
            protocol = self.get_protocol_details(match['filename'])
            if protocol:
                match['conditions'] = self.extract_conditions(protocol)
        
        return results


def recommend_protocol(
    reaction_smiles: str,
    k: int = 5,
    reaction_family: Optional[str] = None,
    tags: Optional[List[str]] = None,
    index_path: Optional[Path] = None
) -> Dict[str, Any]:
    """
    Standalone function to recommend protocols
    
    Args:
        reaction_smiles: Query reaction SMILES
        k: Number of recommendations
        reaction_family: Optional family filter
        tags: Optional tag filter
        index_path: Optional index file path
    
    Returns:
        Recommendation results dictionary
    
    Example:
        from chemtools.protocol.recommend import recommend_protocol
        
        results = recommend_protocol(
            'CCBr.c1ccccc1B(O)O>>CCc1ccccc1',
            k=3,
            tags=['Suzuki']
        )
        
        for match in results['matches']:
            print(f"{match['similarity']:.3f}: {match['source_title']}")
    """
    recommender = ProtocolRecommender(index_path=index_path)
    return recommender.recommend(
        reaction_smiles=reaction_smiles,
        k=k,
        reaction_family=reaction_family,
        tags=tags
    )

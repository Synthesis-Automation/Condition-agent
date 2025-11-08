"""
Unified recommendation engine using DRFP-based similarity matching.

This module provides a unified recommendation system that combines both
protocols and rules from a single indexed source. It uses DRFP (Differentiable
Reaction Fingerprints) for similarity-based matching.

Key Features:
- Single unified index (protocols + rules)
- DRFP-based similarity search
- Combines specific protocols with general rules
- Ranked recommendations with confidence scores
- Source attribution (protocol vs rule)

Architecture:
- UnifiedRecommender: Main class for similarity-based recommendations
- Loads unified index from build output
- Uses DRFP for reaction similarity matching
- Returns ranked results with metadata

Usage:
    >>> from chemtools.recommend.unified import UnifiedRecommender
    >>> recommender = UnifiedRecommender("build/unified_index_complete")
    >>> results = recommender.recommend("CCBr.c1ccccc1B(O)O>>CCc1ccccc1")
    >>> for result in results:
    ...     print(f"{result['name']}: {result['similarity']:.3f}")
"""

from __future__ import annotations
from typing import List, Dict, Any, Optional, Tuple
from pathlib import Path
import json
import numpy as np
from dataclasses import dataclass


# Check DRFP availability
try:
    from drfp import DrfpEncoder
    DRFP_AVAILABLE = True
except ImportError:
    DRFP_AVAILABLE = False
    DrfpEncoder = None


@dataclass
class RecommendationResult:
    """A single recommendation result with metadata and similarity score."""
    
    # Identification
    id: str
    name: str
    source_type: str  # "protocol" or "rule"
    family: str
    
    # Similarity
    similarity: float
    rank: int
    
    # Metadata
    tags: List[str]
    version: str
    source_file: str
    
    # Optional: Full source data (loaded on demand)
    full_data: Optional[Dict[str, Any]] = None
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        return {
            "id": self.id,
            "name": self.name,
            "source_type": self.source_type,
            "family": self.family,
            "similarity": self.similarity,
            "rank": self.rank,
            "tags": self.tags,
            "version": self.version,
            "source_file": self.source_file,
        }


class UnifiedRecommender:
    """
    Unified recommendation engine using DRFP similarity matching.
    
    This class provides a single interface for recommending reaction conditions
    from both specific protocols and general rules, ranked by DRFP similarity.
    
    Attributes:
        index_dir: Path to unified index directory
        index: Loaded index metadata
        fingerprints: NumPy array of DRFP fingerprints
        encoder: DRFP encoder for computing query fingerprints
    
    Example:
        >>> recommender = UnifiedRecommender("build/unified_index_complete")
        >>> results = recommender.recommend(
        ...     "CCBr.c1ccccc1B(O)O>>CCc1ccccc1",
        ...     top_k=5
        ... )
        >>> print(f"Top match: {results[0].name}")
        >>> print(f"Similarity: {results[0].similarity:.3f}")
    """
    
    def __init__(self, index_dir: str | Path | None = None):
        """
        Initialize unified recommender with pre-built index.
        
        Args:
            index_dir: Path to directory containing:
                - index.json: Metadata for all sources
                - fingerprints.npz: DRFP fingerprints
                - stats.json: Build statistics (optional)
                If None, uses default "build/unified_index_complete"
        
        Raises:
            FileNotFoundError: If index files are missing
            ImportError: If DRFP library is not available
        """
        if not DRFP_AVAILABLE:
            raise ImportError(
                "DRFP library is required for UnifiedRecommender. "
                "Install with: pip install drfp"
            )
        
        # Use default index if not specified
        if index_dir is None:
            index_dir = Path(__file__).parent.parent.parent / "build" / "unified_index_complete"
        
        self.index_dir = Path(index_dir)
        
        # Load index metadata
        index_path = self.index_dir / "index.json"
        if not index_path.exists():
            raise FileNotFoundError(f"Index file not found: {index_path}")
        
        with open(index_path, 'r', encoding='utf-8') as f:
            self.index = json.load(f)
        
        # Load DRFP fingerprints
        fp_path = self.index_dir / "fingerprints.npz"
        if not fp_path.exists():
            raise FileNotFoundError(f"Fingerprints file not found: {fp_path}")
        
        fp_data = np.load(fp_path)
        
        # Concatenate protocol and rule fingerprints
        self.protocol_fps = fp_data['protocol_fps']
        self.rule_fps = fp_data['rule_fps']
        self.fingerprints = np.vstack([self.protocol_fps, self.rule_fps])
        
        # Initialize DRFP encoder
        self.encoder = DrfpEncoder()
        
        # Build lookup for fast access
        self._build_source_lookup()
    
    def _build_source_lookup(self):
        """Build internal lookup structures for fast source access."""
        # Load fingerprint indices from NPZ file
        fp_path = self.index_dir / "fingerprints.npz"
        fp_data = np.load(fp_path)
        
        # Build rule ID to fingerprint index range mapping
        rule_fp_mapping = {}
        if 'rule_fp_indices' in fp_data:
            for entry in fp_data['rule_fp_indices']:
                rule_id = str(entry[0])
                start_idx = int(entry[1])
                end_idx = int(entry[2])
                rule_fp_mapping[rule_id] = (start_idx, end_idx)
        
        # Combine protocols and rules into single list with indices
        self.sources = []
        
        # Protocol fingerprints: 1-to-1 mapping
        for idx, protocol in enumerate(self.index.get("protocols", [])):
            protocol['_index'] = idx
            protocol['_fp_index'] = idx  # Direct index into protocol_fps
            protocol['_fp_offset'] = 0  # Offset in concatenated array
            self.sources.append(protocol)
        
        # Rule fingerprints: multiple per rule
        num_protocol_fps = len(self.protocol_fps)
        for idx, rule in enumerate(self.index.get("rules", [])):
            rule['_index'] = len(self.sources)
            rule_id = rule.get('id', 'unknown')
            
            # Get fingerprint range for this rule
            if rule_id in rule_fp_mapping:
                start_idx, end_idx = rule_fp_mapping[rule_id]
                # Indices in concatenated array (protocol fps come first)
                rule['_fp_indices'] = list(range(
                    num_protocol_fps + start_idx,
                    num_protocol_fps + end_idx
                ))
            else:
                rule['_fp_indices'] = []
            
            self.sources.append(rule)
    
    def recommend(
        self,
        reaction_smiles: str,
        top_k: int = 5,
        min_similarity: float = 0.0,
        source_types: Optional[List[str]] = None,
    ) -> List[RecommendationResult]:
        """
        Recommend conditions for a reaction using DRFP similarity.
        
        Args:
            reaction_smiles: Reaction SMILES string (reactants>>products)
            top_k: Number of recommendations to return
            min_similarity: Minimum similarity threshold (0.0-1.0)
            source_types: Filter by source type(s): ["protocol"], ["rule"], or None (both)
        
        Returns:
            List of RecommendationResult objects, sorted by similarity (descending)
        
        Example:
            >>> results = recommender.recommend(
            ...     "CCBr.c1ccccc1B(O)O>>CCc1ccccc1",
            ...     top_k=3,
            ...     min_similarity=0.5
            ... )
            >>> for r in results:
            ...     print(f"{r.name}: {r.similarity:.3f}")
        """
        # Compute DRFP for query reaction
        try:
            query_fp = self.encoder.encode([reaction_smiles])[0]
        except Exception as e:
            # Return empty list for invalid SMILES instead of raising
            return []
        
        # Compute similarities with all fingerprints in index
        similarities = self._compute_similarities(query_fp)
        
        # Map fingerprint similarities back to sources
        source_similarities = self._aggregate_source_similarities(similarities)
        
        # Filter by source type if specified
        if source_types:
            source_similarities = [
                (source, sim) for source, sim in source_similarities
                if source['source_type'] in source_types
            ]
        
        # Filter by minimum similarity
        source_similarities = [
            (source, sim) for source, sim in source_similarities
            if sim >= min_similarity
        ]
        
        # Sort by similarity (descending) and take top_k
        source_similarities.sort(key=lambda x: x[1], reverse=True)
        top_sources = source_similarities[:top_k]
        
        # Build result objects
        results = []
        for rank, (source, similarity) in enumerate(top_sources, start=1):
            result = RecommendationResult(
                id=source['id'],
                name=source['name'],
                source_type=source['source_type'],
                family=source['family'],
                similarity=float(similarity),
                rank=rank,
                tags=source.get('tags', []),
                version=source.get('version', ''),
                source_file=source.get('source_file', ''),
            )
            results.append(result)
        
        return results
    
    def _compute_similarities(self, query_fp: np.ndarray) -> np.ndarray:
        """
        Compute Tanimoto similarities between query and all indexed fingerprints.
        
        Args:
            query_fp: Query DRFP fingerprint (2048-dimensional binary vector)
        
        Returns:
            Array of similarity scores (one per indexed fingerprint)
        """
        # Tanimoto similarity for binary fingerprints
        # sim = (A & B).sum() / (A | B).sum()
        # Vectorized: sim = intersection / (count_a + count_b - intersection)
        
        intersection = np.logical_and(query_fp, self.fingerprints).sum(axis=1)
        count_query = query_fp.sum()
        count_index = self.fingerprints.sum(axis=1)
        union = count_query + count_index - intersection
        
        # Avoid division by zero
        union = np.maximum(union, 1)
        
        similarities = intersection / union
        return similarities
    
    def _aggregate_source_similarities(
        self, 
        fp_similarities: np.ndarray
    ) -> List[Tuple[Dict[str, Any], float]]:
        """
        Aggregate fingerprint similarities back to source level.
        
        For protocols: Direct 1-to-1 mapping
        For rules: Take max similarity across reference reactions
        
        Args:
            fp_similarities: Array of similarities (one per fingerprint)
        
        Returns:
            List of (source_dict, similarity) tuples
        """
        source_similarities = []
        
        for source in self.sources:
            if source['source_type'] == 'protocol':
                # Protocol has single fingerprint
                fp_idx = source['_fp_index']
                similarity = fp_similarities[fp_idx]
                source_similarities.append((source, similarity))
            
            elif source['source_type'] == 'rule':
                # Rule has multiple reference fingerprints - take max
                fp_indices = source.get('_fp_indices', [])
                if fp_indices:
                    max_similarity = fp_similarities[fp_indices].max()
                    source_similarities.append((source, max_similarity))
        
        return source_similarities
    
    def get_source_details(self, source_id: str) -> Optional[Dict[str, Any]]:
        """
        Load full source data from file.
        
        Args:
            source_id: ID of the source to load
        
        Returns:
            Full source data dictionary, or None if not found
        """
        # Find source in index
        source = None
        for s in self.sources:
            if s['id'] == source_id:
                source = s
                break
        
        if not source:
            return None
        
        # Load from file
        source_file = Path(source['source_file'])
        if not source_file.exists():
            return None
        
        try:
            with open(source_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            # Handle array format (protocols can be arrays)
            if isinstance(data, list):
                # Find matching entry by family or take first
                for entry in data:
                    if entry.get('metadata', {}).get('id') == source_id:
                        return entry
                return data[0] if data else None
            
            return data
        
        except Exception as e:
            print(f"Error loading source {source_id}: {e}")
            return None
    
    def get_statistics(self) -> Dict[str, Any]:
        """
        Get statistics about the loaded index.
        
        Returns:
            Dictionary with:
                - num_protocols: Number of protocols
                - num_rules: Number of rules
                - num_fingerprints: Total fingerprints
                - num_families: Unique reaction families
                - version: Index version
        """
        stats_path = self.index_dir / "stats.json"
        
        if stats_path.exists():
            with open(stats_path, 'r', encoding='utf-8') as f:
                stats = json.load(f)
            return stats
        
        # Fallback: compute basic stats from index
        protocol_families = set(p['family'] for p in self.index.get('protocols', []))
        rule_families = set(r['family'] for r in self.index.get('rules', []))
        
        return {
            "num_protocols": len(self.index.get('protocols', [])),
            "num_rules": len(self.index.get('rules', [])),
            "num_fingerprints": len(self.fingerprints),
            "num_families": len(protocol_families | rule_families),
            "version": self.index.get('version', 'unknown'),
        }

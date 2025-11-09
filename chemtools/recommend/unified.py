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

# Import rule-to-protocol converter
from chemtools.formatters.rule_to_protocol import rule_conditions_to_reaction_setup

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
        validate_rules: bool = True,
        format_for_automation: bool = False,
        scale_mmol: float = 1.0,
    ) -> List[RecommendationResult]:
        """
        Recommend conditions for a reaction using DRFP similarity.
        
        Args:
            reaction_smiles: Reaction SMILES string (reactants>>products)
            top_k: Number of recommendations to return
            min_similarity: Minimum similarity threshold (0.0-1.0)
            source_types: Filter by source type(s): ["protocol"], ["rule"], or None (both)
            validate_rules: If True, validates sources post-similarity:
                - Rules: checks applies_if criteria (functional groups)
                - Protocols: checks reaction_SMARTS patterns (exact transformations)
                Default: True (recommended for chemical accuracy)
            format_for_automation: If True, converts rule conditions to protocol-compatible
                reaction_setup format with ordered addition sequences. Protocols already
                have this structure. Default: False
            scale_mmol: Reaction scale in mmol for automated format. Default: 1.0
        
        Returns:
            List of RecommendationResult objects, sorted by similarity (descending)
            
        Note:
            Post-similarity validation provides chemical appropriateness filtering:
            - Rules use applies_if (functional group detection) to verify substrates
            - Protocols use reaction_SMARTS (transformation matching) for exact validation
            Both mechanisms fail-open if validation cannot be performed.
            
            When format_for_automation=True, both rules and protocols will have a
            "reaction_setup" field with ordered chemicals and conditions for automation.
        
        Example:
            >>> results = recommender.recommend(
            ...     "CCBr.c1ccccc1B(O)O>>CCc1ccccc1",
            ...     top_k=3,
            ...     min_similarity=0.5,
            ...     format_for_automation=True,
            ...     scale_mmol=1.0
            ... )
            >>> for r in results:
            ...     print(f"{r.name}: {r.similarity:.3f}")
            ...     if hasattr(r, 'full_data') and 'reaction_setup' in r.full_data:
            ...         setup = r.full_data['reaction_setup'][0]
            ...         print(f"  Chemicals: {[c['name'] for c in setup['chemicals']]}")
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
        
        # Validate rules' applies_if criteria if requested
        if validate_rules:
            source_similarities = self._validate_rule_applicability(
                source_similarities, 
                reaction_smiles
            )
        
        # Validate protocols' reaction_SMARTS patterns if requested
        if validate_rules:
            source_similarities = self._validate_protocol_smarts(
                source_similarities,
                reaction_smiles
            )
        
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
            
            # Load and format full data if automation format requested
            if format_for_automation:
                full_data = self.get_source_data(source['id'])
                if full_data:
                    # For rules: convert to protocol-compatible format
                    if source['source_type'] == 'rule':
                        conditions = full_data.get('default_rule', {}).get('conditions', {})
                        if conditions:
                            formatted = rule_conditions_to_reaction_setup(
                                conditions=conditions,
                                scale_mmol=scale_mmol,
                                reaction_family=source['family']
                            )
                            result.full_data = formatted
                    
                    # For protocols: use existing reaction_setup
                    elif source['source_type'] == 'protocol':
                        # Protocols already have reaction_setup structure
                        result.full_data = {
                            'reaction_setup': full_data.get('reaction_setup', []),
                            'metadata': {
                                'generated_from': 'protocol',
                                'format': 'protocol-native',
                                'scale_mmol': scale_mmol
                            }
                        }
            
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
    
    def _validate_rule_applicability(
        self,
        source_similarities: List[Tuple[Dict[str, Any], float]],
        reaction_smiles: str
    ) -> List[Tuple[Dict[str, Any], float]]:
        """
        Validate that rules meet their applies_if criteria based on detected features.
        
        This acts as a post-similarity filter to ensure recommended rules are
        chemically appropriate for the query reaction.
        
        Args:
            source_similarities: List of (source, similarity) tuples
            reaction_smiles: Query reaction SMILES
        
        Returns:
            Filtered list containing only rules that pass applies_if validation
            (protocols are always included)
        """
        try:
            from ..rule.analyzer import FeatureAnalyzer
            analyzer = FeatureAnalyzer()
            
            # Detect features from the query reaction
            features = analyzer.analyze_reaction(reaction_smiles, combine_method="union")
            
            if not features:
                # If feature detection fails, return all sources (fail open)
                return source_similarities
            
            validated = []
            
            for source, similarity in source_similarities:
                # Protocols don't have applies_if - always include
                if source['source_type'] == 'protocol':
                    validated.append((source, similarity))
                    continue
                
                # For rules, check applies_if criteria
                if source['source_type'] == 'rule':
                    # Load full rule details to get applies_if
                    rule_details = self.get_source_details(source['id'])
                    
                    if not rule_details or 'applies_if' not in rule_details:
                        # If no applies_if, include the rule (permissive)
                        validated.append((source, similarity))
                        continue
                    
                    applies_if = rule_details['applies_if']
                    
                    # Check if applies_if criteria are met
                    if self._check_applies_if(features, applies_if):
                        validated.append((source, similarity))
                    # If not met, silently exclude this rule
            
            return validated
            
        except Exception:
            # If validation fails for any reason, fail open (return all sources)
            return source_similarities
    
    def _check_applies_if(self, features: Dict[str, bool], applies_if: Dict[str, Any]) -> bool:
        """
        Check if detected features satisfy applies_if criteria.
        
        Args:
            features: Dictionary of detected features {feature_name: bool}
            applies_if: applies_if dictionary with 'all' and/or 'any' keys
        
        Returns:
            True if criteria are met, False otherwise
        """
        # Check 'all' conditions (all must be true)
        if 'all' in applies_if:
            all_conditions = applies_if['all']
            if not all(features.get(condition, False) for condition in all_conditions):
                return False
        
        # Check 'any' conditions (at least one must be true)
        if 'any' in applies_if:
            any_conditions = applies_if['any']
            if not any(features.get(condition, False) for condition in any_conditions):
                return False
        
        # If we get here, all criteria are satisfied
        return True
    
    def _validate_protocol_smarts(
        self,
        source_similarities: List[Tuple[Dict[str, Any], float]],
        reaction_smiles: str
    ) -> List[Tuple[Dict[str, Any], float]]:
        """
        Validate protocols against their reaction_SMARTS patterns.
        
        For protocols with reaction_SMARTS field, checks if the pattern matches
        the query reaction. This provides more accurate filtering than DRFP
        similarity alone, as it verifies the exact transformation.
        
        Args:
            source_similarities: List of (source, similarity) tuples
            reaction_smiles: Query reaction SMILES
        
        Returns:
            Filtered list excluding protocols whose reaction_SMARTS don't match
            
        Note:
            - Fails open: if RDKit unavailable or matching fails, includes the protocol
            - Rules are always included (no reaction_SMARTS validation)
            - Protocols without reaction_SMARTS are always included
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
        except ImportError:
            # RDKit not available - fail open
            return source_similarities
        
        try:
            # Parse query reaction
            if '>>' not in reaction_smiles:
                # Invalid reaction SMILES - fail open
                return source_similarities
            
            parts = reaction_smiles.split('>>')
            reactants_smiles = parts[0].split('>')[-1]  # Handle agents if present
            products_smiles = parts[1].split('>')[0]
            
            # Parse reactants
            query_reactants = []
            for r_smi in reactants_smiles.split('.'):
                r_smi = r_smi.strip()
                if r_smi:
                    mol = Chem.MolFromSmiles(r_smi)
                    if mol:
                        query_reactants.append(mol)
            
            if not query_reactants:
                # Can't parse reactants - fail open
                return source_similarities
            
            validated = []
            
            for source, similarity in source_similarities:
                # Rules don't have reaction_SMARTS - always include
                if source['source_type'] == 'rule':
                    validated.append((source, similarity))
                    continue
                
                # For protocols, check reaction_SMARTS if present
                if source['source_type'] == 'protocol':
                    # Load full protocol details to get reaction_SMARTS
                    protocol_details = self.get_source_details(source['id'])
                    
                    if not protocol_details:
                        # Can't load details - fail open
                        validated.append((source, similarity))
                        continue
                    
                    reaction_data = protocol_details.get('reaction', {})
                    reaction_smarts_list = reaction_data.get('reaction_SMARTS', [])
                    
                    if not reaction_smarts_list:
                        # No reaction_SMARTS - include protocol (permissive)
                        validated.append((source, similarity))
                        continue
                    
                    # Try to match any of the reaction_SMARTS patterns
                    matches = False
                    for smarts in reaction_smarts_list:
                        if self._check_reaction_smarts_match(smarts, query_reactants):
                            matches = True
                            break
                    
                    if matches:
                        validated.append((source, similarity))
                    # If no patterns match, silently exclude this protocol
            
            return validated
            
        except Exception:
            # If validation fails for any reason, fail open (return all sources)
            return source_similarities
    
    def _check_reaction_smarts_match(
        self,
        smarts: str,
        query_reactants: List[Any]
    ) -> bool:
        """
        Check if a reaction SMARTS pattern matches query reactants.
        
        Args:
            smarts: Reaction SMARTS string (e.g., "IC=C.CC(O)(C#N)C>>N#CC=C")
            query_reactants: List of RDKit Mol objects for query reactants
        
        Returns:
            True if pattern matches, False otherwise
        """
        try:
            from rdkit.Chem import AllChem
            
            # Parse reaction SMARTS
            rxn = AllChem.ReactionFromSmarts(smarts)
            if rxn is None:
                # Invalid SMARTS - fail open
                return True
            
            # Try to match reactants with pattern
            # RunReactants returns products if reactants match the pattern
            try:
                products = rxn.RunReactants(tuple(query_reactants))
                # If we get any products, the pattern matched
                return len(products) > 0
            except Exception:
                # Matching failed - try reverse order of reactants
                if len(query_reactants) >= 2:
                    try:
                        products = rxn.RunReactants(tuple(reversed(query_reactants)))
                        return len(products) > 0
                    except Exception:
                        pass
                
                # Could not match - fail open (include protocol)
                return True
            
        except Exception:
            # Any error - fail open
            return True
    
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

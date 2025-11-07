"""
Feature Analyzer for Reaction SMILES
=====================================

Extracts calculable features from reaction components for rule matching.
"""

from __future__ import annotations
from typing import Dict, Any, List, Tuple, Optional
import logging

from ..featurizers.calculable import detect_all_features
from ..featurizers.molecular import featurize as molecular_featurize
from ..analysis.smiles import normalize_reaction

logger = logging.getLogger(__name__)


class FeatureAnalyzer:
    """Analyzes reaction SMILES to extract calculable features."""
    
    def __init__(self):
        """Initialize the feature analyzer."""
        pass
    
    def analyze_reaction(
        self,
        reaction_smiles: str,
        combine_method: str = "union"
    ) -> Dict[str, Any]:
        """
        Analyze a reaction SMILES and extract features from all reactants.
        
        Args:
            reaction_smiles: Reaction SMILES string (reactants>>products)
            combine_method: How to combine features from multiple reactants:
                - "union": Feature is True if present in ANY reactant
                - "all": Feature is True only if present in ALL reactants
                - "first": Use only features from first reactant
                - "separate": Return separate feature dicts for each reactant
        
        Returns:
            Dictionary of combined features
        
        Example:
            >>> analyzer = FeatureAnalyzer()
            >>> features = analyzer.analyze_reaction("Brc1ccccc1.OB(O)c1ccccc1>>...")
            >>> features['sp2_halide_present']
            True
        """
        try:
            # Parse reaction SMILES
            reactants, products = self._parse_reaction(reaction_smiles)
            
            if not reactants:
                logger.warning(f"No reactants found in reaction: {reaction_smiles}")
                return {}
            
            # Special case: For 2-component reactions (electrophile + nucleophile),
            # use the full molecular featurizer which includes functional group enrichment
            if len(reactants) == 2:
                try:
                    features = molecular_featurize(reactants[0], reactants[1])
                    logger.debug(f"Used molecular featurizer: {sum(1 for v in features.values() if v)} features detected")
                    return features
                except Exception as e:
                    logger.warning(f"Molecular featurizer failed, falling back to calculable: {e}")
                    # Fall through to calculable featurizer below
            
            # Detect features for each reactant using calculable featurizer
            reactant_features = []
            for i, reactant in enumerate(reactants):
                try:
                    features = detect_all_features(reactant)
                    reactant_features.append(features)
                    logger.debug(f"Reactant {i+1} ({reactant}): {sum(1 for v in features.values() if v)} features detected")
                except Exception as e:
                    logger.error(f"Error detecting features for reactant {reactant}: {e}")
                    reactant_features.append({})
            
            # Combine features based on method
            if combine_method == "separate":
                return {"reactants": reactant_features}
            elif combine_method == "first":
                return reactant_features[0] if reactant_features else {}
            elif combine_method == "all":
                return self._combine_features_all(reactant_features)
            else:  # union (default)
                return self._combine_features_union(reactant_features)
        
        except Exception as e:
            logger.error(f"Error analyzing reaction {reaction_smiles}: {e}")
            return {}
    
    def analyze_reactant(self, smiles: str) -> Dict[str, Any]:
        """
        Analyze a single reactant SMILES.
        
        Args:
            smiles: SMILES string for a single molecule
        
        Returns:
            Dictionary of detected features
        """
        try:
            return detect_all_features(smiles)
        except Exception as e:
            logger.error(f"Error analyzing reactant {smiles}: {e}")
            return {}
    
    def _parse_reaction(self, reaction_smiles: str) -> Tuple[List[str], List[str]]:
        """
        Parse reaction SMILES into reactants and products.
        
        Args:
            reaction_smiles: Full reaction SMILES
        
        Returns:
            (reactants_list, products_list)
        """
        # Handle simple >> splitting
        if ">>" in reaction_smiles:
            parts = reaction_smiles.split(">>")
            reactant_part = parts[0]
            product_part = parts[1] if len(parts) > 1 else ""
        elif ">" in reaction_smiles:
            # Handle single > (reactants>reagents>products)
            parts = reaction_smiles.split(">")
            reactant_part = parts[0]
            product_part = parts[-1] if len(parts) > 2 else ""
        else:
            # No reaction arrow, treat as single molecule
            return [reaction_smiles], []
        
        # Split reactants by period
        reactants = [r.strip() for r in reactant_part.split(".") if r.strip()]
        products = [p.strip() for p in product_part.split(".") if p.strip()]
        
        return reactants, products
    
    def _combine_features_union(self, feature_dicts: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Combine features using union logic (OR).
        
        For boolean features: True if True in ANY reactant
        For integer features: Maximum value across all reactants
        """
        if not feature_dicts:
            return {}
        
        # Start with first dict
        combined = feature_dicts[0].copy()
        
        # Merge remaining dicts
        for feat_dict in feature_dicts[1:]:
            for key, value in feat_dict.items():
                if key not in combined:
                    combined[key] = value
                elif isinstance(value, bool):
                    # OR logic for booleans
                    combined[key] = combined[key] or value
                elif isinstance(value, int):
                    # Max for integers
                    combined[key] = max(combined[key], value)
                elif isinstance(value, (float, str)):
                    # Keep first value for others
                    pass
        
        return combined
    
    def _combine_features_all(self, feature_dicts: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Combine features using ALL logic (AND).
        
        For boolean features: True only if True in ALL reactants
        For integer features: Minimum value across all reactants
        """
        if not feature_dicts:
            return {}
        
        # Start with first dict
        combined = feature_dicts[0].copy()
        
        # Merge remaining dicts
        for feat_dict in feature_dicts[1:]:
            for key in list(combined.keys()):
                if key not in feat_dict:
                    # Feature not in this reactant, set to False/0
                    if isinstance(combined[key], bool):
                        combined[key] = False
                    elif isinstance(combined[key], int):
                        combined[key] = 0
                elif isinstance(combined[key], bool):
                    # AND logic for booleans
                    combined[key] = combined[key] and feat_dict[key]
                elif isinstance(combined[key], int):
                    # Min for integers
                    combined[key] = min(combined[key], feat_dict[key])
        
        return combined
    
    def get_relevant_features(
        self,
        features: Dict[str, Any],
        include_false: bool = False
    ) -> Dict[str, Any]:
        """
        Filter features to show only relevant ones.
        
        Args:
            features: Full feature dictionary
            include_false: If True, include False boolean values
        
        Returns:
            Filtered feature dictionary
        """
        relevant = {}
        
        for key, value in features.items():
            if isinstance(value, bool):
                if value or include_false:
                    relevant[key] = value
            elif isinstance(value, int):
                if value > 0:
                    relevant[key] = value
            else:
                # Include other types (float, str)
                relevant[key] = value
        
        return relevant

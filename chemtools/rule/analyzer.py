"""
Feature Analyzer for Reaction SMILES
=====================================

Extracts calculable features from reaction components for rule matching.
"""

from __future__ import annotations
from typing import Dict, Any, List, Tuple, Optional
import logging

from ..featurizers.unified import (
    featurize_reaction as featurize_unified_reaction,
    featurize_molecule as featurize_unified_molecule,
)

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
        """
        try:
            # Use the unified featurizer bundle
            bundle = featurize_unified_reaction(reaction_smiles)
            reaction_payload = bundle.get("reaction") or bundle
            
            reactant_features = []
            all_motifs = set()
            
            for r_entry in reaction_payload.get("reactants", []):
                r_analysis = self._extract_reactant_analysis(r_entry)
                r_feat = self._flatten_analysis(r_analysis)
                reactant_features.append(r_feat)
                
                for motif in r_analysis.get("motifs", []):
                    all_motifs.add(motif["compound_id"])
            
            # Combine features based on method
            if combine_method == "separate":
                return {"reactants": reactant_features, "motifs": sorted(list(all_motifs))}
            elif combine_method == "first":
                combined = reactant_features[0] if reactant_features else {}
            elif combine_method == "all":
                combined = self._combine_features_all(reactant_features)
            else:  # union (default)
                combined = self._combine_features_union(reactant_features)

            # Add reaction detection results
            detection = reaction_payload.get("detection", {})
            for match in detection.get("matches", []):
                combined[f"rxn_{match['reaction_type']}"] = True
                combined[match['reaction_type']] = True
            
            # Add all motifs to the flat dict
            for motif_id in all_motifs:
                combined[motif_id] = True
                
            return combined
        
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
            bundle = featurize_unified_molecule(smiles)
            analysis = bundle.get("molecule") or bundle
            features = self._flatten_analysis(analysis)
            return features
        except Exception as e:
            logger.error(f"Error analyzing reactant {smiles}: {e}")
            return {}

    def _extract_reactant_analysis(self, entry: Dict[str, Any]) -> Dict[str, Any]:
        """Handle both legacy and unified reactant bundle shapes."""
        if "analysis" in entry and isinstance(entry.get("analysis"), dict):
            return entry["analysis"]
        if "molecule" in entry and isinstance(entry.get("molecule"), dict):
            return entry["molecule"]
        return entry

    def _flatten_analysis(self, analysis: Dict[str, Any]) -> Dict[str, Any]:
        """Flatten V2 analysis into a boolean feature dictionary."""
        features = {}
        
        # 1. Motifs
        for motif in analysis.get("motifs", []):
            motif_id = motif["compound_id"]
            features[motif_id] = True
            
            # Add category as a feature
            if "category" in motif:
                cat_feature = motif["category"].lower().replace(" ", "_").replace("-", "_") + "_present"
                features[cat_feature] = True

        # 2. Sterics (Aryl)
        for entry in analysis.get("steric", {}).get("aryl", []):
            res = entry.get("result", {})
            score = res.get("score_0_10", 0)
            ortho = res.get("ortho", [])
            
            num_ortho_subs = sum(1 for o in ortho if o.get("bulk", 0) > 0)
            if num_ortho_subs >= 1:
                features["ortho_1"] = True
            if num_ortho_subs >= 2:
                features["ortho_2plus"] = True
            
            if score >= 4.0:
                features["ortho_hindered"] = True
            if score >= 7.0:
                features["ortho_very_hindered"] = True
            
            # Add raw score for numeric rules
            features["aryl_steric_score"] = max(features.get("aryl_steric_score", 0), score)

        # 3. Electronics (Aryl)
        for entry in analysis.get("electronics", {}).get("aryl", []):
            res = entry.get("result", {})
            score = res.get("score_0_10", 5.0)
            
            if score > 6.5:
                features["electron_poor_aryl"] = True
            if score > 8.5:
                features["very_electron_poor_aryl"] = True
            if score < 3.5:
                features["electron_rich_aryl"] = True
            
            # Add raw score
            features["aryl_electronic_score"] = max(features.get("aryl_electronic_score", 0), score)

        # 4. Sterics (Alkyl)
        for entry in analysis.get("steric", {}).get("alkyl", []):
            res = entry.get("result", {})
            score = res.get("score_0_10", 0)
            if score >= 5.0:
                features["alkyl_hindered"] = True
            features["alkyl_steric_score"] = max(features.get("alkyl_steric_score", 0), score)

        return features

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

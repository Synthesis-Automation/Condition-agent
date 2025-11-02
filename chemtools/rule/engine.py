"""
Rule-Based Recommendation Engine
=================================

Core orchestration for feature-driven condition recommendations.
"""

from __future__ import annotations
from typing import Dict, Any, List, Optional
from pathlib import Path
import logging

from .database import RuleDatabase
from .analyzer import FeatureAnalyzer
from .models import ConditionRecommendation, AppliedRule, AppliedModifier

logger = logging.getLogger(__name__)


class RuleEngine:
    """
    Orchestrates rule-based condition recommendations.
    
    Workflow:
        1. Parse reaction SMILES
        2. Detect features from reactants
        3. Check if database applies
        4. Find matching base rule
        5. Apply relevant modifiers
        6. Generate recommendation
    """
    
    def __init__(
        self,
        database: Optional[RuleDatabase] = None,
        analyzer: Optional[FeatureAnalyzer] = None
    ):
        """
        Initialize the recommendation engine.
        
        Args:
            database: RuleDatabase instance (can be set later)
            analyzer: FeatureAnalyzer instance (creates default if None)
        """
        self.database = database
        self.analyzer = analyzer or FeatureAnalyzer()
    
    @classmethod
    def from_file(cls, database_path: str | Path) -> RuleEngine:
        """
        Create engine with database loaded from file.
        
        Args:
            database_path: Path to rule database JSON
        
        Returns:
            RuleEngine instance
        
        Example:
            >>> engine = RuleEngine.from_file("data/suzuki.json")
            >>> rec = engine.recommend("Brc1ccccc1.OB(O)c1ccccc1>>...")
        """
        db = RuleDatabase.from_file(database_path)
        return cls(database=db)
    
    def load_database(self, path: str | Path) -> None:
        """
        Load a rule database from file.
        
        Args:
            path: Path to rule database JSON
        """
        self.database = RuleDatabase.from_file(path)
        logger.info(f"Loaded database: {self.database.metadata.get('name', 'Unnamed')}")
    
    def recommend(
        self,
        reaction_smiles: str,
        symptoms: Optional[List[str]] = None,
        combine_method: str = "union",
        include_reasoning: bool = True
    ) -> ConditionRecommendation:
        """
        Generate condition recommendation for a reaction.
        
        Args:
            reaction_smiles: Reaction SMILES string
            symptoms: Optional list of observed symptoms (e.g., ["low_yield"])
            combine_method: How to combine features from multiple reactants
            include_reasoning: Include matched features in recommendation
        
        Returns:
            ConditionRecommendation with conditions and modifiers
        
        Raises:
            ValueError: If no database is loaded
            RuntimeError: If database doesn't apply to this reaction
        
        Example:
            >>> engine = RuleEngine.from_file("data/suzuki.json")
            >>> rec = engine.recommend("Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1")
            >>> print(rec.format_summary())
        """
        if not self.database:
            raise ValueError("No database loaded. Use load_database() or from_file()")
        
        # Step 1: Analyze reaction features
        logger.info(f"Analyzing reaction: {reaction_smiles[:50]}...")
        features = self.analyzer.analyze_reaction(
            reaction_smiles,
            combine_method=combine_method
        )
        
        if not features:
            raise RuntimeError(f"Failed to extract features from reaction: {reaction_smiles}")
        
        logger.debug(f"Detected features: {sum(1 for v in features.values() if v)} active")
        
        # Step 2: Check if database applies
        if not self.database.check_applies(features):
            applies_if = self.database.applies_if
            raise RuntimeError(
                f"Database does not apply to this reaction. "
                f"Required: {applies_if}, "
                f"Got: {self.analyzer.get_relevant_features(features)}"
            )
        
        # Step 3: Find matching base rule
        matching_rule = self.database.find_matching_rule(features)
        
        if not matching_rule:
            raise RuntimeError(
                f"No matching rule found for features: "
                f"{self.analyzer.get_relevant_features(features)}"
            )
        
        # Determine which features matched (for reasoning)
        rule_matches, matched_features = matching_rule.matches(features)
        
        applied_rule = AppliedRule(
            name=matching_rule.name,
            conditions=matching_rule.conditions,
            matched_features=matched_features if include_reasoning else [],
            confidence=1.0  # Can be enhanced with scoring logic
        )
        
        logger.info(f"Matched rule: {matching_rule.name}")
        
        # Step 4: Find matching modifiers
        matching_modifiers = self.database.find_matching_modifiers(
            features,
            symptoms=symptoms
        )
        
        applied_modifiers = []
        for modifier in matching_modifiers:
            # Get matched conditions for reasoning
            matched_conditions = []
            if include_reasoning:
                # Check which specific conditions matched
                for condition in modifier.when:
                    if condition.startswith("symptom:"):
                        symptom_name = condition.replace("symptom:", "").strip()
                        if symptoms and symptom_name in symptoms:
                            matched_conditions.append(condition)
                    elif features.get(condition, False):
                        matched_conditions.append(condition)
            
            applied_modifiers.append(AppliedModifier(
                suggestion=modifier.suggestion,
                rationale=modifier.rationale,
                matched_conditions=matched_conditions if include_reasoning else []
            ))
        
        logger.info(f"Applied {len(applied_modifiers)} modifiers")
        
        # Step 5: Generate recommendation
        recommendation = ConditionRecommendation(
            reaction_smiles=reaction_smiles,
            base_rule=applied_rule,
            modifiers=applied_modifiers,
            detected_features=self.analyzer.get_relevant_features(features) if include_reasoning else {}
        )
        
        return recommendation
    
    def batch_recommend(
        self,
        reactions: List[str],
        symptoms_list: Optional[List[Optional[List[str]]]] = None,
        combine_method: str = "union",
        include_reasoning: bool = True
    ) -> List[ConditionRecommendation]:
        """
        Generate recommendations for multiple reactions.
        
        Args:
            reactions: List of reaction SMILES
            symptoms_list: Optional list of symptom lists (one per reaction)
            combine_method: How to combine features from multiple reactants
            include_reasoning: Include matched features in recommendations
        
        Returns:
            List of ConditionRecommendation objects
        """
        recommendations = []
        
        if symptoms_list is None:
            symptoms_list = [None] * len(reactions)
        
        for i, (reaction, symptoms) in enumerate(zip(reactions, symptoms_list)):
            try:
                rec = self.recommend(
                    reaction,
                    symptoms=symptoms,
                    combine_method=combine_method,
                    include_reasoning=include_reasoning
                )
                recommendations.append(rec)
            except Exception as e:
                logger.error(f"Error processing reaction {i+1}: {e}")
                # Create error recommendation
                recommendations.append(ConditionRecommendation(
                    reaction_smiles=reaction,
                    base_rule=AppliedRule(
                        name="error",
                        conditions={"error": str(e)},
                        matched_features=[],
                        confidence=0.0
                    ),
                    modifiers=[],
                    detected_features={}
                ))
        
        return recommendations
    
    def validate_database(self) -> List[str]:
        """
        Validate the loaded database.
        
        Returns:
            List of validation issues (empty if valid)
        """
        if not self.database:
            return ["No database loaded"]
        
        return self.database.validate()
    
    def get_database_summary(self) -> str:
        """
        Get summary of loaded database.
        
        Returns:
            Formatted summary string
        """
        if not self.database:
            return "No database loaded"
        
        return self.database.get_summary()

"""
Data Models for Rule-Based Recommendation System
=================================================

Structured data classes for representing recommendations, rules, and modifiers.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional


@dataclass
class AppliedRule:
    """Represents a rule that was applied during recommendation."""
    
    id: str
    description: str
    conditions: Dict[str, Any]
    matched_features: List[str] = field(default_factory=list)
    confidence: float = 1.0
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        return {
            "id": self.id,
            "description": self.description,
            "conditions": self.conditions,
            "matched_features": self.matched_features,
            "confidence": self.confidence,
        }


@dataclass
class AppliedModifier:
    """Represents a modifier that was applied during recommendation."""
    
    id: str
    when: List[str]
    suggest: str
    matched_conditions: List[str] = field(default_factory=list)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        return {
            "id": self.id,
            "when": self.when,
            "suggest": self.suggest,
            "matched_conditions": self.matched_conditions,
        }


@dataclass
class ConditionRecommendation:
    """Complete recommendation output with base conditions and modifiers."""
    
    reaction_smiles: str
    reaction_type: str
    base_conditions: Dict[str, Any]
    applied_rule: AppliedRule
    applied_modifiers: List[AppliedModifier] = field(default_factory=list)
    detected_features: Dict[str, Any] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        return {
            "reaction_smiles": self.reaction_smiles,
            "reaction_type": self.reaction_type,
            "base_conditions": self.base_conditions,
            "applied_rule": self.applied_rule.to_dict(),
            "applied_modifiers": [m.to_dict() for m in self.applied_modifiers],
            "detected_features": self.detected_features,
            "warnings": self.warnings,
            "metadata": self.metadata,
        }
    
    def format_summary(self) -> str:
        """Format a human-readable summary."""
        lines = [
            f"{'='*70}",
            f"REACTION CONDITIONS RECOMMENDATION",
            f"{'='*70}",
            f"Reaction Type: {self.reaction_type}",
            f"Reaction SMILES: {self.reaction_smiles}",
            f"",
            f"Applied Rule: {self.applied_rule.id}",
            f"Description: {self.applied_rule.description}",
            f"",
            f"BASE CONDITIONS:",
            f"{'-'*70}",
        ]
        
        for key, value in self.base_conditions.items():
            lines.append(f"  {key:25s}: {value}")
        
        if self.applied_modifiers:
            lines.append(f"")
            lines.append(f"MODIFIERS APPLIED:")
            lines.append(f"{'-'*70}")
            for mod in self.applied_modifiers:
                lines.append(f"  [{mod.id}]")
                lines.append(f"  → {mod.suggest}")
                lines.append(f"")
        
        if self.warnings:
            lines.append(f"")
            lines.append(f"WARNINGS:")
            lines.append(f"{'-'*70}")
            for warning in self.warnings:
                lines.append(f"  ⚠ {warning}")
        
        lines.append(f"")
        lines.append(f"KEY FEATURES DETECTED:")
        lines.append(f"{'-'*70}")
        
        # Show only True boolean features and non-zero integer features
        relevant_features = {
            k: v for k, v in self.detected_features.items()
            if (isinstance(v, bool) and v) or (isinstance(v, int) and v > 0)
        }
        
        for key, value in sorted(relevant_features.items()):
            if isinstance(value, bool):
                lines.append(f"  ✓ {key}")
            else:
                lines.append(f"  • {key} = {value}")
        
        lines.append(f"{'='*70}")
        
        return "\n".join(lines)


@dataclass
class RuleSpec:
    """Specification for a single rule in the database."""
    
    id: str
    description: str
    reactant_features: Optional[Dict[str, Any]] = None
    conditions: Dict[str, Any] = field(default_factory=dict)
    priority: int = 0
    
    def matches(self, features: Dict[str, Any]) -> tuple[bool, List[str]]:
        """
        Check if this rule matches the given features.
        
        Returns:
            (matches, matched_feature_list)
        """
        if not self.reactant_features:
            return True, []
        
        matched = []
        
        # Handle AND logic
        if "and" in self.reactant_features:
            required = self.reactant_features["and"]
            for feat in required:
                if not features.get(feat, False):
                    return False, []
                matched.append(feat)
            return True, matched
        
        # Handle OR logic (any)
        if "any" in self.reactant_features:
            options = self.reactant_features["any"]
            for feat in options:
                if features.get(feat, False):
                    matched.append(feat)
            return len(matched) > 0, matched
        
        # Handle ALL logic (explicit)
        if "all" in self.reactant_features:
            required = self.reactant_features["all"]
            for feat in required:
                if not features.get(feat, False):
                    return False, []
                matched.append(feat)
            return True, matched
        
        # Default: no feature requirements
        return True, []


@dataclass
class ModifierSpec:
    """Specification for a modifier in the database."""
    
    id: str
    when: List[str]
    suggest: str
    priority: int = 0
    
    def matches(self, features: Dict[str, Any], symptoms: List[str] = None) -> tuple[bool, List[str]]:
        """
        Check if this modifier should be applied.
        
        Args:
            features: Detected molecular features
            symptoms: Optional list of observed symptoms (e.g., "hydrodehalogenation")
        
        Returns:
            (should_apply, matched_conditions)
        """
        symptoms = symptoms or []
        matched = []
        
        for condition in self.when:
            # Check for symptom-based conditions
            if condition.startswith("symptom:"):
                symptom_text = condition.replace("symptom:", "").strip()
                if any(symptom_text.lower() in s.lower() for s in symptoms):
                    matched.append(condition)
                    continue
            
            # Check for feature-based conditions
            if features.get(condition, False):
                matched.append(condition)
        
        return len(matched) > 0, matched

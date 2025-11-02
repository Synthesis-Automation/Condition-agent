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
    
    name: str
    conditions: Dict[str, Any]
    matched_features: List[str] = field(default_factory=list)
    confidence: float = 1.0
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        return {
            "name": self.name,
            "conditions": self.conditions,
            "matched_features": self.matched_features,
            "confidence": self.confidence,
        }


@dataclass
class AppliedModifier:
    """Represents a modifier that was applied during recommendation."""
    
    suggestion: str
    rationale: Optional[str] = None
    matched_conditions: List[str] = field(default_factory=list)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        result = {
            "suggestion": self.suggestion,
            "matched_conditions": self.matched_conditions,
        }
        if self.rationale:
            result["rationale"] = self.rationale
        return result


@dataclass
class ConditionRecommendation:
    """Complete recommendation output with base conditions and modifiers."""
    
    reaction_smiles: str
    base_rule: AppliedRule
    modifiers: List[AppliedModifier] = field(default_factory=list)
    detected_features: Dict[str, Any] = field(default_factory=dict)
    applied_modifiers: List[AppliedModifier] = field(default_factory=list)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        return {
            "reaction_smiles": self.reaction_smiles,
            "base_rule": self.base_rule.to_dict(),
            "modifiers": [m.to_dict() for m in self.modifiers],
            "detected_features": self.detected_features,
        }
    
    def format_summary(self) -> str:
        """Format a human-readable summary."""
        lines = [
            f"{'='*70}",
            f"CONDITION RECOMMENDATION",
            f"{'='*70}",
            f"Reaction: {self.reaction_smiles}",
            f"",
            f"Applied Rule: {self.base_rule.name}",
            f"Confidence: {self.base_rule.confidence:.2f}",
            f"",
            f"RECOMMENDED CONDITIONS:",
            f"{'-'*70}",
        ]
        
        for key, value in self.base_rule.conditions.items():
            lines.append(f"  {key:25s}: {value}")
        
        if self.base_rule.matched_features:
            lines.append(f"")
            lines.append(f"Matched features: {', '.join(self.base_rule.matched_features)}")
        
        if self.modifiers:
            lines.append(f"")
            lines.append(f"MODIFIERS APPLIED:")
            lines.append(f"{'-'*70}")
            for mod in self.modifiers:
                lines.append(f"  • {mod.suggestion}")
                if mod.rationale:
                    lines.append(f"    Reason: {mod.rationale}")
                if mod.matched_conditions:
                    lines.append(f"    Triggers: {', '.join(mod.matched_conditions)}")
                lines.append(f"")
        
        if self.detected_features:
            lines.append(f"KEY FEATURES DETECTED:")
            lines.append(f"{'-'*70}")
            
            # Show only True boolean features and non-zero values
            for key, value in sorted(self.detected_features.items()):
                if isinstance(value, bool) and value:
                    lines.append(f"  ✓ {key}")
                elif isinstance(value, int) and value > 0:
                    lines.append(f"  • {key} = {value}")
        
        lines.append(f"{'='*70}")
        
        return "\n".join(lines)


@dataclass
class RuleSpec:
    """Specification for a single rule in the database."""
    
    name: str
    conditions: Dict[str, Any]
    reactant_features: Optional[Dict[str, Any]] = None
    priority: int = 0
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> RuleSpec:
        """Create RuleSpec from dictionary."""
        return cls(
            name=data.get("name", "unnamed"),
            conditions=data.get("conditions", {}),
            reactant_features=data.get("reactant_features"),
            priority=data.get("priority", 0)
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        result = {
            "name": self.name,
            "conditions": self.conditions,
        }
        if self.reactant_features:
            result["reactant_features"] = self.reactant_features
        if self.priority != 0:
            result["priority"] = self.priority
        return result
    
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
    
    when: List[str]
    suggestion: str
    rationale: Optional[str] = None
    priority: int = 0
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> ModifierSpec:
        """Create ModifierSpec from dictionary."""
        return cls(
            when=data.get("when", []),
            suggestion=data.get("suggestion", data.get("suggest", "")),
            rationale=data.get("rationale"),
            priority=data.get("priority", 0)
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        result = {
            "when": self.when,
            "suggestion": self.suggestion,
        }
        if self.rationale:
            result["rationale"] = self.rationale
        if self.priority != 0:
            result["priority"] = self.priority
        return result
    
    def matches(self, features: Dict[str, Any], symptoms: List[str] = None) -> bool:
        """
        Check if this modifier should be applied.
        
        Args:
            features: Detected molecular features
            symptoms: Optional list of observed symptoms
        
        Returns:
            True if any condition matches
        """
        symptoms = symptoms or []
        
        for condition in self.when:
            # Check for symptom-based conditions
            if condition.startswith("symptom:"):
                symptom_text = condition.replace("symptom:", "").strip()
                if any(symptom_text.lower() in s.lower() for s in symptoms):
                    return True
            # Check for feature-based conditions
            elif features.get(condition, False):
                return True
        
        return False

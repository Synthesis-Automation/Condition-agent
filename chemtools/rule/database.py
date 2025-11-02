"""
Rule Database Loader
====================

Loads and validates rule-based recommendation databases from JSON.
"""

from __future__ import annotations
from typing import Dict, Any, List, Optional
from pathlib import Path
import json
import logging

from .models import RuleSpec, ModifierSpec

logger = logging.getLogger(__name__)


class RuleDatabase:
    """Manages rule-based recommendation databases."""
    
    def __init__(
        self,
        applies_if: Optional[Dict[str, Any]] = None,
        default_rule: Optional[Dict[str, Any]] = None,
        base_rules: Optional[List[RuleSpec]] = None,
        modifiers: Optional[List[ModifierSpec]] = None,
        metadata: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize a rule database.
        
        Args:
            applies_if: Global applicability conditions
            default_rule: Fallback conditions if no base_rules match
            base_rules: Priority-ordered list of feature-matched rules
            modifiers: List of conditional modifiers
            metadata: Optional metadata (name, description, version, etc.)
        """
        self.applies_if = applies_if or {}
        self.default_rule = default_rule or {}
        self.base_rules = base_rules or []
        self.modifiers = modifiers or []
        self.metadata = metadata or {}
    
    @classmethod
    def from_file(cls, path: str | Path) -> RuleDatabase:
        """
        Load a rule database from a JSON file.
        
        Args:
            path: Path to JSON file
        
        Returns:
            RuleDatabase instance
        
        Example:
            >>> db = RuleDatabase.from_file("data/suzuki.json")
            >>> db.metadata['name']
            'Suzuki-Miyaura Coupling'
        """
        path = Path(path)
        
        if not path.exists():
            raise FileNotFoundError(f"Rule database not found: {path}")
        
        try:
            with open(path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            logger.info(f"Loaded rule database from {path}")
            return cls.from_dict(data)
        
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in {path}: {e}")
        except Exception as e:
            raise RuntimeError(f"Error loading rule database from {path}: {e}")
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> RuleDatabase:
        """
        Create a RuleDatabase from a dictionary.
        
        Args:
            data: Dictionary with rule database schema
        
        Returns:
            RuleDatabase instance
        """
        # Parse base_rules
        base_rules = []
        for rule_data in data.get("base_rules", []):
            try:
                base_rules.append(RuleSpec.from_dict(rule_data))
            except Exception as e:
                logger.error(f"Error parsing base rule: {e}")
                logger.debug(f"Rule data: {rule_data}")
        
        # Parse modifiers
        modifiers = []
        for mod_data in data.get("modifiers", []):
            try:
                modifiers.append(ModifierSpec.from_dict(mod_data))
            except Exception as e:
                logger.error(f"Error parsing modifier: {e}")
                logger.debug(f"Modifier data: {mod_data}")
        
        # Extract metadata
        metadata = {}
        for key in ["name", "reaction_type", "description", "version", "author", "reference"]:
            if key in data:
                metadata[key] = data[key]
        
        # Use reaction_type as name if name not present
        if "name" not in metadata and "reaction_type" in metadata:
            metadata["name"] = metadata["reaction_type"]
        
        return cls(
            applies_if=data.get("applies_if", {}),
            default_rule=data.get("default_rule", {}),
            base_rules=base_rules,
            modifiers=modifiers,
            metadata=metadata
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert to dictionary representation.
        
        Returns:
            Dictionary with full database schema
        """
        return {
            **self.metadata,
            "applies_if": self.applies_if,
            "default_rule": self.default_rule,
            "base_rules": [rule.to_dict() for rule in self.base_rules],
            "modifiers": [mod.to_dict() for mod in self.modifiers]
        }
    
    def check_applies(self, features: Dict[str, Any]) -> bool:
        """
        Check if this database applies to the given features.
        
        Args:
            features: Detected features from reaction
        
        Returns:
            True if applies_if conditions are satisfied
        
        Example:
            >>> db.applies_if = {"all": ["sp2_halide_present", "sp2_boron_present"]}
            >>> db.check_applies({"sp2_halide_present": True, "sp2_boron_present": True})
            True
        """
        if not self.applies_if:
            # No conditions means always applies
            return True
        
        # Handle "all" logic
        if "all" in self.applies_if:
            required_features = self.applies_if["all"]
            if isinstance(required_features, list):
                return all(features.get(feat, False) for feat in required_features)
        
        # Handle "any" logic
        if "any" in self.applies_if:
            required_features = self.applies_if["any"]
            if isinstance(required_features, list):
                return any(features.get(feat, False) for feat in required_features)
        
        # Handle direct feature requirements
        for key, value in self.applies_if.items():
            if key not in ["all", "any"]:
                if features.get(key) != value:
                    return False
        
        return True
    
    def find_matching_rule(
        self,
        features: Dict[str, Any],
        use_default: bool = True
    ) -> Optional[RuleSpec]:
        """
        Find the best matching base rule for given features.
        
        Args:
            features: Detected features from reaction
            use_default: If True, return default_rule if no base_rules match
        
        Returns:
            Matching RuleSpec or None
        
        Note:
            base_rules are checked in order; first match wins.
        """
        # Try base_rules in priority order
        for rule in self.base_rules:
            matches, matched_features = rule.matches(features)
            if matches:
                logger.debug(f"Matched rule '{rule.name}' with features: {matched_features}")
                return rule
        
        # Fall back to default_rule if available
        if use_default and self.default_rule:
            logger.debug("No base_rules matched, using default_rule")
            # Convert default_rule dict to RuleSpec
            return RuleSpec(
                name="default",
                conditions=self.default_rule,
                reactant_features={}
            )
        
        return None
    
    def find_matching_modifiers(
        self,
        features: Dict[str, Any],
        symptoms: Optional[List[str]] = None
    ) -> List[ModifierSpec]:
        """
        Find all modifiers that match the given features/symptoms.
        
        Args:
            features: Detected features from reaction
            symptoms: Optional list of symptom strings (e.g., ["low_yield", "side_products"])
        
        Returns:
            List of matching ModifierSpec objects
        """
        matching = []
        
        for modifier in self.modifiers:
            if modifier.matches(features, symptoms or []):
                matching.append(modifier)
                logger.debug(f"Matched modifier: {modifier.suggestion}")
        
        return matching
    
    def get_summary(self) -> str:
        """
        Get a human-readable summary of this database.
        
        Returns:
            Formatted summary string
        """
        lines = []
        
        if self.metadata.get("name"):
            lines.append(f"Database: {self.metadata['name']}")
        
        if self.metadata.get("description"):
            lines.append(f"Description: {self.metadata['description']}")
        
        if self.applies_if:
            lines.append(f"Applies when: {self.applies_if}")
        
        lines.append(f"Base rules: {len(self.base_rules)}")
        lines.append(f"Modifiers: {len(self.modifiers)}")
        
        if self.default_rule:
            lines.append(f"Has default rule: Yes")
        
        return "\n".join(lines)
    
    def validate(self) -> List[str]:
        """
        Validate the database schema and return any issues.
        
        Returns:
            List of validation error messages (empty if valid)
        """
        issues = []
        
        # Check for base_rules or default_rule
        if not self.base_rules and not self.default_rule:
            issues.append("No base_rules or default_rule defined")
        
        # Validate base_rules
        for i, rule in enumerate(self.base_rules):
            if not rule.name:
                issues.append(f"base_rules[{i}]: Missing name")
            if not rule.conditions:
                issues.append(f"base_rules[{i}]: Missing conditions")
        
        # Validate modifiers
        for i, mod in enumerate(self.modifiers):
            if not mod.when:
                issues.append(f"modifiers[{i}]: Missing 'when' conditions")
            if not mod.suggestion:
                issues.append(f"modifiers[{i}]: Missing suggestion")
        
        return issues

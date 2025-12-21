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
from ..util.boolean_expr import evaluate as _eval_bool_expr

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
        
        # Extract metadata (supports both legacy flat schema and v2.0 schema).
        metadata: Dict[str, Any] = {}

        meta_block = data.get("metadata")
        if isinstance(meta_block, dict):
            metadata.update({k: v for k, v in meta_block.items() if v is not None})

        reaction_block = data.get("reaction")
        if isinstance(reaction_block, dict):
            family = reaction_block.get("family")
            if isinstance(family, str) and family.strip():
                metadata.setdefault("reaction_family", family.strip())

        for key in ["name", "reaction_type", "description", "version", "author", "reference"]:
            if key in data and data[key] is not None:
                metadata.setdefault(key, data[key])

        # Use a reasonable fallback label if none exists.
        if "name" not in metadata:
            metadata["name"] = metadata.get("reaction_family") or metadata.get("reaction_type") or "Unnamed"
        
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

        # Expression form: applies_if can be a string or {"expr": "..."}
        expr: Optional[str] = None
        if isinstance(self.applies_if, str):
            expr = self.applies_if
        elif isinstance(self.applies_if, dict):
            expr_value = self.applies_if.get("expr")
            if isinstance(expr_value, str):
                expr = expr_value

        if isinstance(expr, str) and expr.strip():
            if not _eval_bool_expr(expr, features):
                return False
        
        if not isinstance(self.applies_if, dict):
            return True

        # Combine "all" and "any" constraints when both exist.
        required_all = self.applies_if.get("all")
        if isinstance(required_all, list) and required_all:
            if not all(features.get(feat, False) for feat in required_all):
                return False
        
        required_any = self.applies_if.get("any")
        if isinstance(required_any, list) and required_any:
            if not any(features.get(feat, False) for feat in required_any):
                return False
        
        # Handle direct feature requirements
        for key, value in self.applies_if.items():
            if key not in ["all", "any", "expr"]:
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
        
        # Fall back to default_rule if available.
        if use_default and self.default_rule:
            logger.debug("No base_rules matched, using default_rule")
            default_rule = self.default_rule

            # Schema v2.0 default_rule: {"id": "...", "description": "...", "conditions": {...}}
            default_conditions: Dict[str, Any] = {}
            default_id: Optional[str] = None
            default_description: Optional[str] = None

            if isinstance(default_rule, dict):
                default_id = default_rule.get("id") if isinstance(default_rule.get("id"), str) else None
                default_description = (
                    default_rule.get("description")
                    if isinstance(default_rule.get("description"), str)
                    else None
                )
                if isinstance(default_rule.get("conditions"), dict):
                    default_conditions = default_rule["conditions"]
                else:
                    default_conditions = default_rule

            name = default_description or default_id or "default"
            return RuleSpec(
                name=name,
                conditions=default_conditions,
                reactant_features={},
                rule_id=default_id,
                description=default_description,
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

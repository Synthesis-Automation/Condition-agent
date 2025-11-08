"""
Tests for Rule-Based Recommendation Engine
===========================================

Test suite for models, database, analyzer, and engine.
"""

import pytest
from pathlib import Path
from typing import Dict, Any

from chemtools.rule.models import RuleSpec, ModifierSpec, AppliedRule, ConditionRecommendation
from chemtools.rule.database import RuleDatabase
from chemtools.rule.analyzer import FeatureAnalyzer
from chemtools.rule.engine import RuleEngine


# ============================================================================
# Test Data
# ============================================================================

SAMPLE_FEATURES = {
    "sp2_halide_present": True,
    "sp2_boron_present": True,
    "sp2_chloride_present": True,
    "aryl_halide_present": True,
    "electron_withdrawing_present": True,
    "halogen_count": 1,
    "ortho_substitution": False
}

SAMPLE_RULE_DICT = {
    "name": "activated_aryl_chloride",
    "reactant_features": {
        "and": ["sp2_chloride_present", "electron_withdrawing_present"]
    },
    "conditions": {
        "pd_source": "Pd(PPh3)4 or Pd(OAc)2",
        "ligand": "XPhos or SPhos",
        "base": "K2CO3 or K3PO4",
        "solvent": "1,4-dioxane or toluene",
        "temp": "100-110°C"
    }
}

SAMPLE_MODIFIER_DICT = {
    "when": ["ortho_substitution", "steric_hindrance"],
    "suggestion": "Use bulky ligand (SPhos, XPhos) and extend reaction time",
    "rationale": "Ortho-substituted substrates require bulky, electron-rich ligands"
}


# ============================================================================
# Test Models
# ============================================================================

class TestRuleSpec:
    """Test RuleSpec model and matching logic."""
    
    def test_from_dict(self):
        """Test creating RuleSpec from dictionary."""
        rule = RuleSpec.from_dict(SAMPLE_RULE_DICT)
        
        assert rule.name == "activated_aryl_chloride"
        assert rule.conditions["pd_source"] == "Pd(PPh3)4 or Pd(OAc)2"
        assert "and" in rule.reactant_features
    
    def test_matches_and_logic(self):
        """Test AND feature matching."""
        rule = RuleSpec.from_dict(SAMPLE_RULE_DICT)
        
        # Should match (has both required features)
        matches, matched = rule.matches(SAMPLE_FEATURES)
        assert matches is True
        assert "sp2_chloride_present" in matched
        assert "electron_withdrawing_present" in matched
        
        # Should not match (missing electron_withdrawing)
        features_no_ew = SAMPLE_FEATURES.copy()
        features_no_ew["electron_withdrawing_present"] = False
        matches, matched = rule.matches(features_no_ew)
        assert matches is False
    
    def test_matches_any_logic(self):
        """Test ANY (OR) feature matching."""
        rule_dict = {
            "name": "any_halide",
            "reactant_features": {
                "any": ["sp2_chloride_present", "sp2_bromide_present", "sp2_iodide_present"]
            },
            "conditions": {"pd_source": "Pd(PPh3)4"}
        }
        rule = RuleSpec.from_dict(rule_dict)
        
        # Should match (has sp2_chloride_present)
        matches, matched = rule.matches(SAMPLE_FEATURES)
        assert matches is True
        assert "sp2_chloride_present" in matched
        
        # Should not match (has none)
        features_no_halide = {k: False for k in SAMPLE_FEATURES}
        matches, matched = rule.matches(features_no_halide)
        assert matches is False
    
    def test_to_dict(self):
        """Test converting RuleSpec to dictionary."""
        rule = RuleSpec.from_dict(SAMPLE_RULE_DICT)
        rule_dict = rule.to_dict()
        
        assert rule_dict["name"] == "activated_aryl_chloride"
        assert "reactant_features" in rule_dict
        assert "conditions" in rule_dict


class TestModifierSpec:
    """Test ModifierSpec model and matching logic."""
    
    def test_from_dict(self):
        """Test creating ModifierSpec from dictionary."""
        modifier = ModifierSpec.from_dict(SAMPLE_MODIFIER_DICT)
        
        assert "ortho_substitution" in modifier.when
        assert "Use bulky ligand" in modifier.suggestion
        assert modifier.rationale is not None
    
    def test_matches_features(self):
        """Test feature-based matching."""
        modifier = ModifierSpec.from_dict(SAMPLE_MODIFIER_DICT)
        
        # Should not match (no ortho_substitution or steric_hindrance)
        assert modifier.matches(SAMPLE_FEATURES, []) is False
        
        # Should match (has ortho_substitution)
        features_ortho = SAMPLE_FEATURES.copy()
        features_ortho["ortho_substitution"] = True
        assert modifier.matches(features_ortho, []) is True
    
    def test_matches_symptoms(self):
        """Test symptom-based matching."""
        modifier_dict = {
            "when": ["symptom:low_yield", "symptom:side_products"],
            "suggestion": "Increase temperature or reaction time"
        }
        modifier = ModifierSpec.from_dict(modifier_dict)
        
        # Should match (has symptom)
        assert modifier.matches(SAMPLE_FEATURES, ["low_yield"]) is True
        
        # Should not match (no symptoms)
        assert modifier.matches(SAMPLE_FEATURES, []) is False


class TestConditionRecommendation:
    """Test ConditionRecommendation model and formatting."""
    
    def test_to_dict(self):
        """Test converting recommendation to dictionary."""
        applied_rule = AppliedRule(
            name="test_rule",
            conditions={"pd_source": "Pd(PPh3)4"},
            matched_features=["sp2_halide_present"],
            confidence=0.95
        )
        
        rec = ConditionRecommendation(
            reaction_smiles="Brc1ccccc1>>c1ccccc1",
            base_rule=applied_rule,
            modifiers=[],
            detected_features=SAMPLE_FEATURES
        )
        
        rec_dict = rec.to_dict()
        assert rec_dict["reaction_smiles"] == "Brc1ccccc1>>c1ccccc1"
        assert rec_dict["base_rule"]["name"] == "test_rule"
        assert "detected_features" in rec_dict
    
    def test_format_summary(self):
        """Test human-readable formatting."""
        applied_rule = AppliedRule(
            name="test_rule",
            conditions={"pd_source": "Pd(PPh3)4", "base": "K2CO3"},
            matched_features=["sp2_halide_present"],
            confidence=0.95
        )
        
        rec = ConditionRecommendation(
            reaction_smiles="Brc1ccccc1>>c1ccccc1",
            base_rule=applied_rule,
            modifiers=[],
            detected_features={}
        )
        
        summary = rec.format_summary()
        assert "CONDITION RECOMMENDATION" in summary
        assert "pd_source" in summary
        assert "Pd(PPh3)4" in summary


# ============================================================================
# Test Database
# ============================================================================

class TestRuleDatabase:
    """Test RuleDatabase loading and validation."""
    
    def test_from_dict_basic(self):
        """Test creating database from dictionary."""
        db_dict = {
            "name": "Test Database",
            "applies_if": {"all": ["sp2_halide_present"]},
            "default_rule": {"pd_source": "Pd(PPh3)4"},
            "base_rules": [SAMPLE_RULE_DICT],
            "modifiers": [SAMPLE_MODIFIER_DICT]
        }
        
        db = RuleDatabase.from_dict(db_dict)
        
        assert db.metadata["name"] == "Test Database"
        assert len(db.base_rules) == 1
        assert len(db.modifiers) == 1
    
    def test_check_applies(self):
        """Test applicability checking."""
        db_dict = {
            "applies_if": {"all": ["sp2_halide_present", "sp2_boron_present"]},
            "default_rule": {}
        }
        db = RuleDatabase.from_dict(db_dict)
        
        # Should apply (has both required features)
        assert db.check_applies(SAMPLE_FEATURES) is True
        
        # Should not apply (missing sp2_boron_present)
        features_no_boron = SAMPLE_FEATURES.copy()
        features_no_boron["sp2_boron_present"] = False
        assert db.check_applies(features_no_boron) is False
    
    def test_find_matching_rule(self):
        """Test finding matching base rule."""
        db_dict = {
            "base_rules": [SAMPLE_RULE_DICT],
            "default_rule": {"pd_source": "default"}
        }
        db = RuleDatabase.from_dict(db_dict)
        
        # Should match the base rule
        rule = db.find_matching_rule(SAMPLE_FEATURES)
        assert rule is not None
        assert rule.name == "activated_aryl_chloride"
        
        # Should fall back to default
        features_no_match = {k: False for k in SAMPLE_FEATURES}
        rule = db.find_matching_rule(features_no_match)
        assert rule is not None
        assert rule.name == "default"
    
    def test_validate(self):
        """Test database validation."""
        # Valid database
        db_dict = {
            "base_rules": [SAMPLE_RULE_DICT],
            "modifiers": [SAMPLE_MODIFIER_DICT]
        }
        db = RuleDatabase.from_dict(db_dict)
        issues = db.validate()
        assert len(issues) == 0
        
        # Invalid database (no rules)
        db_empty = RuleDatabase()
        issues = db_empty.validate()
        assert len(issues) > 0


# ============================================================================
# Test Analyzer
# ============================================================================

class TestFeatureAnalyzer:
    """Test FeatureAnalyzer reaction parsing and feature detection."""
    
    def test_parse_reaction(self):
        """Test reaction SMILES parsing."""
        analyzer = FeatureAnalyzer()
        
        # Standard reaction arrow
        reactants, products = analyzer._parse_reaction("Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1")
        assert len(reactants) == 2
        assert "Br" in reactants[0]
        
        # Single > format
        reactants, products = analyzer._parse_reaction("Brc1ccccc1>solvent>product")
        assert len(reactants) == 1
    
    def test_analyze_reactant(self):
        """Test single reactant analysis."""
        analyzer = FeatureAnalyzer()
        
        # Aryl bromide
        features = analyzer.analyze_reactant("Brc1ccccc1")
        assert features.get("sp2_halide_present") is True
        assert features.get("sp2_bromide_present") is True
    
    def test_analyze_reaction_union(self):
        """Test reaction analysis with union combine method."""
        analyzer = FeatureAnalyzer()
        
        # Suzuki coupling: aryl bromide + boronic acid
        features = analyzer.analyze_reaction(
            "Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1",
            combine_method="union"
        )
        
        # Should detect both halide and boron
        assert features.get("sp2_halide_present") is True
        assert features.get("sp2_boron_present") is True


# ============================================================================
# Test Engine
# ============================================================================

class TestRuleEngine:
    """Test RuleEngine integration."""
    
    @pytest.fixture
    def test_database(self) -> RuleDatabase:
        """Create a test database."""
        db_dict = {
            "name": "Test Suzuki Database",
            "applies_if": {"all": ["sp2_halide_present", "sp2_boron_present"]},
            "default_rule": {
                "pd_source": "Pd(PPh3)4",
                "base": "K2CO3",
                "solvent": "toluene"
            },
            "base_rules": [
                {
                    "name": "activated_chloride",
                    "reactant_features": {
                        "and": ["sp2_chloride_present", "electron_withdrawing_present"]
                    },
                    "conditions": {
                        "pd_source": "Pd(OAc)2",
                        "ligand": "XPhos",
                        "base": "K3PO4"
                    }
                }
            ],
            "modifiers": [
                {
                    "when": ["symptom:low_yield"],
                    "suggestion": "Increase temperature to 110°C"
                }
            ]
        }
        return RuleDatabase.from_dict(db_dict)
    
    def test_recommend_basic(self, test_database):
        """Test basic recommendation workflow."""
        engine = RuleEngine(database=test_database)
        
        # Simple Suzuki coupling (should use default rule)
        rec = engine.recommend("Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1")
        
        assert rec.base_rule.name == "default"
        assert "Pd(PPh3)4" in rec.base_rule.conditions["pd_source"]
    
    def test_recommend_with_symptoms(self, test_database):
        """Test recommendation with symptoms."""
        engine = RuleEngine(database=test_database)
        
        rec = engine.recommend(
            "Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1",
            symptoms=["low_yield"]
        )
        
        # Should apply the low_yield modifier
        assert len(rec.modifiers) > 0
        assert "temperature" in rec.modifiers[0].suggestion.lower()
    
    def test_validate_database(self, test_database):
        """Test database validation."""
        engine = RuleEngine(database=test_database)
        issues = engine.validate_database()
        assert len(issues) == 0


# ============================================================================
# Integration Tests
# ============================================================================

class TestIntegration:
    """Integration tests with real suzuki.json if available."""
    
    def test_load_suzuki_database(self):
        """Test loading the real suzuki.json database."""
        # Try to find suzuki.json
        possible_paths = [
            Path("data/rule_db_v2/suzuki.json"),  # Updated to v2
            Path("data/suzuki.json"),
            Path("data/protocol_db_v2/suzuki.json")  # Updated to v2
        ]
        
        suzuki_path = None
        for path in possible_paths:
            if path.exists():
                suzuki_path = path
                break
        
        if suzuki_path is None:
            pytest.skip("suzuki.json not found")
        
        # Load database
        db = RuleDatabase.from_file(suzuki_path)
        
        assert db.metadata.get("name") is not None
        assert len(db.base_rules) > 0 or db.default_rule
        
        # Validate
        issues = db.validate()
        assert len(issues) == 0, f"Validation issues: {issues}"
    
    def test_recommend_with_suzuki(self):
        """Test full recommendation workflow with suzuki.json."""
        # Try to find suzuki.json
        possible_paths = [
            Path("data/rule_db_v2/suzuki.json"),  # Updated to v2
            Path("data/suzuki.json"),
            Path("data/protocol_db_v2/suzuki.json")  # Updated to v2
        ]
        
        suzuki_path = None
        for path in possible_paths:
            if path.exists():
                suzuki_path = path
                break
        
        if suzuki_path is None:
            pytest.skip("suzuki.json not found")
        
        # Create engine
        engine = RuleEngine.from_file(suzuki_path)
        
        # Test with a simple Suzuki coupling
        rec = engine.recommend("Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1")
        
        assert rec.base_rule is not None
        assert rec.base_rule.conditions is not None
        
        # Should have detected features
        if rec.detected_features:
            assert rec.detected_features.get("sp2_halide_present") is True
            assert rec.detected_features.get("sp2_boron_present") is True


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

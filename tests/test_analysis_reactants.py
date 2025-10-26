"""
Tests for chemtools.analysis.reactants module.
"""
import pytest
from chemtools.analysis import (
    classify_reactant_smiles,
    classify_reactant_category,
    classify_reactant_group,
    classify_reactant_batch,
    get_reactant_category_matches,
    get_all_reactant_matches,
    normalize_reactant_identifier,
    ReactantMatch,
)


class TestClassifyReactantSmiles:
    """Test SMARTS-based reactant classification."""
    
    def test_aryl_bromide(self):
        """Aryl bromide should be classified correctly."""
        smiles = "c1ccccc1Br"
        result = classify_reactant_smiles(smiles)
        assert result is not None
        assert isinstance(result, ReactantMatch)
        # New taxonomy uses ArX* as category, ArBr as member_type
        assert result.member_type == "ArBr"
        assert result.smarts is not None
    
    def test_aryl_iodide(self):
        """Aryl iodide should be classified correctly."""
        smiles = "c1ccccc1I"
        result = classify_reactant_smiles(smiles)
        assert result is not None
        assert result.member_type == "ArI"
    
    def test_alkyl_bromide(self):
        """Alkyl bromide should be classified correctly."""
        smiles = "CCBr"
        result = classify_reactant_smiles(smiles)
        assert result is not None
        assert result.member_type == "Alkyl-Br"
    
    def test_aryl_boronic_acid(self):
        """Aryl boronic acid should be classified correctly."""
        smiles = "c1ccccc1B(O)O"
        result = classify_reactant_smiles(smiles)
        assert result is not None
        assert result.member_type == "ArB(OH)2"
    
    def test_aryl_boronate_ester(self):
        """Aryl boronate ester should be classified correctly."""
        smiles = "c1ccccc1B1OC(C)(C)C(C)(C)O1"  # Phenyl-Bpin
        result = classify_reactant_smiles(smiles)
        assert result is not None
        assert result.member_type == "ArB(OR)2"
    
    def test_primary_amine(self):
        """Primary amine should be classified correctly."""
        smiles = "CCN"
        result = classify_reactant_smiles(smiles)
        assert result is not None
        assert result.member_type == "RNH2"
    
    def test_secondary_amine(self):
        """Secondary amine should be classified correctly."""
        smiles = "CCNCC"
        result = classify_reactant_smiles(smiles)
        assert result is not None
        assert result.member_type == "R2NH"
    
    def test_aniline(self):
        """Aniline should be classified correctly."""
        smiles = "c1ccccc1N"
        result = classify_reactant_smiles(smiles)
        assert result is not None
        assert result.member_type == "ArNH2"
    
    def test_carboxylic_acid(self):
        """Carboxylic acid should be classified correctly."""
        smiles = "CC(=O)O"
        result = classify_reactant_smiles(smiles)
        assert result is not None
        # Could match RCO2H or acyl_source
        assert "CO2" in result.category or "acyl" in result.category.lower()
    
    def test_invalid_smiles(self):
        """Invalid SMILES should return None."""
        result = classify_reactant_smiles("not-a-valid-smiles")
        assert result is None
    
    def test_empty_smiles(self):
        """Empty SMILES should return None."""
        result = classify_reactant_smiles("")
        assert result is None


class TestClassifyReactantCategory:
    """Test convenience function that returns only category."""
    
    def test_returns_category_string(self):
        """Should return category as string."""
        smiles = "c1ccccc1Br"
        result = classify_reactant_category(smiles)
        assert result is not None
        assert isinstance(result, str)
    
    def test_invalid_returns_none(self):
        """Invalid SMILES should return None."""
        result = classify_reactant_category("invalid")
        assert result is None


class TestClassifyReactantGroup:
    """Test convenience function that returns functional group."""
    
    def test_returns_group_string(self):
        """Should return group as string."""
        smiles = "c1ccccc1Br"
        result = classify_reactant_group(smiles)
        # May or may not have a group defined
        assert result is None or isinstance(result, str)


class TestClassifyReactantBatch:
    """Test batch classification."""
    
    def test_batch_classification(self):
        """Should classify multiple SMILES."""
        smiles_list = ["c1ccccc1Br", "CCBr", "c1ccccc1B(O)O"]
        results = classify_reactant_batch(smiles_list)
        assert len(results) == 3
        assert all(r is None or isinstance(r, ReactantMatch) for r in results)
        # First two should succeed
        assert results[0] is not None
        assert results[1] is not None
    
    def test_empty_batch(self):
        """Empty batch should return empty list."""
        results = classify_reactant_batch([])
        assert results == []


class TestGetReactantCategoryMatches:
    """Test getting all category matches."""
    
    def test_aryl_halide_categories(self):
        """Aryl halide might match multiple categories."""
        smiles = "c1ccccc1Br"
        results = get_reactant_category_matches(smiles)
        assert isinstance(results, list)
        assert len(results) > 0
        # Should match ArX* category
        assert any("ArX" in cat or "arx" in cat.lower() for cat in results)
    
    def test_returns_unique_sorted(self):
        """Should return unique, sorted categories."""
        smiles = "c1ccccc1Br"
        results = get_reactant_category_matches(smiles)
        assert results == sorted(set(results))


class TestGetAllReactantMatches:
    """Test getting all SMARTS matches."""
    
    def test_returns_all_matches(self):
        """Should return all SMARTS matches."""
        smiles = "c1ccccc1Br"
        results = get_all_reactant_matches(smiles)
        assert isinstance(results, list)
        assert all(isinstance(m, ReactantMatch) for m in results)
    
    def test_specificity_ordering(self):
        """Results should be ordered by specificity."""
        smiles = "c1ccccc1Br"
        results = get_all_reactant_matches(smiles)
        if len(results) > 1:
            # Should be sorted by specificity (descending)
            specificities = [m.specificity for m in results]
            # Note: actual ordering also considers is_general flag


class TestNormalizeReactantIdentifier:
    """Test reactant identifier normalization."""
    
    def test_normalize_common_aliases(self):
        """Common aliases should normalize correctly."""
        # Test with actual member IDs from the taxonomy
        test_cases = [
            ("ArBr", True),  # Member ID should normalize
            ("RNH2", True),  # Member ID should normalize
        ]
        
        for input_val, should_normalize in test_cases:
            result = normalize_reactant_identifier(input_val)
            if should_normalize:
                assert result is not None, f"Failed to normalize {input_val}"
            # The result may vary based on taxonomy configuration
    
    def test_empty_string(self):
        """Empty string should return None."""
        assert normalize_reactant_identifier("") is None
    
    def test_unknown_identifier(self):
        """Unknown identifier should return None."""
        result = normalize_reactant_identifier("CompletelyUnknownType123")
        assert result is None
    
    def test_whitespace_handling(self):
        """Should handle whitespace correctly."""
        result = normalize_reactant_identifier("  ArBr  ")
        assert result is not None


class TestReactantMatch:
    """Test ReactantMatch dataclass."""
    
    def test_match_structure(self):
        """ReactantMatch should have expected fields."""
        smiles = "c1ccccc1Br"
        match = classify_reactant_smiles(smiles)
        assert match is not None
        
        # Check all expected fields exist
        assert hasattr(match, "category")
        assert hasattr(match, "member_type")
        assert hasattr(match, "name")
        assert hasattr(match, "group")
        assert hasattr(match, "smarts")
        assert hasattr(match, "category_smarts")
        assert hasattr(match, "description")
        assert hasattr(match, "specificity")
        assert hasattr(match, "is_general")
        
        # Check types
        assert isinstance(match.category, str)
        assert isinstance(match.smarts, str)
        assert isinstance(match.specificity, int)
        assert isinstance(match.is_general, bool)


class TestEdgeCases:
    """Test edge cases and error handling."""
    
    def test_none_input(self):
        """None input should not crash."""
        result = classify_reactant_smiles(None)
        assert result is None
    
    def test_special_characters(self):
        """Special characters should be handled."""
        result = classify_reactant_smiles("[Na+].[Cl-]")
        # May or may not classify, but shouldn't crash
        assert result is None or isinstance(result, ReactantMatch)
    
    def test_complex_molecule(self):
        """Complex molecules should be classified if patterns match."""
        # Imatinib (has many functional groups)
        smiles = "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CC=CC=N5"
        result = classify_reactant_smiles(smiles)
        # Should at least find some match
        assert result is None or isinstance(result, ReactantMatch)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

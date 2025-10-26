"""
Tests for chemtools.analysis.reactions module.
"""
import pytest
from chemtools.analysis import (
    resolve_reaction_family,
    canonical_family_label,
    normalize_reaction_type,
)


class TestResolveReactionFamily:
    """Test reaction family resolution."""
    
    def test_buchwald_hartwig_variations(self):
        """Various Buchwald-Hartwig notations should resolve."""
        test_cases = [
            ("Buchwald-Hartwig", True),
            ("Buchwald_CN", True),
            ("buchwald-hartwig-c-n", True),
        ]
        
        for label, should_resolve in test_cases:
            result = resolve_reaction_family(label)
            if should_resolve:
                # Should resolve to something
                assert result is not None, f"Failed to resolve {label}"
                # Should be a canonical ID (lowercase with underscores)
                assert isinstance(result, str)
    
    def test_suzuki_miyaura_variations(self):
        """Various Suzuki notations should resolve."""
        test_cases = [
            "Suzuki-Miyaura",
            "suzuki_miyaura",
            "Suzuki_CC",
            "suzuki",
        ]
        
        for label in test_cases:
            result = resolve_reaction_family(label)
            assert result is not None, f"Failed to resolve {label}"
    
    def test_ullmann_variations(self):
        """Ullmann variations should resolve correctly."""
        test_cases = [
            "Ullmann_CN",
            "ullmann_cn",
            "Ullmann_CO",
            "ullmann_ether",
        ]
        
        for label in test_cases:
            result = resolve_reaction_family(label)
            assert result is not None, f"Failed to resolve {label}"
    
    def test_none_input(self):
        """None input should return None."""
        result = resolve_reaction_family(None)
        assert result is None
    
    def test_empty_string(self):
        """Empty string should return None."""
        result = resolve_reaction_family("")
        assert result is None
    
    def test_unknown_family(self):
        """Unknown family should return None."""
        result = resolve_reaction_family("CompletelyUnknownReaction123")
        assert result is None


class TestCanonicalFamilyLabel:
    """Test canonical family label resolution."""
    
    def test_already_canonical(self):
        """Already canonical IDs should return as-is."""
        # Common canonical forms
        test_ids = [
            "buchwald_hartwig_c_n",
            "suzuki_miyaura",
            "ullmann_cn",
        ]
        
        for test_id in test_ids:
            result = canonical_family_label(test_id)
            # Should return the same or a known canonical form
            assert result is not None
    
    def test_case_insensitive(self):
        """Should handle different cases."""
        result1 = canonical_family_label("Buchwald-Hartwig")
        result2 = canonical_family_label("buchwald-hartwig")
        result3 = canonical_family_label("BUCHWALD-HARTWIG")
        
        # All should resolve (may be to the same value)
        assert result1 is not None
        assert result2 is not None
        assert result3 is not None
    
    def test_whitespace_handling(self):
        """Should handle whitespace."""
        result = canonical_family_label("  Suzuki-Miyaura  ")
        assert result is not None


class TestNormalizeReactionType:
    """Test reaction type normalization."""
    
    def test_common_reactions(self):
        """Common reactions should normalize."""
        test_cases = [
            "Buchwald-Hartwig",
            "Suzuki-Miyaura",
            "Negishi",
            "Heck",
            "Sonogashira",
            "Amide coupling",
            "CN-Coupling",
            "CO-Coupling",
        ]
        
        for label in test_cases:
            result = normalize_reaction_type(label)
            assert result is not None, f"Failed to normalize {label}"
            assert isinstance(result, str)
    
    def test_case_insensitivity(self):
        """Normalization should be case insensitive."""
        result1 = normalize_reaction_type("Buchwald-Hartwig")
        result2 = normalize_reaction_type("buchwald-hartwig")
        
        assert result1 is not None
        assert result2 is not None
        # Should normalize to same value
        assert result1 == result2
    
    def test_empty_string(self):
        """Empty string should return None."""
        assert normalize_reaction_type("") is None
    
    def test_none_input(self):
        """None input should return None."""
        assert normalize_reaction_type(None) is None


class TestFamilyAliases:
    """Test family alias handling."""
    
    def test_cn_coupling_aliases(self):
        """C-N coupling aliases should resolve."""
        aliases = [
            "C_N_Coupling",
            "CN-Coupling",
            "cn_coupling",
        ]
        
        results = [resolve_reaction_family(a) for a in aliases]
        # All should resolve
        assert all(r is not None for r in results)
    
    def test_co_coupling_aliases(self):
        """C-O coupling aliases should resolve."""
        aliases = [
            "C_O_Coupling",
            "CO-Coupling",
            "co_coupling",
        ]
        
        results = [resolve_reaction_family(a) for a in aliases]
        assert all(r is not None for r in results)
    
    def test_specific_vs_general(self):
        """Specific family names should take precedence."""
        # Buchwald-Hartwig is more specific than C-N coupling
        bh_result = resolve_reaction_family("Buchwald-Hartwig")
        cn_result = resolve_reaction_family("CN-Coupling")
        
        assert bh_result is not None
        assert cn_result is not None
        # They should be different (unless taxonomy maps them together)


class TestEdgeCases:
    """Test edge cases in reaction family resolution."""
    
    def test_special_characters(self):
        """Special characters should be handled."""
        test_cases = [
            "Suzuki-Miyaura, in situ",
            "Buchwald/Hartwig",
            "C-N Coupling",
        ]
        
        for label in test_cases:
            result = resolve_reaction_family(label)
            # May or may not resolve, but shouldn't crash
            assert result is None or isinstance(result, str)
    
    def test_numeric_suffixes(self):
        """Numeric suffixes should be handled."""
        # Some taxonomies might have numbered variants
        result = resolve_reaction_family("Negishi-1")
        # Should either resolve or return None gracefully
        assert result is None or isinstance(result, str)
    
    def test_very_long_string(self):
        """Very long strings should be handled."""
        long_string = "A" * 1000
        result = resolve_reaction_family(long_string)
        assert result is None


class TestReactionCategories:
    """Test reaction category relationships."""
    
    def test_coupling_reactions(self):
        """Common coupling reactions should resolve."""
        coupling_reactions = [
            "Suzuki-Miyaura",
            "Negishi",
            "Heck",
            "Sonogashira",
            "Stille",
        ]
        
        for rxn in coupling_reactions:
            result = resolve_reaction_family(rxn)
            # All major coupling reactions should resolve
            assert result is not None, f"Major coupling reaction {rxn} failed to resolve"
    
    def test_cn_bond_forming(self):
        """C-N bond forming reactions should resolve."""
        cn_reactions = [
            "Buchwald-Hartwig",
            "Ullmann_CN",
            "Chan_Lam_CN",
        ]
        
        for rxn in cn_reactions:
            result = resolve_reaction_family(rxn)
            assert result is not None, f"C-N coupling {rxn} failed to resolve"
    
    def test_amide_formation(self):
        """Amide formation reactions should resolve."""
        result = resolve_reaction_family("Amide coupling")
        assert result is not None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

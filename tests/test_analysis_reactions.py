"""
Tests for chemtools.analysis.reactions module (taxonomy v2).
"""
import pytest
from chemtools.analysis import (
    resolve_reaction_family,
    canonical_family_label,
    normalize_reaction_type,
)


class TestResolveReactionFamily:
    """Test reaction family resolution for taxonomy v2."""

    def test_known_aliases(self):
        """Known v2 aliases should resolve to canonical IDs."""
        test_cases = {
            "Buchwald-Hartwig": "c_n_cross_coupling",
            "Suzuki": "suzuki_miyaura",
            "Sonogashira cross-coupling": "sonogashira",
            "SNAr C-N": "snar_cn",
            "SN2": "sn2_substitution",
            "Amide formation": "amide_coupling",
            "Esterification": "esterification",
        }

        for label, expected in test_cases.items():
            assert resolve_reaction_family(label) == expected

    def test_slug_variants(self):
        """Slugged variants should resolve when possible."""
        assert resolve_reaction_family("Suzuki-Miyaura") == "suzuki_miyaura"

    def test_none_input(self):
        """None input should return None."""
        assert resolve_reaction_family(None) is None

    def test_empty_string(self):
        """Empty string should return None."""
        assert resolve_reaction_family("") is None

    def test_unknown_family(self):
        """Unknown family should return None."""
        assert resolve_reaction_family("CompletelyUnknownReaction123") is None


class TestCanonicalFamilyLabel:
    """Test canonical family label resolution."""

    def test_already_canonical(self):
        """Already canonical IDs should return as-is."""
        test_ids = [
            "suzuki_miyaura",
            "c_n_cross_coupling",
            "snar_cn",
        ]
        for test_id in test_ids:
            assert canonical_family_label(test_id) == test_id

    def test_case_insensitive(self):
        """Should handle different cases."""
        assert canonical_family_label("Buchwald-Hartwig") == "c_n_cross_coupling"
        assert canonical_family_label("buchwald-hartwig") == "c_n_cross_coupling"
        assert canonical_family_label("BUCHWALD-HARTWIG") == "c_n_cross_coupling"

    def test_whitespace_handling(self):
        """Should handle whitespace."""
        assert canonical_family_label("  Suzuki  ") == "suzuki_miyaura"


class TestNormalizeReactionType:
    """Test reaction type normalization."""

    def test_common_reactions(self):
        """Common reactions should normalize."""
        test_cases = [
            ("Buchwald-Hartwig", "c_n_cross_coupling"),
            ("Suzuki coupling", "suzuki_miyaura"),
            ("Sonogashira cross-coupling", "sonogashira"),
            ("SNAr C-N", "snar_cn"),
            ("SN2", "sn2_substitution"),
            ("Amide formation", "amide_coupling"),
        ]

        for label, expected in test_cases:
            assert normalize_reaction_type(label) == expected

    def test_case_insensitivity(self):
        """Normalization should be case insensitive."""
        assert normalize_reaction_type("Buchwald-Hartwig") == "c_n_cross_coupling"
        assert normalize_reaction_type("buchwald-hartwig") == "c_n_cross_coupling"

    def test_empty_string(self):
        """Empty string should return None."""
        assert normalize_reaction_type("") is None

    def test_none_input(self):
        """None input should return None."""
        assert normalize_reaction_type(None) is None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

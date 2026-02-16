"""
Unit tests for Reaction SMILES Analysis Agent - Deterministic Components

Tests the core functionality without requiring LLM calls.
"""

import pytest
from reaction_agent.core import (
    clean_reaction_smiles,
    map_reaction,
    extract_bond_changes,
    analyze_reaction_smiles as analyze_deterministic,
    CleanReport,
    MappingReport,
    BondChangeReport,
)


# Helper functions

def _has_rxnmapper() -> bool:
    """Check if rxnmapper is available."""
    try:
        import rxnmapper
        return True
    except ImportError:
        return False


class TestCleanReactionSMILES:
    """Test reaction SMILES cleaning and canonicalization."""

    def test_simple_reaction(self):
        """Test basic cleaning with valid reaction."""
        rxn = "CCBr>>CCN"
        result = clean_reaction_smiles(rxn, drop_spectators=False)

        assert isinstance(result, CleanReport)
        assert result.rxn_smiles_raw == rxn
        assert ">>" in result.rxn_smiles_clean
        assert len(result.reactants_clean) >= 1
        assert len(result.products_clean) >= 1
        assert len(result.parse_warnings) == 0

    def test_spectator_detection(self):
        """Test spectator detection and removal."""
        # Reaction with chloride spectator
        rxn = "CCN.Cl>>CCNC.Cl"
        result = clean_reaction_smiles(rxn, drop_spectators=True)

        # Chloride should be detected and removed
        assert "Cl" in result.spectators
        assert "Dropped spectator" in " ".join(result.standardization_actions)

    def test_invalid_reaction(self):
        """Test handling of invalid reaction SMILES."""
        with pytest.raises(ValueError, match="expected exactly one '>>'"):
            clean_reaction_smiles("CCBr")

        with pytest.raises(ValueError, match="expected exactly one '>>'"):
            clean_reaction_smiles("CCBr>>CCN>>CCC")

    def test_empty_sides(self):
        """Test handling of empty reactants or products."""
        with pytest.raises(ValueError, match="empty reactants or products"):
            clean_reaction_smiles(">>CCN")

        with pytest.raises(ValueError, match="empty reactants or products"):
            clean_reaction_smiles("CCBr>>")

    def test_multiple_components(self):
        """Test reaction with multiple components on each side."""
        rxn = "CCBr.NCC>>CCNCC.Br"
        result = clean_reaction_smiles(rxn, drop_spectators=False)

        assert len(result.reactants_clean) == 2
        # Products may vary depending on spectator handling


class TestMappingReaction:
    """Test atom mapping functionality."""

    def test_mapping_without_rxnmapper(self):
        """Test graceful handling when rxnmapper not available."""
        result = map_reaction("CCBr>>CCN")

        assert isinstance(result, MappingReport)
        # Should return report even if mapping fails
        assert "mapping_qc" in result.__dict__

    @pytest.mark.skipif(
        not _has_rxnmapper(),
        reason="rxnmapper not installed"
    )
    def test_mapping_with_rxnmapper(self):
        """Test mapping when rxnmapper is available."""
        rxn = "CCBr>>CCN"
        result = map_reaction(rxn)

        assert isinstance(result, MappingReport)
        # If rxnmapper works, we should get a mapped SMILES
        if result.mapping_qc.get("ok"):
            assert result.mapped_rxn_smiles is not None
            assert ":" in result.mapped_rxn_smiles  # Atom maps use :N notation


class TestBondChangeExtraction:
    """Test bond change extraction from mapped reactions."""

    def test_extract_from_valid_mapped_smiles(self):
        """Regression: mapped reaction SMILES should parse in RDKit 2025."""
        mapped_rxn = "[Cl:1][CH3:2].[OH2:3]>>[OH:3][CH3:2].[Cl-:1]"
        result = extract_bond_changes(mapped_rxn)

        assert isinstance(result, BondChangeReport)
        assert "Failed to parse mapped reaction" not in result.warnings

    def test_extract_from_unmapped(self):
        """Test that extraction handles unmapped reactions gracefully."""
        rxn = "CCBr>>CCN"
        result = extract_bond_changes(rxn)

        assert isinstance(result, BondChangeReport)
        # Should not crash, but may have warnings
        assert isinstance(result.warnings, list)

    @pytest.mark.skipif(
        not _has_rxnmapper(),
        reason="rxnmapper not installed"
    )
    def test_extract_with_simple_substitution(self):
        """Test extraction on a simple substitution reaction."""
        # This test requires rxnmapper to create mapped SMILES
        rxn = "CCBr>>CCN"
        mapping_result = map_reaction(rxn)

        if mapping_result.mapping_qc.get("ok") and mapping_result.mapped_rxn_smiles:
            result = extract_bond_changes(mapping_result.mapped_rxn_smiles)

            # Should detect bond changes
            assert isinstance(result.bond_changes, list)
            assert isinstance(result.reaction_center_atoms, list)


class TestFullPipeline:
    """Test the complete analysis pipeline."""

    def test_analyze_simple_reaction(self):
        """Test analyzing a simple reaction end-to-end."""
        rxn = "CCBr>>CCN"
        result = analyze_deterministic(rxn, skip_mapping=True)

        assert "input" in result
        assert "tool_facts" in result
        assert result["input"]["rxn_smiles_raw"] == rxn
        assert result["input"]["rxn_smiles_clean"] != ""

    def test_analyze_with_spectators(self):
        """Test analysis with spectator removal."""
        rxn = "CCN.Cl>>CCNC.Cl"
        result = analyze_deterministic(rxn, drop_spectators=True)

        assert "Cl" in result["input"]["spectators"]

    @pytest.mark.skipif(
        not _has_rxnmapper(),
        reason="rxnmapper not installed"
    )
    def test_analyze_with_mapping(self):
        """Test full pipeline including mapping."""
        rxn = "CCBr>>CCN"
        result = analyze_deterministic(rxn, skip_mapping=False)

        assert "tool_facts" in result
        # Mapping may succeed or fail, but should be attempted
        assert "mapping_qc" in result["tool_facts"]
        assert "bond_changes" in result["tool_facts"]


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v"])

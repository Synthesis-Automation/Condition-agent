"""
Integration tests for SMARTS generator CLI using refactored utilities.

Tests that the CLI correctly uses SubstrateClassifier and SmartsBuilder
for chemistry-aware pattern generation.
"""

import pytest
from chemtools.protocol.smarts_generator_cli import SmartsGenerator, ReactionSmartsApplicability


class TestSmartsGeneratorIntegration:
    """Test integration with SubstrateClassifier and SmartsBuilder"""
    
    def test_primary_alkyl_iodide_borylation(self):
        """Test octyl iodide borylation uses chemistry-aware patterns"""
        reaction = "CCCCCCCCI.B1OC(C)(C)C(C)(C)O1>>CCCCCCCCB1OC(C)(C)C(C)(C)O1"
        
        gen = SmartsGenerator(reaction)
        core = gen.generate_core_smarts()
        guards = gen.suggest_guard_patterns()
        
        # Should generate primary alkyl pattern
        assert "[CX4;H2,H3]-[I]" in core
        
        # Should have context-aware guards
        assert len(guards) >= 3
        # Check for secondary exclusion
        assert any("secondary" in g.lower() or "[CX4;H1]-[I]" in g for g in guards)
        # Check for tertiary exclusion
        assert any("tertiary" in g.lower() or "[CX4;H0]-[I]" in g for g in guards)
        # Check for benzylic exclusion
        assert any("benzylic" in g.lower() or "[CH2" in g for g in guards)
    
    def test_aryl_bromide_buchwald_hartwig(self):
        """Test Buchwald-Hartwig uses aryl-specific patterns"""
        reaction = "c1ccccc1Br.c1ccccc1N>>c1ccccc1Nc2ccccc2"
        
        gen = SmartsGenerator(reaction)
        core = gen.generate_core_smarts()
        guards = gen.suggest_guard_patterns()
        
        # Should generate aryl patterns
        assert "c-[Br]" in core
        assert "c-[N" in core or "[NX3" in core
        
        # Should exclude aliphatic versions
        assert any("aliphatic" in g.lower() or "[CX4]" in g for g in guards)
    
    def test_benzylic_chloride_detection(self):
        """Test benzylic chloride gets special pattern"""
        reaction = "c1ccccc1CCl.CCN>>c1ccccc1CNC"
        
        gen = SmartsGenerator(reaction)
        core = gen.generate_core_smarts()
        
        # Should generate benzylic-specific pattern
        assert "[CH2;$([CH2][c])]-[Cl]" in core
    
    def test_secondary_alkyl_bromide(self):
        """Test secondary alkyl bromide pattern"""
        reaction = "CC(C)Br.CCN>>CC(C)NC"
        
        gen = SmartsGenerator(reaction)
        core = gen.generate_core_smarts()
        
        # Should generate secondary alkyl pattern
        assert "[CX4;H1]-[Br]" in core
    
    def test_aniline_vs_aliphatic_amine(self):
        """Test aniline pattern vs aliphatic amine are distinguished"""
        # Test that aniline and aliphatic amines get different patterns
        # Focus on reactant patterns (first molecule in core SMARTS)
        
        # Aniline as first reactant
        reaction1 = "c1ccccc1N.CCCCI>>C1CCCCC1"
        gen1 = SmartsGenerator(reaction1)
        core1 = gen1.generate_core_smarts()
        
        # First part of reaction should use aromatic amine pattern
        assert "c-[NX3" in core1 or "c-[N" in core1
        
        # Aliphatic amine as first reactant
        reaction2 = "CCN.CCCCI>>CCCNCC"
        gen2 = SmartsGenerator(reaction2)
        core2 = gen2.generate_core_smarts()
        
        # First part should use aliphatic amine pattern
        assert "[CX4]-[NX3" in core2 or "[CX4]-[N" in core2
    
    def test_amide_exclusion_in_amine_patterns(self):
        """Test that amine patterns exclude amides"""
        reaction = "c1ccccc1Br.c1ccccc1N>>c1ccccc1Nc2ccccc2"
        
        gen = SmartsGenerator(reaction)
        core = gen.generate_core_smarts()
        
        # Amine pattern should include amide exclusion
        assert "!$(NC=O)" in core or "!$(N" in core
    
    def test_guard_pattern_context_awareness(self):
        """Test that guard patterns are context-aware"""
        # Primary alkyl should exclude secondary, tertiary, benzylic, allylic
        gen1 = SmartsGenerator("CCCCI>>CCCB")
        guards1 = gen1.suggest_guard_patterns()
        
        # Should have multiple context-specific guards
        assert len(guards1) >= 3
        
        # Aryl should exclude aliphatic
        gen2 = SmartsGenerator("c1ccccc1Br>>c1ccccc1B")
        guards2 = gen2.suggest_guard_patterns()
        
        assert len(guards2) >= 1
        assert any("[CX4]" in g for g in guards2)  # Exclude aliphatic
    
    def test_reactants_and_products_parsing(self):
        """Test reaction parsing"""
        reaction = "CCCCI.B>>CCCB.I"
        
        gen = SmartsGenerator(reaction)
        
        assert "CCCCI" in gen.reactants_smiles
        assert "CCCB" in gen.products_smiles
    
    def test_multiple_reactants(self):
        """Test reactions with multiple reactants"""
        reaction = "CCCCI.CCN.B>>CCCB.CCCN.I"
        
        gen = SmartsGenerator(reaction)
        guards = gen.suggest_guard_patterns()
        
        # Should analyze all reactants
        assert len(guards) > 0  # At least some guards from the substrates


class TestReactionSmartsApplicability:
    """Test the data model"""
    
    def test_to_dict(self):
        """Test conversion to dict"""
        app = ReactionSmartsApplicability(
            core="[C:1]-[I:2]>>[C:1]-[B:3]",
            guards_forbid=["[CX4;H0]-[I]"],
            notes="Test pattern"
        )
        
        d = app.to_dict()
        
        assert d["core"] == "[C:1]-[I:2]>>[C:1]-[B:3]"
        assert d["guards_forbid"] == ["[CX4;H0]-[I]"]
        assert d["notes"] == "Test pattern"
    
    def test_from_dict(self):
        """Test creation from dict"""
        data = {
            "core": "[C:1]-[I:2]>>[C:1]-[B:3]",
            "guards_forbid": ["[CX4;H0]-[I]"],
            "notes": "Test"
        }
        
        app = ReactionSmartsApplicability.from_dict(data)
        
        assert app.core == "[C:1]-[I:2]>>[C:1]-[B:3]"
        assert app.guards_forbid == ["[CX4;H0]-[I]"]
        assert app.notes == "Test"


class TestPatternQuality:
    """Test that generated patterns have correct properties"""
    
    def test_patterns_are_specific_not_generic(self):
        """Test that patterns reflect actual chemistry, not just atom counts"""
        # Primary alkyl iodide
        gen1 = SmartsGenerator("CCCCCCCCI>>CCCCCCCCB")
        core1 = gen1.generate_core_smarts()
        
        # Should NOT be generic [C]-[I]
        # Should be chemistry-aware [CX4;H2,H3]-[I]
        assert "H2,H3" in core1 or "CX4" in core1
        
        # Aryl bromide
        gen2 = SmartsGenerator("c1ccccc1Br>>c1ccccc1B")
        core2 = gen2.generate_core_smarts()
        
        # Should use aromatic 'c', not generic [C]
        assert "c-[Br]" in core2 or "c-[B" in core2
    
    def test_patterns_match_substrate_type(self):
        """Test that patterns match the correct substrate type"""
        examples = [
            ("CCCCI", "[CX4;H2,H3]-[I]"),  # Primary alkyl
            ("CC(C)I", "[CX4;H1]-[I]"),  # Secondary alkyl
            ("c1ccccc1I", "c-[I]"),  # Aryl
            ("c1ccccc1CI", "[CH2;$([CH2][c])]-[I]"),  # Benzylic
        ]
        
        for smiles, expected_pattern in examples:
            gen = SmartsGenerator(f"{smiles}>>{smiles.replace('I', 'B')}")
            core = gen.generate_core_smarts()
            assert expected_pattern in core, f"Expected {expected_pattern} for {smiles}, got {core}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

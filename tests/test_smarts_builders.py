"""
Tests for smarts_builders module.

Test context-aware SMARTS pattern generation including:
- Halide patterns (alkyl, aryl, benzylic, allylic)
- Amine patterns (aniline, aliphatic)
- Amide patterns (distinguished from amines)
- Alcohol patterns (phenol, benzylic, aliphatic)
- Carbonyl patterns (acid, ester, aldehyde, ketone)
- Boron patterns (boronic acid, pinacol ester)
- Guard pattern generation
- Pattern matching
"""

import pytest
from chemtools.util.smarts_builders import (
    SmartsBuilder,
    SmartsPatternMatcher,
    build_smarts,
    build_smarts_with_guards,
    match_smarts,
)
from chemtools.util.rdkit_helpers import rdkit_available


class TestSmartsBuilder:
    """Test SmartsBuilder class"""
    
    def setup_method(self):
        """Setup test fixtures"""
        self.builder = SmartsBuilder()
    
    def test_init(self):
        """Test builder initialization"""
        assert self.builder is not None
        assert hasattr(self.builder, 'build_from_smiles')
    
    def test_build_from_smiles(self):
        """Test building SMARTS from SMILES"""
        smarts = self.builder.build_from_smiles("CCCCI")
        assert smarts is not None
        assert 'I' in smarts or 'X' in smarts


class TestHalidePatterns:
    """Test halide SMARTS pattern generation"""
    
    def setup_method(self):
        self.builder = SmartsBuilder()
    
    def test_primary_alkyl_iodide(self):
        """Test primary alkyl iodide pattern"""
        smarts = self.builder.build_from_smiles("CCCCI")
        assert "[CX4;H2,H3]-[I]" == smarts
    
    def test_methyl_iodide(self):
        """Test methyl iodide (simplest primary)"""
        smarts = self.builder.build_from_smiles("CI")
        assert "[CX4;H2,H3]-[I]" == smarts
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_secondary_alkyl_bromide(self):
        """Test secondary alkyl bromide pattern"""
        smarts = self.builder.build_from_smiles("CC(C)Br")
        assert "[CX4;H1]-[Br]" == smarts
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_tertiary_alkyl_chloride(self):
        """Test tertiary alkyl chloride pattern"""
        smarts = self.builder.build_from_smiles("CC(C)(C)Cl")
        assert "[CX4;H0]-[Cl]" == smarts
    
    def test_aryl_iodide(self):
        """Test aryl iodide pattern"""
        smarts = self.builder.build_from_smiles("c1ccccc1I")
        assert "c-[I]" == smarts
    
    def test_aryl_bromide(self):
        """Test aryl bromide pattern"""
        smarts = self.builder.build_from_smiles("c1ccccc1Br")
        assert "c-[Br]" == smarts
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_heteroaryl_chloride(self):
        """Test heteroaryl chloride pattern"""
        smarts = self.builder.build_from_smiles("c1cnc(Cl)cc1")
        # Should use generic aromatic 'a' for heteroaromatic
        assert "[a]-[Cl]" == smarts or "c-[Cl]" == smarts
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_benzyl_chloride(self):
        """Test benzyl chloride (benzylic) pattern"""
        smarts = self.builder.build_from_smiles("c1ccccc1CCl")
        assert "[CH2;$([CH2][c])]-[Cl]" == smarts
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_allyl_bromide(self):
        """Test allyl bromide (allylic) pattern"""
        smarts = self.builder.build_from_smiles("C=CCBr")
        assert "[CH2;$([CH2]C=C)]-[Br]" == smarts
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_propargyl_iodide(self):
        """Test propargyl iodide (propargylic) pattern"""
        smarts = self.builder.build_from_smiles("C#CCI")
        assert "[CH2;$([CH2]C#C)]-[I]" == smarts


class TestAminePatterns:
    """Test amine SMARTS pattern generation"""
    
    def setup_method(self):
        self.builder = SmartsBuilder()
    
    def test_aniline(self):
        """Test aniline pattern"""
        smarts = self.builder.build_from_smiles("c1ccccc1N")
        assert "c-[NX3;H2;!$(NC=O)]" == smarts
    
    def test_primary_amine(self):
        """Test primary aliphatic amine pattern"""
        smarts = self.builder.build_from_smiles("CCN")
        assert "[CX4]-[NX3;H2;!$(NC=O)]" == smarts
    
    def test_secondary_amine(self):
        """Test secondary amine pattern"""
        smarts = self.builder.build_from_smiles("CCNC")
        # Should detect secondary amine
        assert "[NX3;H1;!$(NC=O)]" == smarts


class TestAmidePatterns:
    """Test amide SMARTS pattern generation (distinct from amines)"""
    
    def setup_method(self):
        self.builder = SmartsBuilder()
    
    def test_primary_amide(self):
        """Test primary amide pattern"""
        smarts = self.builder.build_from_smiles("CC(=O)N")
        assert "[NX3;H2]-[CX3](=O)" == smarts or "[NX3" in smarts
    
    def test_secondary_amide(self):
        """Test secondary amide pattern"""
        smarts = self.builder.build_from_smiles("CC(=O)NC")
        assert "[NX3;H1]-[CX3](=O)" == smarts or "[NX3" in smarts


class TestAlcoholPatterns:
    """Test alcohol SMARTS pattern generation"""
    
    def setup_method(self):
        self.builder = SmartsBuilder()
    
    def test_phenol(self):
        """Test phenol pattern"""
        smarts = self.builder.build_from_smiles("c1ccccc1O")
        assert "c-[OX2H]" == smarts
    
    def test_aliphatic_alcohol(self):
        """Test aliphatic alcohol pattern"""
        smarts = self.builder.build_from_smiles("CCCO")
        assert "[CX4]-[OX2H]" == smarts
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_benzylic_alcohol(self):
        """Test benzylic alcohol pattern"""
        smarts = self.builder.build_from_smiles("c1ccccc1CO")
        assert "[CH2;$([CH2][c])]-[OX2H]" == smarts
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_allylic_alcohol(self):
        """Test allylic alcohol pattern"""
        smarts = self.builder.build_from_smiles("C=CCO")
        assert "[CH2;$([CH2]C=C)]-[OX2H]" == smarts


class TestCarbonylPatterns:
    """Test carbonyl SMARTS pattern generation"""
    
    def setup_method(self):
        self.builder = SmartsBuilder()
    
    def test_carboxylic_acid(self):
        """Test carboxylic acid pattern"""
        smarts = self.builder.build_from_smiles("CC(=O)O")
        assert "[CX3](=O)-[OX2H]" == smarts
    
    def test_ester(self):
        """Test ester pattern"""
        smarts = self.builder.build_from_smiles("CC(=O)OC")
        assert "[CX3](=O)-[OX2]-[C]" == smarts
    
    def test_aldehyde(self):
        """Test aldehyde pattern"""
        smarts = self.builder.build_from_smiles("CCC=O")
        assert "[CX3;H1](=O)" == smarts
    
    def test_ketone(self):
        """Test ketone pattern"""
        smarts = self.builder.build_from_smiles("CC(=O)C")
        assert "[C]-[CX3](=O)-[C]" == smarts


class TestBoronPatterns:
    """Test boron compound SMARTS pattern generation"""
    
    def setup_method(self):
        self.builder = SmartsBuilder()
    
    def test_boronic_acid(self):
        """Test boronic acid pattern"""
        smarts = self.builder.build_from_smiles("c1ccccc1B(O)O")
        assert "[B]([OH])([OH])" == smarts
    
    def test_boronic_ester_pinacol(self):
        """Test boronic ester pinacol pattern"""
        smarts = self.builder.build_from_smiles("c1ccccc1B1OC(C)(C)C(C)(C)O1")
        assert "[B]1OC(C)(C)C(C)(C)O1" == smarts


class TestGuardPatterns:
    """Test guard pattern generation"""
    
    def setup_method(self):
        self.builder = SmartsBuilder()
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_primary_alkyl_iodide_guards(self):
        """Test guard patterns for primary alkyl iodide"""
        from chemtools.util.substrate_classifier import classify_substrate
        
        info = classify_substrate("CCCCI")
        guards = self.builder.build_guard_patterns(info)
        
        # Should exclude secondary, tertiary, benzylic, allylic
        assert len(guards) >= 3
        assert any('[H1]' in g or 'secondary' in g.lower() for g in guards)
        assert any('[H0]' in g or 'tertiary' in g.lower() for g in guards)
        assert any('benzylic' in g.lower() for g in guards)
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_aryl_halide_guards(self):
        """Test guard patterns for aryl halide"""
        from chemtools.util.substrate_classifier import classify_substrate
        
        info = classify_substrate("c1ccccc1Br")
        guards = self.builder.build_guard_patterns(info)
        
        # Should exclude aliphatic halides
        assert len(guards) >= 1
        assert any('CX4' in g or 'aliphatic' in g.lower() for g in guards)
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_aniline_guards(self):
        """Test guard patterns for aniline"""
        from chemtools.util.substrate_classifier import classify_substrate
        
        info = classify_substrate("c1ccccc1N")
        guards = self.builder.build_guard_patterns(info)
        
        # Should exclude aliphatic amines
        assert len(guards) >= 1
        assert any('CX4' in g or 'aliphatic' in g.lower() for g in guards)


class TestSmartsPatternMatcher:
    """Test SMARTS pattern matching"""
    
    def setup_method(self):
        self.matcher = SmartsPatternMatcher()
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_match_primary_alkyl_iodide(self):
        """Test matching primary alkyl iodides"""
        pattern = "[CX4;H2,H3]-[I]"
        
        # Should match
        assert self.matcher.match("CI", pattern)
        assert self.matcher.match("CCI", pattern)
        assert self.matcher.match("CCCCI", pattern)
        
        # Should not match secondary
        assert not self.matcher.match("CC(C)I", pattern)
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_match_aryl_bromide(self):
        """Test matching aryl bromides"""
        pattern = "c-[Br]"
        
        # Should match
        assert self.matcher.match("c1ccccc1Br", pattern)
        assert self.matcher.match("c1ccc(Br)cc1", pattern)
        
        # Should not match aliphatic
        assert not self.matcher.match("CCBr", pattern)
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_match_aniline(self):
        """Test matching aniline"""
        pattern = "c-[NX3;H2;!$(NC=O)]"
        
        # Should match
        assert self.matcher.match("c1ccccc1N", pattern)
        
        # Should not match aliphatic amine
        assert not self.matcher.match("CCN", pattern)
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_explain_match(self):
        """Test match explanation"""
        pattern = "[CX4;H2,H3]-[I]"
        
        explanation = self.matcher.explain_match("CCCCI", pattern)
        assert "match" in explanation.lower()
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_find_matching_atoms(self):
        """Test finding matching atom indices"""
        pattern = "[I]"
        
        atoms = self.matcher.find_matching_atoms("CCCCI", pattern)
        assert len(atoms) >= 1


class TestConvenienceFunctions:
    """Test convenience functions"""
    
    def test_build_smarts(self):
        """Test build_smarts convenience function"""
        smarts = build_smarts("CCCCI")
        assert smarts == "[CX4;H2,H3]-[I]"
    
    def test_build_smarts_with_guards(self):
        """Test build_smarts_with_guards"""
        result = build_smarts_with_guards("CCCCI")
        
        assert 'core' in result
        assert 'guards_forbid' in result
        assert 'substrate_class' in result
        assert result['core'] == "[CX4;H2,H3]-[I]"
        assert len(result['guards_forbid']) >= 3
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_match_smarts(self):
        """Test match_smarts convenience function"""
        assert match_smarts("CCCCI", "[CX4;H2,H3]-[I]")
        assert not match_smarts("CC(C)I", "[CX4;H2,H3]-[I]")


class TestRealWorldExamples:
    """Test with real-world reaction substrates"""
    
    def setup_method(self):
        self.builder = SmartsBuilder()
    
    def test_octyl_iodide(self):
        """Test octyl iodide (from protocol example)"""
        smarts = self.builder.build_from_smiles("CCCCCCCCI")
        assert smarts == "[CX4;H2,H3]-[I]"
    
    def test_phenyl_bromide(self):
        """Test phenyl bromide"""
        smarts = self.builder.build_from_smiles("c1ccccc1Br")
        assert smarts == "c-[Br]"
    
    def test_4_bromobiphenyl(self):
        """Test 4-bromobiphenyl"""
        smarts = self.builder.build_from_smiles("c1ccc(Br)cc1-c2ccccc2")
        assert smarts == "c-[Br]"
    
    def test_benzylamine(self):
        """Test benzylamine"""
        smarts = self.builder.build_from_smiles("c1ccccc1CN")
        # Primary amine
        assert "[NX3;H2" in smarts or "[CX4]-[NX3" in smarts
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_phenylboronic_acid(self):
        """Test phenylboronic acid"""
        smarts = self.builder.build_from_smiles("c1ccccc1B(O)O")
        assert smarts == "[B]([OH])([OH])"


class TestEdgeCases:
    """Test edge cases and error handling"""
    
    def setup_method(self):
        self.builder = SmartsBuilder()
    
    def test_empty_smiles(self):
        """Test with empty SMILES"""
        smarts = self.builder.build_from_smiles("")
        # Should not crash, return generic pattern
        assert smarts is not None
    
    def test_invalid_smiles(self):
        """Test with invalid SMILES"""
        smarts = self.builder.build_from_smiles("INVALID")
        # Should not crash
        assert smarts is not None


class TestPatternConsistency:
    """Test that patterns are consistent and chemically meaningful"""
    
    def setup_method(self):
        self.builder = SmartsBuilder()
        self.matcher = SmartsPatternMatcher()
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_primary_pattern_matches_primaries(self):
        """Test that primary pattern matches all primary halides"""
        pattern = "[CX4;H2,H3]-[I]"
        
        primaries = ["CI", "CCI", "CCCCI", "CCCCCCCCI", "ICCCCCCCCCCCC"]
        for smiles in primaries:
            assert self.matcher.match(smiles, pattern), f"{smiles} should match primary pattern"
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_primary_pattern_rejects_secondary(self):
        """Test that primary pattern with negative constraints rejects secondary halides
        
        Note: The pattern [CX4;H2,H3]-[I] alone will match CC(C)CI because it HAS
        a primary CH2-I group. We need guard patterns to exclude the full molecule.
        This test verifies secondary halides match the pattern (they have primary carbons)
        but would be excluded by guard patterns in real usage.
        """
        pattern = "[CX4;H2,H3]-[I]"
        
        # These molecules DO have primary CH2 groups bonded to I
        # (That's chemically correct - they contain primary carbons)
        secondary_with_primary_I = ["CC(C)CI"]  # has primary CH2-I 
        for smiles in secondary_with_primary_I:
            # This should match because there IS a primary carbon bonded to I
            assert self.matcher.match(smiles, pattern), f"{smiles} HAS primary C-I"
        
        # Pure secondary iodides (no primary carbons)
        pure_secondary = ["CC(C)I", "CCC(C)I"]  # C bonded to I has only 1 H
        for smiles in pure_secondary:
            assert not self.matcher.match(smiles, pattern), f"{smiles} should NOT match primary pattern"
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_aryl_pattern_matches_aryls(self):
        """Test that aryl pattern matches aryl halides"""
        pattern = "c-[Br]"
        
        aryls = ["c1ccccc1Br", "c1ccc(Br)cc1", "c1ccc(Br)cc1C"]
        for smiles in aryls:
            assert self.matcher.match(smiles, pattern), f"{smiles} should match aryl pattern"
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_aryl_pattern_rejects_alkyl(self):
        """Test that aryl pattern rejects alkyl halides"""
        pattern = "c-[Br]"
        
        alkyls = ["CCBr", "CC(C)Br", "CCCCBr"]
        for smiles in alkyls:
            assert not self.matcher.match(smiles, pattern), f"{smiles} should NOT match aryl pattern"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

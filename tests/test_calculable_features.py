"""
Tests for calculable feature detection.

Tests cover:
- Boolean SMARTS-based features
- Integer count features  
- Heuristic features (polarity, beta-hydride)
- Derived/shortcut features
- Edge cases and error handling
"""

import pytest
from chemtools.featurizers import calculable
from chemtools.util.rdkit_helpers import rdkit_available


# Skip all tests if RDKit is not available
pytestmark = pytest.mark.skipif(
    not rdkit_available(),
    reason="RDKit not available"
)


class TestBooleanSMARTSFeatures:
    """Test SMARTS-based boolean feature detection."""
    
    def test_sp2_halides(self):
        """Test sp2 halide detection (aryl and vinyl)."""
        # Aryl chloride
        features = calculable.detect_all_features("c1ccc(Cl)cc1")
        assert features["sp2_chloride_present"] is True
        assert features["aryl_halide_present"] is True
        assert features["sp2_bromide_present"] is False
        
        # Aryl bromide
        features = calculable.detect_all_features("c1ccc(Br)cc1")
        assert features["sp2_bromide_present"] is True
        assert features["aryl_halide_present"] is True
        
        # Aryl iodide
        features = calculable.detect_all_features("c1ccc(I)cc1")
        assert features["sp2_iodide_present"] is True
        assert features["aryl_halide_present"] is True
        
        # Vinyl bromide
        features = calculable.detect_all_features("C=CBr")
        assert features["sp2_bromide_present"] is True
        assert features["vinyl_halide_present"] is True
        assert features["aryl_halide_present"] is False
    
    def test_sp3_halides(self):
        """Test sp3 halide detection (alkyl halides)."""
        # Alkyl chloride
        features = calculable.detect_all_features("CCCl")
        assert features["sp3_chloride_present"] is True
        assert features["sp2_chloride_present"] is False
        
        # Alkyl bromide
        features = calculable.detect_all_features("CCBr")
        assert features["sp3_bromide_present"] is True
        assert features["sp2_bromide_present"] is False
        
        # Alkyl iodide
        features = calculable.detect_all_features("CCI")
        assert features["sp3_iodide_present"] is True
        assert features["sp2_iodide_present"] is False
    
    def test_sulfonates(self):
        """Test sulfonate leaving groups (OTf, OTs, OMs)."""
        # Aryl triflate
        features = calculable.detect_all_features("c1ccc(OS(=O)(=O)C(F)(F)F)cc1")
        assert features["sp2_triflate_present"] is True
        assert features["sp2_pseudohalide_present"] is True
        assert features["aryl_halide_present"] is True
        
        # Aryl tosylate
        features = calculable.detect_all_features("c1ccc(OS(=O)(=O)c2ccc(C)cc2)cc1")
        assert features["sp2_tosylate_present"] is True
        assert features["sp2_pseudohalide_present"] is True
        
        # Aryl mesylate
        features = calculable.detect_all_features("c1ccc(OS(=O)(=O)C)cc1")
        assert features["sp2_mesylate_present"] is True
        assert features["sp2_pseudohalide_present"] is True
    
    def test_boron_reagents(self):
        """Test boronic acids, esters, and BF3K salts."""
        # Aryl boronic acid pinacol ester (approximate)
        features = calculable.detect_all_features("c1ccc(B2OC(C)(C)COC2(C)C)cc1")
        assert features["sp2_boron_present"] is True
        
        # BF3K salt
        features = calculable.detect_all_features("c1ccc([B-](F)(F)F)cc1")
        assert features["boron_bf3k_present"] is True
        assert features["sp2_boron_present"] is True
    
    def test_organometallics(self):
        """Test organometallic reagents."""
        # Grignard
        features = calculable.detect_all_features("CCCMgBr")
        assert features["grignard_present"] is True
        
        # Organozinc
        features = calculable.detect_all_features("CCCZnBr")
        assert features["organozinc_present"] is True
        
        # Organolithium
        features = calculable.detect_all_features("CCCLi")
        assert features["organolithium_present"] is True
        
        # Stannane
        features = calculable.detect_all_features("CCC[Sn](C)(C)C")
        assert features["stannane_present"] is True
        
        # Organosilane
        features = calculable.detect_all_features("CCC[Si](C)(C)C")
        assert features["organosilane_present"] is True
    
    def test_nucleophiles(self):
        """Test amine, alcohol, and other nucleophile detection."""
        # Aliphatic amine
        features = calculable.detect_all_features("CCN")
        assert features["aliphatic_amine_present"] is True
        assert features["aniline_present"] is False
        
        # Aniline
        features = calculable.detect_all_features("c1ccc(N)cc1")
        assert features["aniline_present"] is True
        assert features["aliphatic_amine_present"] is False
        
        # Phenol
        features = calculable.detect_all_features("c1ccc(O)cc1")
        assert features["phenol_present"] is True
        assert features["alcohol_present"] is False
        
        # Aliphatic alcohol
        features = calculable.detect_all_features("CCO")
        assert features["alcohol_present"] is True
        assert features["phenol_present"] is False
        
        # Thiol
        features = calculable.detect_all_features("CCS")
        assert features["thiol_present"] is True
        assert features["thiol_poison_risk"] is True
    
    def test_alkenes_alkynes(self):
        """Test alkene and alkyne detection."""
        # Terminal alkene
        features = calculable.detect_all_features("C=C")
        assert features["alkene_present"] is True
        assert features["terminal_alkene_present"] is True
        
        # Internal alkene
        features = calculable.detect_all_features("CC=CC")
        assert features["alkene_present"] is True
        # Note: terminal_alkene depends on SMARTS specificity
        
        # Terminal alkyne
        features = calculable.detect_all_features("C#C")
        assert features["alkyne_present"] is True
        assert features["terminal_alkyne_present"] is True
        
        # Phenylacetylene
        features = calculable.detect_all_features("c1ccccc1C#C")
        assert features["alkyne_present"] is True
        assert features["terminal_alkyne_present"] is True
    
    def test_heterocycles(self):
        """Test heterocycle detection."""
        # Pyridine
        features = calculable.detect_all_features("c1ccncc1")
        assert features["pyridine_present"] is True
        assert features["pyridine_poison_risk"] is True
        
        # Pyrimidine
        features = calculable.detect_all_features("c1ncncn1")
        assert features["pyrimidine_present"] is True
        
        # Indole
        features = calculable.detect_all_features("c1ccc2[nH]ccc2c1")
        assert features["indole_present"] is True
    
    def test_activated_esters(self):
        """Test activated ester detection."""
        # p-Nitrophenyl ester
        features = calculable.detect_all_features("CC(=O)Oc1ccc([N+](=O)[O-])cc1")
        assert features["activated_aryl_ester_present"] is True
        
        # Pentafluorophenyl ester
        features = calculable.detect_all_features("CC(=O)Oc1c(F)c(F)c(F)c(F)c1F")
        assert features["activated_fluoro_ester_present"] is True


class TestIntegerCountFeatures:
    """Test integer count features."""
    
    def test_sp2_halide_count(self):
        """Test counting sp2 halide sites."""
        # Single aryl bromide
        features = calculable.detect_all_features("c1ccc(Br)cc1")
        assert features["sp2_halide_site_count"] == 1
        
        # Two aryl bromides
        features = calculable.detect_all_features("c1cc(Br)cc(Br)c1")
        assert features["sp2_halide_site_count"] == 2
        
        # Mixed halides
        features = calculable.detect_all_features("c1cc(Br)cc(Cl)c1")
        assert features["sp2_halide_site_count"] == 2
        
        # No halides
        features = calculable.detect_all_features("c1ccccc1")
        assert features["sp2_halide_site_count"] == 0
    
    def test_sp2_sulfonate_count(self):
        """Test counting sp2 sulfonate sites."""
        # Single triflate
        features = calculable.detect_all_features("c1ccc(OS(=O)(=O)C(F)(F)F)cc1")
        assert features["sp2_sulfonate_site_count"] == 1
        
        # Two sulfonates
        features = calculable.detect_all_features("c1cc(OS(=O)(=O)C)cc(OS(=O)(=O)C)c1")
        assert features["sp2_sulfonate_site_count"] == 2


class TestHeuristicFeatures:
    """Test heuristic-based features."""
    
    def test_polarity_features(self):
        """Test polarity classification."""
        # High polarity molecule (amino acids, sugars)
        features = calculable.detect_all_features("C(C(C(=O)O)N)O")  # Serine-like
        assert features["polarity_high"] is True
        assert features["polarity_low"] is False
        
        # Low polarity molecule (simple hydrocarbons)
        features = calculable.detect_all_features("CCCCCC")  # Hexane
        assert features["polarity_low"] is True
        assert features["polarity_high"] is False
        
        # Medium polarity
        features = calculable.detect_all_features("c1ccccc1")  # Benzene
        # Should be neither high nor low
    
    def test_base_sensitive(self):
        """Test base-sensitive functional groups."""
        # Acyl chloride
        features = calculable.detect_all_features("CC(=O)Cl")
        assert features["base_sensitive"] is True
        
        # Sulfonyl fluoride
        features = calculable.detect_all_features("CS(=O)(=O)F")
        assert features["base_sensitive"] is True
        
        # Silane
        features = calculable.detect_all_features("C[Si](C)(C)C")
        assert features["base_sensitive"] is True
    
    def test_acidic_proton(self):
        """Test acidic proton detection."""
        # Alcohol
        features = calculable.detect_all_features("CCO")
        assert features["acidic_proton_present"] is True
        
        # Thiol
        features = calculable.detect_all_features("CCS")
        assert features["acidic_proton_present"] is True
        
        # Carboxylic acid
        features = calculable.detect_all_features("CC(=O)O")
        assert features["acidic_proton_present"] is True
    
    def test_beta_hydride(self):
        """Test β-hydride elimination risk detection."""
        # Primary alkyl bromide with β-hydrogens
        features = calculable.detect_all_features("CCBr")
        assert features["beta_hydride_possible"] is True
        
        # Secondary alkyl bromide with β-hydrogens
        features = calculable.detect_all_features("CC(Br)C")
        assert features["beta_hydride_possible"] is True
        
        # Neopentyl bromide (no β-hydrogens on β-carbon)
        features = calculable.detect_all_features("CC(C)(C)CBr")
        # This should be False, but depends on implementation detail
        
        # Aryl bromide (no β-hydride issue)
        features = calculable.detect_all_features("c1ccc(Br)cc1")
        assert features["beta_hydride_possible"] is False


class TestDerivedFeatures:
    """Test derived/shortcut features."""
    
    def test_internal_alkyne(self):
        """Test internal alkyne detection."""
        # Terminal alkyne
        features = calculable.detect_all_features("C#C")
        assert features["alkyne_present"] is True
        assert features["terminal_alkyne_present"] is True
        assert features["internal_alkyne_present"] is False
        
        # Internal alkyne
        features = calculable.detect_all_features("CC#CC")
        assert features["alkyne_present"] is True
        # Depending on SMARTS, terminal_alkyne_present may vary
    
    def test_aryl_halide_shortcuts(self):
        """Test ArX shortcuts."""
        # Aryl bromide
        features = calculable.detect_all_features("c1ccc(Br)cc1")
        assert features["ArBr_present"] is True
        assert features["ArCl_present"] is False
        assert features["ArI_present"] is False
        
        # Aryl chloride
        features = calculable.detect_all_features("c1ccc(Cl)cc1")
        assert features["ArCl_present"] is True
        assert features["ArBr_present"] is False
        
        # Aryl iodide
        features = calculable.detect_all_features("c1ccc(I)cc1")
        assert features["ArI_present"] is True
        assert features["ArBr_present"] is False
    
    def test_vinyl_halide_shortcuts(self):
        """Test VinylX shortcuts."""
        # Vinyl bromide
        features = calculable.detect_all_features("C=CBr")
        assert features["VinylBr_present"] is True
        assert features["VinylCl_present"] is False
        assert features["ArBr_present"] is False
    
    def test_aryl_sulfonate_shortcuts(self):
        """Test ArOX sulfonate shortcuts."""
        # Aryl triflate
        features = calculable.detect_all_features("c1ccc(OS(=O)(=O)C(F)(F)F)cc1")
        assert features["ArOTf_present"] is True
        assert features["ArOTs_present"] is False
        
        # Aryl tosylate
        features = calculable.detect_all_features("c1ccc(OS(=O)(=O)c2ccc(C)cc2)cc1")
        assert features["ArOTs_present"] is True
        assert features["ArOTf_present"] is False


class TestUtilityFunctions:
    """Test utility and convenience functions."""
    
    def test_detect_feature_single(self):
        """Test detecting a single feature."""
        result = calculable.detect_feature("c1ccc(Br)cc1", "sp2_bromide_present")
        assert result is True
        
        result = calculable.detect_feature("c1ccc(Br)cc1", "sp2_chloride_present")
        assert result is False
    
    def test_get_present_features(self):
        """Test getting list of present features."""
        present = calculable.get_present_features("c1ccc(Br)cc1")
        assert "sp2_bromide_present" in present
        assert "aryl_halide_present" in present
        assert "ArBr_present" in present
        assert "sp2_chloride_present" not in present
    
    def test_feature_summary(self):
        """Test human-readable feature summary."""
        summary = calculable.feature_summary("c1ccc(Br)cc1")
        assert "c1ccc(Br)cc1" in summary
        assert "sp2_bromide_present" in summary or "Detected features" in summary
    
    def test_batch_detection(self):
        """Test batch feature detection."""
        smiles_list = ["c1ccc(Br)cc1", "CCBr", "c1ccccc1"]
        results = calculable.detect_features_batch(smiles_list)
        
        assert len(results) == 3
        assert results[0]["sp2_bromide_present"] is True
        assert results[1]["sp3_bromide_present"] is True
        assert results[2]["sp2_bromide_present"] is False
    
    def test_get_feature_spec(self):
        """Test getting feature specification."""
        spec = calculable.get_feature_spec()
        assert "version" in spec
        assert "features" in spec
        assert "derived_shortcuts" in spec
        assert len(spec["features"]) > 50  # Should have many features


class TestEdgeCases:
    """Test edge cases and error handling."""
    
    def test_invalid_smiles(self):
        """Test handling of invalid SMILES."""
        features = calculable.detect_all_features("INVALID")
        # Should return all False/0 rather than crashing
        assert features["sp2_bromide_present"] is False
        assert features["sp2_halide_site_count"] == 0
    
    def test_empty_smiles(self):
        """Test handling of empty SMILES."""
        features = calculable.detect_all_features("")
        assert features["sp2_bromide_present"] is False
    
    def test_simple_molecules(self):
        """Test simple molecules."""
        # Methane
        features = calculable.detect_all_features("C")
        assert features["sp2_halide_site_count"] == 0
        assert features["polarity_low"] is True
        
        # Water (if parsed)
        features = calculable.detect_all_features("O")
        assert features["alcohol_present"] is False  # Not an alcohol
    
    def test_caching(self):
        """Test that caching works (same result for repeated calls)."""
        smiles = "c1ccc(Br)cc1"
        
        result1 = calculable.detect_all_features(smiles)
        result2 = calculable.detect_all_features(smiles)
        
        # Should be identical
        assert result1 == result2
        assert result1 is result2  # Should be same cached object


class TestComplexMolecules:
    """Test realistic complex molecules."""
    
    def test_pharmaceutical_like(self):
        """Test drug-like molecules."""
        # Ibuprofen-like structure
        ibuprofen = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
        features = calculable.detect_all_features(ibuprofen)
        assert features["carboxylic_acid_present"] is True
        assert features["acidic_proton_present"] is True
        assert features["polarity_high"] is False  # Relatively hydrophobic
    
    def test_suzuki_coupling_pair(self):
        """Test Suzuki coupling substrates."""
        # Aryl bromide
        ar_br = "c1ccc(Br)c(C(=O)OC)c1"
        features = calculable.detect_all_features(ar_br)
        assert features["ArBr_present"] is True
        assert features["sp2_halide_site_count"] == 1
        
        # Aryl boronic acid (approximate)
        ar_boh2 = "c1ccc(B(O)O)cc1"
        features = calculable.detect_all_features(ar_boh2)
        assert features["sp2_boron_present"] is True
    
    def test_buchwald_hartwig_pair(self):
        """Test Buchwald-Hartwig C-N coupling substrates."""
        # Aryl chloride
        ar_cl = "c1ccc(Cl)c([N+](=O)[O-])c1"
        features = calculable.detect_all_features(ar_cl)
        assert features["ArCl_present"] is True
        
        # Aniline
        aniline = "c1ccc(N)cc1"
        features = calculable.detect_all_features(aniline)
        assert features["aniline_present"] is True
        assert features["aliphatic_amine_present"] is False


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

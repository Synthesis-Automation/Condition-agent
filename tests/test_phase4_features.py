"""
Test suite for Phase 4 calculable features expansion.

Tests:
- EWG/EDG classification (strong_ewg_present, strong_edg_present)
- Chelating group detection (bidentate_chelator_present, phosphine_present)
- Molecular weight categories (low_molecular_weight, high_molecular_weight)
- Ring complexity (fused_ring_system, spirocyclic_present)
- Chirality indicators (chiral_center_present, chiral_center_count)
"""

import pytest
from chemtools.featurizers.calculable import detect_all_features


class TestEWGEDGDetection:
    """Test electron-withdrawing and electron-donating group detection."""
    
    def test_strong_ewg_nitro(self):
        """Nitrobenzene has strong EWG."""
        result = detect_all_features("c1ccc(cc1)[N+](=O)[O-]")
        assert result["strong_ewg_present"] is True
    
    def test_strong_ewg_nitrile(self):
        """Benzonitrile has strong EWG."""
        result = detect_all_features("N#Cc1ccccc1")
        assert result["strong_ewg_present"] is True
    
    def test_strong_ewg_carboxylic_acid(self):
        """Benzoic acid has strong EWG."""
        result = detect_all_features("O=C(O)c1ccccc1")
        assert result["strong_ewg_present"] is True
    
    def test_strong_edg_phenol(self):
        """Phenol has strong EDG."""
        result = detect_all_features("Oc1ccccc1")
        assert result["strong_edg_present"] is True
    
    def test_strong_edg_methoxy(self):
        """Anisole (methoxybenzene) has strong EDG."""
        result = detect_all_features("COc1ccccc1")
        assert result["strong_edg_present"] is True
    
    def test_strong_edg_aniline(self):
        """Aniline has strong EDG."""
        result = detect_all_features("Nc1ccccc1")
        assert result["strong_edg_present"] is True
    
    def test_no_ewg_edg(self):
        """Benzene has neither strong EWG nor EDG."""
        result = detect_all_features("c1ccccc1")
        assert result["strong_ewg_present"] is False
        assert result["strong_edg_present"] is False


class TestChelatingGroups:
    """Test chelating group detection."""
    
    def test_bidentate_chelator_2_phenylpyridine(self):
        """2-Phenylpyridine is a bidentate chelator."""
        result = detect_all_features("c1ccccc1-c1ccccn1")
        assert result["bidentate_chelator_present"] is True
    
    def test_bidentate_chelator_amino_acid(self):
        """Amino acid (glycine) is a bidentate chelator."""
        result = detect_all_features("NCC(=O)O")
        assert result["bidentate_chelator_present"] is True
    
    def test_phosphine_present(self):
        """Triphenylphosphine."""
        result = detect_all_features("c1ccccc1P(c2ccccc2)c3ccccc3")
        assert result["phosphine_present"] is True
    
    def test_phosphine_absent(self):
        """Benzene has no phosphine."""
        result = detect_all_features("c1ccccc1")
        assert result["phosphine_present"] is False


class TestMolecularWeight:
    """Test molecular weight categories."""
    
    def test_low_molecular_weight(self):
        """Benzene (MW=78) is low MW."""
        result = detect_all_features("c1ccccc1")
        assert result["low_molecular_weight"] is True
        assert result["high_molecular_weight"] is False
    
    def test_medium_molecular_weight(self):
        """Biphenyl (MW=154) is still considered low MW (boundary case)."""
        result = detect_all_features("c1ccc(cc1)c2ccccc2")
        # MW is 154, which is < 200, so still low_molecular_weight
        assert result["low_molecular_weight"] is True
        assert result["high_molecular_weight"] is False
    
    def test_high_molecular_weight(self):
        """Large peptide-like molecule is high MW."""
        # Approximate tripeptide
        result = detect_all_features("CC(C)CC(NC(=O)C(Cc1ccccc1)NC(=O)C(Cc2ccccc2)NC(=O)C(Cc3ccccc3)N)C(=O)O")
        assert result["high_molecular_weight"] is True
        assert result["low_molecular_weight"] is False


class TestRingComplexity:
    """Test ring system complexity detection."""
    
    def test_fused_ring_naphthalene(self):
        """Naphthalene has fused ring system."""
        result = detect_all_features("c1ccc2ccccc2c1")
        assert result["fused_ring_system"] is True
    
    def test_fused_ring_quinoline(self):
        """Quinoline has fused ring system."""
        result = detect_all_features("c1ccc2ncccc2c1")
        assert result["fused_ring_system"] is True
    
    def test_fused_ring_absent(self):
        """Biphenyl has separate rings, not fused."""
        result = detect_all_features("c1ccc(cc1)c2ccccc2")
        assert result["fused_ring_system"] is False
    
    def test_spirocyclic_present(self):
        """Spirocyclic compound."""
        result = detect_all_features("C1CCC2(CC1)CCCC2")
        assert result["spirocyclic_present"] is True
    
    def test_spirocyclic_absent(self):
        """Cyclohexane is not spirocyclic."""
        result = detect_all_features("C1CCCCC1")
        assert result["spirocyclic_present"] is False


class TestChirality:
    """Test chiral center detection."""
    
    def test_chiral_center_present_single(self):
        """(S)-2-Butanol has one chiral center."""
        result = detect_all_features("CC[C@H](C)O")
        assert result["chiral_center_present"] is True
        assert result["chiral_center_count"] >= 1
    
    def test_chiral_center_present_multiple(self):
        """Molecule with multiple chiral centers."""
        result = detect_all_features("C[C@H](O)[C@@H](C)O")
        assert result["chiral_center_present"] is True
        assert result["chiral_center_count"] >= 2
    
    def test_chiral_center_absent(self):
        """Achiral molecule (benzene)."""
        result = detect_all_features("c1ccccc1")
        assert result["chiral_center_present"] is False
        assert result["chiral_center_count"] == 0
    
    def test_chiral_center_amino_acid(self):
        """L-Alanine has one chiral center."""
        result = detect_all_features("C[C@H](N)C(=O)O")
        assert result["chiral_center_present"] is True
        assert result["chiral_center_count"] >= 1


class TestPhase4Integration:
    """Integration tests combining Phase 4 features."""
    
    def test_ewg_with_fused_rings(self):
        """Nitronaphthalene combines EWG and fused rings."""
        result = detect_all_features("c1ccc2cc([N+](=O)[O-])ccc2c1")
        assert result["strong_ewg_present"] is True
        assert result["fused_ring_system"] is True
    
    def test_edg_with_low_mw(self):
        """Phenol is low MW with strong EDG."""
        result = detect_all_features("Oc1ccccc1")
        assert result["strong_edg_present"] is True
        assert result["low_molecular_weight"] is True
    
    def test_chiral_chelator(self):
        """Chiral amino acid is also a chelator."""
        result = detect_all_features("C[C@H](N)C(=O)O")
        assert result["chiral_center_present"] is True
        assert result["bidentate_chelator_present"] is True
    
    def test_complex_pharmaceutical(self):
        """Complex molecule with multiple Phase 4 features."""
        # Chiral, fused rings, high MW
        result = detect_all_features("CC1=C2[C@@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H](C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C")
        # This is paclitaxel (Taxol) - high MW, multiple chiral centers, fused rings
        assert result["high_molecular_weight"] is True
        assert result["chiral_center_present"] is True
        assert result["fused_ring_system"] is True


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

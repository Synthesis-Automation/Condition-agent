"""
Tests for substrate_classifier module.

Test comprehensive substrate classification including:
- Halides (alkyl, aryl, benzylic, allylic)
- Amines (aniline, aliphatic, primary/secondary/tertiary)
- Alcohols (phenol, benzylic, aliphatic)
- Carbonyls (acid, ester, amide, ketone)
- Boron compounds
- Special positions (benzylic, allylic, propargylic)
"""

import pytest
from chemtools.util.substrate_classifier import (
    SubstrateClassifier,
    classify_substrate,
    get_substrate_class,
    get_substrate_family,
    SpecialPositions,
)
from chemtools.util.rdkit_helpers import rdkit_available


class TestSubstrateClassifier:
    """Test SubstrateClassifier class"""
    
    def setup_method(self):
        """Setup test fixtures"""
        self.classifier = SubstrateClassifier()
    
    def test_init(self):
        """Test classifier initialization"""
        assert self.classifier is not None
        assert hasattr(self.classifier, 'classify')
    
    def test_empty_input(self):
        """Test with empty input"""
        info = self.classifier.classify("")
        assert info.substrate_class == "unknown"
        assert info.substrate_family == "unknown"


class TestAlkylHalides:
    """Test alkyl halide classification"""
    
    def setup_method(self):
        self.classifier = SubstrateClassifier()
    
    def test_primary_alkyl_iodide(self):
        """Test primary alkyl iodide"""
        info = self.classifier.classify("CCCCI")
        assert info.substrate_family == "halide"
        assert "iodide" in info.substrate_class
        assert "alkyl_iodide" in info.functional_groups
    
    def test_methyl_iodide(self):
        """Test methyl iodide (simplest primary)"""
        info = self.classifier.classify("CI")
        assert info.substrate_family == "halide"
        assert "iodide" in info.substrate_class
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_secondary_alkyl_bromide(self):
        """Test secondary alkyl bromide"""
        info = self.classifier.classify("CC(C)Br")
        assert info.substrate_family == "halide"
        assert "bromide" in info.substrate_class
        # With RDKit, should detect secondary
        if rdkit_available():
            assert "secondary" in info.substrate_class
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_tertiary_alkyl_chloride(self):
        """Test tertiary alkyl chloride"""
        info = self.classifier.classify("CC(C)(C)Cl")
        assert info.substrate_family == "halide"
        assert "chloride" in info.substrate_class
        if rdkit_available():
            assert "tertiary" in info.substrate_class
    
    def test_long_chain_iodide(self):
        """Test long chain alkyl iodide"""
        info = self.classifier.classify("CCCCCCCCI")
        assert info.substrate_family == "halide"
        assert "iodide" in info.substrate_class


class TestArylHalides:
    """Test aryl halide classification"""
    
    def setup_method(self):
        self.classifier = SubstrateClassifier()
    
    def test_aryl_iodide(self):
        """Test aryl iodide"""
        info = self.classifier.classify("c1ccccc1I")
        assert info.substrate_family == "halide"
        assert "aryl" in info.substrate_class
        assert "iodide" in info.substrate_class
        assert "aryl_iodide" in info.functional_groups
    
    def test_aryl_bromide(self):
        """Test aryl bromide"""
        info = self.classifier.classify("c1ccccc1Br")
        assert info.substrate_family == "halide"
        assert "bromide" in info.substrate_class
    
    def test_aryl_chloride(self):
        """Test aryl chloride"""
        info = self.classifier.classify("c1ccc(Cl)cc1")
        assert info.substrate_family == "halide"
        assert "chloride" in info.substrate_class
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_heteroaryl_bromide(self):
        """Test heteroaryl bromide (pyridine)"""
        info = self.classifier.classify("c1cnc(Br)cc1")
        assert info.substrate_family == "halide"
        assert "bromide" in info.substrate_class
        if rdkit_available():
            assert info.has_heteroaromatic


class TestBenzylicPositions:
    """Test benzylic position detection"""
    
    def setup_method(self):
        self.classifier = SubstrateClassifier()
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_benzyl_chloride(self):
        """Test benzyl chloride (benzylic)"""
        info = self.classifier.classify("c1ccccc1CCl")
        assert info.substrate_family == "halide"
        if rdkit_available():
            assert len(info.special_positions.benzylic) > 0
            assert "benzylic" in info.substrate_class
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_benzyl_iodide(self):
        """Test benzyl iodide"""
        info = self.classifier.classify("c1ccccc1CI")
        assert info.substrate_family == "halide"
        if rdkit_available():
            assert len(info.special_positions.benzylic) > 0
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_benzylic_alcohol(self):
        """Test benzylic alcohol"""
        info = self.classifier.classify("c1ccccc1CO")
        assert info.substrate_family == "alcohol"
        if rdkit_available():
            assert len(info.special_positions.benzylic) > 0


class TestAllylicPositions:
    """Test allylic position detection"""
    
    def setup_method(self):
        self.classifier = SubstrateClassifier()
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_allylic_bromide(self):
        """Test allylic bromide"""
        info = self.classifier.classify("C=CCBr")
        assert info.substrate_family == "halide"
        if rdkit_available():
            assert len(info.special_positions.allylic) > 0
            assert "allylic" in info.substrate_class
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_allylic_alcohol(self):
        """Test allylic alcohol"""
        info = self.classifier.classify("C=CCO")
        assert info.substrate_family == "alcohol"
        if rdkit_available():
            assert len(info.special_positions.allylic) > 0


class TestAmines:
    """Test amine classification"""
    
    def setup_method(self):
        self.classifier = SubstrateClassifier()
    
    def test_aniline(self):
        """Test aniline (aromatic amine)"""
        info = self.classifier.classify("c1ccccc1N")
        assert info.substrate_family == "amine"
        assert info.substrate_class == "aniline"
        assert "aniline" in info.functional_groups
    
    def test_primary_amine(self):
        """Test primary aliphatic amine"""
        info = self.classifier.classify("CCN")
        assert info.substrate_family == "amine"
        assert "primary" in info.substrate_class
        assert "amine_primary" in info.functional_groups
    
    def test_secondary_amine(self):
        """Test secondary amine"""
        info = self.classifier.classify("CC(C)N")
        # Note: This SMILES might parse as primary, use explicit
        info2 = self.classifier.classify("CCNC")
        assert info2.substrate_family == "amine"
    
    def test_amide_not_amine(self):
        """Test that amide is classified separately from amine"""
        info = self.classifier.classify("CC(=O)N")
        assert info.substrate_family == "amide"
        assert info.substrate_class in ["primary_amide", "secondary_amide", "tertiary_amide"]
        # Should NOT be classified as amine
        assert info.substrate_family != "amine"


class TestAlcohols:
    """Test alcohol classification"""
    
    def setup_method(self):
        self.classifier = SubstrateClassifier()
    
    def test_phenol(self):
        """Test phenol"""
        info = self.classifier.classify("c1ccccc1O")
        assert info.substrate_family == "alcohol"
        assert info.substrate_class == "phenol"
        assert "phenol" in info.functional_groups
    
    def test_aliphatic_alcohol(self):
        """Test aliphatic alcohol"""
        info = self.classifier.classify("CCCO")
        assert info.substrate_family == "alcohol"
        assert "alcohol" in info.substrate_class
        assert "alcohol" in info.functional_groups


class TestCarbonyls:
    """Test carbonyl compound classification"""
    
    def setup_method(self):
        self.classifier = SubstrateClassifier()
    
    def test_carboxylic_acid(self):
        """Test carboxylic acid"""
        info = self.classifier.classify("CC(=O)O")
        assert info.substrate_family == "carbonyl"
        assert info.substrate_class == "carboxylic_acid"
        assert "carboxylic_acid" in info.functional_groups
    
    def test_ester(self):
        """Test ester"""
        info = self.classifier.classify("CC(=O)OC")
        assert info.substrate_family == "carbonyl"
        assert info.substrate_class == "ester"
        assert "ester" in info.functional_groups
    
    def test_aldehyde(self):
        """Test aldehyde"""
        info = self.classifier.classify("CCC=O")
        assert info.substrate_family == "carbonyl"
        assert info.substrate_class == "aldehyde"
        assert "aldehyde" in info.functional_groups
    
    def test_ketone(self):
        """Test ketone"""
        info = self.classifier.classify("CC(=O)C")
        assert info.substrate_family == "carbonyl"
        assert info.substrate_class == "ketone"
        assert "ketone" in info.functional_groups


class TestBoronCompounds:
    """Test boron compound classification"""
    
    def setup_method(self):
        self.classifier = SubstrateClassifier()
    
    def test_boronic_acid(self):
        """Test boronic acid"""
        info = self.classifier.classify("c1ccccc1B(O)O")
        assert info.substrate_family == "boron"
        assert "boronic" in info.substrate_class
    
    def test_boronic_ester_pinacol(self):
        """Test boronic ester (pinacol)"""
        info = self.classifier.classify("c1ccccc1B1OC(C)(C)C(C)(C)O1")
        assert info.substrate_family == "boron"
        assert "boronic_ester" in info.substrate_class or "pinacol" in info.substrate_class


class TestCarbonTypes:
    """Test carbon hybridization detection"""
    
    def setup_method(self):
        self.classifier = SubstrateClassifier()
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_sp3_carbons(self):
        """Test sp3 carbon detection"""
        info = self.classifier.classify("CCCC")
        if rdkit_available():
            assert len(info.carbon_types) > 0
            assert all(ct == 'sp3' for ct in info.carbon_types.values())
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_sp2_carbons(self):
        """Test sp2 carbon detection"""
        info = self.classifier.classify("C=C")
        if rdkit_available():
            assert len(info.carbon_types) > 0
            assert any(ct == 'sp2' for ct in info.carbon_types.values())
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_sp_carbons(self):
        """Test sp carbon detection"""
        info = self.classifier.classify("C#C")
        if rdkit_available():
            assert len(info.carbon_types) > 0
            assert any(ct == 'sp' for ct in info.carbon_types.values())
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_aromatic_carbons(self):
        """Test aromatic carbon detection"""
        info = self.classifier.classify("c1ccccc1")
        if rdkit_available():
            assert len(info.carbon_types) > 0
            assert any(ct == 'aromatic' for ct in info.carbon_types.values())


class TestReactiveCenters:
    """Test reactive center detection"""
    
    def setup_method(self):
        self.classifier = SubstrateClassifier()
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_halide_reactive_center(self):
        """Test halide as reactive center"""
        info = self.classifier.classify("CCCCI")
        if rdkit_available():
            assert len(info.reactive_centers) > 0
            assert len(info.reactive_center_types) > 0
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_amine_reactive_center(self):
        """Test amine as reactive center"""
        info = self.classifier.classify("CCN")
        if rdkit_available():
            assert len(info.reactive_centers) > 0


class TestConvenienceFunctions:
    """Test convenience functions"""
    
    def test_get_substrate_class(self):
        """Test get_substrate_class function"""
        cls = get_substrate_class("c1ccccc1I")
        assert "aryl" in cls
        assert "iodide" in cls
    
    def test_get_substrate_family(self):
        """Test get_substrate_family function"""
        family = get_substrate_family("c1ccccc1I")
        assert family == "halide"
    
    def test_classify_substrate(self):
        """Test classify_substrate function"""
        info = classify_substrate("CCCCI")
        assert info.substrate_family == "halide"
        assert "iodide" in info.substrate_class


class TestRealWorldExamples:
    """Test real-world substrate examples"""
    
    def setup_method(self):
        self.classifier = SubstrateClassifier()
    
    def test_octyl_iodide(self):
        """Test octyl iodide (from protocol example)"""
        info = self.classifier.classify("CCCCCCCCI")
        assert info.substrate_family == "halide"
        assert "iodide" in info.substrate_class
        assert info.has_aromatic is False
    
    def test_phenyl_bromide(self):
        """Test phenyl bromide"""
        info = self.classifier.classify("c1ccccc1Br")
        assert info.substrate_family == "halide"
        assert "aryl" in info.substrate_class
        assert info.has_aromatic is True
    
    def test_4_bromobiphenyl(self):
        """Test 4-bromobiphenyl"""
        info = self.classifier.classify("c1ccc(Br)cc1-c2ccccc2")
        assert info.substrate_family == "halide"
        assert "aryl" in info.substrate_class
        assert info.has_aromatic is True
    
    def test_benzylamine(self):
        """Test benzylamine"""
        info = self.classifier.classify("c1ccccc1CN")
        assert info.substrate_family == "amine"
        assert "primary" in info.substrate_class
    
    @pytest.mark.skipif(not rdkit_available(), reason="RDKit required")
    def test_complex_molecule(self):
        """Test complex molecule with multiple functional groups"""
        # 4-iodoaniline
        info = self.classifier.classify("c1ccc(I)cc1N")
        assert info.substrate_family == "halide" or info.substrate_family == "amine"
        if rdkit_available():
            assert len(info.functional_groups) >= 2  # Has both halide and amine


class TestTextBasedFallback:
    """Test text-based classification when RDKit unavailable"""
    
    def test_aromatic_detection_without_rdkit(self):
        """Test aromatic detection works without RDKit"""
        classifier = SubstrateClassifier()
        info = classifier.classify("c1ccccc1Br")
        # Should still detect aromatic even without RDKit
        assert info.has_aromatic is True
    
    def test_functional_group_detection_without_rdkit(self):
        """Test functional group detection works without RDKit"""
        classifier = SubstrateClassifier()
        info = classifier.classify("CCCCI")
        # Should still detect functional groups
        assert len(info.functional_groups) > 0


class TestEdgeCases:
    """Test edge cases and error handling"""
    
    def setup_method(self):
        self.classifier = SubstrateClassifier()
    
    def test_invalid_smiles(self):
        """Test with invalid SMILES"""
        info = self.classifier.classify("INVALID")
        # Should not crash, return unknown
        assert info.substrate_class in ["unknown", "alkane", "hydrocarbon"]
    
    def test_none_input(self):
        """Test with None input"""
        info = self.classifier.classify(None)
        assert info.substrate_class == "unknown"
    
    def test_simple_alkane(self):
        """Test simple alkane (no functional groups)"""
        info = self.classifier.classify("CCCC")
        assert info.substrate_family == "hydrocarbon"
        assert info.substrate_class == "alkane"
    
    def test_benzene(self):
        """Test benzene (simple aromatic)"""
        info = self.classifier.classify("c1ccccc1")
        assert info.substrate_family == "hydrocarbon"
        assert info.substrate_class == "aromatic"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

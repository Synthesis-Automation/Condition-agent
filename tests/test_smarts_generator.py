"""Tests for reaction SMARTS applicability pattern generator."""

import pytest
from pathlib import Path
import json
import tempfile

from chemtools.protocol.smarts_generator_cli import (
    SmartsGenerator,
    ReactionSmartsApplicability,
    visualize_smarts_pattern,
    visualize_reaction_smarts,
    visualize_pattern_with_examples,
    main,
    RDKIT_AVAILABLE
)


class TestReactionSmartsApplicability:
    """Test data model for SMARTS applicability patterns."""
    
    def test_to_dict_minimal(self):
        """Test conversion to dict with minimal fields."""
        pattern = ReactionSmartsApplicability(
            core="[C:1]-[I:2]>>[C:1]-[B:3]"
        )
        result = pattern.to_dict()
        assert result == {"core": "[C:1]-[I:2]>>[C:1]-[B:3]"}
    
    def test_to_dict_full(self):
        """Test conversion to dict with all fields."""
        pattern = ReactionSmartsApplicability(
            core="[C:1]-[I:2]>>[C:1]-[B:3]",
            guards_forbid=["[CH]-I", "[C;H0]-I"],
            guards_require=["[C;H2]"],
            notes="Primary alkyl iodides only"
        )
        result = pattern.to_dict()
        assert result["core"] == "[C:1]-[I:2]>>[C:1]-[B:3]"
        assert result["guards_forbid"] == ["[CH]-I", "[C;H0]-I"]
        assert result["guards_require"] == ["[C;H2]"]
        assert result["notes"] == "Primary alkyl iodides only"
    
    def test_from_dict(self):
        """Test creation from dict."""
        data = {
            "core": "[C:1]-[I:2]>>[C:1]-[B:3]",
            "guards_forbid": ["[CH]-I"],
            "notes": "Test pattern"
        }
        pattern = ReactionSmartsApplicability.from_dict(data)
        assert pattern.core == "[C:1]-[I:2]>>[C:1]-[B:3]"
        assert pattern.guards_forbid == ["[CH]-I"]
        assert pattern.notes == "Test pattern"
        assert pattern.guards_require == []


class TestSmartsGenerator:
    """Test SMARTS pattern generator."""
    
    def test_parse_reaction_double_arrow(self):
        """Test parsing reaction with >> separator."""
        gen = SmartsGenerator("CCCCI>>CCCB")
        assert gen.reactants_smiles == "CCCCI"
        assert gen.products_smiles == "CCCB"
        assert gen.reagents_smiles == ""
    
    def test_parse_reaction_single_arrow(self):
        """Test parsing reaction with > separator."""
        gen = SmartsGenerator("CCCCI>CCCB")
        assert gen.reactants_smiles == "CCCCI"
        assert gen.products_smiles == "CCCB"
    
    def test_parse_reaction_with_reagents(self):
        """Test parsing reaction with reagents."""
        gen = SmartsGenerator("CCCCI>B2pin2>CCCBpin")
        assert gen.reactants_smiles == "CCCCI"
        assert gen.reagents_smiles == "B2pin2"
        assert gen.products_smiles == "CCCBpin"
    
    def test_parse_reaction_invalid(self):
        """Test parsing invalid reaction."""
        with pytest.raises(ValueError, match="Invalid reaction SMILES"):
            SmartsGenerator("CCCCI>>>CCCB")
    
    def test_parse_reaction_no_arrow(self):
        """Test parsing reaction without arrow."""
        with pytest.raises(ValueError, match="No reaction arrow"):
            SmartsGenerator("CCCCI CCCB")
    
    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not available")
    def test_generate_core_smarts_simple(self):
        """Test basic SMARTS generation."""
        gen = SmartsGenerator("CCCCI>>CCCB")
        core = gen.generate_core_smarts()
        assert isinstance(core, str)
        assert ">>" in core
    
    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not available")
    def test_generate_core_smarts_custom(self):
        """Test SMARTS generation with custom patterns."""
        gen = SmartsGenerator("CCCCI>>CCCB")
        core = gen.generate_core_smarts(
            reactant_pattern="[C:1]-[I:2]",
            product_pattern="[C:1]-[B:3]"
        )
        assert core == "[C:1]-[I:2]>>[C:1]-[B:3]"
    
    def test_generate_core_smarts_no_rdkit(self):
        """Test SMARTS generation without RDKit."""
        gen = SmartsGenerator("CCCCI>>CCCB")
        if not RDKIT_AVAILABLE:
            with pytest.raises(RuntimeError, match="RDKit required"):
                gen.generate_core_smarts()
    
    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not available")
    def test_suggest_guard_patterns_iodide(self):
        """Test guard pattern suggestions for iodide reactions."""
        gen = SmartsGenerator("CCCCI>>CCCB")
        guards = gen.suggest_guard_patterns()
        assert isinstance(guards, list)
        # Should suggest exclusions for tertiary, secondary, benzylic, allylic
        assert any("H0" in g for g in guards)  # tertiary
        assert any("CH" in g for g in guards)   # secondary
    
    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not available")
    def test_suggest_guard_patterns_bromide(self):
        """Test guard pattern suggestions for bromide reactions."""
        gen = SmartsGenerator("CCCBr>>CCCOH")
        guards = gen.suggest_guard_patterns()
        assert any("Br" in g for g in guards)


class TestCLI:
    """Test CLI interface."""
    
    def test_main_no_args(self):
        """Test CLI with no arguments prints help."""
        result = main([])
        assert result == 1
    
    def test_main_check_rdkit(self):
        """Test --check-rdkit flag."""
        result = main(['--check-rdkit'])
        assert result in (0, 1)  # Depends on whether RDKit is installed
    
    def test_main_single_reaction(self):
        """Test processing single reaction."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            output_path = Path(f.name)
        
        try:
            result = main([
                '--reaction', 'CCCCI>>CCCB',
                '--output', str(output_path)
            ])
            
            if RDKIT_AVAILABLE:
                assert result == 0
                assert output_path.exists()
                
                # Verify output format
                data = json.loads(output_path.read_text())
                assert "core" in data
                assert ">>" in data["core"]
            else:
                assert result == 1
        finally:
            if output_path.exists():
                output_path.unlink()
    
    def test_main_batch_processing(self):
        """Test batch processing from file."""
        # Create input file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("CCCCI>>CCCB\n")
            f.write("CCBr>>CCOH\n")
            input_path = Path(f.name)
        
        output_path = input_path.with_suffix('.json')
        
        try:
            result = main([
                '--batch', str(input_path),
                '--output', str(output_path)
            ])
            
            if RDKIT_AVAILABLE:
                assert result == 0
                assert output_path.exists()
                
                # Verify output format
                data = json.loads(output_path.read_text())
                assert isinstance(data, list)
                assert len(data) == 2
                
                for item in data:
                    assert "reaction_smiles" in item
                    if "smarts_applicability" in item:
                        assert "core" in item["smarts_applicability"]
            else:
                assert result == 0  # Should still work but with errors
        finally:
            if input_path.exists():
                input_path.unlink()
            if output_path.exists():
                output_path.unlink()


class TestRealWorldExample:
    """Test with real protocol example."""
    
    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not available")
    def test_alkyl_iodide_borylation(self):
        """Test with the alkyl iodide borylation example from the protocol DB."""
        reaction_smiles = "CCCCCCCCI.B(B1OC(C)(C)C(C)(C)O1)B2OC(C)(C)C(C)(C)O2>>CCCCCCCCB1OC(C)(C)C(C)(C)O1"
        
        gen = SmartsGenerator(reaction_smiles)
        
        # Should parse correctly
        assert "CCCCCCCCI" in gen.reactants_smiles
        assert "B" in gen.reactants_smiles
        assert "CCCCCCCCB" in gen.products_smiles
        
        # Generate pattern
        core = gen.generate_core_smarts()
        assert isinstance(core, str)
        assert ">>" in core
        
        # Should suggest relevant guards
        guards = gen.suggest_guard_patterns()
        assert len(guards) > 0
        assert any("I" in g for g in guards)
    
    def test_expected_output_format(self):
        """Test that output matches expected protocol JSON format."""
        expected_format = {
            "core": "[C:1;H2,H3;X4;!$(C-[a]);!$(C-[C]=[C])]-[I:2]>>[C:1]-[B:3;X3](-[O;H0])(-[O;H0])",
            "guards_forbid": ["[CH]-I", "[C;H0]-I", "[CH2]-[a]-I", "[CH2]-[C]=[C]-I"]
        }
        
        # Verify structure can be created
        pattern = ReactionSmartsApplicability.from_dict(expected_format)
        output = pattern.to_dict()
        
        assert output["core"] == expected_format["core"]
        assert output["guards_forbid"] == expected_format["guards_forbid"]


class TestVisualization:
    """Test visualization functions."""
    
    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not available")
    def test_visualize_smarts_pattern(self):
        """Test SMARTS pattern visualization."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test_pattern.png"
            result = visualize_smarts_pattern("[C;H2,H3]-I", output_path)
            assert result is True
            assert output_path.exists()
            assert output_path.stat().st_size > 0
    
    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not available")
    def test_visualize_smarts_pattern_invalid(self):
        """Test visualization with invalid SMARTS."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test_invalid.png"
            result = visualize_smarts_pattern("[C;INVALID", output_path)
            assert result is False
    
    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not available")
    def test_visualize_reaction_smarts(self):
        """Test reaction SMARTS visualization."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test_reaction.png"
            result = visualize_reaction_smarts("[C:1]-[I:2]>>[C:1]-[B:3]", output_path)
            assert result is True
            assert output_path.exists()
            assert output_path.stat().st_size > 0
    
    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not available")
    def test_visualize_reaction_smarts_invalid(self):
        """Test reaction visualization with invalid SMARTS."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test_invalid_rxn.png"
            result = visualize_reaction_smarts("[C;INVALID>>[C]", output_path)
            assert result is False
    
    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not available")
    def test_visualize_pattern_with_examples(self):
        """Test pattern visualization with test molecules."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "viz"
            
            pattern = ReactionSmartsApplicability(
                core="[C:1]-[I:2]>>[C:1]-[B:3]",
                guards_forbid=["[CH]-I", "[C;H0]-I"]
            )
            
            test_smiles = ["CCCCI", "CC(C)I", "CC(C)(C)I"]
            results = visualize_pattern_with_examples(pattern, test_smiles, output_dir)
            
            assert output_dir.exists()
            assert (output_dir / "core_transformation.png").exists()
            assert (output_dir / "guard_forbid_1.png").exists()
            assert (output_dir / "guard_forbid_2.png").exists()
            
            # Check test results
            assert len(results) == 3
            assert results["CCCCI"] is True  # Primary iodide - should match
            assert results["CC(C)I"] is False  # Secondary iodide - forbidden
            assert results["CC(C)(C)I"] is False  # Tertiary iodide - forbidden
    
    def test_visualize_no_rdkit(self):
        """Test visualization when RDKit not available."""
        if RDKIT_AVAILABLE:
            pytest.skip("RDKit is available")
        
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test.png"
            result = visualize_smarts_pattern("[C]-[I]", output_path)
            assert result is False


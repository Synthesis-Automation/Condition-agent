"""
Integration tests for chemtools.analysis module.
Tests the analyze_reaction function and end-to-end workflows.
"""
import pytest
from chemtools.analysis import analyze_reaction, normalize_reaction


class TestAnalyzeReaction:
    """Test the analyze_reaction function."""
    
    def test_buchwald_hartwig_reaction(self):
        """Analyze a Buchwald-Hartwig C-N coupling."""
        # Aryl bromide + amine -> aryl amine
        rxn_smiles = "c1ccccc1Br.CCN>>c1ccccc1NCC"
        
        result = analyze_reaction(rxn_smiles)
        
        # Check top-level structure
        assert "input" in result
        assert "normalized" in result
        assert "reactants" in result
        assert "agents" in result
        assert "products" in result
        assert "family" in result
        
        # Check input preserved
        assert result["input"] == rxn_smiles
        
        # Check reactants analyzed
        assert len(result["reactants"]) == 2
        
        # First reactant should be aryl bromide
        reactant1 = result["reactants"][0]
        assert "normalized" in reactant1
        assert "taxonomy" in reactant1
        taxonomy1 = reactant1["taxonomy"]
        assert "best_match" in taxonomy1
        
        # Second reactant should be amine
        reactant2 = result["reactants"][1]
        taxonomy2 = reactant2["taxonomy"]
        assert "best_match" in taxonomy2
        
        # Check family detection
        family = result["family"]
        assert "detected" in family
        assert "canonical_id" in family
        
        # Should detect as C-N coupling or Buchwald-Hartwig
        canonical = family.get("canonical_id")
        if canonical:
            assert "cn" in canonical.lower() or "buchwald" in canonical.lower()
    
    def test_suzuki_reaction(self):
        """Analyze a Suzuki-Miyaura coupling."""
        # Aryl bromide + aryl boronic acid -> biaryl
        rxn_smiles = "c1ccccc1Br.c1ccccc1B(O)O>>c1ccccc1c1ccccc1"
        
        result = analyze_reaction(rxn_smiles)
        
        # Check structure
        assert len(result["reactants"]) == 2
        
        # Should detect as Suzuki
        family = result["family"]
        canonical = family.get("canonical_id")
        if canonical:
            assert "suzuki" in canonical.lower()
    
    def test_reaction_with_agents(self):
        """Analyze reaction with explicit agents/catalysts."""
        # Include catalyst in the middle section
        rxn_smiles = "c1ccccc1Br.CCN>[Pd]>c1ccccc1NCC"
        
        result = analyze_reaction(rxn_smiles)
        
        # Check agents present
        assert "agents" in result
        agents = result["agents"]
        assert len(agents) >= 1
    
    def test_simple_reaction(self):
        """Test a simple 2-component reaction."""
        rxn_smiles = "CCBr.CCN>>CCNCC"
        
        result = analyze_reaction(rxn_smiles)
        
        # Should complete without errors
        assert result is not None
        assert len(result["reactants"]) == 2
    
    def test_reaction_type_metadata(self):
        """Check that reaction type metadata is populated when available."""
        rxn_smiles = "c1ccccc1Br.c1ccccc1B(O)O>>c1ccccc1c1ccccc1"
        
        result = analyze_reaction(rxn_smiles)
        
        family = result["family"]
        
        # If a canonical family is detected, check metadata
        if family.get("canonical_id"):
            # Should have reaction type info
            assert "reaction_type" in family
            # May be None if not in taxonomy
            
            # Should have required roles info
            assert "required_roles" in family
            
            # Should have reactant requirements
            assert "reactant_requirements" in family
    
    def test_use_rxn_insight_flag(self):
        """Test that use_rxn_insight flag is passed through."""
        rxn_smiles = "c1ccccc1Br.CCN>>c1ccccc1NCC"
        
        # Should work with both True and False
        result_with = analyze_reaction(rxn_smiles, use_rxn_insight=True)
        result_without = analyze_reaction(rxn_smiles, use_rxn_insight=False)
        
        assert result_with is not None
        assert result_without is not None
    
    def test_invalid_reaction(self):
        """Invalid reaction SMILES should be handled gracefully."""
        rxn_smiles = "not>>valid"
        
        result = analyze_reaction(rxn_smiles)
        
        # Should return result but possibly with errors
        assert result is not None
        assert "input" in result


class TestReactantTaxonomyMatches:
    """Test reactant taxonomy matching within analyze_reaction."""
    
    def test_all_matches_included(self):
        """Check that all_matches field is populated."""
        rxn_smiles = "c1ccccc1Br.CCN>>c1ccccc1NCC"
        
        result = analyze_reaction(rxn_smiles)
        
        for reactant in result["reactants"]:
            taxonomy = reactant["taxonomy"]
            assert "all_matches" in taxonomy
            assert isinstance(taxonomy["all_matches"], list)
    
    def test_category_matches_included(self):
        """Check that category_matches field is populated."""
        rxn_smiles = "c1ccccc1Br.CCN>>c1ccccc1NCC"
        
        result = analyze_reaction(rxn_smiles)
        
        for reactant in result["reactants"]:
            taxonomy = reactant["taxonomy"]
            assert "category_matches" in taxonomy
            assert isinstance(taxonomy["category_matches"], list)
    
    def test_best_match_structure(self):
        """Check best_match has expected structure."""
        rxn_smiles = "c1ccccc1Br>>c1ccccc1"
        
        result = analyze_reaction(rxn_smiles)
        
        reactant = result["reactants"][0]
        best_match = reactant["taxonomy"]["best_match"]
        
        if best_match is not None:
            # Should be a dict (converted from ReactantMatch)
            assert isinstance(best_match, dict)
            assert "category" in best_match
            assert "smarts" in best_match
            assert "specificity" in best_match


class TestNormalizedComponents:
    """Test that normalized components are correct."""
    
    def test_reactants_normalized(self):
        """Reactants should be normalized."""
        rxn_smiles = "c1ccccc1Br.CCN>>c1ccccc1NCC"
        
        result = analyze_reaction(rxn_smiles)
        
        for reactant in result["reactants"]:
            assert "normalized" in reactant
            norm = reactant["normalized"]
            assert "smiles_norm" in norm or "largest_smiles" in norm
    
    def test_products_normalized(self):
        """Products should be normalized."""
        rxn_smiles = "c1ccccc1Br.CCN>>c1ccccc1NCC"
        
        result = analyze_reaction(rxn_smiles)
        
        assert "products" in result
        assert len(result["products"]) >= 1
    
    def test_normalized_reaction_smiles(self):
        """Normalized reaction SMILES should be generated."""
        rxn_smiles = "c1ccccc1Br.CCN>>c1ccccc1NCC"
        
        result = analyze_reaction(rxn_smiles)
        
        assert "normalized" in result
        norm_dict = result["normalized"]
        assert "normalized" in norm_dict  # Nested normalized SMILES string


class TestEdgeCases:
    """Test edge cases and error handling."""
    
    def test_single_reactant(self):
        """Reaction with single reactant should work."""
        rxn_smiles = "CC(O)>>CC(=O)"
        
        result = analyze_reaction(rxn_smiles)
        
        assert result is not None
        assert len(result["reactants"]) == 1
    
    def test_multiple_products(self):
        """Reaction with multiple products should work."""
        rxn_smiles = "CC(O)>>CC(=O).O"
        
        result = analyze_reaction(rxn_smiles)
        
        assert result is not None
        assert len(result["products"]) >= 2
    
    def test_empty_agents(self):
        """Reaction with no agents should work."""
        rxn_smiles = "CCBr.CCN>>CCNCC"
        
        result = analyze_reaction(rxn_smiles)
        
        assert "agents" in result
        # Should be empty or minimal
        assert len(result["agents"]) == 0
    
    def test_complex_reaction(self):
        """Complex multi-step-like reaction should be handled."""
        rxn_smiles = "c1ccccc1Br.c1ccccc1B(O)O.CCN>[Pd].Base>c1ccccc1c1ccccc1.c1ccccc1NCC"
        
        result = analyze_reaction(rxn_smiles)
        
        # Should complete without crashing
        assert result is not None
        assert "reactants" in result
        assert "agents" in result
        assert "products" in result


class TestRequiredRoles:
    """Test required roles metadata when available."""
    
    def test_required_roles_structure(self):
        """Required roles should have expected structure."""
        rxn_smiles = "c1ccccc1Br.c1ccccc1B(O)O>>c1ccccc1c1ccccc1"
        
        result = analyze_reaction(rxn_smiles)
        
        family = result["family"]
        required_roles = family.get("required_roles")
        
        if required_roles:
            assert isinstance(required_roles, list)
            for role in required_roles:
                assert "role_id" in role
                assert "required" in role


class TestReactantRequirements:
    """Test reactant requirements metadata when available."""
    
    def test_reactant_requirements_structure(self):
        """Reactant requirements should have expected structure."""
        rxn_smiles = "c1ccccc1Br.c1ccccc1B(O)O>>c1ccccc1c1ccccc1"
        
        result = analyze_reaction(rxn_smiles)
        
        family = result["family"]
        reactant_reqs = family.get("reactant_requirements")
        
        if reactant_reqs:
            assert isinstance(reactant_reqs, list)
            for req in reactant_reqs:
                assert "reactant_type_id" in req


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

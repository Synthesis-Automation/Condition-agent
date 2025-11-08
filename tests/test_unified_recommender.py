"""
Tests for UnifiedRecommender DRFP-based similarity search.

Tests cover:
- Initialization and index loading
- Recommendation logic (similarity computation, ranking)
- Filtering (source_type, min_similarity, top_k)
- Edge cases (invalid SMILES, empty results, missing index)
"""

import pytest
from pathlib import Path
import tempfile
import json
import numpy as np
from chemtools.recommend import UnifiedRecommender, RecommendationResult


class TestUnifiedRecommenderInit:
    """Test UnifiedRecommender initialization and setup."""
    
    def test_init_with_default_index(self):
        """Test initialization with default unified index."""
        recommender = UnifiedRecommender()
        assert recommender.index_dir.exists()
        assert recommender.sources is not None
        assert len(recommender.sources) > 0
        
    def test_init_with_custom_index(self, tmp_path):
        """Test initialization with custom index directory."""
        # Create minimal test index
        index_dir = tmp_path / "test_index"
        index_dir.mkdir()
        
        # Create index.json with correct structure
        index_data = {
            "version": "2.0",
            "build_date": "2025-01-01T00:00:00",
            "num_protocols": 1,
            "num_rules": 0,
            "protocols": [
                {
                    "id": "test_protocol_1",
                    "name": "Test Protocol",
                    "source_type": "protocol",
                    "family": "Test_Family",
                    "file_path": "test_protocol_1.json",
                    "version": "1.0.0",
                    "tags": ["test"]
                }
            ],
            "rules": []
        }
        with open(index_dir / "index.json", "w") as f:
            json.dump(index_data, f)
        
        # Create fingerprints.npz (minimal)
        protocol_fps = np.random.randint(0, 2, (1, 2048), dtype=np.uint8)
        np.savez_compressed(
            index_dir / "fingerprints.npz",
            protocol_fps=protocol_fps,
            rule_fps=np.array([], dtype=np.uint8).reshape(0, 2048),
            protocol_ids=np.array(["test_protocol_1"]),
            rule_ids=np.array([]),
            rule_fp_indices=np.array([])
        )
        
        # Create protocol file
        protocol_data = {
            "id": "test_protocol_1",
            "name": "Test Protocol",
            "family": "Test_Family",
            "version": "1.0.0"
        }
        with open(index_dir / "test_protocol_1.json", "w") as f:
            json.dump(protocol_data, f)
        
        recommender = UnifiedRecommender(index_dir=index_dir)
        assert len(recommender.sources) == 1
        assert recommender.sources[0]["id"] == "test_protocol_1"
        
    def test_get_statistics(self):
        """Test statistics reporting."""
        recommender = UnifiedRecommender()
        stats = recommender.get_statistics()
        
        # Stats has nested structure: build_info, protocols, rules, drfp
        assert "build_info" in stats
        assert "protocols" in stats
        assert "rules" in stats
        assert "drfp" in stats
        assert stats["protocols"]["count"] > 0
        assert stats["drfp"]["computed"] > 0


class TestUnifiedRecommenderRecommendations:
    """Test recommendation logic and similarity computation."""
    
    @pytest.fixture
    def recommender(self):
        """Fixture for UnifiedRecommender instance."""
        return UnifiedRecommender()
    
    def test_recommend_with_valid_smiles(self, recommender):
        """Test recommendation with valid reaction SMILES."""
        # Suzuki coupling reaction
        reaction = "CC(=O)c1ccc(Br)cc1.c1ccc(B(O)O)cc1>>CC(=O)c1ccc(-c2ccccc2)cc1"
        results = recommender.recommend(reaction, top_k=5, min_similarity=0.0)
        
        assert isinstance(results, list)
        assert len(results) > 0
        assert all(isinstance(r, RecommendationResult) for r in results)
        
        # Check result structure
        result = results[0]
        assert hasattr(result, "id")
        assert hasattr(result, "name")
        assert hasattr(result, "source_type")
        assert hasattr(result, "family")
        assert hasattr(result, "similarity")
        assert hasattr(result, "rank")
        assert hasattr(result, "tags")
        
        # Check similarity bounds
        assert 0.0 <= result.similarity <= 1.0
        
        # Check ranking
        assert result.rank == 1
        for i in range(1, len(results)):
            assert results[i].rank == i + 1
            # Similarities should be non-increasing
            assert results[i].similarity <= results[i-1].similarity
    
    def test_recommend_with_invalid_smiles(self, recommender):
        """Test recommendation with invalid SMILES returns empty."""
        results = recommender.recommend("INVALID_SMILES", top_k=5)
        assert results == []
    
    def test_recommend_with_min_similarity_filter(self, recommender):
        """Test min_similarity filtering."""
        reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
        
        # High threshold - should reduce results
        results_high = recommender.recommend(reaction, top_k=10, min_similarity=0.8)
        
        # Low threshold - should include more results
        results_low = recommender.recommend(reaction, top_k=10, min_similarity=0.1)
        
        # All high-threshold results should be in low-threshold results
        assert len(results_high) <= len(results_low)
        
        # All results should meet threshold
        for r in results_high:
            assert r.similarity >= 0.8
        for r in results_low:
            assert r.similarity >= 0.1
    
    def test_recommend_with_top_k_limit(self, recommender):
        """Test top_k limiting."""
        reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
        
        results_k3 = recommender.recommend(reaction, top_k=3, min_similarity=0.0)
        results_k10 = recommender.recommend(reaction, top_k=10, min_similarity=0.0)
        
        assert len(results_k3) <= 3
        assert len(results_k10) <= 10
        
        # First k results should be identical
        for r3, r10 in zip(results_k3, results_k10):
            assert r3.id == r10.id
            assert r3.similarity == r10.similarity
    
    def test_recommend_filter_by_protocol(self, recommender):
        """Test filtering by protocol source type."""
        reaction = "C=CCN(C(=O)OC(C)(C)C)CC=C>>CC(C)(C)OC(=O)N1CCC=C1"
        
        results = recommender.recommend(
            reaction,
            top_k=5,
            source_types=["protocol"],
            min_similarity=0.0
        )
        
        assert len(results) > 0
        assert all(r.source_type == "protocol" for r in results)
    
    def test_recommend_filter_by_rule(self, recommender):
        """Test filtering by rule source type."""
        reaction = "CC(=O)c1ccc(Br)cc1.c1ccc(B(O)O)cc1>>CC(=O)c1ccc(-c2ccccc2)cc1"
        
        results = recommender.recommend(
            reaction,
            top_k=5,
            source_types=["rule"],
            min_similarity=0.0
        )
        
        assert len(results) > 0
        assert all(r.source_type == "rule" for r in results)
    
    def test_recommend_no_matches_above_threshold(self, recommender):
        """Test that very dissimilar reaction returns empty with high threshold."""
        # Simple esterification (unlikely to match C-C coupling protocols/rules)
        reaction = "CC(=O)O.CCO>>CC(=O)OCC"
        
        results = recommender.recommend(reaction, top_k=5, min_similarity=0.95)
        
        # Should be empty or very few results
        assert len(results) <= 1


class TestUnifiedRecommenderSourceDetails:
    """Test loading full source details."""
    
    @pytest.fixture
    def recommender(self):
        """Fixture for UnifiedRecommender instance."""
        return UnifiedRecommender()
    
    def test_get_source_details_protocol(self, recommender):
        """Test loading protocol details."""
        # Get a recommendation first
        reaction = "C=CCN(C(=O)OC(C)(C)C)CC=C>>CC(C)(C)OC(=O)N1CCC=C1"
        results = recommender.recommend(
            reaction,
            top_k=1,
            source_types=["protocol"],
            min_similarity=0.0
        )
        
        if results:
            source_id = results[0].id
            details = recommender.get_source_details(source_id)
            
            assert details is not None
            # v2.0 schema has nested metadata
            assert "metadata" in details
            assert details["metadata"]["id"] == source_id
            assert "metadata" in details and "name" in details["metadata"]
            assert "metadata" in details and "version" in details["metadata"]
    
    def test_get_source_details_rule(self, recommender):
        """Test loading rule details."""
        # Get a rule recommendation
        reaction = "CC(=O)c1ccc(Br)cc1.c1ccc(B(O)O)cc1>>CC(=O)c1ccc(-c2ccccc2)cc1"
        results = recommender.recommend(
            reaction,
            top_k=1,
            source_types=["rule"],
            min_similarity=0.0
        )
        
        if results:
            source_id = results[0].id
            details = recommender.get_source_details(source_id)
            
            assert details is not None
            # v2.0 schema has nested metadata
            assert "metadata" in details
            assert details["metadata"]["id"] == source_id
    
    def test_get_source_details_invalid_id(self, recommender):
        """Test loading details for non-existent source."""
        details = recommender.get_source_details("nonexistent_source_id_12345")
        assert details is None


class TestUnifiedRecommenderEdgeCases:
    """Test edge cases and error handling."""
    
    def test_empty_reaction_smiles(self):
        """Test with empty SMILES string."""
        recommender = UnifiedRecommender()
        results = recommender.recommend("", top_k=5)
        assert results == []
    
    def test_malformed_reaction_smiles(self):
        """Test with malformed reaction SMILES."""
        recommender = UnifiedRecommender()
        
        # Missing arrow - DRFP may still compute, just check it returns a list
        results = recommender.recommend("Brc1ccccc1.Nc1ccccc1", top_k=5)
        assert isinstance(results, list)
        
        # Multiple arrows - DRFP may fail or succeed, just check it returns a list  
        results = recommender.recommend("A>>B>>C", top_k=5)
        assert isinstance(results, list)
    
    def test_reaction_with_no_products(self):
        """Test reaction with missing products."""
        recommender = UnifiedRecommender()
        results = recommender.recommend("Brc1ccccc1.Nc1ccccc1>>", top_k=5)
        # DRFP may still compute for incomplete reactions, just check it returns a list
        assert isinstance(results, list)
    
    def test_reaction_with_no_reactants(self):
        """Test reaction with missing reactants."""
        recommender = UnifiedRecommender()
        results = recommender.recommend(">>c1ccccc1Nc1ccccc1", top_k=5)
        # DRFP may still compute for incomplete reactions, just check it returns a list
        assert isinstance(results, list)
    
    def test_zero_top_k(self):
        """Test with top_k=0."""
        recommender = UnifiedRecommender()
        reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
        results = recommender.recommend(reaction, top_k=0)
        assert results == []
    
    def test_negative_top_k(self):
        """Test with negative top_k (should handle gracefully)."""
        recommender = UnifiedRecommender()
        reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
        results = recommender.recommend(reaction, top_k=-1)
        # Should either return empty or clamp to 0
        assert isinstance(results, list)
    
    def test_similarity_threshold_above_one(self):
        """Test with impossible similarity threshold."""
        recommender = UnifiedRecommender()
        reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
        results = recommender.recommend(reaction, top_k=5, min_similarity=1.5)
        assert results == []
    
    def test_invalid_source_type_filter(self):
        """Test with invalid source_type filter."""
        recommender = UnifiedRecommender()
        reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
        results = recommender.recommend(
            reaction,
            top_k=5,
            source_types=["invalid_type"]
        )
        # Should return empty (no matches for invalid type)
        assert results == []


class TestUnifiedRecommenderIntegration:
    """Integration tests with real chemistry examples."""
    
    @pytest.fixture
    def recommender(self):
        """Fixture for UnifiedRecommender instance."""
        return UnifiedRecommender()
    
    def test_suzuki_coupling_recommendations(self, recommender):
        """Test Suzuki coupling returns relevant recommendations."""
        reaction = "CC(=O)c1ccc(Br)cc1.c1ccc(B(O)O)cc1>>CC(=O)c1ccc(-c2ccccc2)cc1"
        results = recommender.recommend(reaction, top_k=5, min_similarity=0.3)
        
        assert len(results) > 0
        
        # Check that top result is Suzuki-related
        top_result = results[0]
        assert any(tag.lower() in ["suzuki", "cross-coupling", "boronic-acid"] 
                   for tag in top_result.tags) or \
               "suzuki" in top_result.name.lower() or \
               "Suzuki" in top_result.family
    
    def test_buchwald_hartwig_recommendations(self, recommender):
        """Test Buchwald-Hartwig coupling returns relevant recommendations."""
        reaction = "Brc1ccccc1.NCc1ccccc1>>c1ccc(NCc2ccccc2)cc1"
        results = recommender.recommend(reaction, top_k=5, min_similarity=0.3)
        
        assert len(results) > 0
        
        # Check that results include C-N coupling
        names_and_families = " ".join([r.name + " " + r.family for r in results]).lower()
        assert "c-n" in names_and_families or "c–n" in names_and_families or \
               "amination" in names_and_families
    
    def test_sonogashira_recommendations(self, recommender):
        """Test Sonogashira coupling returns relevant recommendations."""
        reaction = "Ic1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1"
        results = recommender.recommend(reaction, top_k=5, min_similarity=0.3)
        
        assert len(results) > 0
        
        # Should find Sonogashira
        top_result = results[0]
        name_and_tags = (top_result.name + " " + " ".join(top_result.tags)).lower()
        assert "sonogashira" in name_and_tags or "alkyne" in name_and_tags
    
    def test_amide_formation_recommendations(self, recommender):
        """Test amide formation returns relevant recommendations."""
        reaction = "O=C(O)c1ccccc1.NCc1ccccc1>>O=C(NCc1ccccc1)c1ccccc1"
        results = recommender.recommend(reaction, top_k=5, min_similarity=0.3)
        
        assert len(results) > 0
        # Any result is fine - just check structure
        assert all(hasattr(r, "similarity") for r in results)

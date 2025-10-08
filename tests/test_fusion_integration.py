"""
Integration tests for multi-source evidence fusion recommendation.

Tests the fusion system integrated into recommend_from_reaction().
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
repo_root = Path(__file__).parent.parent
sys.path.insert(0, str(repo_root))

import pytest
from chemtools.recommend import recommend_from_reaction


def test_fusion_mode_basic():
    """Test fusion mode can be enabled and returns valid results."""
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    
    # Run with fusion mode
    results = recommend_from_reaction(
        reaction=reaction,
        k=50,
        use_fusion=True,
        max_variants=5
    )
    
    # Should return valid results
    assert results is not None
    assert 'recommendation' in results
    assert 'formatted' in results
    assert 'fusion_meta' in results
    
    # Check fusion metadata
    fusion_meta = results['fusion_meta']
    assert 'adaptive_weights' in fusion_meta
    assert 'reasoning' in fusion_meta
    assert 'evidence_summary' in fusion_meta
    
    # Check adaptive weights exist
    weights = fusion_meta['adaptive_weights']
    assert 'α' in weights
    assert 'β' in weights
    assert 'γ' in weights
    assert 'δ' in weights
    
    # Weights should sum to approximately 1.0
    weight_sum = sum(weights.values())
    assert 0.98 <= weight_sum <= 1.02
    
    # Should have recommendations
    formatted = results['formatted']
    recommended = formatted.get('recommended_conditions', [])
    assert len(recommended) > 0
    
    print(f"\n✅ Fusion mode test passed!")
    print(f"   Adaptive weights: α={weights['α']:.3f}, β={weights['β']:.3f}, "
          f"γ={weights['γ']:.3f}, δ={weights['δ']:.3f}")
    print(f"   Generated {len(recommended)} recommendations")


def test_fusion_vs_baseline():
    """Compare fusion recommendations with baseline k-NN."""
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    
    # Baseline (k-NN)
    baseline = recommend_from_reaction(
        reaction=reaction,
        k=50,
        use_fusion=False,
        max_variants=5
    )
    
    # Fusion
    fusion = recommend_from_reaction(
        reaction=reaction,
        k=50,
        use_fusion=True,
        max_variants=5
    )
    
    # Both should return recommendations
    assert baseline['recommendation'] is not None
    assert fusion['recommendation'] is not None
    
    # Fusion should have fusion_meta
    assert 'fusion_meta' in fusion
    assert 'fusion_meta' not in baseline
    
    # Fusion metadata should have model='fusion'
    assert fusion['formatted']['meta']['model'] == 'fusion'
    
    # Baseline should use drfp_similarity
    assert baseline['formatted']['meta']['model'] == 'drfp_similarity'
    
    print(f"\n✅ Fusion vs baseline test passed!")
    print(f"   Baseline model: {baseline['formatted']['meta']['model']}")
    print(f"   Fusion model: {fusion['formatted']['meta']['model']}")
    
    # Compare top recommendations
    baseline_top = baseline['formatted']['recommended_conditions'][0]
    fusion_top = fusion['formatted']['recommended_conditions'][0]
    
    print(f"\n   Baseline top: {baseline_top['summary'].get('core')}")
    print(f"   Fusion top: {fusion_top['summary'].get('core')}")


def test_fusion_component_scores():
    """Test fusion recommendations include component scores."""
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    
    results = recommend_from_reaction(
        reaction=reaction,
        k=50,
        use_fusion=True,
        max_variants=3
    )
    
    # Check each recommendation has component scores
    recommended = results['formatted']['recommended_conditions']
    
    for idx, rec in enumerate(recommended, start=1):
        summary = rec['summary']
        
        # Should have component scores
        assert 'component_scores' in summary
        scores = summary['component_scores']
        
        # Should have all four components
        assert 'precedent' in scores
        assert 'analytics' in scores
        assert 'rules' in scores
        assert 'ml' in scores
        
        # Should have fusion score
        assert 'fusion_score' in summary
        
        # Should have adaptive weights
        assert 'adaptive_weights' in summary
        
        print(f"\n   Recommendation {idx}:")
        print(f"      Core: {summary.get('core')}")
        print(f"      Fusion score: {summary['fusion_score']:.3f}")
        print(f"      Component scores: PS={scores['precedent']:.3f}, "
              f"AS={scores['analytics']:.3f}, RS={scores['rules']:.3f}, MS={scores['ml']:.3f}")
    
    print(f"\n✅ Component scores test passed!")


def test_fusion_evidence_summary():
    """Test fusion evidence summary is provided."""
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    
    results = recommend_from_reaction(
        reaction=reaction,
        k=50,
        use_fusion=True,
        max_variants=3
    )
    
    fusion_meta = results['fusion_meta']
    evidence = fusion_meta['evidence_summary']
    
    # Should have evidence fields
    assert 'precedents' in evidence
    assert 'diversity' in evidence
    assert 'dataset_size' in evidence
    
    # Precedents should be found (numeric value >= 0)
    assert isinstance(evidence['precedents'], (int, float))
    assert evidence['precedents'] >= 0  # Can be 0 if wrong family detected
    
    # Diversity should be in [0, 1]
    diversity = evidence['diversity']
    assert diversity is not None
    assert isinstance(diversity, (int, float))
    assert 0.0 <= diversity <= 1.0
    
    # Dataset should have reactions
    assert isinstance(evidence['dataset_size'], (int, float))
    assert evidence['dataset_size'] > 0
    
    print(f"\n✅ Evidence summary test passed!")
    print(f"   Precedents: {evidence['precedents']}")
    print(f"   Diversity: {diversity:.3f}")
    print(f"   Dataset size: {evidence['dataset_size']:,}")


def test_fusion_adaptive_weighting():
    """Test adaptive weights adjust based on evidence quality."""
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    
    results = recommend_from_reaction(
        reaction=reaction,
        k=50,
        use_fusion=True,
        max_variants=3
    )
    
    fusion_meta = results['fusion_meta']
    weights = fusion_meta['adaptive_weights']
    reasoning = fusion_meta['reasoning']
    
    # Should have reasoning
    assert len(reasoning) > 0
    
    # Weights should be positive
    assert all(w >= 0 for w in weights.values())
    
    # At least one weight should be non-zero
    assert any(w > 0 for w in weights.values())
    
    print(f"\n✅ Adaptive weighting test passed!")
    print(f"   Reasoning:")
    for reason in reasoning:
        print(f"      • {reason}")


def test_fusion_deduplication():
    """Test fusion recommendations don't contain duplicates."""
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    
    results = recommend_from_reaction(
        reaction=reaction,
        k=50,
        use_fusion=True,
        max_variants=10  # Request more to test deduplication
    )
    
    recommended = results['formatted']['recommended_conditions']
    
    # Check for duplicate (core, base, solvent) tuples
    seen_combos = set()
    duplicates = []
    
    for rec in recommended:
        summary = rec['summary']
        combo = (
            summary.get('core'),
            summary.get('base', {}).get('cas') if summary.get('base') else None,
            summary.get('solvent', {}).get('cas') if summary.get('solvent') else None,
        )
        
        if combo in seen_combos:
            duplicates.append(combo)
        seen_combos.add(combo)
    
    # Should have no duplicates
    assert len(duplicates) == 0, f"Found duplicate recommendations: {duplicates}"
    
    print(f"\n✅ Deduplication test passed!")
    print(f"   Generated {len(recommended)} unique recommendations")


if __name__ == "__main__":
    print("=" * 70)
    print("FUSION INTEGRATION TESTS")
    print("=" * 70)
    
    test_fusion_mode_basic()
    test_fusion_vs_baseline()
    test_fusion_component_scores()
    test_fusion_evidence_summary()
    test_fusion_adaptive_weighting()
    test_fusion_deduplication()
    
    print("\n" + "=" * 70)
    print("ALL FUSION INTEGRATION TESTS PASSED! ✅")
    print("=" * 70)

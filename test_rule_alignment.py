"""
Test Rule-Alignment Scoring for ML Recommendations

This script demonstrates the new rule-alignment scoring system that reranks
ML-based recommendations based on their similarity to rule-based results.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from chemtools.ml.rule_alignment import (
    rerank_ml_by_rule_alignment,
    compute_rule_alignment_score,
    explain_alignment,
    DEFAULT_ALIGNMENT_WEIGHTS
)


def create_mock_ml_recommendations():
    """Create mock ML recommendations for testing."""
    return [
        {
            'rank': 1,
            'reaction': {'smiles': 'Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1'},
            'chemicals': [
                {'name': 'Bromobenzene', 'cas': None, 'smiles': 'Brc1ccccc1', 'role': 'starting_material', 'equivalents': None},
                {'name': 'Aniline', 'cas': None, 'smiles': 'Nc1ccccc1', 'role': 'starting_material', 'equivalents': None},
                {'name': 'Palladium(II) acetate', 'cas': '3375-31-3', 'role': 'metal_precursor', 'smiles': None, 'equivalents': None},
                {'name': 'XPhos', 'cas': '564483-18-7', 'role': 'ligand', 'smiles': None, 'equivalents': None},
                {'name': 'Sodium tert-butoxide', 'cas': '865-48-5', 'role': 'base', 'smiles': None, 'equivalents': None},
                {'name': 'Toluene', 'cas': '108-88-3', 'role': 'solvent', 'smiles': None, 'equivalents': None},
            ],
            'conditions': {
                'temperature': {'value': 100, 'unit': '°C'},
                'time': {'value': 12, 'unit': 'hours'}
            },
            'summary': {
                'rank': 1,
                'core': 'Pd/XPhos',
                'base': {'name': 'Sodium tert-butoxide', 'cas': '865-48-5', 'role': 'base'},
                'solvent': {'name': 'Toluene', 'cas': '108-88-3', 'role': 'solvent'},
                'confidence': 0.85
            }
        },
        {
            'rank': 2,
            'reaction': {'smiles': 'Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1'},
            'chemicals': [
                {'name': 'Bromobenzene', 'cas': None, 'smiles': 'Brc1ccccc1', 'role': 'starting_material', 'equivalents': None},
                {'name': 'Aniline', 'cas': None, 'smiles': 'Nc1ccccc1', 'role': 'starting_material', 'equivalents': None},
                {'name': 'Palladium(II) acetate', 'cas': '3375-31-3', 'role': 'metal_precursor', 'smiles': None, 'equivalents': None},
                {'name': 'BINAP', 'cas': '98327-87-8', 'role': 'ligand', 'smiles': None, 'equivalents': None},
                {'name': 'Potassium carbonate', 'cas': '584-08-7', 'role': 'base', 'smiles': None, 'equivalents': None},
                {'name': 'Dioxane', 'cas': '123-91-1', 'role': 'solvent', 'smiles': None, 'equivalents': None},
            ],
            'conditions': {
                'temperature': {'value': 110, 'unit': '°C'},
                'time': {'value': 16, 'unit': 'hours'}
            },
            'summary': {
                'rank': 2,
                'core': 'Pd/BINAP',
                'base': {'name': 'Potassium carbonate', 'cas': '584-08-7', 'role': 'base'},
                'solvent': {'name': 'Dioxane', 'cas': '123-91-1', 'role': 'solvent'},
                'confidence': 0.78
            }
        },
        {
            'rank': 3,
            'reaction': {'smiles': 'Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1'},
            'chemicals': [
                {'name': 'Bromobenzene', 'cas': None, 'smiles': 'Brc1ccccc1', 'role': 'starting_material', 'equivalents': None},
                {'name': 'Aniline', 'cas': None, 'smiles': 'Nc1ccccc1', 'role': 'starting_material', 'equivalents': None},
                {'name': 'Tris(dibenzylideneacetone)dipalladium(0)', 'cas': '51364-51-3', 'role': 'metal_precursor', 'smiles': None, 'equivalents': None},
                {'name': 'SPhos', 'cas': '657408-07-6', 'role': 'ligand', 'smiles': None, 'equivalents': None},
                {'name': 'Cesium carbonate', 'cas': '534-17-8', 'role': 'base', 'smiles': None, 'equivalents': None},
                {'name': 'Toluene', 'cas': '108-88-3', 'role': 'solvent', 'smiles': None, 'equivalents': None},
            ],
            'conditions': {
                'temperature': {'value': 90, 'unit': '°C'},
                'time': {'value': 8, 'unit': 'hours'}
            },
            'summary': {
                'rank': 3,
                'core': 'Pd/SPhos',
                'base': {'name': 'Cesium carbonate', 'cas': '534-17-8', 'role': 'base'},
                'solvent': {'name': 'Toluene', 'cas': '108-88-3', 'role': 'solvent'},
                'confidence': 0.72
            }
        }
    ]


def create_mock_rule_recommendations():
    """Create mock rule-based recommendations for testing."""
    return [
        {
            'rank': 1,
            'chemicals': [
                {'name': 'Palladium(II) acetate', 'cas': '3375-31-3', 'role': 'catalyst', 'smiles': None, 'equivalents': 0.01},
                {'name': 'XPhos', 'cas': '564483-18-7', 'role': 'ligand', 'smiles': None, 'equivalents': 0.02},
                {'name': 'Sodium tert-butoxide', 'cas': '865-48-5', 'role': 'base', 'smiles': None, 'equivalents': 1.5},
                {'name': 'Toluene', 'cas': '108-88-3', 'role': 'solvent', 'smiles': None, 'equivalents': None},
            ],
            'conditions': {
                'temperature': {'value': 100, 'unit': '°C'},
                'time': {'value': 12, 'unit': 'hours'}
            },
            'core': 'Pd/XPhos',
            'confidence': 0.95
        },
        {
            'rank': 2,
            'chemicals': [
                {'name': 'Tris(dibenzylideneacetone)dipalladium(0)', 'cas': '51364-51-3', 'role': 'catalyst', 'smiles': None, 'equivalents': 0.005},
                {'name': 'SPhos', 'cas': '657408-07-6', 'role': 'ligand', 'smiles': None, 'equivalents': 0.01},
                {'name': 'Potassium tert-butoxide', 'cas': '865-47-4', 'role': 'base', 'smiles': None, 'equivalents': 1.4},
                {'name': 'Toluene', 'cas': '108-88-3', 'role': 'solvent', 'smiles': None, 'equivalents': None},
            ],
            'conditions': {
                'temperature': {'value': 90, 'unit': '°C'},
                'time': {'value': 8, 'unit': 'hours'}
            },
            'core': 'Pd/SPhos',
            'confidence': 0.88
        }
    ]


def test_rule_alignment():
    """Test the rule-alignment scoring system."""
    print("=" * 80)
    print("Testing Rule-Alignment Scoring for ML Recommendations")
    print("=" * 80)
    print()
    
    # Create mock data
    ml_recs = create_mock_ml_recommendations()
    rule_recs = create_mock_rule_recommendations()
    
    print("Original ML Rankings:")
    print("-" * 80)
    for rec in ml_recs:
        summary = rec['summary']
        print(f"Rank {rec['rank']}: {summary['core']} | "
              f"Base: {summary['base']['name']} | "
              f"Solvent: {summary['solvent']['name']} | "
              f"ML Confidence: {summary['confidence']:.2f}")
    print()
    
    print("Rule-Based Recommendations:")
    print("-" * 80)
    for rec in rule_recs:
        chems = {c['role']: c['name'] for c in rec['chemicals']}
        print(f"Rank {rec['rank']}: {rec['core']} | "
              f"Base: {chems.get('base', 'N/A')} | "
              f"Solvent: {chems.get('solvent', 'N/A')} | "
              f"Rule Confidence: {rec['confidence']:.2f}")
    print()
    
    # Show alignment weights
    print("Alignment Weights:")
    print("-" * 80)
    for component, weight in DEFAULT_ALIGNMENT_WEIGHTS.items():
        print(f"  {component:15s}: {weight:.2f}")
    print()
    
    # Compute individual alignment scores
    print("Individual Alignment Scores:")
    print("-" * 80)
    for ml_rec in ml_recs:
        print(f"\nML Rank {ml_rec['rank']} ({ml_rec['summary']['core']}):")
        for rule_rec in rule_recs:
            score, components = compute_rule_alignment_score(ml_rec, rule_rec)
            print(f"  vs Rule {rule_rec['rank']} ({rule_rec['core']}): {score:.3f}")
            # Show top component matches
            sorted_comps = sorted(components.items(), key=lambda x: x[1], reverse=True)
            print(f"    Top matches: {', '.join([f'{k}={v:.2f}' for k, v in sorted_comps[:3]])}")
    print()
    
    # Rerank ML recommendations
    print("Reranking ML Recommendations by Rule Alignment...")
    print("-" * 80)
    reranked = rerank_ml_by_rule_alignment(
        ml_recs,
        rule_recs,
        weights=DEFAULT_ALIGNMENT_WEIGHTS,
        ml_weight=0.6,
        alignment_weight=0.4
    )
    
    print("\nReranked ML Recommendations:")
    print("-" * 80)
    for rec in reranked:
        summary = rec['summary']
        alignment = rec['alignment_meta']
        
        print(f"\nRank {rec['rank']} (was Rank {alignment['original_rank']}): {summary['core']}")
        print(f"  Base: {summary['base']['name']}")
        print(f"  Solvent: {summary['solvent']['name']}")
        print(f"  Original ML Score: {alignment['original_ml_score']:.3f}")
        print(f"  Rule Alignment Score: {alignment['rule_alignment']['alignment_score']:.3f}")
        print(f"  Combined Score: {alignment['combined_score']:.3f}")
        print(f"  Best Rule Match: #{alignment['rule_alignment']['best_match_index'] + 1}")
        print(f"  Reasoning:")
        for reason in alignment['rule_alignment']['reasoning'][:3]:
            print(f"    - {reason}")
    
    print()
    print("=" * 80)
    print("Detailed Alignment Explanation for Top Recommendation:")
    print("=" * 80)
    
    top_rec = reranked[0]
    best_rule_idx = top_rec['alignment_meta']['rule_alignment']['best_match_index']
    explanation = explain_alignment(top_rec, rule_recs[best_rule_idx])
    
    print(f"\nML Rec: {top_rec['summary']['core']} vs Rule Rec: {rule_recs[best_rule_idx]['core']}")
    print(f"Total Alignment Score: {explanation['total_alignment_score']:.3f}")
    print()
    print("Component Breakdown:")
    print("-" * 80)
    print(f"{'Component':<15} {'Score':<8} {'Weight':<8} {'Contrib.':<10} {'Quality':<12}")
    print("-" * 80)
    for item in explanation['breakdown']:
        print(f"{item['component']:<15} {item['score']:<8.3f} {item['weight']:<8.2f} "
              f"{item['contribution']:<10.3f} {item['match_quality']:<12}")
    
    print()
    print("✅ Test Complete!")
    print()


if __name__ == "__main__":
    test_rule_alignment()

"""
Demo: Fusion Recommendation System for Suzuki Reactions
========================================================

This script demonstrates the new multi-source evidence fusion recommendation system
on Suzuki coupling reactions from the sample dataset. It compares:

1. Baseline k-NN recommendations (use_fusion=False)
2. Fusion recommendations (use_fusion=True) with adaptive weighting

The fusion system balances:
- Precedent Evidence (PE): DRFP-based k-NN with diversity scoring
- Dataset Statistics (DS): Analytics from common catalytic systems, bases, solvents
- Rule Evidence (RE): SMARTS-based scheme matching
- ML Predictions (ML): DRFP yield predictor (currently δ=0)

Key Features Demonstrated:
- Adaptive weight adjustments based on evidence quality
- Diversity scoring to detect precedent bias
- Component score breakdown (PS, AS, RS, MS)
- Evidence quality assessment
- Detailed reasoning for recommendations
"""

import sys
from pathlib import Path

# Add repo root to path
repo_root = Path(__file__).parent
sys.path.insert(0, str(repo_root))

from chemtools.recommend.core import recommend_from_reaction
import json

# Note: We don't actually need SAMPLE_REACTIONS import - reactions defined inline below


def print_section(title):
    """Print a formatted section header"""
    print("\n" + "=" * 80)
    print(f"  {title}")
    print("=" * 80)


def print_recommendations(recs, mode="Baseline"):
    """Print recommendations in a readable format"""
    print(f"\n{mode} Recommendations ({len(recs)} total):")
    print("-" * 80)
    
    for i, rec in enumerate(recs[:5], 1):  # Show top 5
        # Handle different possible structures
        if 'summary' in rec:
            # Structure from Suzuki formatted output
            summary = rec['summary']
            core = summary.get('core', 'N/A')
            base = summary.get('base', {})
            solvent = summary.get('solvent', {})
            confidence = summary.get('confidence', 'UNKNOWN')
            
            base_name = base.get('name', 'N/A') if isinstance(base, dict) else str(base)
            solvent_name = solvent.get('name', 'N/A') if isinstance(solvent, dict) else str(solvent)
            
            # Score from component_scores if available, otherwise use rank
            score = rec.get('score', rec.get('rank', 0.0))
            
            print(f"\n{i}. Core: {core} | Confidence: {confidence}")
            print(f"   Base: {base_name}")
            print(f"   Solvent: {solvent_name}")
            
            # Component scores (fusion only)
            if 'component_scores' in rec:
                cs = rec['component_scores']
                print(f"   Component Scores: PS={cs.get('PS', 0):.3f}, AS={cs.get('AS', 0):.3f}, " + 
                      f"RS={cs.get('RS', 0):.3f}, MS={cs.get('MS', 0):.3f}")
            
            # Precedents
            precedents = summary.get('precedents', [])
            if precedents:
                print(f"   Precedents: {len(precedents)} supporting examples")
                
        else:
            # Original structure
            catalyst = rec.get('catalyst', {})
            base = rec.get('base', {})
            solvent = rec.get('solvent', {})
            
            cat_core = catalyst.get('core', 'N/A')
            cat_ligand = catalyst.get('ligand', 'N/A')
            base_name = base.get('name', 'N/A')
            solvent_name = solvent.get('name', 'N/A')
            score = rec.get('score', 0.0)
            confidence = rec.get('confidence', 'UNKNOWN')
            
            print(f"\n{i}. Score: {score:.3f} | Confidence: {confidence}")
            print(f"   Catalyst: {cat_core} + {cat_ligand}")
            print(f"   Base: {base_name}")
            print(f"   Solvent: {solvent_name}")
            
            # Component scores (fusion only)
            if 'component_scores' in rec:
                cs = rec['component_scores']
                print(f"   Component Scores: PS={cs.get('PS', 0):.3f}, AS={cs.get('AS', 0):.3f}, " + 
                      f"RS={cs.get('RS', 0):.3f}, MS={cs.get('MS', 0):.3f}")
            
            # Evidence count
            if 'evidence_count' in rec:
                print(f"   Evidence: {rec['evidence_count']} precedents")


def print_fusion_metadata(result):
    """Print fusion-specific metadata"""
    # Check for fusion_meta key
    if 'fusion_meta' not in result:
        print("\nNo fusion metadata available")
        return
        
    metadata = result.get('fusion_meta', {})
    
    print("\nFusion System Metadata:")
    print("-" * 80)
    
    # Adaptive weights
    weights = metadata.get('adaptive_weights', {})
    if weights:
        print(f"Adaptive Weights: α={weights.get('alpha', 0):.3f} (precedents), " + 
              f"β={weights.get('beta', 0):.3f} (analytics), " + 
              f"γ={weights.get('gamma', 0):.3f} (rules), " + 
              f"δ={weights.get('delta', 0):.3f} (ML)")
    
    # Evidence summary
    evidence = metadata.get('evidence_summary', {})
    if evidence:
        prec_count = evidence.get('precedent_count', evidence.get('precedents', 0))
        diversity = evidence.get('diversity_score', evidence.get('diversity', 0.0))
        avg_sim = evidence.get('avg_similarity', 0.0)
        
        print(f"\nEvidence Quality:")
        print(f"  - Precedent count: {prec_count}")
        print(f"  - Diversity score: {diversity:.3f} " + 
              ("(LOW - may indicate batch effect)" if diversity < 0.3 else "(OK)"))
        if avg_sim > 0:
            print(f"  - Avg similarity: {avg_sim:.3f}")
    
    # Reasoning
    reasoning = metadata.get('reasoning', [])
    if reasoning:
        print(f"\nAdaptive Weight Reasoning:")
        for reason in reasoning:
            print(f"  - {reason}")


def test_suzuki_reaction(smiles_str, description):
    """Test a single Suzuki reaction with both baseline and fusion"""
    print_section(f"Testing: {description}")
    print(f"SMILES: {smiles_str}\n")
    
    # Baseline recommendations
    print("\n[1] BASELINE k-NN RECOMMENDATIONS")
    print("=" * 80)
    try:
        baseline = recommend_from_reaction(
            reaction=smiles_str,
            use_fusion=False,
            k=5
        )
        
        print(f"DEBUG - Baseline result keys: {list(baseline.keys())}")
        print(f"DEBUG - Baseline result type: {type(baseline)}")
        
        if 'recommended_conditions' in baseline and baseline['recommended_conditions']:
            print_recommendations(baseline['recommended_conditions'], mode="Baseline")
        elif 'recommendations' in baseline and baseline['recommendations']:
            print_recommendations(baseline['recommendations'], mode="Baseline")
        else:
            print(f"No recommendations returned. Full result: {baseline}")
            
    except Exception as e:
        print(f"ERROR in baseline: {e}")
    
    # Fusion recommendations
    print("\n\n[2] FUSION MULTI-SOURCE RECOMMENDATIONS")
    print("=" * 80)
    try:
        fusion = recommend_from_reaction(
            reaction=smiles_str,
            use_fusion=True,
            k=5
        )
        
        print(f"DEBUG - Fusion result keys: {list(fusion.keys())}")
        print(f"DEBUG - Fusion 'formatted' keys: {list(fusion.get('formatted', {}).keys())}")
        
        # Check multiple possible keys
        if 'formatted' in fusion and 'recommended_conditions' in fusion['formatted']:
            print_recommendations(fusion['formatted']['recommended_conditions'], mode="Fusion")
            print_fusion_metadata(fusion)
        elif 'recommended_conditions' in fusion:
            print_recommendations(fusion['recommended_conditions'], mode="Fusion")
            print_fusion_metadata(fusion)
        else:
            print(f"No recommendations returned. Keys available: {list(fusion.keys())}")
            if 'formatted' in fusion:
                print(f"Formatted keys: {list(fusion['formatted'].keys())}")
            
    except Exception as e:
        print(f"ERROR in fusion: {e}")
    
    print("\n" + "=" * 80)


def main():
    """Run demo on selected Suzuki reactions"""
    print_section("FUSION RECOMMENDATION SYSTEM DEMO - SUZUKI REACTIONS")
    
    print("""
This demo compares baseline k-NN recommendations with the new fusion system.

The fusion system addresses the concern:
  "ML based recommendation should be based on both precedents and dataset analytics
   and could consider the rule-based result. How to use these data is critical,
   we do not want a random precedents (could be top-k precedent) take too much weight."

Key improvements:
  1. Diversity scoring: Detects when top-k precedents are biased (batch effects)
  2. Dataset analytics: Uses statistics to validate precedent recommendations
  3. Adaptive weighting: Reduces precedent weight when diversity is low
  4. Rule-based evidence: Incorporates chemistry knowledge (SMARTS patterns)
  5. Transparency: Shows component scores and reasoning for each recommendation
""")
    
    # Select diverse Suzuki reactions
    suzuki_reactions = [
        ("Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1", 
         "Suzuki - Simple Ph-Ph (benchmark)"),
        
        ("Clc1ccc(C#N)cc1.c1ccc(B(O)O)cc1>>N#Cc1ccc(-c2ccccc2)cc1", 
         "Suzuki - Electron-poor ArCl (challenging)"),
        
        ("Brc1ccc(OC)cc1.c1ccc(B(O)O)cc1>>COc1ccc(-c2ccccc2)cc1", 
         "Suzuki - Electron-rich ArBr (donor groups)"),
        
        ("Ic1ccncc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccncc2)cc1", 
         "Suzuki - Heteroaryl pyridine (coordination issues)"),
        
        ("Brc1ccccc1OCC.c1ccc(B(O)O)cc1C>>CCOc1ccccc1-c1ccccc1C", 
         "Suzuki - Ortho-substituted (sterically hindered)"),
    ]
    
    # Test each reaction
    for i, (smiles, desc) in enumerate(suzuki_reactions, 1):
        if i > 1:
            input("\nPress Enter to continue to next reaction...")
        test_suzuki_reaction(smiles, desc)
    
    print_section("DEMO COMPLETE")
    print("""
Summary:
  - Tested 5 diverse Suzuki coupling reactions
  - Compared baseline k-NN vs multi-source fusion
  - Demonstrated adaptive weight adjustments
  - Showed component score breakdowns
  - Highlighted diversity detection

The fusion system prevents precedent bias by:
  1. Measuring diversity of top-k precedents
  2. Using dataset statistics for validation
  3. Adjusting weights when evidence quality is low
  4. Providing transparent reasoning for recommendations
""")


if __name__ == "__main__":
    main()

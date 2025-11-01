#!/usr/bin/env python3
"""
Comprehensive Cross-Family Recommendation Test

Tests the enhanced cross-family workflow using diverse sample reactions 
to demonstrate marking and ranking functionality across different reaction families.
"""

import sys
import json
import time
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))

# Import not needed for this test - reactions defined inline


def run_comprehensive_cross_family_test():
    """Run comprehensive tests using sample reactions from different families."""
    
    print("🧪 Comprehensive Cross-Family Recommendation Test")
    print("=" * 70)
    print()
    
    # Select representative reactions from different families
    test_reactions = {
        "Suzuki (C-C)": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
        "C-N Coupling (Pd)": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1", 
        "C-N Coupling (Cu)": "Brc1ccc(C(F)(F)F)cc1.Nc1ccccc1>>FC(F)(F)c1ccc(Nc2ccccc2)cc1",
        "Heck": "Brc1ccccc1.C=C>>C(=Cc1ccccc1)",
        "Sonogashira": "Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1",
        "Amidation": "CC(=O)O.Nc1ccccc1>>CC(=O)Nc1ccccc1",
    }
    
    try:
        from chemtools import chem
        
        for reaction_name, reaction_smiles in test_reactions.items():
            print(f"🔬 Testing: {reaction_name}")
            print(f"📝 Reaction: {reaction_smiles}")
            print()
            
            # Test cross-family search with enhanced marking
            start_time = time.time()
            result = chem.recommend.conditions(
                reaction=reaction_smiles,
                k=30,
                search_all_families=True,
                relax={
                    'reaction_type_threshold': 0.12,  # 12% threshold
                    'mechanism_similarity_threshold': 0.35,  # Keep more for analysis
                    'mechanism_weight': 0.25,
                    'debug_timing': False
                },
                rerank_strategy='none'
            )
            
            processing_time = time.time() - start_time
            
            # Analyze results
            analyze_cross_family_results(reaction_name, result, processing_time)
            print("-" * 50)
            print()
            
    except Exception as e:
        print(f"❌ Test failed: {e}")
        import traceback
        traceback.print_exc()


def analyze_cross_family_results(reaction_name: str, result: dict, processing_time: float):
    """Analyze and display cross-family recommendation results."""
    
    recommendations = result.get('recommendations', [])
    detection = result.get('detection', {})
    
    print(f"⚡ Processing time: {processing_time:.2f}s")
    print(f"🎯 Detected family: {detection.get('canonical_family', 'Unknown')}")
    print(f"📊 Total recommendations: {len(recommendations)}")
    
    if not recommendations:
        print("❌ No recommendations found")
        return
    
    # Analyze compatibility breakdown
    compatibility_stats = {
        'compatible': 0,
        'moderate': 0,
        'incompatible': 0,
        'unknown': 0
    }
    
    representation_stats = {
        'well_represented': 0,
        'underrepresented': 0,
        'unknown': 0
    }
    
    family_distribution = {}
    similarity_scores = []
    mechanism_scores = []
    
    for rec in recommendations:
        conditions = rec.get('conditions', {})
        cross_meta = conditions.get('cross_family_metadata', {})
        
        # Compatibility stats
        mech_status = cross_meta.get('mechanism_status', 'unknown')
        compatibility_stats[mech_status] = compatibility_stats.get(mech_status, 0) + 1
        
        # Representation stats
        type_status = cross_meta.get('reaction_type_status', 'unknown')
        representation_stats[type_status] = representation_stats.get(type_status, 0) + 1
        
        # Family distribution
        prec_family = cross_meta.get('precedent_family', 'Unknown')
        family_distribution[prec_family] = family_distribution.get(prec_family, 0) + 1
        
        # Similarity scores
        sim = conditions.get('similarity', 0)
        if sim > 0:
            similarity_scores.append(sim)
        
        mech_sim = cross_meta.get('mechanism_similarity', 0)
        if mech_sim > 0:
            mechanism_scores.append(mech_sim)
    
    # Display compatibility breakdown
    print("🏷️ Compatibility Breakdown:")
    print(f"   • Compatible: {compatibility_stats['compatible']} ✅")
    print(f"   • Moderate: {compatibility_stats['moderate']} ⚠️")
    print(f"   • Incompatible: {compatibility_stats['incompatible']} ❌")
    
    print("📈 Representation Status:")
    print(f"   • Well represented: {representation_stats['well_represented']} 📊")
    print(f"   • Underrepresented: {representation_stats['underrepresented']} 📉")
    
    # Display family distribution (top 5)
    if family_distribution:
        sorted_families = sorted(family_distribution.items(), key=lambda x: x[1], reverse=True)[:5]
        print("🎯 Top Families Found:")
        for family, count in sorted_families:
            percentage = (count / len(recommendations)) * 100
            print(f"   • {family}: {count} ({percentage:.1f}%)")
    
    # Display similarity statistics
    if similarity_scores:
        avg_sim = sum(similarity_scores) / len(similarity_scores)
        min_sim = min(similarity_scores)
        max_sim = max(similarity_scores)
        print(f"📊 Similarity Stats: Avg={avg_sim:.3f}, Min={min_sim:.3f}, Max={max_sim:.3f}")
    
    if mechanism_scores:
        avg_mech = sum(mechanism_scores) / len(mechanism_scores)
        min_mech = min(mechanism_scores)
        max_mech = max(mechanism_scores)
        print(f"🔧 Mechanism Stats: Avg={avg_mech:.3f}, Min={min_mech:.3f}, Max={max_mech:.3f}")
    
    # Show top 3 recommendations with detailed status
    print("🏆 Top 3 Recommendations:")
    for i, rec in enumerate(recommendations[:3], 1):
        conditions = rec.get('conditions', {})
        cross_meta = conditions.get('cross_family_metadata', {})
        
        core = rec.get('core', 'Unknown')
        family = cross_meta.get('precedent_family', 'Unknown')
        similarity = conditions.get('similarity', 0)
        mech_sim = cross_meta.get('mechanism_similarity', 0)
        mech_status = cross_meta.get('mechanism_status', 'unknown')
        type_status = cross_meta.get('reaction_type_status', 'unknown')
        
        # Status icons
        mech_icon = {'compatible': '✅', 'moderate': '⚠️', 'incompatible': '❌'}.get(mech_status, '❓')
        type_icon = {'well_represented': '📊', 'underrepresented': '📉'}.get(type_status, '❓')
        
        print(f"   {i}. {core[:40]}")
        print(f"      Family: {family} | Sim: {similarity:.3f} | Mech: {mech_sim:.3f}")
        print(f"      Status: {mech_icon} {mech_status}, {type_icon} {type_status.replace('_', ' ')}")


def test_family_specific_vs_cross_family():
    """Compare family-specific vs cross-family search results."""
    
    print("\n🔍 Family-Specific vs Cross-Family Comparison")
    print("=" * 70)
    
    # Use a C-N coupling reaction for comparison
    test_reaction = "Brc1ccc(C(F)(F)F)cc1.Nc1ccccc1>>FC(F)(F)c1ccc(Nc2ccccc2)cc1"
    
    try:
        from chemtools import chem
        
        print(f"📝 Test Reaction: C-N coupling with CF3 substrate")
        print()
        
        # Family-specific search
        print("1️⃣ Family-Specific Search:")
        family_result = chem.recommend.conditions(
            reaction=test_reaction,
            k=25,
            search_all_families=False
        )
        
        family_recs = family_result.get('recommendations', [])
        detected_family = family_result.get('detection', {}).get('canonical_family', 'Unknown')
        print(f"   Detected family: {detected_family}")
        print(f"   Found: {len(family_recs)} recommendations")
        
        if family_recs:
            family_cores = [rec.get('core', 'Unknown') for rec in family_recs[:3]]
            print(f"   Top cores: {', '.join(family_cores)}")
        print()
        
        # Cross-family search
        print("2️⃣ Cross-Family Search:")
        cross_result = chem.recommend.conditions(
            reaction=test_reaction,
            k=40,
            search_all_families=True,
            relax={
                'reaction_type_threshold': 0.10,
                'mechanism_similarity_threshold': 0.30,
                'mechanism_weight': 0.25
            }
        )
        
        cross_recs = cross_result.get('recommendations', [])
        print(f"   Found: {len(cross_recs)} recommendations")
        
        if cross_recs:
            # Count compatible vs incompatible
            compatible = sum(1 for rec in cross_recs 
                           if rec.get('conditions', {}).get('cross_family_metadata', {}).get('mechanism_status') == 'compatible')
            moderate = sum(1 for rec in cross_recs 
                         if rec.get('conditions', {}).get('cross_family_metadata', {}).get('mechanism_status') == 'moderate')
            incompatible = sum(1 for rec in cross_recs 
                             if rec.get('conditions', {}).get('cross_family_metadata', {}).get('mechanism_status') == 'incompatible')
            
            print(f"   Compatible: {compatible} ✅ | Moderate: {moderate} ⚠️ | Incompatible: {incompatible} ❌")
            
            # Show family diversity
            families = {}
            for rec in cross_recs:
                family = rec.get('conditions', {}).get('cross_family_metadata', {}).get('precedent_family', 'Unknown')
                families[family] = families.get(family, 0) + 1
            
            sorted_families = sorted(families.items(), key=lambda x: x[1], reverse=True)[:5]
            print(f"   Family diversity: {len(families)} families")
            for family, count in sorted_families:
                print(f"     • {family}: {count}")
        
        print()
        print("📊 Summary:")
        print(f"   Family-specific: {len(family_recs)} recommendations")
        print(f"   Cross-family: {len(cross_recs)} recommendations")
        if cross_recs:
            coverage_increase = ((len(cross_recs) - len(family_recs)) / max(len(family_recs), 1)) * 100
            print(f"   Coverage increase: {coverage_increase:.1f}%")
        
    except Exception as e:
        print(f"❌ Comparison failed: {e}")


def test_parameter_sensitivity():
    """Test how different parameter settings affect results."""
    
    print("\n🎛️ Parameter Sensitivity Test")
    print("=" * 70)
    
    test_reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"  # Simple C-N coupling
    
    parameter_sets = [
        {
            'name': 'Conservative',
            'params': {
                'reaction_type_threshold': 0.20,
                'mechanism_similarity_threshold': 0.60,
                'mechanism_weight': 0.40
            }
        },
        {
            'name': 'Balanced',
            'params': {
                'reaction_type_threshold': 0.15,
                'mechanism_similarity_threshold': 0.40,
                'mechanism_weight': 0.30
            }
        },
        {
            'name': 'Exploratory',
            'params': {
                'reaction_type_threshold': 0.05,
                'mechanism_similarity_threshold': 0.25,
                'mechanism_weight': 0.20
            }
        }
    ]
    
    try:
        from chemtools import chem
        
        print(f"📝 Test Reaction: Simple C-N coupling")
        print()
        
        for param_set in parameter_sets:
            print(f"🔧 {param_set['name']} Parameters:")
            
            result = chem.recommend.conditions(
                reaction=test_reaction,
                k=30,
                search_all_families=True,
                relax=param_set['params']
            )
            
            recs = result.get('recommendations', [])
            
            if recs:
                # Analyze compatibility distribution
                compatible = sum(1 for rec in recs 
                               if rec.get('conditions', {}).get('cross_family_metadata', {}).get('mechanism_status') == 'compatible')
                moderate = sum(1 for rec in recs 
                             if rec.get('conditions', {}).get('cross_family_metadata', {}).get('mechanism_status') == 'moderate')
                incompatible = sum(1 for rec in recs 
                                 if rec.get('conditions', {}).get('cross_family_metadata', {}).get('mechanism_status') == 'incompatible')
                
                total = len(recs)
                print(f"   Total: {total} | Compatible: {compatible} ({compatible/total:.1%}) | Moderate: {moderate} ({moderate/total:.1%}) | Incompatible: {incompatible} ({incompatible/total:.1%})")
                
                # Show top mechanism similarity
                mech_sims = [rec.get('conditions', {}).get('cross_family_metadata', {}).get('mechanism_similarity', 0) 
                            for rec in recs[:5]]
                avg_top_mech = sum(mech_sims) / len(mech_sims) if mech_sims else 0
                print(f"   Avg top-5 mechanism similarity: {avg_top_mech:.3f}")
            else:
                print(f"   No recommendations found")
            
            print()
        
    except Exception as e:
        print(f"❌ Parameter test failed: {e}")


def main():
    """Run all test suites."""
    
    print("🚀 Enhanced Cross-Family Recommendation Test Suite")
    print("=" * 80)
    print()
    
    # Run comprehensive tests
    run_comprehensive_cross_family_test()
    
    # Run comparison test  
    test_family_specific_vs_cross_family()
    
    # Run parameter sensitivity test
    test_parameter_sensitivity()
    
    print("\n✅ Test suite completed!")
    print("=" * 80)


if __name__ == "__main__":
    main()
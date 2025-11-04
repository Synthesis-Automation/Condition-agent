#!/usr/bin/env python3
"""
Enhanced Cross-Family Recommendation CLI with Mechanism-Aware Filtering

This enhanced version includes:
- Reaction type threshold filtering
- Mechanism similarity scoring
- Configurable filtering parameters
- Quality metrics for cross-family results
"""

import sys
import json
from pathlib import Path

# Add project root to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from chemtools import chem
from chemtools.output_formatter import ensure_standard_output


def enhanced_cross_family_recommend(
    reaction_smiles: str,
    k: int = 50,
    reaction_type_threshold: float = 0.15,
    mechanism_similarity_threshold: float = 0.4,
    mechanism_weight: float = 0.3,
    show_debug: bool = False
) -> dict:
    """
    Get cross-family recommendations with enhanced filtering.
    
    Args:
        reaction_smiles: Input reaction SMILES
        k: Number of precedents to retrieve
        reaction_type_threshold: Minimum representation for reaction types (0.0-1.0)
        mechanism_similarity_threshold: Minimum mechanism similarity (0.0-1.0)
        mechanism_weight: Weight for mechanism-enhanced similarity (0.0-1.0)
        show_debug: Show debug information about filtering
        
    Returns:
        Formatted recommendation results
    """
    print(f"🔍 Enhanced Cross-Family Search for: {reaction_smiles}")
    print(f"📊 Parameters:")
    print(f"   - k: {k}")
    print(f"   - reaction_type_threshold: {reaction_type_threshold:.1%}")
    print(f"   - mechanism_similarity_threshold: {mechanism_similarity_threshold:.1%}")
    print(f"   - mechanism_weight: {mechanism_weight:.1f}")
    print()
    
    # Configure relaxation parameters for cross-family search with filtering
    relax_config = {
        'reaction_type_threshold': reaction_type_threshold,
        'mechanism_similarity_threshold': mechanism_similarity_threshold,
        'mechanism_weight': mechanism_weight,
        'debug_timing': show_debug,
        'use_drfp': True,  # Enable DRFP for better similarity
    }
    
    try:
        # Get recommendations with cross-family search enabled
        result = chem.recommend.conditions(
            reaction=reaction_smiles,
            k=k,
            search_all_families=True,
            relax=relax_config,
            rerank_strategy='none',  # Cross-family search uses mechanism reranking instead
            filter_unknown_reagents=True  # Remove precedents with unknown reagents
        )
        
        # Format output for readability
        formatted_result = ensure_standard_output(
            result,
            default_model="Enhanced-Cross-Family-Search",
            fallback_reaction_smiles=reaction_smiles
        )
        
        # Add cross-family quality metrics
        _add_cross_family_metrics(formatted_result, result)
        
        return formatted_result
        
    except Exception as e:
        return {
            'error': str(e),
            'input_reaction': reaction_smiles,
            'search_type': 'enhanced_cross_family'
        }


def _add_cross_family_metrics(formatted_result: dict, raw_result: dict) -> None:
    """Add cross-family search quality metrics to the formatted result."""
    
    # Extract precedents and metadata
    recommendations = formatted_result.get('recommendations', [])
    meta = formatted_result.setdefault('meta', {})
    
    if not recommendations:
        return
    
    # Analyze reaction family diversity and compatibility
    family_counts = {}
    mechanism_scores = []
    compatibility_stats = {
        'compatible': 0,
        'moderate': 0, 
        'incompatible': 0,
        'well_represented': 0,
        'underrepresented': 0
    }
    
    for rec in recommendations:
        # Extract family from conditions or metadata
        family = 'Unknown'
        cross_family_meta = {}
        
        if 'conditions' in rec:
            cond = rec['conditions']
            if isinstance(cond, dict):
                family = cond.get('reaction_family', cond.get('rxn_type', 'Unknown'))
                cross_family_meta = cond.get('cross_family_metadata', {})
        
        family_counts[family] = family_counts.get(family, 0) + 1
        
        # Extract mechanism similarity if available
        mech_sim = cross_family_meta.get('mechanism_similarity')
        if mech_sim is not None:
            mechanism_scores.append(mech_sim)
        
        # Count compatibility statuses
        mech_status = cross_family_meta.get('mechanism_status', 'unknown')
        if mech_status in compatibility_stats:
            compatibility_stats[mech_status] += 1
            
        type_status = cross_family_meta.get('reaction_type_status', 'unknown') 
        if type_status in compatibility_stats:
            compatibility_stats[type_status] += 1
    
    # Calculate diversity metrics
    total_recs = len(recommendations)
    family_diversity = len(family_counts) / max(total_recs, 1)
    
    # Add metrics to meta section
    cross_family_metrics = {
        'total_recommendations': total_recs,
        'family_diversity_score': round(family_diversity, 3),
        'unique_families': len(family_counts),
        'family_distribution': family_counts,
        'compatibility_breakdown': {
            'mechanism_compatibility': {
                'compatible': compatibility_stats['compatible'],
                'moderate': compatibility_stats['moderate'],
                'incompatible': compatibility_stats['incompatible']
            },
            'reaction_type_representation': {
                'well_represented': compatibility_stats['well_represented'],
                'underrepresented': compatibility_stats['underrepresented']
            }
        }
    }
    
    if mechanism_scores:
        cross_family_metrics.update({
            'avg_mechanism_similarity': round(sum(mechanism_scores) / len(mechanism_scores), 3),
            'min_mechanism_similarity': round(min(mechanism_scores), 3),
            'max_mechanism_similarity': round(max(mechanism_scores), 3)
        })
    
    meta['cross_family_metrics'] = cross_family_metrics


def main():
    """Main CLI interface."""
    
    if len(sys.argv) < 2:
        print("Usage: python enhanced_cross_family_cli.py <reaction_smiles> [options]")
        print()
        print("Options:")
        print("  --k <int>                              Number of precedents (default: 50)")
        print("  --reaction-type-threshold <float>      Min reaction type representation (default: 0.15)")
        print("  --mechanism-threshold <float>          Min mechanism similarity (default: 0.4)")
        print("  --mechanism-weight <float>             Mechanism similarity weight (default: 0.3)")
        print("  --debug                                Show debug information")
        print("  --json                                 Output JSON format")
        print()
        print("Example:")
        print("  python enhanced_cross_family_cli.py 'CC(C)c1ccc(Br)cc1.Nc2ccccc2>>CC(C)c1ccc(Nc2ccccc2)cc1' --k 25 --mechanism-threshold 0.5")
        sys.exit(1)
    
    reaction_smiles = sys.argv[1]
    
    # Parse command line arguments
    k = 50
    reaction_type_threshold = 0.15
    mechanism_similarity_threshold = 0.4
    mechanism_weight = 0.3
    show_debug = False
    json_output = False
    
    i = 2
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '--k' and i + 1 < len(sys.argv):
            k = int(sys.argv[i + 1])
            i += 2
        elif arg == '--reaction-type-threshold' and i + 1 < len(sys.argv):
            reaction_type_threshold = float(sys.argv[i + 1])
            i += 2
        elif arg == '--mechanism-threshold' and i + 1 < len(sys.argv):
            mechanism_similarity_threshold = float(sys.argv[i + 1])
            i += 2
        elif arg == '--mechanism-weight' and i + 1 < len(sys.argv):
            mechanism_weight = float(sys.argv[i + 1])
            i += 2
        elif arg == '--debug':
            show_debug = True
            i += 1
        elif arg == '--json':
            json_output = True
            i += 1
        else:
            print(f"Unknown argument: {arg}")
            sys.exit(1)
    
    # Validate parameters
    if not (0.0 <= reaction_type_threshold <= 1.0):
        print(f"Error: reaction_type_threshold must be between 0.0 and 1.0, got {reaction_type_threshold}")
        sys.exit(1)
    
    if not (0.0 <= mechanism_similarity_threshold <= 1.0):
        print(f"Error: mechanism_similarity_threshold must be between 0.0 and 1.0, got {mechanism_similarity_threshold}")
        sys.exit(1)
    
    if not (0.0 <= mechanism_weight <= 1.0):
        print(f"Error: mechanism_weight must be between 0.0 and 1.0, got {mechanism_weight}")
        sys.exit(1)
    
    # Get enhanced recommendations
    result = enhanced_cross_family_recommend(
        reaction_smiles=reaction_smiles,
        k=k,
        reaction_type_threshold=reaction_type_threshold,
        mechanism_similarity_threshold=mechanism_similarity_threshold,
        mechanism_weight=mechanism_weight,
        show_debug=show_debug
    )
    
    # Output results
    if json_output:
        print(json.dumps(result, indent=2))
    else:
        # Pretty print key information
        if 'error' in result:
            print(f"❌ Error: {result['error']}")
        else:
            print("✅ Enhanced Cross-Family Recommendations:")
            print()
            
            # Show cross-family metrics
            metrics = result.get('meta', {}).get('cross_family_metrics', {})
            if metrics:
                print(f"📈 Quality Metrics:")
                print(f"   - Total recommendations: {metrics.get('total_recommendations', 0)}")
                print(f"   - Family diversity: {metrics.get('family_diversity_score', 0):.1%}")
                print(f"   - Unique families: {metrics.get('unique_families', 0)}")
                
                if 'avg_mechanism_similarity' in metrics:
                    print(f"   - Avg mechanism similarity: {metrics['avg_mechanism_similarity']:.1%}")
                
                # Show compatibility breakdown
                compat = metrics.get('compatibility_breakdown', {})
                if compat:
                    mech_compat = compat.get('mechanism_compatibility', {})
                    type_repr = compat.get('reaction_type_representation', {})
                    
                    print(f"   - Mechanism compatibility:")
                    print(f"     • Compatible: {mech_compat.get('compatible', 0)}")
                    print(f"     • Moderate: {mech_compat.get('moderate', 0)}")
                    print(f"     • Incompatible: {mech_compat.get('incompatible', 0)}")
                    
                    print(f"   - Reaction type representation:")
                    print(f"     • Well represented: {type_repr.get('well_represented', 0)}")
                    print(f"     • Underrepresented: {type_repr.get('underrepresented', 0)}")
                
                family_dist = metrics.get('family_distribution', {})
                if family_dist:
                    print(f"   - Family distribution: {dict(sorted(family_dist.items(), key=lambda x: x[1], reverse=True))}")
                print()
            
            # Show top recommendations with compatibility status
            recommendations = result.get('recommendations', [])[:5]  # Show top 5
            for i, rec in enumerate(recommendations, 1):
                print(f"{i}. {rec.get('core', 'Unknown')}")
                if 'conditions' in rec and isinstance(rec['conditions'], dict):
                    cond = rec['conditions']
                    family = cond.get('reaction_family', cond.get('rxn_type', 'Unknown'))
                    sim = cond.get('similarity', 'N/A')
                    
                    # Get cross-family metadata
                    cross_meta = cond.get('cross_family_metadata', {})
                    mech_sim = cross_meta.get('mechanism_similarity', 'N/A')
                    mech_status = cross_meta.get('mechanism_status', 'unknown')
                    type_status = cross_meta.get('reaction_type_status', 'unknown')
                    
                    print(f"   Family: {family}")
                    print(f"   Similarity: {sim}")
                    if mech_sim != 'N/A':
                        print(f"   Mechanism Similarity: {mech_sim:.1%}")
                    
                    # Add compatibility status indicators
                    status_icons = {
                        'compatible': '✅',
                        'moderate': '⚠️',
                        'incompatible': '❌',
                        'well_represented': '📊',
                        'underrepresented': '📉'
                    }
                    
                    mech_icon = status_icons.get(mech_status, '❓')
                    type_icon = status_icons.get(type_status, '❓')
                    
                    print(f"   Status: {mech_icon} {mech_status.title()} mechanism, {type_icon} {type_status.replace('_', ' ')}")
                print()


if __name__ == "__main__":
    main()
"""
Compare V1 vs V2 Prompt Performance

Tests the optimized V2 prompt against the original V1 prompt to measure:
- Latency reduction (target: 30%)
- Token usage reduction
- Quality maintenance (correctness, explanation quality)
- Confidence accuracy improvement

Usage:
    python test_v1_vs_v2_comparison.py

Author: Condition-Agent Team
Date: October 12, 2025
"""

import sys
import json
import time
from pathlib import Path
from typing import Dict, List
from datetime import datetime

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from llmtools.clients import LLMClient
from llmtools.recommendation_llm import synthesize_recommendations_llm


# Use subset of benchmark reactions for faster testing
TEST_REACTIONS = [
    {
        "id": 1,
        "name": "Simple Suzuki - High Consensus",
        "reaction_smiles": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
        "ml_result": {
            "meta": {"model": "ml_drfp"},
            "recommended_conditions": [{
                "catalyst": "Pd(PPh3)4",
                "solvent": "toluene",
                "temperature": "80°C",
                "base": "K3PO4",
                "similarity": 0.90
            }]
        },
        "rule_result": {
            "meta": {"model": "rule_scdb"},
            "recommended_conditions": [{
                "catalyst": "Pd(PPh3)4",
                "solvent": "toluene",
                "temperature": "80°C",
                "base": "K3PO4"
            }]
        },
        "protocol_result": {
            "meta": {"model": "protocol_literature"},
            "recommended_conditions": [{
                "catalyst": "Pd(PPh3)4",
                "solvent": "toluene",
                "temperature": "80°C",
                "base": "K3PO4",
                "time": "12-16h"
            }]
        },
        "expected_confidence": "high"
    },
    {
        "id": 2,
        "name": "Electron-Poor Suzuki - Medium Consensus",
        "reaction_smiles": "Brc1ccc([N+](=O)[O-])cc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccc([N+](=O)[O-])cc2)cc1",
        "ml_result": {
            "meta": {"model": "ml_drfp"},
            "recommended_conditions": [{
                "catalyst": "Pd(PPh3)4",
                "solvent": "toluene",
                "temperature": "80°C",
                "base": "K3PO4",
                "similarity": 0.75
            }]
        },
        "rule_result": {
            "meta": {"model": "rule_scdb"},
            "recommended_conditions": [{
                "catalyst": "Pd(PPh3)4",
                "solvent": "toluene",
                "temperature": "80°C",
                "base": "K3PO4"
            }]
        },
        "protocol_result": {
            "meta": {"model": "protocol_literature"},
            "recommended_conditions": [{
                "catalyst": "Pd(dppf)Cl2",
                "solvent": "DMF",
                "temperature": "100°C",
                "base": "K3PO4",
                "time": "12-16h",
                "note": "Electron-rich ligand for electron-poor substrate"
            }]
        },
        "constraints": {"scale": "multigram", "cost": "low"},
        "expected_confidence": "medium"
    },
    {
        "id": 3,
        "name": "Buchwald-Hartwig - Medium Consensus",
        "reaction_smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
        "ml_result": {
            "meta": {"model": "ml_drfp"},
            "recommended_conditions": [{
                "catalyst": "Pd2(dba)3",
                "ligand": "BINAP",
                "solvent": "toluene",
                "temperature": "100°C",
                "base": "NaOtBu",
                "similarity": 0.78
            }]
        },
        "rule_result": {
            "meta": {"model": "rule_scdb"},
            "recommended_conditions": [{
                "catalyst": "Pd(OAc)2",
                "ligand": "XPhos",
                "solvent": "dioxane",
                "temperature": "90°C",
                "base": "Cs2CO3"
            }]
        },
        "protocol_result": {
            "meta": {"model": "protocol_literature"},
            "recommended_conditions": [{
                "catalyst": "Pd2(dba)3",
                "ligand": "BINAP",
                "solvent": "toluene",
                "temperature": "110°C",
                "base": "NaOtBu",
                "time": "16-24h"
            }]
        },
        "expected_confidence": "medium"
    },
    {
        "id": 4,
        "name": "Novel Substrate - Low Consensus Expected",
        "reaction_smiles": "Brc1c(OC)cc(OC)c(OC)c1.c1ccc(B(O)O)cc1>>c1ccc(-c2c(OC)cc(OC)c(OC)c2)cc1",
        "ml_result": {
            "meta": {"model": "ml_drfp"},
            "recommended_conditions": [{
                "catalyst": "Pd(PPh3)4",
                "solvent": "toluene",
                "temperature": "80°C",
                "base": "K3PO4",
                "similarity": 0.55  # Low similarity
            }]
        },
        "rule_result": {
            "meta": {"model": "rule_scdb"},
            "recommended_conditions": [{
                "catalyst": "Pd(PPh3)4",
                "solvent": "toluene",
                "temperature": "80°C",
                "base": "K3PO4"
            }]
        },
        "protocol_result": {
            "meta": {"model": "protocol_literature"},
            "recommended_conditions": [{
                "catalyst": "Pd(dppf)Cl2",
                "solvent": "DME",
                "temperature": "100°C",
                "base": "K2CO3",
                "time": "16-24h",
                "note": "Electron-rich substrate adaptation"
            }]
        },
        "expected_confidence": "low"
    },
    {
        "id": 5,
        "name": "Heteroaryl - Chemistry Guidelines Test",
        "reaction_smiles": "Brc1ccccn1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccn2)cc1",
        "ml_result": {
            "meta": {"model": "ml_drfp"},
            "recommended_conditions": [{
                "catalyst": "Pd(PPh3)4",
                "solvent": "DME",
                "temperature": "90°C",
                "base": "Na2CO3",
                "similarity": 0.82
            }]
        },
        "rule_result": {
            "meta": {"model": "rule_scdb"},
            "recommended_conditions": [{
                "catalyst": "Pd(PPh3)4",
                "solvent": "toluene",
                "temperature": "80°C",
                "base": "K3PO4"
            }]
        },
        "protocol_result": {
            "meta": {"model": "protocol_literature"},
            "recommended_conditions": [{
                "catalyst": "Pd(PPh3)4",
                "solvent": "DME",
                "temperature": "100°C",
                "base": "Na2CO3",
                "time": "12-16h"
            }]
        },
        "expected_confidence": "medium"
    }
]


def run_synthesis_with_version(reaction: Dict, llm_client: LLMClient, version: str) -> Dict:
    """Run synthesis with specified prompt version."""
    
    result = synthesize_recommendations_llm(
        reaction_smiles=reaction['reaction_smiles'],
        ml_results=reaction['ml_result'],
        rule_results=reaction['rule_result'],
        protocol_results=reaction['protocol_result'],
        constraints=reaction.get('constraints'),
        llm_client=llm_client,
        prompt_version=version
    )
    
    return result


def evaluate_result(result: Dict, expected_confidence: str) -> Dict:
    """Evaluate synthesis result quality."""
    
    scores = {
        "has_synthesis": False,
        "has_consensus": False,
        "has_confidence": False,
        "confidence_match": False,
        "has_backups": False,
        "has_warnings": False,
        "has_source_comparison": False,
        "explanation_quality": 0
    }
    
    if result.get('status') != 'success':
        return scores
    
    synthesis = result.get('synthesis', {})
    
    scores['has_synthesis'] = bool(synthesis)
    scores['has_consensus'] = bool(synthesis.get('consensus_analysis'))
    scores['has_confidence'] = 'confidence_level' in synthesis
    scores['confidence_match'] = synthesis.get('confidence_level') == expected_confidence
    scores['has_backups'] = bool(synthesis.get('backup_conditions'))
    scores['has_warnings'] = bool(synthesis.get('warnings'))
    scores['has_source_comparison'] = bool(synthesis.get('source_comparison'))
    
    # Evaluate explanation quality (1-5)
    if synthesis.get('recommended_condition', {}).get('rationale'):
        scores['explanation_quality'] = 4
    if synthesis.get('backup_conditions'):
        scores['explanation_quality'] += 0.5
    if synthesis.get('warnings'):
        scores['explanation_quality'] += 0.5
    scores['explanation_quality'] = min(5, scores['explanation_quality'])
    
    return scores


def compare_versions():
    """Run comparison between V1 and V2 prompts."""
    
    print("="*80)
    print("V1 vs V2 PROMPT COMPARISON")
    print("="*80)
    print(f"Testing {len(TEST_REACTIONS)} reactions with both prompt versions\n")
    
    # Initialize LLM client
    try:
        llm_client = LLMClient(provider="aliyun", model="deepseek-v3.2-exp")
        print(f"✓ LLM client initialized: {llm_client.model}\n")
    except Exception as e:
        print(f"❌ Failed to initialize LLM client: {e}")
        return
    
    results = {
        'v1': [],
        'v2': []
    }
    
    # Run tests for each version
    for version in ['v1', 'v2']:
        print(f"\n{'='*80}")
        print(f"TESTING PROMPT {version.upper()}")
        print(f"{'='*80}\n")
        
        for i, reaction in enumerate(TEST_REACTIONS, 1):
            print(f"[{i}/{len(TEST_REACTIONS)}] {reaction['name']}")
            print(f"  Expected confidence: {reaction['expected_confidence']}")
            
            try:
                result = run_synthesis_with_version(reaction, llm_client, version)
                
                if result.get('status') == 'success':
                    synthesis = result['synthesis']
                    metadata = result.get('llm_metadata', {})
                    
                    print(f"  ✓ Success")
                    print(f"    Confidence: {synthesis.get('confidence_level')}")
                    print(f"    Latency: {metadata.get('latency_ms', 0)/1000:.1f}s")
                    print(f"    Tokens: {metadata.get('tokens', 0)}")
                    
                    # Evaluate quality
                    scores = evaluate_result(result, reaction['expected_confidence'])
                    print(f"    Quality: {scores['explanation_quality']:.1f}/5")
                    
                    results[version].append({
                        'reaction_id': reaction['id'],
                        'reaction_name': reaction['name'],
                        'result': result,
                        'scores': scores,
                        'latency_ms': metadata.get('latency_ms', 0),
                        'tokens': metadata.get('tokens', 0)
                    })
                else:
                    print(f"  ❌ Error: {result.get('error', 'Unknown')}")
                    results[version].append({
                        'reaction_id': reaction['id'],
                        'reaction_name': reaction['name'],
                        'result': result,
                        'scores': evaluate_result(result, reaction['expected_confidence']),
                        'latency_ms': 0,
                        'tokens': 0
                    })
                    
            except Exception as e:
                print(f"  ❌ Exception: {e}")
                results[version].append({
                    'reaction_id': reaction['id'],
                    'reaction_name': reaction['name'],
                    'error': str(e),
                    'scores': {},
                    'latency_ms': 0,
                    'tokens': 0
                })
    
    # Analyze results
    print(f"\n\n{'='*80}")
    print("ANALYSIS")
    print(f"{'='*80}\n")
    
    # Calculate averages
    for version in ['v1', 'v2']:
        version_results = results[version]
        successful = [r for r in version_results if r.get('result', {}).get('status') == 'success']
        
        if successful:
            avg_latency = sum(r['latency_ms'] for r in successful) / len(successful)
            avg_tokens = sum(r['tokens'] for r in successful) / len(successful)
            avg_quality = sum(r['scores'].get('explanation_quality', 0) for r in successful) / len(successful)
            confidence_matches = sum(1 for r in successful if r['scores'].get('confidence_match', False))
            
            print(f"{version.upper()} Results:")
            print(f"  Success rate: {len(successful)}/{len(version_results)} ({len(successful)/len(version_results)*100:.0f}%)")
            print(f"  Avg latency: {avg_latency/1000:.1f}s")
            print(f"  Avg tokens: {avg_tokens:.0f}")
            print(f"  Avg quality: {avg_quality:.2f}/5.00")
            print(f"  Confidence accuracy: {confidence_matches}/{len(successful)} ({confidence_matches/len(successful)*100:.0f}%)")
            print()
    
    # Calculate improvement
    v1_successful = [r for r in results['v1'] if r.get('result', {}).get('status') == 'success']
    v2_successful = [r for r in results['v2'] if r.get('result', {}).get('status') == 'success']
    
    if v1_successful and v2_successful:
        v1_avg_latency = sum(r['latency_ms'] for r in v1_successful) / len(v1_successful)
        v2_avg_latency = sum(r['latency_ms'] for r in v2_successful) / len(v2_successful)
        
        v1_avg_tokens = sum(r['tokens'] for r in v1_successful) / len(v1_successful)
        v2_avg_tokens = sum(r['tokens'] for r in v2_successful) / len(v2_successful)
        
        latency_reduction = (v1_avg_latency - v2_avg_latency) / v1_avg_latency * 100
        token_reduction = (v1_avg_tokens - v2_avg_tokens) / v1_avg_tokens * 100
        
        print(f"{'='*80}")
        print(f"IMPROVEMENTS (V2 vs V1)")
        print(f"{'='*80}")
        print(f"Latency reduction: {latency_reduction:+.1f}%")
        print(f"Token reduction: {token_reduction:+.1f}%")
        
        if latency_reduction >= 25:
            print(f"\n✅ TARGET MET: Latency reduced by ≥25% (actual: {latency_reduction:.1f}%)")
        else:
            print(f"\n⚠️  TARGET MISSED: Latency reduced by {latency_reduction:.1f}% (target: ≥25%)")
    
    # Save detailed results
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = f"v1_vs_v2_comparison_{timestamp}.json"
    
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    print(f"\n✅ Detailed results saved to: {output_file}")


if __name__ == "__main__":
    compare_versions()

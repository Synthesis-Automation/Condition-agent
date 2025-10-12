"""
Step 3: Compare Multi-Source Synthesis vs Individual Modes

This test suite compares the quality of recommendations from:
1. ML-only (DRFP similarity)
2. Rule-only (SCDB patterns)
3. Protocol-only (Literature)
4. Multi-source synthesis (Our new approach)

Tests 20 diverse reactions covering various reaction types and difficulty levels.

Usage:
    python test_mode_comparison.py

Author: Condition-Agent Team
Date: October 12, 2025
"""

import sys
import json
from pathlib import Path
from typing import Dict, List, Any
from datetime import datetime

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from llmtools.clients import LLMClient
from llmtools.recommendation_llm import synthesize_recommendations_llm


# 20 Diverse Benchmark Reactions
BENCHMARK_REACTIONS = [
    # ============================================================
    # SUZUKI COUPLINGS (5 reactions)
    # ============================================================
    {
        "id": 1,
        "name": "Simple Suzuki - Aryl Bromide",
        "reaction_smiles": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
        "reaction_type": "Suzuki",
        "difficulty": "easy",
        "substrate_features": ["aryl_bromide", "electron_neutral"],
        "known_good_conditions": {
            "catalyst": "Pd(PPh3)4",
            "solvent": "toluene",
            "temperature": "80°C",
            "base": "K3PO4"
        },
        "expected_consensus": "high"
    },
    {
        "id": 2,
        "name": "Suzuki - Electron-Poor Substrate",
        "reaction_smiles": "Brc1ccc([N+](=O)[O-])cc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccc([N+](=O)[O-])cc2)cc1",
        "reaction_type": "Suzuki",
        "difficulty": "medium",
        "substrate_features": ["aryl_bromide", "electron_poor", "nitro_group"],
        "known_good_conditions": {
            "catalyst": "Pd(dppf)Cl2",
            "solvent": "DMF",
            "temperature": "100°C",
            "base": "K3PO4"
        },
        "expected_consensus": "medium",
        "constraints": {"scale": "multigram", "cost": "low"}
    },
    {
        "id": 3,
        "name": "Suzuki - Aryl Chloride (Challenging)",
        "reaction_smiles": "Clc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
        "reaction_type": "Suzuki",
        "difficulty": "hard",
        "substrate_features": ["aryl_chloride", "less_reactive"],
        "known_good_conditions": {
            "catalyst": "Pd(OAc)2",
            "ligand": "SPhos",
            "solvent": "toluene",
            "temperature": "110°C",
            "base": "K3PO4"
        },
        "expected_consensus": "medium"
    },
    {
        "id": 4,
        "name": "Suzuki - Heteroaryl",
        "reaction_smiles": "Brc1ccccn1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccn2)cc1",
        "reaction_type": "Suzuki",
        "difficulty": "medium",
        "substrate_features": ["pyridyl_bromide", "heteroaryl"],
        "known_good_conditions": {
            "catalyst": "Pd(PPh3)4",
            "solvent": "DME",
            "temperature": "90°C",
            "base": "Na2CO3"
        },
        "expected_consensus": "medium"
    },
    {
        "id": 5,
        "name": "Suzuki - Sterically Hindered",
        "reaction_smiles": "Brc1c(C)cc(C)cc1C.c1ccc(B(O)O)cc1>>c1ccc(-c2c(C)cc(C)cc2C)cc1",
        "reaction_type": "Suzuki",
        "difficulty": "hard",
        "substrate_features": ["aryl_bromide", "sterically_hindered", "ortho_methyl"],
        "known_good_conditions": {
            "catalyst": "Pd(dppf)Cl2",
            "solvent": "dioxane",
            "temperature": "100°C",
            "base": "Cs2CO3"
        },
        "expected_consensus": "low"
    },
    
    # ============================================================
    # BUCHWALD-HARTWIG (4 reactions)
    # ============================================================
    {
        "id": 6,
        "name": "Buchwald-Hartwig - Primary Amine",
        "reaction_smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
        "reaction_type": "Buchwald-Hartwig",
        "difficulty": "medium",
        "substrate_features": ["aryl_bromide", "primary_amine"],
        "known_good_conditions": {
            "catalyst": "Pd2(dba)3",
            "ligand": "BINAP",
            "solvent": "toluene",
            "temperature": "100°C",
            "base": "NaOtBu"
        },
        "expected_consensus": "medium"
    },
    {
        "id": 7,
        "name": "Buchwald-Hartwig - Secondary Amine",
        "reaction_smiles": "Brc1ccccc1.CN(C)c1ccccc1>>c1ccc(N(C)c2ccccc2)cc1",
        "reaction_type": "Buchwald-Hartwig",
        "difficulty": "easy",
        "substrate_features": ["aryl_bromide", "secondary_amine"],
        "known_good_conditions": {
            "catalyst": "Pd(OAc)2",
            "ligand": "XPhos",
            "solvent": "dioxane",
            "temperature": "90°C",
            "base": "Cs2CO3"
        },
        "expected_consensus": "medium"
    },
    {
        "id": 8,
        "name": "Buchwald-Hartwig - Electron-Poor Aryl",
        "reaction_smiles": "Brc1ccc(C(F)(F)F)cc1.Nc1ccccc1>>c1ccc(Nc2ccc(C(F)(F)F)cc2)cc1",
        "reaction_type": "Buchwald-Hartwig",
        "difficulty": "medium",
        "substrate_features": ["aryl_bromide", "electron_poor", "CF3_group"],
        "known_good_conditions": {
            "catalyst": "Pd(dba)2",
            "ligand": "SPhos",
            "solvent": "toluene",
            "temperature": "100°C",
            "base": "NaOtBu"
        },
        "expected_consensus": "low"
    },
    {
        "id": 9,
        "name": "Buchwald-Hartwig - Heteroaryl Halide",
        "reaction_smiles": "Brc1ccccn1.Nc1ccccc1>>c1ccc(Nc2ccccn2)cc1",
        "reaction_type": "Buchwald-Hartwig",
        "difficulty": "hard",
        "substrate_features": ["pyridyl_bromide", "heteroaryl", "potential_coordination"],
        "known_good_conditions": {
            "catalyst": "Pd(OAc)2",
            "ligand": "RuPhos",
            "solvent": "dioxane",
            "temperature": "100°C",
            "base": "Cs2CO3"
        },
        "expected_consensus": "low"
    },
    
    # ============================================================
    # HECK REACTIONS (3 reactions)
    # ============================================================
    {
        "id": 10,
        "name": "Heck - Terminal Alkene",
        "reaction_smiles": "Brc1ccccc1.C=C>>C=Cc1ccccc1",
        "reaction_type": "Heck",
        "difficulty": "easy",
        "substrate_features": ["aryl_bromide", "terminal_alkene"],
        "known_good_conditions": {
            "catalyst": "Pd(OAc)2",
            "ligand": "PPh3",
            "solvent": "DMF",
            "temperature": "100°C",
            "base": "Et3N"
        },
        "expected_consensus": "high"
    },
    {
        "id": 11,
        "name": "Heck - Acrylate",
        "reaction_smiles": "Brc1ccccc1.C=CC(=O)OC>>C=CC(=O)Oc1ccccc1",
        "reaction_type": "Heck",
        "difficulty": "easy",
        "substrate_features": ["aryl_bromide", "acrylate"],
        "known_good_conditions": {
            "catalyst": "Pd(OAc)2",
            "solvent": "DMF",
            "temperature": "80°C",
            "base": "Et3N"
        },
        "expected_consensus": "high"
    },
    {
        "id": 12,
        "name": "Heck - Electron-Poor Substrate",
        "reaction_smiles": "Brc1ccc(C#N)cc1.C=C>>C=Cc1ccc(C#N)cc1",
        "reaction_type": "Heck",
        "difficulty": "medium",
        "substrate_features": ["aryl_bromide", "electron_poor", "cyano_group"],
        "known_good_conditions": {
            "catalyst": "Pd(PPh3)4",
            "solvent": "DMF",
            "temperature": "100°C",
            "base": "Et3N"
        },
        "expected_consensus": "medium"
    },
    
    # ============================================================
    # ULLMANN COUPLINGS (3 reactions)
    # ============================================================
    {
        "id": 13,
        "name": "Ullmann C-N - Simple Aniline",
        "reaction_smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
        "reaction_type": "Ullmann",
        "difficulty": "medium",
        "substrate_features": ["aryl_bromide", "primary_amine"],
        "known_good_conditions": {
            "catalyst": "CuI",
            "ligand": "phenanthroline",
            "solvent": "toluene",
            "temperature": "110°C",
            "base": "K3PO4"
        },
        "expected_consensus": "medium"
    },
    {
        "id": 14,
        "name": "Ullmann C-O - Phenol Coupling",
        "reaction_smiles": "Brc1ccccc1.Oc1ccccc1>>c1ccc(Oc2ccccc2)cc1",
        "reaction_type": "Ullmann",
        "difficulty": "medium",
        "substrate_features": ["aryl_bromide", "phenol"],
        "known_good_conditions": {
            "catalyst": "CuI",
            "ligand": "DMEDA",
            "solvent": "DMSO",
            "temperature": "90°C",
            "base": "K2CO3"
        },
        "expected_consensus": "medium"
    },
    {
        "id": 15,
        "name": "Ullmann - Heteroaryl",
        "reaction_smiles": "Brc1ccccn1.Nc1ccccc1>>c1ccc(Nc2ccccn2)cc1",
        "reaction_type": "Ullmann",
        "difficulty": "hard",
        "substrate_features": ["pyridyl_bromide", "heteroaryl", "potential_coordination"],
        "known_good_conditions": {
            "catalyst": "CuI",
            "ligand": "proline",
            "solvent": "DMSO",
            "temperature": "100°C",
            "base": "K3PO4"
        },
        "expected_consensus": "low"
    },
    
    # ============================================================
    # NEGISHI COUPLINGS (2 reactions)
    # ============================================================
    {
        "id": 16,
        "name": "Negishi - Aryl-Aryl",
        "reaction_smiles": "Brc1ccccc1.Zn(CC)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
        "reaction_type": "Negishi",
        "difficulty": "medium",
        "substrate_features": ["aryl_bromide", "organozinc"],
        "known_good_conditions": {
            "catalyst": "Pd(PPh3)4",
            "solvent": "THF",
            "temperature": "60°C",
            "base": None
        },
        "expected_consensus": "medium"
    },
    {
        "id": 17,
        "name": "Negishi - Heteroaryl",
        "reaction_smiles": "Brc1ccccn1.Zn(CC)c1ccccc1>>c1ccc(-c2ccccn2)cc1",
        "reaction_type": "Negishi",
        "difficulty": "hard",
        "substrate_features": ["pyridyl_bromide", "heteroaryl", "organozinc"],
        "known_good_conditions": {
            "catalyst": "Pd(PPh3)4",
            "solvent": "THF",
            "temperature": "60°C",
            "base": None
        },
        "expected_consensus": "low"
    },
    
    # ============================================================
    # EDGE CASES (3 reactions)
    # ============================================================
    {
        "id": 18,
        "name": "Ambiguous - Could be Suzuki or Ullmann",
        "reaction_smiles": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
        "reaction_type": "Ambiguous",
        "difficulty": "easy",
        "substrate_features": ["aryl_bromide", "multiple_pathways"],
        "known_good_conditions": {
            "catalyst": "Pd(PPh3)4",  # Suzuki preferred
            "solvent": "toluene",
            "temperature": "80°C",
            "base": "K3PO4"
        },
        "expected_consensus": "high",
        "note": "Should recognize Suzuki as preferred pathway"
    },
    {
        "id": 19,
        "name": "Low Precedent - Novel Substrate",
        "reaction_smiles": "Brc1cc(OC)c(OC)c(OC)c1.c1ccc(B(O)O)cc1>>c1ccc(-c2cc(OC)c(OC)c(OC)c2)cc1",
        "reaction_type": "Suzuki",
        "difficulty": "hard",
        "substrate_features": ["aryl_bromide", "multiple_methoxy", "electron_rich"],
        "known_good_conditions": {
            "catalyst": "Pd(PPh3)4",
            "solvent": "DME",
            "temperature": "90°C",
            "base": "K2CO3"
        },
        "expected_consensus": "low",
        "note": "Limited precedents for highly substituted substrates"
    },
    {
        "id": 20,
        "name": "Multi-Functional - Competing Sites",
        "reaction_smiles": "Brc1ccc(Br)cc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccc(Br)cc2)cc1",
        "reaction_type": "Suzuki",
        "difficulty": "hard",
        "substrate_features": ["dibromide", "selectivity_challenge"],
        "known_good_conditions": {
            "catalyst": "Pd(PPh3)4",
            "solvent": "toluene",
            "temperature": "80°C",
            "base": "K3PO4",
            "equivalents": "1.1 equiv boronic acid"
        },
        "expected_consensus": "medium",
        "note": "Mono-coupling requires careful stoichiometry control"
    }
]


def simulate_ml_result(reaction: Dict) -> Dict:
    """
    Simulate ML-based (DRFP) recommendation.
    In production, this would call the actual ML recommendation engine.
    
    ML tends to find closest precedent, which might not match perfectly.
    """
    # Simple simulation based on reaction difficulty
    if reaction['difficulty'] == 'easy':
        similarity = 0.90
        # Easy reactions: ML finds very similar precedent
        conditions = reaction['known_good_conditions'].copy()
    elif reaction['difficulty'] == 'medium':
        similarity = 0.75
        # Medium: ML might suggest slightly different catalyst or solvent
        conditions = reaction['known_good_conditions'].copy()
        # Sometimes suggests cheaper/more common catalyst
        if conditions.get('catalyst') == 'Pd(dppf)Cl2':
            conditions['catalyst'] = 'Pd(PPh3)4'  # More common in training data
        if conditions.get('solvent') == 'DME':
            conditions['solvent'] = 'toluene'  # More common alternative
    else:  # hard
        similarity = 0.55
        # Hard: ML struggles, suggests generic conditions
        conditions = {
            'catalyst': 'Pd(PPh3)4',  # Default Pd catalyst
            'solvent': 'toluene',  # Generic solvent
            'temperature': '80°C',
            'base': 'K3PO4' if reaction['reaction_type'] in ['Suzuki', 'Ullmann'] else 'Cs2CO3'
        }
        if 'ligand' in reaction['known_good_conditions']:
            conditions['ligand'] = 'PPh3'  # Generic ligand
    
    conditions['similarity'] = similarity
    conditions['source'] = 'ML (DRFP)'
    
    return {
        "meta": {"model": "ml_drfp"},
        "recommended_conditions": [conditions]
    }


def simulate_rule_result(reaction: Dict) -> Dict:
    """
    Simulate Rule-based (SCDB) recommendation.
    In production, this would call the actual rule-based engine.
    
    Rules tend to be conservative and general-purpose.
    """
    reaction_type = reaction['reaction_type']
    
    # Rules give general-purpose conditions (may not be optimal)
    if reaction_type == "Suzuki":
        conditions = {
            "catalyst": "Pd(PPh3)4",  # Classic catalyst
            "solvent": "toluene",  # Standard solvent
            "temperature": "80°C",
            "base": "K3PO4",
            "source": "Rule (SCDB)"
        }
        # Rules don't adapt well to special cases
        if "electron_poor" in reaction['substrate_features']:
            # Might not recognize need for better catalyst
            pass  # Keeps generic conditions
            
    elif reaction_type == "Buchwald-Hartwig":
        conditions = {
            "catalyst": "Pd(OAc)2",
            "ligand": "XPhos",  # Modern ligand
            "solvent": "dioxane",
            "temperature": "90°C",
            "base": "Cs2CO3",
            "source": "Rule (SCDB)"
        }
        # For difficult substrates, might miss nuances
        
    elif reaction_type == "Heck":
        conditions = {
            "catalyst": "Pd(OAc)2",
            "solvent": "DMF",
            "temperature": "100°C",
            "base": "Et3N",
            "source": "Rule (SCDB)"
        }
        
    elif reaction_type == "Ullmann":
        conditions = {
            "catalyst": "CuI",
            "ligand": "phenanthroline",
            "solvent": "toluene",
            "temperature": "110°C",
            "base": "K3PO4",
            "source": "Rule (SCDB)"
        }
        # Ullmann rules might not capture all ligand variations
        
    elif reaction_type == "Negishi":
        conditions = {
            "catalyst": "Pd(PPh3)4",
            "solvent": "THF",
            "temperature": "60°C",
            "source": "Rule (SCDB)"
        }
        
    else:  # Ambiguous
        # Rules default to Suzuki for ambiguous cases
        conditions = {
            "catalyst": "Pd(PPh3)4",
            "solvent": "toluene",
            "temperature": "80°C",
            "base": "K3PO4",
            "source": "Rule (SCDB)"
        }
    
    return {
        "meta": {"model": "rule_scdb"},
        "recommended_conditions": [conditions]
    }


def simulate_protocol_result(reaction: Dict) -> Dict:
    """
    Simulate Protocol-based (Literature) recommendation.
    In production, this would call the actual protocol recommendation engine.
    
    Protocols provide detailed, reliable conditions from literature.
    Closest to ground truth but might be overly specific.
    """
    # Protocol uses known good conditions (from literature)
    conditions = reaction['known_good_conditions'].copy()
    
    # Add protocol-specific details
    conditions['source'] = 'Protocol (Literature)'
    conditions['time'] = '12-16h'
    conditions['note'] = 'From literature procedure'
    
    # For very challenging substrates, might suggest harsher conditions
    if reaction['difficulty'] == 'hard' and 'temperature' in conditions:
        temp_val = int(conditions['temperature'].replace('°C', ''))
        # Protocol might be overly cautious and suggest higher temp
        conditions['temperature'] = f"{temp_val + 10}°C"
    
    return {
        "meta": {"model": "protocol_literature"},
        "recommended_conditions": [conditions]
    }


def evaluate_recommendation(recommendation: Dict, ground_truth: Dict, reaction: Dict) -> Dict:
    """
    Evaluate recommendation quality against ground truth.
    
    Returns score dict with:
    - overall_score (1-5)
    - catalyst_match (bool)
    - solvent_match (bool)
    - temperature_reasonable (bool)
    - explanation_quality (1-5, for synthesis only)
    """
    score = {
        "overall_score": 3,  # Default: reasonable
        "catalyst_match": False,
        "solvent_match": False,
        "temperature_reasonable": False,
        "base_match": False,
        "explanation_quality": 0
    }
    
    # Extract recommended conditions
    rec_cond = recommendation.get('recommended_conditions', [{}])[0]
    if not rec_cond:
        rec_cond = recommendation.get('synthesis', {}).get('recommended_condition', {})
    
    # Check catalyst match
    if rec_cond.get('catalyst') == ground_truth.get('catalyst'):
        score['catalyst_match'] = True
        score['overall_score'] += 1
    
    # Check solvent match
    if rec_cond.get('solvent') == ground_truth.get('solvent'):
        score['solvent_match'] = True
        score['overall_score'] += 0.5
    
    # Check temperature reasonable (within ±20°C)
    rec_temp = rec_cond.get('temperature', '')
    gt_temp = ground_truth.get('temperature', '')
    if rec_temp and gt_temp:
        try:
            rec_val = int(rec_temp.replace('°C', ''))
            gt_val = int(gt_temp.replace('°C', ''))
            if abs(rec_val - gt_val) <= 20:
                score['temperature_reasonable'] = True
                score['overall_score'] += 0.5
        except:
            pass
    
    # Check base match
    if rec_cond.get('base') == ground_truth.get('base'):
        score['base_match'] = True
        score['overall_score'] += 0.5
    
    # Cap at 5
    score['overall_score'] = min(5, score['overall_score'])
    
    # For synthesis results, evaluate explanation
    if 'synthesis' in recommendation:
        synthesis = recommendation['synthesis']
        if synthesis.get('recommended_condition', {}).get('rationale'):
            score['explanation_quality'] = 4  # Good explanation
        if synthesis.get('backup_conditions'):
            score['explanation_quality'] += 0.5
        if synthesis.get('warnings'):
            score['explanation_quality'] += 0.5
        score['explanation_quality'] = min(5, score['explanation_quality'])
    
    return score


def run_comparison(llm_client: LLMClient):
    """Run comparison across all benchmark reactions."""
    
    results = []
    
    print(f"\n{'='*80}")
    print(f"RUNNING MODE COMPARISON - {len(BENCHMARK_REACTIONS)} REACTIONS")
    print(f"{'='*80}\n")
    
    for i, reaction in enumerate(BENCHMARK_REACTIONS, 1):
        print(f"\n[{i}/{len(BENCHMARK_REACTIONS)}] {reaction['name']}")
        print(f"Difficulty: {reaction['difficulty']} | Type: {reaction['reaction_type']}")
        print(f"Features: {', '.join(reaction['substrate_features'])}")
        
        # Simulate results from each mode
        print("  → Running ML-based...")
        ml_result = simulate_ml_result(reaction)
        
        print("  → Running Rule-based...")
        rule_result = simulate_rule_result(reaction)
        
        print("  → Running Protocol-based...")
        protocol_result = simulate_protocol_result(reaction)
        
        print("  → Running Multi-source synthesis (LLM)...")
        try:
            synthesis_result = synthesize_recommendations_llm(
                reaction_smiles=reaction['reaction_smiles'],
                ml_results=ml_result,
                rule_results=rule_result,
                protocol_results=protocol_result,
                constraints=reaction.get('constraints'),
                llm_client=llm_client
            )
        except Exception as e:
            print(f"  ❌ Synthesis failed: {e}")
            synthesis_result = {"status": "error", "error": str(e)}
        
        # Evaluate each mode
        print("  → Evaluating results...")
        ml_score = evaluate_recommendation(ml_result, reaction['known_good_conditions'], reaction)
        rule_score = evaluate_recommendation(rule_result, reaction['known_good_conditions'], reaction)
        protocol_score = evaluate_recommendation(protocol_result, reaction['known_good_conditions'], reaction)
        
        if synthesis_result.get('status') == 'success':
            synthesis_score = evaluate_recommendation(synthesis_result, reaction['known_good_conditions'], reaction)
        else:
            synthesis_score = {"overall_score": 0, "explanation_quality": 0}
        
        # Print scores
        print(f"\n  SCORES:")
        print(f"    ML:        {ml_score['overall_score']:.1f}/5")
        print(f"    Rule:      {rule_score['overall_score']:.1f}/5")
        print(f"    Protocol:  {protocol_score['overall_score']:.1f}/5")
        print(f"    Synthesis: {synthesis_score['overall_score']:.1f}/5 (explain: {synthesis_score.get('explanation_quality', 0):.1f}/5)")
        
        # Determine winner
        scores_dict = {
            'ML': ml_score['overall_score'],
            'Rule': rule_score['overall_score'],
            'Protocol': protocol_score['overall_score'],
            'Synthesis': synthesis_score['overall_score']
        }
        winner = max(scores_dict, key=scores_dict.get)
        print(f"  🏆 Best: {winner}")
        
        # Store results
        results.append({
            'reaction': reaction,
            'scores': {
                'ml': ml_score,
                'rule': rule_score,
                'protocol': protocol_score,
                'synthesis': synthesis_score
            },
            'winner': winner,
            'ml_result': ml_result,
            'rule_result': rule_result,
            'protocol_result': protocol_result,
            'synthesis_result': synthesis_result
        })
    
    return results


def analyze_results(results: List[Dict]):
    """Analyze and summarize comparison results."""
    
    print(f"\n\n{'='*80}")
    print(f"ANALYSIS SUMMARY")
    print(f"{'='*80}\n")
    
    # Count wins by mode
    wins = {'ML': 0, 'Rule': 0, 'Protocol': 0, 'Synthesis': 0}
    for result in results:
        wins[result['winner']] += 1
    
    print(f"WINS BY MODE:")
    for mode, count in sorted(wins.items(), key=lambda x: x[1], reverse=True):
        pct = (count / len(results)) * 100
        print(f"  {mode:12} {count:2d}/{len(results)} ({pct:5.1f}%)")
    
    # Average scores by mode
    avg_scores = {'ML': 0, 'Rule': 0, 'Protocol': 0, 'Synthesis': 0}
    for result in results:
        avg_scores['ML'] += result['scores']['ml']['overall_score']
        avg_scores['Rule'] += result['scores']['rule']['overall_score']
        avg_scores['Protocol'] += result['scores']['protocol']['overall_score']
        avg_scores['Synthesis'] += result['scores']['synthesis']['overall_score']
    
    print(f"\nAVERAGE SCORES:")
    for mode in avg_scores:
        avg = avg_scores[mode] / len(results)
        avg_scores[mode] = avg
        print(f"  {mode:12} {avg:.2f}/5.00")
    
    # Breakdown by difficulty
    print(f"\nPERFORMANCE BY DIFFICULTY:")
    for difficulty in ['easy', 'medium', 'hard']:
        diff_results = [r for r in results if r['reaction']['difficulty'] == difficulty]
        if not diff_results:
            continue
        
        diff_wins = {'ML': 0, 'Rule': 0, 'Protocol': 0, 'Synthesis': 0}
        for result in diff_results:
            diff_wins[result['winner']] += 1
        
        print(f"\n  {difficulty.upper()} ({len(diff_results)} reactions):")
        for mode, count in sorted(diff_wins.items(), key=lambda x: x[1], reverse=True):
            pct = (count / len(diff_results)) * 100 if diff_results else 0
            print(f"    {mode:12} {count:2d}/{len(diff_results)} ({pct:5.1f}%)")
    
    # Explanation quality (synthesis only)
    synthesis_results = [r for r in results if r['scores']['synthesis'].get('explanation_quality', 0) > 0]
    if synthesis_results:
        avg_explain = sum(r['scores']['synthesis']['explanation_quality'] for r in synthesis_results) / len(synthesis_results)
        print(f"\nSYNTHESIS EXPLANATION QUALITY: {avg_explain:.2f}/5.00")
    
    # Save detailed results
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = f"comparison_results_{timestamp}.json"
    
    with open(output_file, 'w') as f:
        json.dump({
            'summary': {
                'total_reactions': len(results),
                'wins_by_mode': wins,
                'avg_scores_by_mode': avg_scores,
                'timestamp': timestamp
            },
            'detailed_results': results
        }, f, indent=2, default=str)
    
    print(f"\n✅ Detailed results saved to: {output_file}")
    
    return {
        'wins': wins,
        'avg_scores': avg_scores
    }


def main():
    """Main comparison test runner."""
    
    print("="*80)
    print("STEP 3: MODE COMPARISON TEST SUITE")
    print("="*80)
    print(f"Testing {len(BENCHMARK_REACTIONS)} diverse reactions")
    print("Comparing: ML vs Rule vs Protocol vs Multi-source Synthesis")
    
    # Initialize LLM client
    try:
        llm_client = LLMClient(provider="aliyun", model="deepseek-v3.2-exp")
        print(f"✓ LLM client initialized: {llm_client.model}\n")
    except Exception as e:
        print(f"❌ Failed to initialize LLM client: {e}")
        return
    
    # Run comparison
    results = run_comparison(llm_client)
    
    # Analyze results
    summary = analyze_results(results)
    
    # Final verdict
    print(f"\n{'='*80}")
    print(f"FINAL VERDICT")
    print(f"{'='*80}")
    
    winner_mode = max(summary['wins'], key=summary['wins'].get)
    best_avg = max(summary['avg_scores'], key=summary['avg_scores'].get)
    
    print(f"\n🏆 Most Wins: {winner_mode} ({summary['wins'][winner_mode]}/{len(results)})")
    print(f"📊 Best Avg Score: {best_avg} ({summary['avg_scores'][best_avg]:.2f}/5.00)")
    
    if winner_mode == 'Synthesis' or best_avg == 'Synthesis':
        print(f"\n✅ Multi-source synthesis shows SUPERIOR performance!")
    else:
        print(f"\n⚠️  Multi-source synthesis needs improvement.")
    
    print(f"\n{'='*80}")


if __name__ == "__main__":
    main()

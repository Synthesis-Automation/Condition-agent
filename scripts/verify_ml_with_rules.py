"""
ML vs Rule-Based Recommendation Comparison
===============================def get_rule_based_recommendations(reaction_smiles: str,
                                   dataset_path: str = "data/reaction_dataset/C_N_Coupling_Pd.jsonl",
                                   n_precedents: int = 5) -> List[Dict[str, Any]]:========

This script compares ML-based yield predictions with rule-based precedent
recommendations to flag potential issues and validate predictions.

For each reaction:
1. Get ML prediction (DRFP + LightGBM)
2. Get rule-based recommendation (precedent k-NN)
3. Compare conditions and flag discrepancies

Flags raised when:
- Different catalyst cores (e.g., XPhos vs RuPhos)
- Different bases (e.g., NaOtBu vs K3PO4)
- Different solvents (e.g., Toluene vs Dioxane)
- ML yield differs from precedent yield by >20%
- ML recommends conditions with no precedents

Usage:
    python scripts/verify_ml_with_rules.py
"""

import sys
import json
import pickle
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from chemtools.ml.drfp_predictor import DRFPYieldPredictor
from chemtools.smiles import _split_reaction_smiles


def load_model(model_path: str) -> DRFPYieldPredictor:
    """Load trained DRFP model."""
    with open(model_path, 'rb') as f:
        model_dict = pickle.load(f)
    
    predictor = DRFPYieldPredictor(
        n_bits=model_dict['n_bits'],
        radius=model_dict['radius']
    )
    predictor.drfp_encoder = model_dict['drfp_encoder']
    predictor.model = model_dict['model']
    predictor.core_vocab = model_dict['core_vocab']
    predictor.base_vocab = model_dict['base_vocab']
    predictor.solvent_vocab = model_dict['solvent_vocab']
    
    return predictor


def get_ml_prediction(predictor: DRFPYieldPredictor,
                     reaction_smiles: str,
                     core: str,
                     base_uid: str,
                     solvent_uid: str,
                     T_C: float = 100.0,
                     time_h: float = 12.0) -> Optional[float]:
    """Get ML yield prediction."""
    import pandas as pd
    try:
        df = pd.DataFrame([{
            'reaction_smiles': reaction_smiles,
            'core': core,
            'base_uid': base_uid,
            'solvent_uid': solvent_uid,
            'T_C': T_C,
            'time_h': time_h
        }])
        predictions = predictor.predict(df)
        return predictions[0] if len(predictions) > 0 else None
    except Exception as e:
        print(f"  ML prediction failed: {e}")
        return None


def get_rule_based_recommendations(reaction_smiles: str,
                                   dataset_path: str = "data/reaction_dataset/Buchwald2021-2024.jsonl",
                                   n_precedents: int = 10) -> List[Dict]:
    """Get rule-based precedent recommendations.
    
    Returns list of precedents with their conditions and yields.
    """
    try:
        # Load precedent database
        precedents = []
        with open(dataset_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    rec = json.loads(line)
                    precedents.append(rec)
                except:
                    continue
        
        if not precedents:
            return []
        
        # Parse query reaction
        parts = _split_reaction_smiles(reaction_smiles)
        if len(parts) < 2:
            return []
        
        reactants_smiles = parts[0]
        products_smiles = parts[1] if len(parts) > 1 else ""
        
        # For each precedent, extract key info
        precedent_list = []
        for prec in precedents:
            rxn_smiles = prec.get('smiles', {})
            prec_reactants = rxn_smiles.get('reactants', '')
            
            if not prec_reactants:
                continue
            
            precedent_list.append({
                'reaction_id': prec.get('reaction_id', ''),
                'reactants': prec_reactants,
                'products': rxn_smiles.get('products', ''),
                'core': prec.get('condition_core', ''),
                'base_uid': _extract_base_uid(prec),
                'solvent_uid': _extract_solvent_uid(prec),
                'T_C': prec.get('conditions', {}).get('temperature_C', 100),
                'time_h': prec.get('conditions', {}).get('time_h', 12),
                'yield_pct': prec.get('conditions', {}).get('yield_pct', 0),
                'smiles': f"{prec_reactants}>>{rxn_smiles.get('products', '')}"
            })
        
        # Simple SMILES-based similarity (just check if reactants are similar)
        # In a real system, you'd use Tanimoto similarity
        similar_precedents = []
        for prec in precedent_list:
            # Very simple: check if any reactant SMILES overlaps
            if any(r in reaction_smiles for r in prec['reactants'].split('.')):
                similar_precedents.append(prec)
        
        # Return top N by yield
        similar_precedents.sort(key=lambda x: x['yield_pct'], reverse=True)
        return similar_precedents[:n_precedents]
        
    except Exception as e:
        print(f"  Rule-based recommendation failed: {e}")
        return []


def _extract_base_uid(record: Dict) -> str:
    """Extract base CAS number from reaction record."""
    reagents = record.get('reagents', [])
    for r in reagents:
        if r.get('role', '').upper() == 'BASE':
            # Use 'cas' not 'uid' - that's the actual field name!
            return r.get('cas', '')
    return ''


def _extract_solvent_uid(record: Dict) -> str:
    """Extract solvent CAS number from reaction record."""
    solvents = record.get('solvents', [])
    if solvents and isinstance(solvents, list):
        # Use 'cas' not 'uid' - that's the actual field name!
        return solvents[0].get('cas', '') if isinstance(solvents[0], dict) else ''
    return ''


def compare_conditions(ml_cond: Dict, rule_cond: Dict) -> Dict[str, bool]:
    """Compare ML and rule-based conditions, return flags for discrepancies."""
    # Safe extraction with defaults
    ml_temp = float(ml_cond.get('T_C', 100))
    rule_temp = float(rule_cond.get('T_C') or 100)
    ml_time = float(ml_cond.get('time_h', 12))
    rule_time = float(rule_cond.get('time_h') or 12)
    
    flags = {
        'core_mismatch': ml_cond.get('core', '') != rule_cond.get('core', ''),
        'base_mismatch': ml_cond.get('base_uid', '') != rule_cond.get('base_uid', ''),
        'solvent_mismatch': ml_cond.get('solvent_uid', '') != rule_cond.get('solvent_uid', ''),
        'temp_large_diff': abs(ml_temp - rule_temp) > 20,
        'time_large_diff': abs(ml_time - rule_time) > 6,
    }
    
    # Calculate yield difference if both available
    if 'predicted_yield' in ml_cond and 'yield_pct' in rule_cond:
        yield_diff = abs(ml_cond['predicted_yield'] - rule_cond['yield_pct'])
        flags['yield_large_diff'] = yield_diff > 20
    else:
        flags['yield_large_diff'] = False
    
    return flags


def verify_predictions(ml_results_path: str,
                      model_path: str,
                      dataset_path: str,
                      output_path: str):
    """Compare ML predictions with rule-based recommendations."""
    
    print("=" * 80)
    print("ML vs RULE-BASED VERIFICATION")
    print("=" * 80)
    
    # Load ML results
    print(f"\n[1/4] Loading ML predictions from {ml_results_path}...")
    with open(ml_results_path, 'r') as f:
        ml_results = json.load(f)
    print(f"  ✓ Loaded {len(ml_results)} ML predictions")
    
    # Load model
    print(f"\n[2/4] Loading trained model...")
    predictor = load_model(model_path)
    print(f"  ✓ Model loaded")
    
    # Filter to Buchwald-Hartwig and C-N coupling only (exclude Chan-Lam)
    buchwald_results = [r for r in ml_results if r['type'] in ['Buchwald-Hartwig', 'C-N Coupling']]
    print(f"\n[3/4] Filtering to Buchwald-Hartwig reactions...")
    print(f"  ✓ {len(buchwald_results)} reactions (excluded {len(ml_results) - len(buchwald_results)} Chan-Lam)")
    
    # Compare each reaction
    print(f"\n[4/4] Comparing ML vs rule-based recommendations...")
    comparisons = []
    
    for i, ml_res in enumerate(buchwald_results, 1):
        print(f"\n  [{i}/{len(buchwald_results)}] {ml_res['description'][:60]}...")
        
        # Get ML best condition
        if not ml_res.get('best_condition'):
            print("    ⚠ No ML prediction available")
            continue
        
        ml_best = ml_res['best_condition']
        
        # Get rule-based recommendations
        precedents = get_rule_based_recommendations(ml_res['smiles'], dataset_path, n_precedents=5)
        
        if not precedents:
            print("    ⚠ No precedents found - ML recommendation cannot be verified")
            comparisons.append({
                'reaction': ml_res['description'],
                'smiles': ml_res['smiles'],
                'ml_yield': ml_best['predicted_yield'],
                'ml_core': ml_best['core'],
                'precedent_count': 0,
                'flags': {'no_precedents': True},
                'status': 'UNVERIFIABLE'
            })
            continue
        
        # Compare with top precedent
        top_precedent = precedents[0]
        flags = compare_conditions(ml_best, top_precedent)
        
        # Count total flags
        flag_count = sum(1 for v in flags.values() if v)
        
        # Determine status
        if flag_count == 0:
            status = 'VERIFIED ✓'
        elif flag_count <= 2:
            status = 'MINOR DISCREPANCY ⚠'
        else:
            status = 'MAJOR DISCREPANCY ⚠⚠'
        
        print(f"    Status: {status}")
        print(f"    ML: {ml_best['core']}, {ml_best['predicted_yield']:.1f}%")
        print(f"    Top Precedent: {top_precedent['core']}, {top_precedent['yield_pct']:.1f}%")
        
        if flags['core_mismatch']:
            print(f"    🚩 Core mismatch: ML={ml_best['core']} vs Precedent={top_precedent['core']}")
        if flags['base_mismatch']:
            print(f"    🚩 Base mismatch: ML={ml_best['base_uid']} vs Precedent={top_precedent['base_uid']}")
        if flags['solvent_mismatch']:
            print(f"    🚩 Solvent mismatch: ML={ml_best['solvent_uid']} vs Precedent={top_precedent['solvent_uid']}")
        if flags['yield_large_diff']:
            print(f"    🚩 Large yield difference: {abs(ml_best['predicted_yield'] - top_precedent['yield_pct']):.1f}%")
        
        comparisons.append({
            'reaction': ml_res['description'],
            'smiles': ml_res['smiles'],
            'ml_yield': ml_best['predicted_yield'],
            'ml_core': ml_best['core'],
            'ml_base': ml_best['base_uid'],
            'ml_solvent': ml_best['solvent_uid'],
            'precedent_yield': top_precedent['yield_pct'],
            'precedent_core': top_precedent['core'],
            'precedent_base': top_precedent['base_uid'],
            'precedent_solvent': top_precedent['solvent_uid'],
            'precedent_count': len(precedents),
            'flags': flags,
            'flag_count': flag_count,
            'status': status
        })
    
    # Generate summary statistics
    print("\n" + "=" * 80)
    print("VERIFICATION SUMMARY")
    print("=" * 80)
    
    verified = [c for c in comparisons if 'VERIFIED' in c['status']]
    minor = [c for c in comparisons if 'MINOR' in c['status']]
    major = [c for c in comparisons if 'MAJOR' in c['status']]
    unverifiable = [c for c in comparisons if 'UNVERIFIABLE' in c['status']]
    
    print(f"\nTotal Reactions: {len(comparisons)}")
    print(f"  ✓ Verified (exact match): {len(verified)} ({len(verified)/len(comparisons)*100:.1f}%)")
    print(f"  ⚠ Minor Discrepancy (1-2 flags): {len(minor)} ({len(minor)/len(comparisons)*100:.1f}%)")
    print(f"  ⚠⚠ Major Discrepancy (3+ flags): {len(major)} ({len(major)/len(comparisons)*100:.1f}%)")
    print(f"  ? Unverifiable (no precedents): {len(unverifiable)} ({len(unverifiable)/len(comparisons)*100:.1f}%)")
    
    # Flag statistics
    print(f"\nMost Common Discrepancies:")
    flag_counts = defaultdict(int)
    for c in comparisons:
        for flag_name, flag_value in c.get('flags', {}).items():
            if flag_value:
                flag_counts[flag_name] += 1
    
    for flag_name, count in sorted(flag_counts.items(), key=lambda x: x[1], reverse=True):
        print(f"  - {flag_name}: {count} reactions")
    
    # Save detailed results
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w') as f:
        json.dump({
            'summary': {
                'total': len(comparisons),
                'verified': len(verified),
                'minor_discrepancy': len(minor),
                'major_discrepancy': len(major),
                'unverifiable': len(unverifiable)
            },
            'flag_counts': dict(flag_counts),
            'comparisons': comparisons
        }, f, indent=2)
    
    print(f"\n✓ Detailed results saved to: {output_file}")
    
    # Generate markdown report
    report_path = output_file.parent / "Verification_Report.md"
    generate_markdown_report(comparisons, verified, minor, major, unverifiable, flag_counts, report_path)
    print(f"✓ Summary report saved to: {report_path}")
    
    print("\n" + "=" * 80)


def generate_markdown_report(comparisons, verified, minor, major, unverifiable, flag_counts, output_path):
    """Generate markdown verification report."""
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("# ML vs Rule-Based Verification Report\n\n")
        f.write("**Purpose:** Validate ML predictions against precedent-based recommendations\n\n")
        f.write("---\n\n")
        
        # Executive Summary
        f.write("## Executive Summary\n\n")
        total = len(comparisons)
        f.write(f"- **Total Reactions Verified:** {total}\n")
        f.write(f"- **Fully Verified:** {len(verified)} ({len(verified)/total*100:.1f}%)\n")
        f.write(f"- **Minor Discrepancies:** {len(minor)} ({len(minor)/total*100:.1f}%)\n")
        f.write(f"- **Major Discrepancies:** {len(major)} ({len(major)/total*100:.1f}%)\n")
        f.write(f"- **Unverifiable:** {len(unverifiable)} ({len(unverifiable)/total*100:.1f}%)\n\n")
        
        # Flag Statistics
        f.write("## Common Discrepancy Types\n\n")
        f.write("| Discrepancy Type | Count | % of Total |\n")
        f.write("|-----------------|-------|------------|\n")
        for flag_name, count in sorted(flag_counts.items(), key=lambda x: x[1], reverse=True):
            pct = count / total * 100
            f.write(f"| {flag_name.replace('_', ' ').title()} | {count} | {pct:.1f}% |\n")
        f.write("\n")
        
        # Major Discrepancies (need attention)
        if major:
            f.write("## ⚠⚠ Major Discrepancies (Requires Review)\n\n")
            for c in major:
                f.write(f"### {c['reaction']}\n\n")
                f.write(f"**SMILES:** `{c['smiles']}`\n\n")
                f.write(f"**ML Prediction:**\n")
                f.write(f"- Core: {c['ml_core']}\n")
                f.write(f"- Predicted Yield: {c['ml_yield']:.1f}%\n\n")
                f.write(f"**Top Precedent:**\n")
                f.write(f"- Core: {c['precedent_core']}\n")
                f.write(f"- Actual Yield: {c['precedent_yield']:.1f}%\n")
                f.write(f"- Precedent Count: {c['precedent_count']}\n\n")
                f.write(f"**Flags:**\n")
                for flag_name, flag_value in c['flags'].items():
                    if flag_value:
                        f.write(f"- 🚩 {flag_name.replace('_', ' ').title()}\n")
                f.write("\n---\n\n")
        
        # Minor Discrepancies
        if minor:
            f.write("## ⚠ Minor Discrepancies\n\n")
            f.write("| Reaction | ML Yield | Precedent Yield | Flags |\n")
            f.write("|----------|----------|----------------|-------|\n")
            for c in minor:
                desc = c['reaction'][:40] + "..." if len(c['reaction']) > 40 else c['reaction']
                flag_str = ", ".join(k.replace('_', ' ') for k, v in c['flags'].items() if v)
                f.write(f"| {desc} | {c['ml_yield']:.1f}% | {c['precedent_yield']:.1f}% | {flag_str} |\n")
            f.write("\n")
        
        # Verified (agreement)
        if verified:
            f.write("## ✓ Verified Predictions (ML-Precedent Agreement)\n\n")
            f.write(f"**{len(verified)} reactions** where ML and precedent recommendations agree\n\n")
            f.write("| Reaction | Core | Yield Comparison |\n")
            f.write("|----------|------|------------------|\n")
            for c in verified[:20]:  # Show first 20
                desc = c['reaction'][:40] + "..." if len(c['reaction']) > 40 else c['reaction']
                f.write(f"| {desc} | {c['ml_core']} | ML: {c['ml_yield']:.1f}% / Prec: {c['precedent_yield']:.1f}% |\n")
            if len(verified) > 20:
                f.write(f"\n*...and {len(verified) - 20} more*\n")
            f.write("\n")
        
        # Recommendations
        f.write("## Recommendations\n\n")
        if len(major) > 0:
            f.write(f"1. **Review {len(major)} major discrepancies** - These predictions may be unreliable\n")
        if len(unverifiable) > 0:
            f.write(f"2. **{len(unverifiable)} predictions lack precedents** - Consider running these experimentally to validate\n")
        if len(verified) + len(minor) > total * 0.7:
            f.write(f"3. **{(len(verified) + len(minor))/total*100:.0f}% agreement** - ML model shows good alignment with precedents\n")
        f.write("\n")


if __name__ == "__main__":
    # Paths
    ml_results_path = "results/cn_coupling_test/cn_coupling_results.json"
    model_path = "models/cn_coupling_pd_buchwald_v1.pkl"
    dataset_path = "data/reaction_dataset/C_N_Coupling_Pd.jsonl"
    output_path = "results/cn_coupling_test/verification_results.json"
    
    # Run verification
    verify_predictions(ml_results_path, model_path, dataset_path, output_path)

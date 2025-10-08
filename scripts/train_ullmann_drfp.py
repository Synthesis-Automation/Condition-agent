"""
Train DRFP-based yield prediction model for Ullmann C-N coupling reactions.

Similar to Buchwald model, but trained on Cu-catalyzed reactions.
"""

import json
import pickle
import sys
from pathlib import Path
from collections import Counter
from typing import Dict, List, Any

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, r2_score
import matplotlib.pyplot as plt

# Add parent to path
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from chemtools.ml.drfp_predictor import DRFPYieldPredictor


def load_ullmann_dataset(jsonl_path: str) -> pd.DataFrame:
    """Load and parse Ullmann dataset from JSONL."""
    records = []
    
    with open(jsonl_path, 'r', encoding='utf-8') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            
            try:
                rec = json.loads(line)
            except json.JSONDecodeError as e:
                print(f"Warning: Failed to parse line {line_num}: {e}")
                continue
            
            # Extract yield
            yield_pct = rec.get('conditions', {}).get('yield_pct')
            if yield_pct is None or not isinstance(yield_pct, (int, float)):
                continue
            
            # Extract SMILES
            smiles_data = rec.get('smiles', {})
            reactants = smiles_data.get('reactants', '')
            products = smiles_data.get('products', '')
            
            if not reactants or not products:
                continue
            
            reaction_smiles = f"{reactants}>>{products}"
            
            # Extract core (catalyst)
            core = rec.get('condition_core', '')
            if not core:
                continue
            
            # Extract base (from reagents with role='BASE')
            base_uid = ''
            for reagent in rec.get('reagents', []):
                if reagent.get('role') == 'BASE':
                    base_uid = reagent.get('cas', '')
                    break
            
            # Extract solvent (first solvent)
            solvents = rec.get('solvents', [])
            solvent_uid = ''
            if solvents and isinstance(solvents, list):
                first_solv = solvents[0]
                if isinstance(first_solv, dict):
                    solvent_uid = first_solv.get('cas', '')
            
            # Extract temperature and time (often null in Ullmann dataset)
            T_C = rec.get('conditions', {}).get('temperature_c')
            time_h = rec.get('conditions', {}).get('time_h')
            
            # Use defaults if missing
            if T_C is None or not isinstance(T_C, (int, float)):
                T_C = 110.0  # Typical Ullmann temperature
            if time_h is None or not isinstance(time_h, (int, float)):
                time_h = 12.0  # Typical Ullmann time
            
            records.append({
                'reaction_id': rec.get('reaction_id', ''),
                'reaction_smiles': reaction_smiles,
                'core': core,
                'base_uid': base_uid,
                'solvent_uid': solvent_uid,
                'T_C': float(T_C),
                'time_h': float(time_h),
                'yield_pct': float(yield_pct),
                'reference': rec.get('reference', ''),
            })
    
    df = pd.DataFrame(records)
    print(f"\n✓ Loaded {len(df)} Ullmann reactions from {jsonl_path}")
    return df


def filter_and_clean_data(df: pd.DataFrame, min_support: int = 10) -> pd.DataFrame:
    """Filter reactions with sufficient representation."""
    print(f"\n[1/4] Filtering dataset...")
    print(f"  Initial: {len(df)} reactions")
    
    # Remove reactions with empty core, base, or solvent
    df_clean = df[
        (df['core'] != '') & 
        (df['base_uid'] != '') & 
        (df['solvent_uid'] != '')
    ].copy()
    print(f"  After removing empty fields: {len(df_clean)} reactions")
    
    # Count occurrences
    core_counts = Counter(df_clean['core'])
    base_counts = Counter(df_clean['base_uid'])
    solvent_counts = Counter(df_clean['solvent_uid'])
    
    # Filter to cores with min_support
    valid_cores = {c for c, n in core_counts.items() if n >= min_support}
    df_filtered = df_clean[df_clean['core'].isin(valid_cores)].copy()
    print(f"  After filtering cores (min {min_support}): {len(df_filtered)} reactions")
    
    # Report vocabulary sizes
    print(f"\n  Vocabulary sizes:")
    print(f"    Cores: {df_filtered['core'].nunique()}")
    print(f"    Bases: {df_filtered['base_uid'].nunique()}")
    print(f"    Solvents: {df_filtered['solvent_uid'].nunique()}")
    
    return df_filtered


def main():
    """Train Ullmann DRFP yield prediction model."""
    
    print("=" * 80)
    print("ULLMANN C-N COUPLING - DRFP YIELD PREDICTION MODEL")
    print("=" * 80)
    
    # Paths
    dataset_path = ROOT / "data" / "reaction_dataset" / "C_N_Coupling_Cu.jsonl"
    output_dir = ROOT / "models"
    output_dir.mkdir(exist_ok=True)
    
    model_path = output_dir / "cn_coupling_cu_ullmann_v1.pkl"
    report_path = output_dir / "cn_coupling_cu_ullmann_training_report.txt"
    
    # Load data
    df = load_ullmann_dataset(str(dataset_path))
    
    # Filter and clean
    df_clean = filter_and_clean_data(df, min_support=10)
    
    if len(df_clean) < 100:
        print(f"\n❌ Error: Only {len(df_clean)} reactions after filtering. Need at least 100.")
        return 1
    
    # Split data
    print(f"\n[2/4] Splitting data...")
    train_df, test_df = train_test_split(
        df_clean, 
        test_size=0.15, 
        random_state=42,
        stratify=pd.cut(df_clean['yield_pct'], bins=[0, 60, 80, 100])
    )
    
    train_df, val_df = train_test_split(
        train_df,
        test_size=0.15,
        random_state=42,
        stratify=pd.cut(train_df['yield_pct'], bins=[0, 60, 80, 100])
    )
    
    print(f"  Train: {len(train_df)} reactions")
    print(f"  Val:   {len(val_df)} reactions")
    print(f"  Test:  {len(test_df)} reactions")
    
    # Train model
    print(f"\n[3/4] Training DRFP model...")
    print(f"  DRFP parameters: n_bits=2048, radius=3")
    print(f"  Model: LightGBM with 500 estimators")
    
    predictor = DRFPYieldPredictor(
        n_bits=2048, 
        radius=3,
        n_estimators=500,
        learning_rate=0.05,
        max_depth=6
    )
    
    predictor.fit(
        train_df,
        val_df=val_df,
        verbose=True
    )
    
    # Evaluate
    print(f"\n[4/4] Evaluating model...")
    
    train_preds = predictor.predict(train_df)
    val_preds = predictor.predict(val_df)
    test_preds = predictor.predict(test_df)
    
    train_mae = mean_absolute_error(train_df['yield_pct'], train_preds)
    val_mae = mean_absolute_error(val_df['yield_pct'], val_preds)
    test_mae = mean_absolute_error(test_df['yield_pct'], test_preds)
    
    train_r2 = r2_score(train_df['yield_pct'], train_preds)
    val_r2 = r2_score(val_df['yield_pct'], val_preds)
    test_r2 = r2_score(test_df['yield_pct'], test_preds)
    
    print(f"\n  Performance Metrics:")
    print(f"  {'Split':<10} {'MAE':<10} {'R²':<10} {'Samples':<10}")
    print(f"  {'-'*40}")
    print(f"  {'Train':<10} {train_mae:>6.2f}%    {train_r2:>6.3f}    {len(train_df):<10}")
    print(f"  {'Val':<10} {val_mae:>6.2f}%    {val_r2:>6.3f}    {len(val_df):<10}")
    print(f"  {'Test':<10} {test_mae:>6.2f}%    {test_r2:>6.3f}    {len(test_df):<10}")
    
    # Save model
    print(f"\n[5/5] Saving model...")
    
    model_data = {
        'model': predictor.model,
        'drfp_encoder': predictor.drfp_encoder,
        'core_vocab': predictor.core_vocab,
        'base_vocab': predictor.base_vocab,
        'solvent_vocab': predictor.solvent_vocab,
        'n_bits': predictor.n_bits,
        'radius': predictor.radius,
        'train_mae': train_mae,
        'val_mae': val_mae,
        'test_mae': test_mae,
        'train_r2': train_r2,
        'val_r2': val_r2,
        'test_r2': test_r2,
        'n_train': len(train_df),
        'n_val': len(val_df),
        'n_test': len(test_df),
    }
    
    with open(model_path, 'wb') as f:
        pickle.dump(model_data, f)
    
    print(f"  ✓ Model saved to: {model_path}")
    
    # Generate report
    report = f"""ULLMANN C-N COUPLING - DRFP YIELD PREDICTION MODEL
{'='*80}

DATASET SUMMARY
{'-'*80}
Total reactions:     {len(df):>6}
After filtering:     {len(df_clean):>6}
Training set:        {len(train_df):>6}
Validation set:      {len(val_df):>6}
Test set:            {len(test_df):>6}

VOCABULARY SIZES
{'-'*80}
Cores (catalysts):   {len(predictor.core_vocab):>6}
Bases:               {len(predictor.base_vocab):>6}
Solvents:            {len(predictor.solvent_vocab):>6}

TOP 10 CORES (CATALYSTS)
{'-'*80}
"""
    
    core_counts = Counter(df_clean['core'])
    for core, count in core_counts.most_common(10):
        report += f"{core:<50} {count:>6} reactions\n"
    
    report += f"""
TOP 10 BASES
{'-'*80}
"""
    base_counts = Counter(df_clean['base_uid'])
    for base, count in base_counts.most_common(10):
        report += f"{base:<50} {count:>6} reactions\n"
    
    report += f"""
TOP 10 SOLVENTS
{'-'*80}
"""
    solvent_counts = Counter(df_clean['solvent_uid'])
    for solvent, count in solvent_counts.most_common(10):
        report += f"{solvent:<50} {count:>6} reactions\n"
    
    report += f"""
MODEL PERFORMANCE
{'-'*80}
Split          MAE        R²         Samples
{'-'*80}
Train      {train_mae:>7.2f}%   {train_r2:>6.3f}    {len(train_df):>6}
Val        {val_mae:>7.2f}%   {val_r2:>6.3f}    {len(val_df):>6}
Test       {test_mae:>7.2f}%   {test_r2:>6.3f}    {len(test_df):>6}

YIELD DISTRIBUTION
{'-'*80}
High (≥80%):     {len(df_clean[df_clean['yield_pct'] >= 80]):>6} ({len(df_clean[df_clean['yield_pct'] >= 80])/len(df_clean)*100:>5.1f}%)
Medium (60-79%): {len(df_clean[(df_clean['yield_pct'] >= 60) & (df_clean['yield_pct'] < 80)]):>6} ({len(df_clean[(df_clean['yield_pct'] >= 60) & (df_clean['yield_pct'] < 80)])/len(df_clean)*100:>5.1f}%)
Low (<60%):      {len(df_clean[df_clean['yield_pct'] < 60]):>6} ({len(df_clean[df_clean['yield_pct'] < 60])/len(df_clean)*100:>5.1f}%)

Average yield:   {df_clean['yield_pct'].mean():>6.1f}%
Median yield:    {df_clean['yield_pct'].median():>6.1f}%
Min yield:       {df_clean['yield_pct'].min():>6.1f}%
Max yield:       {df_clean['yield_pct'].max():>6.1f}%

NOTES
{'-'*80}
- Dataset: C_N_Coupling_Cu.jsonl (5,552 Cu-catalyzed C-N coupling reactions)
- DRFP fingerprints: 2048-bit, radius=3
- Model: LightGBM (500 estimators, learning_rate=0.05)
- Missing T/time values filled with defaults (110°C, 12h)
- Filtered to cores with ≥10 examples for robust training

COMPARISON: Ullmann (Cu) vs Buchwald (Pd)
{'-'*80}
                    Ullmann (Cu)     Buchwald (Pd)
Dataset size:       ~4,700           ~1,343
Avg yield:          ~74%             ~73%
Catalysts:          Cu-based         Pd/ligand complexes
Typical temp:       110°C            100°C
Ligands:            Simple/none      Complex (XPhos, RuPhos, etc.)

Cu catalysis is simpler but comparable yields to Pd!
"""
    
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(report)
    
    print(f"  ✓ Report saved to: {report_path}")
    
    print(f"\n{'='*80}")
    print(f"✓ Training complete!")
    print(f"{'='*80}")
    print(f"\nModel: {model_path}")
    print(f"Report: {report_path}")
    print(f"\nNext steps:")
    print(f"  1. Test on sample reactions: python scripts/test_ullmann_reactions.py")
    print(f"  2. Verify predictions: python scripts/verify_ullmann_ml_with_rules.py")
    print(f"  3. Compare Cu vs Pd: python scripts/compare_ullmann_buchwald.py")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python
"""
Prepare dataset for ML training.

Extracts reaction data from data/reaction_dataset/*.jsonl and creates
train/val/test splits stratified by reaction family.

Usage:
    python scripts/ml/prepare_dataset.py
    python scripts/ml/prepare_dataset.py --source data/reaction_dataset --output data/ml_splits
"""

from __future__ import annotations

import json
import argparse
from pathlib import Path
from typing import List, Dict, Any
import sys

try:
    import pandas as pd
    from sklearn.model_selection import train_test_split
except ImportError:
    print("ERROR: Required dependencies not installed.")
    print("Please install: pip install pandas scikit-learn")
    sys.exit(1)


def load_reaction_dataset(dataset_dir: str) -> List[Dict[str, Any]]:
    """Load all JSONL files from dataset directory."""
    rows = []
    dataset_path = Path(dataset_dir)
    
    if not dataset_path.exists():
        print(f"WARNING: Dataset directory not found: {dataset_dir}")
        return rows
    
    jsonl_files = list(dataset_path.glob('*.jsonl'))
    
    if not jsonl_files:
        print(f"WARNING: No .jsonl files found in {dataset_dir}")
        return rows
    
    print(f"Found {len(jsonl_files)} JSONL files")
    
    for jsonl_file in jsonl_files:
        print(f"  Loading: {jsonl_file.name}")
        with open(jsonl_file, 'r', encoding='utf-8') as f:
            file_count = 0
            for line in f:
                if line.strip():
                    try:
                        row = json.loads(line)
                        rows.append(row)
                        file_count += 1
                    except json.JSONDecodeError as e:
                        print(f"    WARNING: JSON decode error: {e}")
                        continue
            print(f"    Loaded {file_count} reactions")
    
    return rows


def extract_features(row: Dict[str, Any]) -> Dict[str, Any]:
    """Extract features from dataset row for ML training."""
    # Get reaction SMILES
    smiles_block = row.get('smiles', {})
    reactants = smiles_block.get('reactants', '')
    products = smiles_block.get('products', '')
    
    # Normalize empty strings
    if not reactants or not products:
        return {}
    
    reaction_smiles = f"{reactants}>>{products}"
    
    # Get conditions
    conditions = row.get('conditions', {})
    core = row.get('condition_core', 'Unknown')
    
    # Get base (first base reagent)
    base_uid = None
    for reagent in row.get('reagents', []):
        if reagent.get('role', '').upper() == 'BASE':
            base_uid = reagent.get('cas') or reagent.get('uid') or reagent.get('name')
            break
    
    # Get solvent (first solvent)
    solvents = row.get('solvents', [])
    solvent_uid = None
    if solvents and isinstance(solvents, list):
        solvent_uid = solvents[0].get('cas') or solvents[0].get('uid') or solvents[0].get('name')
    
    # Get yield
    yield_pct = conditions.get('yield_pct')
    
    # Get T and time
    temp_c = conditions.get('temperature_c')
    time_h = conditions.get('time_h')
    
    # Only return if we have minimum required data
    if not all([reaction_smiles, yield_pct is not None]):
        return {}
    
    return {
        'reaction_id': row.get('reaction_id'),
        'reaction_smiles': reaction_smiles,
        'reaction_type': row.get('reaction_type', 'Unknown'),
        'core': core,
        'base_uid': base_uid or 'Unknown',
        'solvent_uid': solvent_uid or 'Unknown',
        'T_C': temp_c if temp_c is not None else 80.0,  # Default
        'time_h': time_h if time_h is not None else 12.0,  # Default
        'yield_pct': yield_pct,
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Prepare ML training dataset from reaction JSONL files"
    )
    parser.add_argument(
        '--source',
        default='data/reaction_dataset',
        help='Source directory containing .jsonl files (default: data/reaction_dataset)'
    )
    parser.add_argument(
        '--output-dir',
        default='data',
        help='Output directory for train/val/test splits (default: data/)'
    )
    parser.add_argument(
        '--train-size',
        type=float,
        default=0.7,
        help='Training set fraction (default: 0.7)'
    )
    parser.add_argument(
        '--val-size',
        type=float,
        default=0.15,
        help='Validation set fraction (default: 0.15)'
    )
    parser.add_argument(
        '--random-state',
        type=int,
        default=42,
        help='Random seed for reproducibility (default: 42)'
    )
    parser.add_argument(
        '--min-yield',
        type=float,
        default=0.0,
        help='Minimum yield to include (default: 0.0)'
    )
    parser.add_argument(
        '--max-yield',
        type=float,
        default=100.0,
        help='Maximum yield to include (default: 100.0)'
    )
    
    args = parser.parse_args(argv)
    
    # Validate split fractions
    test_size = 1.0 - args.train_size - args.val_size
    if test_size <= 0 or test_size >= 1:
        print(f"ERROR: Invalid split fractions (train={args.train_size}, val={args.val_size})")
        return 1
    
    print("=" * 70)
    print("ML Dataset Preparation")
    print("=" * 70)
    print(f"Source: {args.source}")
    print(f"Output: {args.output_dir}")
    print(f"Split: {args.train_size:.0%} train / {args.val_size:.0%} val / {test_size:.0%} test")
    print()
    
    # Load data
    print("Loading dataset...")
    rows = load_reaction_dataset(args.source)
    
    if not rows:
        print("ERROR: No reactions loaded. Check source directory.")
        return 1
    
    print(f"Loaded {len(rows)} total reactions")
    print()
    
    # Extract features
    print("Extracting features...")
    features_list = []
    for row in rows:
        feat = extract_features(row)
        if feat:  # Only add if valid
            features_list.append(feat)
    
    if not features_list:
        print("ERROR: No valid reactions after feature extraction")
        return 1
    
    df = pd.DataFrame(features_list)
    print(f"Extracted features from {len(df)} reactions")
    print()
    
    # Filter by yield range
    df = df[
        (df['yield_pct'] >= args.min_yield) &
        (df['yield_pct'] <= args.max_yield)
    ]
    print(f"After yield filtering ({args.min_yield}-{args.max_yield}%): {len(df)} reactions")
    print()
    
    if len(df) < 10:
        print("ERROR: Too few reactions after filtering (need at least 10)")
        return 1
    
    # Statistics
    print("Dataset Statistics:")
    print("-" * 70)
    print(f"  Reactions:      {len(df)}")
    print(f"  Reaction types: {df['reaction_type'].nunique()}")
    print(f"  Unique cores:   {df['core'].nunique()}")
    print(f"  Unique bases:   {df['base_uid'].nunique()}")
    print(f"  Unique solvents: {df['solvent_uid'].nunique()}")
    print()
    print(f"  Yield (mean):   {df['yield_pct'].mean():.1f}%")
    print(f"  Yield (median): {df['yield_pct'].median():.1f}%")
    print(f"  Yield (std):    {df['yield_pct'].std():.1f}%")
    print(f"  Yield (min):    {df['yield_pct'].min():.1f}%")
    print(f"  Yield (max):    {df['yield_pct'].max():.1f}%")
    print()
    
    # Reaction type distribution
    print("Reaction type distribution:")
    for rxn_type, count in df['reaction_type'].value_counts().items():
        print(f"  {rxn_type:30s} {count:5d} ({count/len(df)*100:.1f}%)")
    print()
    
    # Train/val/test split (stratified by reaction_type)
    # Check if we have enough samples per class for stratification
    min_class_size = df['reaction_type'].value_counts().min()
    
    if min_class_size < 3:
        print("WARNING: Some reaction types have <3 samples, stratification disabled")
        stratify_col = None
    else:
        stratify_col = df['reaction_type']
    
    # First split: train vs (val + test)
    train_df, temp_df = train_test_split(
        df,
        test_size=(args.val_size + test_size),
        random_state=args.random_state,
        stratify=stratify_col
    )
    
    # Second split: val vs test
    if stratify_col is not None:
        temp_stratify = temp_df['reaction_type']
    else:
        temp_stratify = None
    
    val_df, test_df = train_test_split(
        temp_df,
        test_size=(test_size / (args.val_size + test_size)),
        random_state=args.random_state,
        stratify=temp_stratify
    )
    
    print("Split Sizes:")
    print(f"  Train: {len(train_df):5d} ({len(train_df)/len(df)*100:.1f}%)")
    print(f"  Val:   {len(val_df):5d} ({len(val_df)/len(df)*100:.1f}%)")
    print(f"  Test:  {len(test_df):5d} ({len(test_df)/len(df)*100:.1f}%)")
    print()
    
    # Save splits
    output_path = Path(args.output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    train_file = output_path / 'ml_train.jsonl'
    val_file = output_path / 'ml_val.jsonl'
    test_file = output_path / 'ml_test.jsonl'
    
    print("Saving splits...")
    train_df.to_json(train_file, orient='records', lines=True, force_ascii=False)
    val_df.to_json(val_file, orient='records', lines=True, force_ascii=False)
    test_df.to_json(test_file, orient='records', lines=True, force_ascii=False)
    
    print(f"  Train: {train_file}")
    print(f"  Val:   {val_file}")
    print(f"  Test:  {test_file}")
    print()
    
    print("=" * 70)
    print("✓ Dataset preparation complete!")
    print("=" * 70)
    print()
    print("Next steps:")
    print("  1. Install DRFP: pip install drfp")
    print("  2. Train model: python scripts/ml/train_drfp_model.py")
    print()
    
    return 0


if __name__ == '__main__':
    sys.exit(main())

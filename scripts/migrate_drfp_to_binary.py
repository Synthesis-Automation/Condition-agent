#!/usr/bin/env python3
"""
Migrate DRFP fingerprints from JSONL to binary NPZ format.

This script extracts the 4096-element drfp_fp arrays from JSONL files
and saves them to compressed NPZ files, dramatically reducing file size.

Usage:
    # Migrate single family
    python scripts/migrate_drfp_to_binary.py --family C_N_Coupling_Cu
    
    # Migrate all families
    python scripts/migrate_drfp_to_binary.py --all
    
    # Migrate and remove drfp_fp from JSONL (saves space)
    python scripts/migrate_drfp_to_binary.py --all --remove-from-jsonl
    
    # Just analyze without migrating
    python scripts/migrate_drfp_to_binary.py --all --dry-run

File Size Comparison:
    Before: Suzuki.jsonl = 670 MB (with embedded drfp_fp arrays)
    After:  Suzuki.jsonl = ~70 MB + Suzuki_drfp.npz = ~12 MB
    Savings: ~90% total storage
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Optional, List, Tuple

# Ensure project root on path
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from chemtools.util.drfp_storage import save_drfp_index


def get_file_size_mb(path: str) -> float:
    """Get file size in MB."""
    if not os.path.exists(path):
        return 0.0
    return os.path.getsize(path) / (1024 * 1024)


def migrate_family(
    family: str,
    data_dir: str = "data/reaction_dataset",
    remove_from_jsonl: bool = False,
    dry_run: bool = False
) -> Tuple[int, float, float]:
    """
    Migrate DRFP fingerprints for a single reaction family.
    
    Args:
        family: Reaction family name (e.g., "C_N_Coupling_Cu")
        data_dir: Directory containing JSONL files
        remove_from_jsonl: If True, remove drfp_fp from JSONL after migration
        dry_run: If True, only analyze without writing files
    
    Returns:
        Tuple of (num_fingerprints, original_size_mb, new_size_mb)
    """
    jsonl_path = os.path.join(data_dir, f"{family}.jsonl")
    npz_path = os.path.join(data_dir, f"{family}_drfp.npz")
    
    if not os.path.exists(jsonl_path):
        print(f"❌ File not found: {jsonl_path}")
        return 0, 0.0, 0.0
    
    original_size = get_file_size_mb(jsonl_path)
    
    print(f"\n{'='*70}")
    print(f"Processing: {family}")
    print(f"{'='*70}")
    print(f"Input:  {jsonl_path} ({original_size:.2f} MB)")
    print(f"Output: {npz_path}")
    
    # Extract fingerprints
    fingerprints = []
    reaction_ids = []
    n_bits = 4096
    radius = 3
    
    print("\nExtracting DRFP fingerprints from JSONL...")
    with open(jsonl_path, 'r', encoding='utf-8') as f:
        for line_num, line in enumerate(f, 1):
            if line_num % 1000 == 0:
                print(f"  Processed {line_num} lines, found {len(fingerprints)} fingerprints...", end='\r')
            
            line = line.strip()
            if not line:
                continue
            
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                continue
            
            # Get reaction_id
            reaction_id = record.get('reaction_id')
            if not reaction_id:
                continue
            
            # Extract DRFP from precomputed field
            precomputed = record.get('precomputed', {})
            if not isinstance(precomputed, dict):
                continue
            
            drfp_fp = precomputed.get('drfp_fp')
            if drfp_fp is None:
                continue
            
            # Get metadata
            if 'drfp_n_bits' in precomputed:
                n_bits = precomputed['drfp_n_bits']
            if 'drfp_radius' in precomputed:
                radius = precomputed['drfp_radius']
            
            fingerprints.append(drfp_fp)
            reaction_ids.append(reaction_id)
    
    print(f"\n✓ Extracted {len(fingerprints)} DRFP fingerprints")
    
    if not fingerprints:
        print("⚠️  No fingerprints found - skipping this family")
        return 0, original_size, original_size
    
    # Save to NPZ
    if not dry_run:
        print(f"\nSaving to NPZ format...")
        save_drfp_index(fingerprints, reaction_ids, npz_path, n_bits, radius)
        npz_size = get_file_size_mb(npz_path)
        print(f"✓ Created {npz_path} ({npz_size:.2f} MB)")
    else:
        print(f"\n[DRY RUN] Would save to {npz_path}")
        # Estimate NPZ size (typically ~10-15% of uncompressed)
        est_uncompressed = len(fingerprints) * 4096  # bytes
        npz_size = est_uncompressed / (1024 * 1024) * 0.12  # Estimate 12% compression
    
    # Remove drfp_fp from JSONL if requested
    new_jsonl_size = original_size
    if remove_from_jsonl and not dry_run:
        print(f"\nRemoving drfp_fp from JSONL to save space...")
        temp_path = jsonl_path + ".tmp"
        
        with open(jsonl_path, 'r', encoding='utf-8') as fin:
            with open(temp_path, 'w', encoding='utf-8') as fout:
                for line_num, line in enumerate(fin, 1):
                    if line_num % 1000 == 0:
                        print(f"  Cleaned {line_num} lines...", end='\r')
                    
                    line = line.strip()
                    if not line:
                        fout.write('\n')
                        continue
                    
                    try:
                        record = json.loads(line)
                        
                        # Remove drfp_fp from precomputed
                        if 'precomputed' in record and isinstance(record['precomputed'], dict):
                            record['precomputed'].pop('drfp_fp', None)
                        
                        fout.write(json.dumps(record, ensure_ascii=False))
                        fout.write('\n')
                    except Exception as e:
                        # Keep original line if error
                        fout.write(line)
                        fout.write('\n')
        
        # Replace original with cleaned version
        os.replace(temp_path, jsonl_path)
        new_jsonl_size = get_file_size_mb(jsonl_path)
        print(f"\n✓ Cleaned JSONL: {jsonl_path} ({new_jsonl_size:.2f} MB)")
    elif remove_from_jsonl and dry_run:
        print(f"\n[DRY RUN] Would remove drfp_fp from JSONL")
        # Estimate ~90% size reduction
        new_jsonl_size = original_size * 0.10
    
    # Summary
    print(f"\n{'-'*70}")
    print(f"SUMMARY: {family}")
    print(f"{'-'*70}")
    print(f"Fingerprints migrated: {len(fingerprints)}")
    print(f"Original JSONL size:   {original_size:.2f} MB")
    if remove_from_jsonl:
        print(f"New JSONL size:        {new_jsonl_size:.2f} MB")
        print(f"NPZ file size:         {npz_size:.2f} MB")
        print(f"Total new size:        {new_jsonl_size + npz_size:.2f} MB")
        savings_pct = (1 - (new_jsonl_size + npz_size) / original_size) * 100
        print(f"Space saved:           {original_size - new_jsonl_size - npz_size:.2f} MB ({savings_pct:.1f}%)")
    else:
        print(f"NPZ file size:         {npz_size:.2f} MB")
        print(f"Space saved:           0 MB (drfp_fp still in JSONL)")
    
    return len(fingerprints), original_size, new_jsonl_size + npz_size


def main(argv: Optional[List[str]] = None) -> int:
    """Main entry point."""
    ap = argparse.ArgumentParser(
        description="Migrate DRFP fingerprints from JSONL to binary NPZ format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Migrate single family
  python scripts/migrate_drfp_to_binary.py --family C_N_Coupling_Cu
  
  # Migrate all families
  python scripts/migrate_drfp_to_binary.py --all
  
  # Migrate and clean JSONL (saves ~90% space)
  python scripts/migrate_drfp_to_binary.py --all --remove-from-jsonl
  
  # Dry run to see what would happen
  python scripts/migrate_drfp_to_binary.py --all --remove-from-jsonl --dry-run
        """
    )
    ap.add_argument(
        '--family',
        help='Reaction family name (e.g., C_N_Coupling_Cu, Suzuki, Amide_formation)'
    )
    ap.add_argument(
        '--all',
        action='store_true',
        help='Migrate all families in data/reaction_dataset/'
    )
    ap.add_argument(
        '--data-dir',
        default='data/reaction_dataset',
        help='Directory containing JSONL files (default: data/reaction_dataset)'
    )
    ap.add_argument(
        '--remove-from-jsonl',
        action='store_true',
        help='Remove drfp_fp from JSONL after migration (saves ~90%% space)'
    )
    ap.add_argument(
        '--dry-run',
        action='store_true',
        help='Analyze without writing files'
    )
    
    args = ap.parse_args(argv)
    
    if not args.family and not args.all:
        ap.error("Must specify --family <name> or --all")
    
    # Determine families to process
    families = []
    if args.all:
        # Find all JSONL files
        if not os.path.exists(args.data_dir):
            print(f"❌ Data directory not found: {args.data_dir}")
            return 1
        
        for filename in os.listdir(args.data_dir):
            if filename.endswith('.jsonl'):
                family = filename[:-6]  # Remove .jsonl
                families.append(family)
        
        families.sort()
        print(f"Found {len(families)} families: {', '.join(families)}")
    else:
        families = [args.family]
    
    if not families:
        print("❌ No families found to process")
        return 1
    
    # Process each family
    total_fps = 0
    total_original = 0.0
    total_new = 0.0
    
    for family in families:
        num_fps, orig_size, new_size = migrate_family(
            family,
            data_dir=args.data_dir,
            remove_from_jsonl=args.remove_from_jsonl,
            dry_run=args.dry_run
        )
        total_fps += num_fps
        total_original += orig_size
        total_new += new_size
    
    # Overall summary
    print(f"\n{'='*70}")
    print(f"OVERALL SUMMARY")
    print(f"{'='*70}")
    print(f"Families processed:    {len(families)}")
    print(f"Total fingerprints:    {total_fps}")
    print(f"Total original size:   {total_original:.2f} MB")
    print(f"Total new size:        {total_new:.2f} MB")
    if args.remove_from_jsonl:
        savings_pct = (1 - total_new / total_original) * 100 if total_original > 0 else 0
        print(f"Total space saved:     {total_original - total_new:.2f} MB ({savings_pct:.1f}%)")
    
    if args.dry_run:
        print(f"\n⚠️  DRY RUN - No files were modified")
        print(f"Remove --dry-run to perform actual migration")
    
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

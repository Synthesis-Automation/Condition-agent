#!/usr/bin/env python3
"""
Build Unified DRFP Index for Cross-Family Search

This script creates a single NPZ file containing DRFP fingerprints from all
reaction families. This enables fast cross-family search with DRFP similarity.

Usage:
    python scripts/build_unified_drfp_index.py
    
    # Custom output path
    python scripts/build_unified_drfp_index.py --output data/reaction_dataset/unified_drfp.npz
    
    # Specify families to include
    python scripts/build_unified_drfp_index.py --families C_N_Coupling Suzuki Amide_formation
    
    # Rebuild (overwrite existing)
    python scripts/build_unified_drfp_index.py --force

Features:
    - Loads all family NPZ files
    - Combines into single unified index
    - Preserves reaction_id → fingerprint mapping
    - Validates no duplicate reaction_ids
    - Shows progress and statistics
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

# Add project root to path
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

import numpy as np


def get_available_family_npz_files(data_dir: str = "data/reaction_dataset") -> List[str]:
    """Find all family-specific DRFP NPZ files."""
    data_path = Path(data_dir)
    if not data_path.exists():
        return []
    
    # Find all *_drfp.npz files except unified
    npz_files = []
    for npz_file in data_path.glob("*_drfp.npz"):
        # Skip unified file if it exists
        if "unified" not in npz_file.name.lower() and "all_families" not in npz_file.name.lower():
            npz_files.append(str(npz_file))
    
    return sorted(npz_files)


def load_family_npz(npz_path: str) -> Dict[str, np.ndarray]:
    """Load fingerprints from a family NPZ file."""
    data = np.load(npz_path, allow_pickle=True)
    
    # Extract reaction_ids and fingerprints
    # Try both 'fingerprints' and 'fps' keys (different versions use different names)
    reaction_ids = data.get('reaction_ids', [])
    fingerprints = data.get('fingerprints', data.get('fps', []))
    
    if len(reaction_ids) == 0:
        print(f"  ⚠️  Warning: No reaction_ids found in {npz_path}")
        return {}
    
    if len(fingerprints) == 0:
        print(f"  ⚠️  Warning: No fingerprints found in {npz_path}")
        return {}
    
    # Convert to dict
    fingerprint_dict = {}
    for rid, fp in zip(reaction_ids, fingerprints):
        fingerprint_dict[str(rid)] = np.array(fp, dtype='uint8')
    
    return fingerprint_dict


def build_unified_index(
    family_npz_files: List[str],
    output_path: str,
    verbose: bool = True
) -> Dict[str, any]:
    """
    Build unified DRFP index from multiple family NPZ files.
    
    Args:
        family_npz_files: List of paths to family NPZ files
        output_path: Path to save unified NPZ file
        verbose: Print progress messages
        
    Returns:
        Dict with statistics (total_reactions, total_families, file_size, etc.)
    """
    if verbose:
        print("=" * 80)
        print("Building Unified DRFP Index")
        print("=" * 80)
        print()
    
    # Collect all fingerprints
    all_fingerprints: Dict[str, np.ndarray] = {}
    family_stats = {}
    duplicates = []
    
    for npz_path in family_npz_files:
        family_name = Path(npz_path).stem.replace("_drfp", "")
        
        if verbose:
            print(f"Loading {family_name}...")
        
        # Load fingerprints from this family
        fp_dict = load_family_npz(npz_path)
        
        # Check for duplicates
        for rid in fp_dict.keys():
            if rid in all_fingerprints:
                duplicates.append((rid, family_name))
        
        # Add to unified dict
        all_fingerprints.update(fp_dict)
        
        family_stats[family_name] = len(fp_dict)
        
        if verbose:
            print(f"  ✓ Loaded {len(fp_dict)} reactions from {family_name}")
    
    if verbose:
        print()
        print(f"Total reactions: {len(all_fingerprints)}")
        
        if duplicates:
            print(f"⚠️  Found {len(duplicates)} duplicate reaction_ids:")
            for rid, family in duplicates[:10]:  # Show first 10
                print(f"  - {rid} in {family}")
            if len(duplicates) > 10:
                print(f"  ... and {len(duplicates) - 10} more")
        print()
    
    # Prepare arrays for NPZ storage
    reaction_ids = list(all_fingerprints.keys())
    fingerprints = [all_fingerprints[rid] for rid in reaction_ids]
    
    # Get fingerprint parameters (from first fingerprint)
    if fingerprints:
        first_fp = fingerprints[0]
        n_bits = len(first_fp) * 8  # Convert bytes to bits
        radius = 3  # Default DRFP radius
    else:
        n_bits = 4096
        radius = 3
    
    if verbose:
        print(f"Saving unified index to {output_path}...")
    
    # Save unified NPZ file (use 'fps' key to match existing format)
    np.savez_compressed(
        output_path,
        reaction_ids=np.array(reaction_ids, dtype=object),
        fps=np.array(fingerprints, dtype=object),  # Use 'fps' key for compatibility
        n_bits=n_bits,
        radius=radius,
        families=np.array(list(family_stats.keys()), dtype=object),
        family_counts=np.array(list(family_stats.values()), dtype=int),
    )
    
    # Get file size
    file_size = Path(output_path).stat().st_size / (1024 * 1024)  # MB
    
    if verbose:
        print(f"  ✓ Saved {len(reaction_ids)} fingerprints")
        print(f"  ✓ File size: {file_size:.2f} MB")
        print()
        
        print("Family Statistics:")
        print("-" * 80)
        for family, count in family_stats.items():
            print(f"  {family:30s} {count:5d} reactions")
        print("-" * 80)
        print(f"  {'TOTAL':30s} {len(all_fingerprints):5d} reactions")
        print()
        print("=" * 80)
        print("✓ Unified DRFP index built successfully!")
        print("=" * 80)
    
    return {
        "total_reactions": len(all_fingerprints),
        "total_families": len(family_stats),
        "file_size_mb": file_size,
        "output_path": output_path,
        "family_stats": family_stats,
        "duplicates": len(duplicates),
    }


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Build unified DRFP index for cross-family search",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Build unified index from all families
  python scripts/build_unified_drfp_index.py
  
  # Specify custom output path
  python scripts/build_unified_drfp_index.py --output data/my_unified_drfp.npz
  
  # Select specific families
  python scripts/build_unified_drfp_index.py --families C_N_Coupling Suzuki
  
  # Force rebuild (overwrite existing)
  python scripts/build_unified_drfp_index.py --force
  
  # Quiet mode (less output)
  python scripts/build_unified_drfp_index.py --quiet
        """
    )
    
    parser.add_argument(
        "--data-dir",
        type=str,
        default="data/reaction_dataset",
        help="Directory containing reaction dataset NPZ files (default: data/reaction_dataset)"
    )
    
    parser.add_argument(
        "--output",
        type=str,
        default="data/reaction_dataset/ALL_FAMILIES_drfp.npz",
        help="Output path for unified NPZ file (default: data/reaction_dataset/ALL_FAMILIES_drfp.npz)"
    )
    
    parser.add_argument(
        "--families",
        type=str,
        nargs="+",
        default=None,
        help="Specific families to include (default: all available). "
             "Example: --families C_N_Coupling Suzuki Amide_formation"
    )
    
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force rebuild even if output file exists"
    )
    
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Reduce output verbosity"
    )
    
    args = parser.parse_args()
    
    # Check if output already exists
    if Path(args.output).exists() and not args.force:
        print(f"Error: Output file already exists: {args.output}")
        print("Use --force to overwrite, or specify a different output path with --output")
        return 1
    
    # Get available NPZ files
    available_npz = get_available_family_npz_files(args.data_dir)
    
    if not available_npz:
        print(f"Error: No family DRFP NPZ files found in {args.data_dir}")
        print("Expected files like: C_N_Coupling_drfp.npz, Suzuki_drfp.npz, etc.")
        print()
        print("To generate family NPZ files, run:")
        print("  python scripts/build_family_drfp_index.py")
        return 1
    
    # Filter by specified families if provided
    if args.families:
        filtered_npz = []
        for family in args.families:
            # Look for NPZ file matching this family
            matching = [npz for npz in available_npz if family in npz]
            if matching:
                filtered_npz.extend(matching)
            else:
                print(f"Warning: No NPZ file found for family '{family}'")
        
        if not filtered_npz:
            print("Error: No matching NPZ files found for specified families")
            return 1
        
        available_npz = filtered_npz
    
    if not args.quiet:
        print(f"Found {len(available_npz)} family NPZ files:")
        for npz in available_npz:
            print(f"  - {Path(npz).name}")
        print()
    
    # Build unified index
    try:
        stats = build_unified_index(
            family_npz_files=available_npz,
            output_path=args.output,
            verbose=not args.quiet
        )
        
        if not args.quiet:
            print()
            print("Next steps:")
            print(f"  1. The unified index is ready at: {args.output}")
            print(f"  2. Cross-family search will automatically use this file")
            print(f"  3. Test with: python app/local_recommendation_cli.py --search-all-families --rxn \"...\"")
        
        return 0
        
    except Exception as e:
        print(f"Error building unified index: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())

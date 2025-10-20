#!/usr/bin/env python3
"""
Demonstrate DRFP storage optimization for protocol indexer

This script shows the size difference between:
1. OLD: Embedding DRFP fingerprints (4096 floats) in JSON
2. NEW: Storing DRFP fingerprints in separate NPZ file

Run this to see the space savings!
"""

import json
import tempfile
from pathlib import Path
import sys
import os

# Add project root to path
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from chemtools.protocol.indexer import ProtocolIndexer


def create_test_protocols(protocol_dir: Path, count: int = 20):
    """Create dummy protocol files for testing"""
    print(f"Creating {count} test protocol files...")
    
    for i in range(count):
        protocol = {
            "reaction": {
                "reaction_smiles": f"Br{'C'*i}.N>>{'C'*(i+1)}N",
                "family": "C_N_Coupling",
                "tags": "test;benchmark",
                "notes": f"Test protocol {i+1}",
                "reaction_SMARTS": ["[Br:1].[N:2]>>[C:1][N:2]"]
            },
            "source": {
                "title": f"Test Protocol {i+1}",
                "journal": "J. Test Chemistry",
                "year": 2025,
                "doi": f"10.1234/test.{i+1}"
            },
            "reaction_setup": [
                {
                    "chemicals": [],
                    "conditions": []
                }
            ]
        }
        
        filename = protocol_dir / f"test_protocol_{i+1:03d}.json"
        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(protocol, f, indent=2)
    
    print(f"✓ Created {count} protocol files\n")


def demonstrate_old_approach(protocol_dir: Path):
    """Simulate OLD approach: embedding DRFP in JSON"""
    print("=" * 70)
    print("OLD APPROACH: Embedding DRFP fingerprints in JSON")
    print("=" * 70)
    
    # Build index with DRFP
    indexer = ProtocolIndexer(protocol_dir=protocol_dir)
    indexer.build_index(compute_drfp=True, force_rebuild=True)
    
    # Simulate embedding DRFP in JSON (OLD way)
    old_index_data = {
        'metadata': indexer.metadata,
        'records': {},
        'family_index': indexer.family_index,
        'tag_index': indexer.tag_index
    }
    
    # Add records WITH embedded fingerprints (simulated)
    for filename, record in indexer.records.items():
        record_dict = record.to_dict()
        
        # Simulate embedding 4096 float values (OLD approach)
        if record.has_drfp:
            # Create fake fingerprint data (4096 floats)
            record_dict['drfp_fingerprint'] = [0.0] * 4096
        
        old_index_data['records'][filename] = record_dict
    
    # Save simulated old index
    old_path = protocol_dir / ".protocol_index_OLD.json"
    with open(old_path, 'w', encoding='utf-8') as f:
        json.dump(old_index_data, f, indent=2)
    
    old_size = old_path.stat().st_size
    print(f"Index file size: {old_size:,} bytes ({old_size / 1024:.1f} KB)")
    print(f"Storage format: Single JSON file with embedded DRFP arrays")
    print(f"DRFP per protocol: 4096 floats (~32 KB in JSON)")
    print()
    
    return old_size


def demonstrate_new_approach(protocol_dir: Path):
    """Demonstrate NEW approach: separate NPZ storage"""
    print("=" * 70)
    print("NEW APPROACH: Separate NPZ storage for DRFP")
    print("=" * 70)
    
    # Build index with DRFP (uses new approach automatically)
    indexer = ProtocolIndexer(protocol_dir=protocol_dir)
    indexer.build_index(compute_drfp=True, force_rebuild=True)
    indexer.save()
    
    index_path = protocol_dir / ".protocol_index.json"
    drfp_path = protocol_dir / ".protocol_drfp.npz"
    
    index_size = index_path.stat().st_size if index_path.exists() else 0
    drfp_size = drfp_path.stat().st_size if drfp_path.exists() else 0
    total_size = index_size + drfp_size
    
    print(f"Index JSON size: {index_size:,} bytes ({index_size / 1024:.1f} KB)")
    print(f"DRFP NPZ size:   {drfp_size:,} bytes ({drfp_size / 1024:.1f} KB)")
    print(f"Total size:      {total_size:,} bytes ({total_size / 1024:.1f} KB)")
    print(f"Storage format: JSON (metadata) + NPZ (binary DRFP)")
    print(f"DRFP per protocol: ~{drfp_size / len(indexer.records):.0f} bytes in compressed NPZ")
    print()
    
    return total_size, index_size, drfp_size


def main():
    """Run the demonstration"""
    print("\n" + "=" * 70)
    print("PROTOCOL DRFP STORAGE OPTIMIZATION DEMO")
    print("=" * 70)
    print()
    
    # Create temporary directory
    with tempfile.TemporaryDirectory() as tmpdir:
        protocol_dir = Path(tmpdir) / "protocol_db"
        protocol_dir.mkdir()
        
        # Create test protocols
        NUM_PROTOCOLS = 50
        create_test_protocols(protocol_dir, count=NUM_PROTOCOLS)
        
        # Demonstrate old approach
        old_size = demonstrate_old_approach(protocol_dir)
        
        # Demonstrate new approach
        new_total, new_index, new_drfp = demonstrate_new_approach(protocol_dir)
        
        # Show comparison
        print("=" * 70)
        print("COMPARISON SUMMARY")
        print("=" * 70)
        print(f"Number of protocols: {NUM_PROTOCOLS}")
        print()
        print(f"OLD approach (embedded JSON): {old_size:,} bytes ({old_size / 1024:.1f} KB)")
        print(f"NEW approach (JSON + NPZ):    {new_total:,} bytes ({new_total / 1024:.1f} KB)")
        print()
        
        savings = old_size - new_total
        savings_pct = (savings / old_size) * 100 if old_size > 0 else 0
        
        print(f"Space saved: {savings:,} bytes ({savings / 1024:.1f} KB)")
        print(f"Reduction:   {savings_pct:.1f}%")
        print()
        
        if savings_pct > 80:
            print("✓ Excellent space savings! (>80% reduction)")
        elif savings_pct > 50:
            print("✓ Good space savings! (>50% reduction)")
        else:
            print("✓ Some space savings achieved")
        
        print()
        print("Additional benefits of NEW approach:")
        print("  • Faster JSON parsing (no large arrays)")
        print("  • Lazy loading (DRFP loaded only when needed)")
        print("  • Binary format (faster operations)")
        print("  • Compressed storage (NPZ compression)")
        print()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nDemo interrupted.")
        sys.exit(1)
    except Exception as e:
        print(f"\n\nError: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

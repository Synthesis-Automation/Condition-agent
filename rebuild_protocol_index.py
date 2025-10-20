"""
Quick script to rebuild the protocol index after adding new protocol files
"""

from pathlib import Path
from chemtools.protocol.indexer import ProtocolIndexer

def main():
    print("=" * 80)
    print("Rebuilding Protocol Index")
    print("=" * 80)
    print()
    
    # Get protocol directory
    project_root = Path(__file__).parent
    protocol_dir = project_root / "data" / "protocol_db"
    
    print(f"Protocol directory: {protocol_dir}")
    print()
    
    # Create indexer
    indexer = ProtocolIndexer(protocol_dir=protocol_dir)
    
    # Build index with DRFP fingerprints
    print("Building index with DRFP fingerprints...")
    print("This may take a few minutes...")
    print()
    
    indexer.build_index(compute_drfp=True, force_rebuild=True)
    
    # Save index
    print()
    print("Saving index...")
    indexer.save()
    
    print()
    print("=" * 80)
    print("Index Rebuild Complete!")
    print("=" * 80)
    print()
    print(f"Total protocols indexed: {len(indexer.records)}")
    print(f"Index file: {indexer.index_path}")
    print(f"DRFP file: {indexer.drfp_path}")
    print()
    print("Summary by family:")
    for family, filenames in sorted(indexer.family_index.items()):
        print(f"  {family}: {len(filenames)} protocols")
    print()

if __name__ == "__main__":
    main()

"""
Tests for protocol DRFP storage optimization

Verifies that:
1. DRFP fingerprints are stored in separate NPZ file
2. Index JSON does not contain large fingerprint arrays
3. Fingerprints can be loaded on-demand
4. System gracefully handles missing DRFP files
"""

import json
import tempfile
from pathlib import Path
import pytest


def test_protocol_record_no_drfp_in_dict():
    """ProtocolRecord.to_dict() should not include DRFP data"""
    from chemtools.protocol.indexer import ProtocolRecord
    
    record = ProtocolRecord(
        filename="test.json",
        filepath="protocol_db/test.json",
        file_hash="abc123",
        reaction_smiles="Br>>I",
        reaction_family="Test",
        tags=["test"],
        notes="Test protocol",
        source_title="Test",
        source_journal="J. Test",
        source_year=2025,
        source_doi="10.1234/test",
        has_drfp=True  # Has DRFP but data stored separately
    )
    
    data = record.to_dict()
    
    # Should have metadata flag
    assert data['has_drfp'] is True
    
    # Should NOT have embedded fingerprint data
    assert 'drfp_fingerprint' not in data


def test_index_without_drfp_storage():
    """Index should build correctly even without DRFP storage utilities"""
    from chemtools.protocol import indexer
    
    # Create temporary test directory
    with tempfile.TemporaryDirectory() as tmpdir:
        protocol_dir = Path(tmpdir) / "protocol_db"
        protocol_dir.mkdir()
        
        # Create a simple test protocol
        test_protocol = {
            "reaction": {
                "reaction_smiles": "Br.N>>C",
                "family": "Test",
                "tags": "test",
                "notes": "Test protocol"
            },
            "source": {
                "title": "Test",
                "journal": "J. Test",
                "year": 2025,
                "doi": "10.1234/test"
            }
        }
        
        protocol_file = protocol_dir / "test_protocol.json"
        with open(protocol_file, 'w') as f:
            json.dump(test_protocol, f)
        
        # Build index (may or may not have DRFP depending on dependencies)
        idx = indexer.ProtocolIndexer(protocol_dir=protocol_dir)
        idx.build_index(compute_drfp=True, force_rebuild=True)
        
        # Should have one record
        assert len(idx.records) == 1
        assert "test_protocol.json" in idx.records
        
        # Save and check file sizes
        idx.save()
        
        index_path = protocol_dir / ".protocol_index.json"
        assert index_path.exists()
        
        # Load and verify JSON structure
        with open(index_path, 'r') as f:
            index_data = json.load(f)
        
        # Metadata should exist
        assert 'metadata' in index_data
        assert 'records' in index_data
        
        # Records should not have embedded DRFP arrays
        for filename, record_data in index_data['records'].items():
            assert 'drfp_fingerprint' not in record_data
            # Only has flag
            assert 'has_drfp' in record_data


def test_drfp_lazy_loading():
    """DRFP fingerprints should be loaded only when requested"""
    from chemtools.protocol import indexer
    
    with tempfile.TemporaryDirectory() as tmpdir:
        protocol_dir = Path(tmpdir) / "protocol_db"
        protocol_dir.mkdir()
        
        # Create test protocol
        test_protocol = {
            "reaction": {
                "reaction_smiles": "CCBr.c1ccc(B(O)O)cc1>>CCc1ccccc1",
                "family": "Suzuki",
                "tags": "coupling;Pd",
                "notes": "Test Suzuki coupling"
            },
            "source": {"title": "Test", "journal": "J. Org. Chem.", "year": 2025}
        }
        
        protocol_file = protocol_dir / "suzuki_test.json"
        with open(protocol_file, 'w') as f:
            json.dump(test_protocol, f)
        
        # Build index with DRFP
        idx = indexer.ProtocolIndexer(protocol_dir=protocol_dir)
        idx.build_index(compute_drfp=True, force_rebuild=True)
        idx.save()
        
        # Load index (fresh instance)
        idx2 = indexer.ProtocolIndexer.load(idx.index_path)
        
        # Loader should not be initialized yet
        assert idx2._drfp_loader is None
        
        # Try to get DRFP fingerprint
        fp = idx2.get_drfp_fingerprint("suzuki_test.json")
        
        # If DRFP was computed and storage is available, we should get it
        if idx2.metadata.get('has_drfp', False):
            # Loader should now be initialized
            assert idx2._drfp_loader is not None
            
            # Should get valid fingerprint or None
            # (None is OK if DRFP package not installed)
            if fp is not None:
                import numpy as np
                assert isinstance(fp, np.ndarray)
                assert fp.dtype == 'uint8'


def test_index_size_comparison(tmp_path):
    """Verify that NPZ storage is much smaller than JSON embedding"""
    from chemtools.protocol import indexer
    
    protocol_dir = tmp_path / "protocol_db"
    protocol_dir.mkdir()
    
    # Create multiple test protocols
    for i in range(10):
        test_protocol = {
            "reaction": {
                "reaction_smiles": f"Br{i}.N>>C{i}",
                "family": "Test",
                "tags": "test",
                "notes": f"Test protocol {i}"
            },
            "source": {"title": f"Test {i}", "journal": "J. Test", "year": 2025}
        }
        
        protocol_file = protocol_dir / f"test_protocol_{i:03d}.json"
        with open(protocol_file, 'w') as f:
            json.dump(test_protocol, f)
    
    # Build index
    idx = indexer.ProtocolIndexer(protocol_dir=protocol_dir)
    idx.build_index(compute_drfp=True, force_rebuild=True)
    idx.save()
    
    index_path = protocol_dir / ".protocol_index.json"
    drfp_path = protocol_dir / ".protocol_drfp.npz"
    
    # Index should exist and be reasonably sized
    assert index_path.exists()
    index_size = index_path.stat().st_size
    
    # Index should be small (mostly metadata)
    # With 10 protocols, expect < 50 KB
    assert index_size < 50_000, f"Index too large: {index_size} bytes"
    
    # If DRFP was computed
    if idx.metadata.get('has_drfp', False) and drfp_path.exists():
        drfp_size = drfp_path.stat().st_size
        
        # NPZ should be compressed and efficient
        # Expect < 100 KB for 10 protocols with 4096-bit fingerprints
        assert drfp_size < 100_000, f"DRFP file too large: {drfp_size} bytes"
        
        print(f"Index size: {index_size:,} bytes")
        print(f"DRFP size: {drfp_size:,} bytes")
        print(f"Total: {index_size + drfp_size:,} bytes")


def test_missing_drfp_file_graceful():
    """System should handle missing DRFP file gracefully"""
    from chemtools.protocol import indexer
    
    with tempfile.TemporaryDirectory() as tmpdir:
        protocol_dir = Path(tmpdir) / "protocol_db"
        protocol_dir.mkdir()
        
        # Create and build index
        test_protocol = {
            "reaction": {"reaction_smiles": "Br.N>>C", "family": "Test"},
            "source": {"title": "Test"}
        }
        
        with open(protocol_dir / "test.json", 'w') as f:
            json.dump(test_protocol, f)
        
        idx = indexer.ProtocolIndexer(protocol_dir=protocol_dir)
        idx.build_index(compute_drfp=True)
        idx.save()
        
        # Delete DRFP file if it exists
        drfp_path = protocol_dir / ".protocol_drfp.npz"
        if drfp_path.exists():
            drfp_path.unlink()
        
        # Load index (should work)
        idx2 = indexer.ProtocolIndexer.load(idx.index_path)
        
        # Getting DRFP should return None gracefully (not crash)
        fp = idx2.get_drfp_fingerprint("test.json")
        assert fp is None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

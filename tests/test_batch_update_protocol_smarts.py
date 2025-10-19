"""
Tests for batch protocol SMARTS updater

Tests the batch_update_protocol_smarts.py tool that generates and updates
reaction_smarts_applicability in protocol JSON files.
"""

import json
import pytest
import tempfile
from pathlib import Path

from chemtools.protocol.batch_update_protocol_smarts import (
    ProtocolSmartsUpdater,
    ProcessingResult
)


@pytest.fixture
def temp_protocol_dir(tmp_path):
    """Create temporary protocol directory with test files"""
    protocol_dir = tmp_path / "protocols"
    protocol_dir.mkdir()
    return protocol_dir


@pytest.fixture
def sample_protocol():
    """Sample protocol JSON structure"""
    return {
        "source": {
            "title": "Test Protocol",
            "journal": "Test Journal",
            "year": 2024
        },
        "reaction": {
            "reaction_smiles": "CCCCI.B(B1OC(C)(C)C(C)(C)O1)B2OC(C)(C)C(C)(C)O2>>CCCB1OC(C)(C)C(C)(C)O1",
            "family": "Test_Borylation",
            "notes": "Test reaction"
        }
    }


def test_process_single_protocol_file(temp_protocol_dir, sample_protocol):
    """Test processing a single protocol file"""
    # Write sample protocol
    protocol_file = temp_protocol_dir / "test_protocol.json"
    with open(protocol_file, 'w', encoding='utf-8') as f:
        json.dump(sample_protocol, f)
    
    # Process it
    updater = ProtocolSmartsUpdater(temp_protocol_dir, dry_run=False)
    result = updater.process_protocol_file(protocol_file)
    
    # Check result
    assert result.success
    assert result.filename == "test_protocol.json"
    assert "Pattern generated and applied" in result.message
    
    # Check file was updated
    with open(protocol_file, 'r', encoding='utf-8') as f:
        updated = json.load(f)
    
    assert "reaction_smarts_applicability" in updated["reaction"]
    assert "core" in updated["reaction"]["reaction_smarts_applicability"]
    assert "guards_forbid" in updated["reaction"]["reaction_smarts_applicability"]
    
    # Pattern should be chemistry-aware
    pattern = updated["reaction"]["reaction_smarts_applicability"]
    assert "[CX4;H2,H3]-[I]" in pattern["core"]  # Primary alkyl iodide pattern
    assert len(pattern["guards_forbid"]) > 0  # Should have guard patterns


def test_dry_run_mode(temp_protocol_dir, sample_protocol):
    """Test dry run mode doesn't modify files"""
    # Write sample protocol
    protocol_file = temp_protocol_dir / "test_protocol.json"
    with open(protocol_file, 'w', encoding='utf-8') as f:
        json.dump(sample_protocol, f)
    
    # Process in dry run mode
    updater = ProtocolSmartsUpdater(temp_protocol_dir, dry_run=True)
    result = updater.process_protocol_file(protocol_file)
    
    # Check result is successful
    assert result.success
    
    # Check file was NOT modified
    with open(protocol_file, 'r', encoding='utf-8') as f:
        unchanged = json.load(f)
    
    assert "reaction_smarts_applicability" not in unchanged["reaction"]


def test_overwrite_existing_pattern(temp_protocol_dir, sample_protocol):
    """Test that existing patterns are overwritten"""
    # Add existing pattern
    sample_protocol["reaction"]["reaction_smarts_applicability"] = {
        "core": "[C]-[I]>>old_pattern",
        "guards_forbid": ["old_guard"]
    }
    
    # Write protocol
    protocol_file = temp_protocol_dir / "test_protocol.json"
    with open(protocol_file, 'w', encoding='utf-8') as f:
        json.dump(sample_protocol, f)
    
    # Process it
    updater = ProtocolSmartsUpdater(temp_protocol_dir, dry_run=False)
    result = updater.process_protocol_file(protocol_file)
    
    # Check result
    assert result.success
    assert result.old_pattern is not None  # Should detect old pattern
    
    # Check pattern was updated
    with open(protocol_file, 'r', encoding='utf-8') as f:
        updated = json.load(f)
    
    new_pattern = updated["reaction"]["reaction_smarts_applicability"]
    assert new_pattern["core"] != "[C]-[I]>>old_pattern"  # Should be different
    assert "old_guard" not in new_pattern["guards_forbid"]  # Should have new guards


def test_missing_reaction_smiles(temp_protocol_dir):
    """Test handling of protocol without reaction_smiles"""
    protocol = {
        "source": {"title": "Test"},
        "reaction": {
            "family": "Test"
            # Missing reaction_smiles
        }
    }
    
    protocol_file = temp_protocol_dir / "test_protocol.json"
    with open(protocol_file, 'w', encoding='utf-8') as f:
        json.dump(protocol, f)
    
    # Process it
    updater = ProtocolSmartsUpdater(temp_protocol_dir, dry_run=False)
    result = updater.process_protocol_file(protocol_file)
    
    # Check result
    assert not result.success
    assert "No reaction_smiles found" in result.message


def test_invalid_json_file(temp_protocol_dir):
    """Test handling of invalid JSON"""
    protocol_file = temp_protocol_dir / "invalid.json"
    with open(protocol_file, 'w') as f:
        f.write("{invalid json content")
    
    # Process it
    updater = ProtocolSmartsUpdater(temp_protocol_dir, dry_run=False)
    result = updater.process_protocol_file(protocol_file)
    
    # Check result
    assert not result.success
    assert "JSON parsing error" in result.message


def test_batch_processing(temp_protocol_dir):
    """Test processing multiple protocol files"""
    # Create 3 test protocols
    protocols = [
        {
            "reaction": {
                "reaction_smiles": "CCCCI.B(B1OC(C)(C)C(C)(C)O1)B2OC(C)(C)C(C)(C)O2>>CCCB1OC(C)(C)C(C)(C)O1"
            }
        },
        {
            "reaction": {
                "reaction_smiles": "Brc1ccccc1.CNc2ccccc2>>CNc1ccccc1c2ccccc2"
            }
        },
        {
            "reaction": {
                "family": "No SMILES"
                # Missing reaction_smiles - should fail
            }
        }
    ]
    
    for i, protocol in enumerate(protocols):
        with open(temp_protocol_dir / f"protocol_{i}.json", 'w') as f:
            json.dump(protocol, f)
    
    # Process all
    updater = ProtocolSmartsUpdater(temp_protocol_dir, dry_run=False)
    results = updater.process_all_protocols()
    
    # Check results
    assert len(results) == 3
    successful = [r for r in results if r.success]
    failed = [r for r in results if not r.success]
    
    assert len(successful) == 2
    assert len(failed) == 1


def test_chemistry_awareness(temp_protocol_dir):
    """Test that generated patterns are chemistry-aware"""
    # Test different substrate types
    test_cases = [
        {
            "name": "primary_alkyl_iodide",
            "smiles": "CCCCI.B>>CCCB",
            "expected_in_core": "[CX4;H2,H3]-[I]",
            "expected_in_guards": "[CX4;H1]-[I]"  # Exclude secondary
        },
        {
            "name": "aryl_bromide",
            "smiles": "Brc1ccccc1.N>>Nc1ccccc1",
            "expected_in_core": "c-[Br]",
            "expected_in_guards": "[CX4]-[Br]"  # Should exclude aliphatic halides
        }
    ]
    
    for tc in test_cases:
        protocol = {
            "reaction": {
                "reaction_smiles": tc["smiles"],
                "family": tc["name"]
            }
        }
        
        protocol_file = temp_protocol_dir / f"{tc['name']}.json"
        with open(protocol_file, 'w') as f:
            json.dump(protocol, f)
        
        updater = ProtocolSmartsUpdater(temp_protocol_dir, dry_run=False)
        result = updater.process_protocol_file(protocol_file)
        
        assert result.success
        
        # Check pattern
        with open(protocol_file, 'r') as f:
            updated = json.load(f)
        
        pattern = updated["reaction"]["reaction_smarts_applicability"]
        assert tc["expected_in_core"] in pattern["core"]
        
        if "expected_in_guards" in tc:
            guards_str = " ".join(pattern["guards_forbid"])
            assert tc["expected_in_guards"] in guards_str


def test_buchwald_hartwig_pattern(temp_protocol_dir):
    """Test Buchwald-Hartwig coupling pattern generation"""
    protocol = {
        "reaction": {
            "reaction_smiles": "Brc1ccccc1.CNc2ccccc2>>CNc1ccccc1c2ccccc2",
            "family": "Buchwald_Hartwig"
        }
    }
    
    protocol_file = temp_protocol_dir / "buchwald_hartwig.json"
    with open(protocol_file, 'w') as f:
        json.dump(protocol, f)
    
    updater = ProtocolSmartsUpdater(temp_protocol_dir, dry_run=False)
    result = updater.process_protocol_file(protocol_file)
    
    assert result.success
    
    # Check pattern distinguishes aryl vs alkyl and aniline vs aliphatic
    with open(protocol_file, 'r') as f:
        updated = json.load(f)
    
    pattern = updated["reaction"]["reaction_smarts_applicability"]
    core = pattern["core"]
    guards = " ".join(pattern["guards_forbid"])
    
    # Should have aryl bromide pattern
    assert "c-[Br]" in core
    
    # Should have aniline pattern with amide exclusion
    assert "[NX3" in core
    assert "!$(NC=O)" in core or "!$(NC=O)" in guards


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

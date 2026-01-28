"""Test simplified featurizer system without legacy support."""

from chemtools.featurizers.formatters.molecule import featurize_molecule
from chemtools.featurizers.formatters.reaction import featurize_reaction


def test_molecule_core():
    """Test core molecule format."""
    result = featurize_molecule('c1ccccc1Br')
    print("✓ Core molecule format:")
    print(f"  Schema: {result.get('schema_version')}")
    print(f"  Fields ({len(result)}): {sorted(result.keys())}")
    assert result.get('schema_version') == 'v2'
    assert result.get('kind') == 'molecule'
    assert 'motifs' in result
    assert 'properties' in result
    print()


def test_molecule_extended():
    """Test extended molecule format."""
    result = featurize_molecule('c1ccccc1Br', options={'detailed': True})
    print("✓ Extended molecule format:")
    print(f"  Schema: {result.get('schema_version')}")
    print(f"  Fields ({len(result)}): {sorted(result.keys())}")
    print(f"  Has extended: {'extended' in result}")
    assert result.get('schema_version') == 'v2'
    assert result.get('kind') == 'molecule'
    assert 'extended' in result
    print()


def test_reaction_core():
    """Test core reaction format."""
    result = featurize_reaction('c1ccc(Br)cc1.CC(=O)O>>CC(=O)c1ccccc1')
    print("✓ Core reaction format:")
    print(f"  Schema: {result.get('schema_version')}")
    print(f"  Fields ({len(result)}): {sorted(result.keys())}")
    print(f"  Reaction type: {result.get('reaction_type')}")
    assert result.get('schema_version') == 'v2'
    assert result.get('kind') == 'reaction'
    assert 'reactants' in result
    assert 'products' in result
    print()


def test_reaction_extended():
    """Test extended reaction format."""
    result = featurize_reaction('c1ccc(Br)cc1.CC(=O)O>>CC(=O)c1ccccc1', options={'detailed': True})
    print("✓ Extended reaction format:")
    print(f"  Schema: {result.get('schema_version')}")
    print(f"  Fields ({len(result)}): {sorted(result.keys())}")
    print(f"  Has extended: {'extended' in result}")
    assert result.get('schema_version') == 'v2'
    assert result.get('kind') == 'reaction'
    assert 'extended' in result
    print()


def test_no_legacy_option():
    """Verify legacy option is ignored."""
    result = featurize_molecule('c1ccccc1Br', options={'legacy': True})
    print("✓ Legacy option ignored:")
    print(f"  Schema: {result.get('schema_version')} (always v2)")
    assert result.get('schema_version') == 'v2'
    print()


if __name__ == '__main__':
    print("=" * 72)
    print("Testing Clean Two-Tier System (No Legacy Support)")
    print("=" * 72)
    print()
    
    test_molecule_core()
    test_molecule_extended()
    test_reaction_core()
    test_reaction_extended()
    test_no_legacy_option()
    
    print("=" * 72)
    print("✅ ALL TESTS PASSED - Clean system working!")
    print("=" * 72)

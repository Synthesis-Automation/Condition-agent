"""Test script to verify taxonomy-driven validation is working."""

from chemtools.featurizers.formatters.detection_validation import (
    validate_detection_with_reacted_motifs,
)


def test_suzuki():
    """Test Suzuki detection with HeteroAr-I."""
    result = validate_detection_with_reacted_motifs(
        "Unknown", 
        0.3, 
        ["Ar-B(OH)2", "HeteroAr-I"], 
        ["Ar-Ar"]
    )
    assert result["reaction_type"] == "Suzuki_miyaura"
    assert result["confidence"] == 0.95
    print("✓ Suzuki with HeteroAr-I: PASSED")


def test_suzuki_heteroar_boron():
    """Test Suzuki with HeteroAr boron."""
    result = validate_detection_with_reacted_motifs(
        "Unknown",
        0.3,
        ["HeteroAr-B(OH)2", "Ar-Br"],
        ["Ar-Ar"]
    )
    assert result["reaction_type"] == "Suzuki_miyaura"
    assert result["confidence"] == 0.95
    print("✓ Suzuki with HeteroAr-B(OH)2: PASSED")


def test_cn_coupling():
    """Test C-N coupling detection."""
    result = validate_detection_with_reacted_motifs(
        "Unknown",
        0.3,
        ["Ar-Br", "Ar-NH2"],
        ["Ar-NHR"]
    )
    assert result["reaction_type"] == "C_N_Coupling"
    assert result["confidence"] == 0.95
    print("✓ C-N coupling: PASSED")


def test_cn_coupling_heteroar():
    """Test C-N coupling with HeteroAr."""
    result = validate_detection_with_reacted_motifs(
        "Unknown",
        0.3,
        ["HeteroAr-Br", "Alkyl-NH2"],
        ["Ar-NHR"]
    )
    assert result["reaction_type"] == "C_N_Coupling"
    assert result["confidence"] == 0.95
    print("✓ C-N coupling with HeteroAr-Br: PASSED")


def test_no_correction_needed():
    """Test that correct detection is not changed."""
    result = validate_detection_with_reacted_motifs(
        "Suzuki_miyaura",
        0.95,
        ["Ar-B(OH)2", "Ar-Br"],
        ["Ar-Ar"]
    )
    assert result["reaction_type"] == "Suzuki_miyaura"
    assert result["validation_method"] == "slot_based_confirmed"
    print("✓ No correction needed: PASSED")


def test_organometallic_exclusion():
    """Test Arylation_Ar_H exclusion with organometallics."""
    result = validate_detection_with_reacted_motifs(
        "Arylation_Ar_H",
        0.85,
        ["Ar-B(OH)2", "Ar-H"],
        ["Ar-Ar"]
    )
    assert result["reaction_type"] == "Unknown"
    assert result["confidence"] == 0.3
    print("✓ Organometallic exclusion: PASSED")


def test_sonogashira():
    """Test Sonogashira detection."""
    result = validate_detection_with_reacted_motifs(
        "Unknown",
        0.3,
        ["Ar-I", "Alkyl-Alkynyl_terminal"],
        ["Ar-Alkynyl"]
    )
    assert result["reaction_type"] == "Sonogashira"
    assert result["confidence"] == 0.95
    print("✓ Sonogashira: PASSED")


def test_sonogashira_heteroar():
    """Test Sonogashira with HeteroAr."""
    result = validate_detection_with_reacted_motifs(
        "Unknown",
        0.3,
        ["HeteroAr-Br", "Ar-Alkynyl_terminal"],
        ["HeteroAr-Alkynyl"]
    )
    assert result["reaction_type"] == "Sonogashira"
    assert result["confidence"] == 0.95
    print("✓ Sonogashira with HeteroAr: PASSED")


if __name__ == "__main__":
    print("Testing taxonomy-driven validation...")
    print("=" * 60)
    
    test_suzuki()
    test_suzuki_heteroar_boron()
    test_cn_coupling()
    test_cn_coupling_heteroar()
    test_no_correction_needed()
    test_organometallic_exclusion()
    test_sonogashira()
    test_sonogashira_heteroar()
    
    print("=" * 60)
    print("✅ All taxonomy-driven validation tests PASSED!")

import pytest

from chemtools.featurizers.formatters import reaction as reaction_fmt
from chemtools.util import rdkit_helpers


def test_alt_key_promoted_when_primary_unmatched() -> None:
    primary = "Ar-Br -> Ar-XYZ || []"
    alts = [
        "Alkyl-Bpin|Ar-Br -> Ar-Alkyl || []",
        "Ar-Br|Ar-Bpin -> Ar-Ar || []",
    ]

    new_primary, new_alts = reaction_fmt._select_primary_reaction_key(primary, alts)

    assert new_primary == alts[0]
    assert new_alts[0] == primary
    assert primary in new_alts


def test_alt_key_not_promoted_when_primary_matches() -> None:
    primary = "Alkyl-Bpin|Ar-Br -> Ar-Alkyl || []"
    alts = ["Ar-Br -> Ar-XYZ || []"]

    new_primary, new_alts = reaction_fmt._select_primary_reaction_key(primary, alts)

    assert new_primary == primary
    assert new_alts == alts


def test_multi_event_reaction_promotes_suzuki_alt() -> None:
    if not rdkit_helpers.rdkit_available():
        pytest.skip("rdkit not available")
    rxn = "COC(OC)c1ccccc1Br.C=CC(CCCC)B1OC(C)(C)C(C)(C)O1>>C=CC(CCCC)c1ccccc1C=O"
    result = reaction_fmt.featurize_reaction(rxn)

    assert result.get("reaction_key")
    alt_keys = result.get("reaction_keys_alt") or []
    assert alt_keys

    from chemtools.reaction_key_matcher_v2 import detect_from_reaction_key_v2
    keys = [result["reaction_key"]] + alt_keys
    # At least one key should match Suzuki (primary or alt).
    suzuki_found = False
    for key in keys:
        reacted, formed, spectators = reaction_fmt._parse_reaction_key_parts(key)
        match, _ = detect_from_reaction_key_v2(reacted, formed, spectators)
        if match and match.reaction_type == "Suzuki_miyaura":
            suzuki_found = True
            break
    assert suzuki_found

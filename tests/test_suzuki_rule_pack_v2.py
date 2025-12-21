import pytest

from chemtools.rule.analyzer import FeatureAnalyzer
from chemtools.rule.database import RuleDatabase
from chemtools.util.rdkit_helpers import rdkit_available


@pytest.mark.skipif(not rdkit_available(), reason="RDKit not available")
def test_suzuki_rule_pack_applies_and_matches():
    analyzer = FeatureAnalyzer()
    features = analyzer.analyze_reaction("Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1")

    db = RuleDatabase.from_file("data/rule_db_v2/Suzuki_db.json")
    assert db.check_applies(features) is True
    assert db.find_matching_rule(features) is not None


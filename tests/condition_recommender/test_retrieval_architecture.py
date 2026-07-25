"""Architecture regressions for condition-retrieval entry points."""

from importlib.util import find_spec

import condition_recommender


def test_suzuki_specific_structural_retrieval_is_removed() -> None:
    """Only the generic structural API remains publicly available."""
    assert "recommend_generic_conditions" in condition_recommender.__all__
    assert "recommend_conditions" not in condition_recommender.__all__
    assert not hasattr(condition_recommender, "recommend_conditions")
    assert find_spec("condition_recommender.api") is None
    assert find_spec("condition_recommender.indexing") is None
    assert find_spec("condition_recommender.retrieval") is None

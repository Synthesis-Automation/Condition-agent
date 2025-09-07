from chemtools.precedent import knn


def test_knn_no_precedents_fast():
    # Minimal, fast negative case using the tiny demo dataset
    out = knn(
        family="Suzuki_CC",  # not present in sample dataset
        features={"bin": "LG:Br|NUC:aniline"},
        k=5,
        relax={},
    )
    assert out.get("error") == "NO_PRECEDENTS"
    assert out.get("support") == 0
    assert out.get("precedents") == []

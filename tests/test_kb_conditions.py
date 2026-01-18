from chem_assistant.kb_conditions import search_condition_summaries


def test_kb_condition_search_returns_rows() -> None:
    payload = search_condition_summaries(
        "5-membered heterocyclic boronates", top_k=3
    )
    results = payload.get("results") or []
    assert results
    first = results[0]
    condition = first.get("condition") or {}
    assert "catalyst" in condition
    assert "solvent" in condition

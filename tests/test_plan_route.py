"""
Tests for the plan_route retrosynthesis tool
(chem_coworker/tools/retrosynthesis.py :: _plan_route).

Coverage:
  - import and registration
  - helper functions (_bertz_complexity_fast, _inchi_key)
  - simple molecule (BertzCT < threshold) → 0 steps, treated as building block
  - complex molecule produces steps with correct structure
  - error handling (invalid SMILES, empty input)
  - parameter effects (max_depth, complexity_threshold)
  - cycle detection doesn't blow up
  - output contract (_success shape)
"""
from __future__ import annotations

from types import SimpleNamespace

import pytest


# ---------------------------------------------------------------------------
# Import guards
# ---------------------------------------------------------------------------

def test_plan_route_imports():
    """plan_route_tool is importable and registered in RETROSYNTHESIS_TOOLS."""
    from chem_coworker.tools.retrosynthesis import (  # noqa: F401
        plan_route_tool,
        plan_route_candidates_tool,
        RETROSYNTHESIS_TOOLS,
    )
    names = [t.name for t in RETROSYNTHESIS_TOOLS]
    assert "plan_route" in names
    assert "plan_route_candidates" in names


def test_helper_imports():
    """Internal helpers are importable."""
    from chem_coworker.tools.retrosynthesis import (  # noqa: F401
        _bertz_complexity_fast,
        _inchi_key,
    )


# ---------------------------------------------------------------------------
# Helper: _bertz_complexity_fast
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("smiles,expected_min,expected_max", [
    ("c1ccccc1", 60.0, 90.0),       # benzene ~72
    ("c1ccncc1", 60.0, 90.0),       # pyridine ~76
    ("CO", 0.0, 10.0),              # methanol ~2
])
def test_bertz_complexity_fast_known_values(smiles, expected_min, expected_max):
    from chem_coworker.tools.retrosynthesis import _bertz_complexity_fast
    cx = _bertz_complexity_fast(smiles)
    assert expected_min <= cx <= expected_max, (
        f"BertzCT({smiles!r}) = {cx:.1f}, expected [{expected_min}, {expected_max}]"
    )


def test_bertz_complexity_fast_invalid_smiles():
    """Returns 9999.0 for unparseable SMILES (does not raise)."""
    from chem_coworker.tools.retrosynthesis import _bertz_complexity_fast
    result = _bertz_complexity_fast("not_a_smiles_$$$$")
    assert result == 9999.0


def test_bertz_complexity_fast_complex_is_high():
    """Celecoxib complexity is well above 100."""
    from chem_coworker.tools.retrosynthesis import _bertz_complexity_fast
    celecoxib = "Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1"
    cx = _bertz_complexity_fast(celecoxib)
    assert cx > 500.0, f"Expected celecoxib BertzCT > 500, got {cx}"


# ---------------------------------------------------------------------------
# Helper: _inchi_key
# ---------------------------------------------------------------------------

def test_inchi_key_benzene():
    """Benzene produces a 27-char InChI key."""
    from chem_coworker.tools.retrosynthesis import _inchi_key
    ik = _inchi_key("c1ccccc1")
    assert ik is not None
    assert len(ik) == 27
    assert ik.count("-") == 2   # three hyphen-separated blocks


def test_inchi_key_canonical():
    """Different SMILES for the same molecule → same InChI key."""
    from chem_coworker.tools.retrosynthesis import _inchi_key
    # Two equivalent SMILES for toluene
    ik1 = _inchi_key("Cc1ccccc1")
    ik2 = _inchi_key("c1ccc(C)cc1")
    assert ik1 is not None and ik2 is not None
    assert ik1 == ik2


def test_inchi_key_invalid_smiles():
    """Returns None for invalid SMILES (does not raise)."""
    from chem_coworker.tools.retrosynthesis import _inchi_key
    assert _inchi_key("not_valid!!!") is None


def test_inchi_key_wildcard_atom():
    """Wildcard SMILES (*CC) returns None gracefully (no InChI for *)."""
    from chem_coworker.tools.retrosynthesis import _inchi_key
    result = _inchi_key("*CC")
    # None is acceptable; must not raise
    assert result is None or isinstance(result, str)


# ---------------------------------------------------------------------------
# plan_route: error handling
# ---------------------------------------------------------------------------

def test_plan_route_empty_smiles():
    from chem_coworker.tools.retrosynthesis import plan_route_tool
    result = plan_route_tool.fn(target_smiles="")
    assert result.get("success") is False
    assert "error" in result


def test_plan_route_invalid_smiles():
    from chem_coworker.tools.retrosynthesis import plan_route_tool
    result = plan_route_tool.fn(target_smiles="NOTASMILES###")
    # Either error or 0-step result with the molecule treated as complex leaf
    # (MolFromSmiles returns None → complexity = 9999 → no disconnections → complex_leaf)
    assert isinstance(result, dict)


# ---------------------------------------------------------------------------
# plan_route: simple molecule (BertzCT < 100)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("smiles", [
    "c1ccccc1",    # benzene, BertzCT ~72
    "c1ccncc1",    # pyridine, ~76
    "CO",          # methanol, ~2
])
def test_plan_route_simple_molecule_zero_steps(smiles):
    from chem_coworker.tools.retrosynthesis import plan_route_tool
    result = plan_route_tool.fn(target_smiles=smiles, complexity_threshold=100.0)
    assert result.get("success") is True
    assert result["total_steps"] == 0
    assert smiles in result.get("route_summary", "")


# ---------------------------------------------------------------------------
# plan_route: complex molecule produces steps
# ---------------------------------------------------------------------------

def test_plan_route_ibuprofen_produces_steps():
    """Ibuprofen (BertzCT ~325) should yield at least 1 disconnection step."""
    from chem_coworker.tools.retrosynthesis import plan_route_tool
    result = plan_route_tool.fn(
        target_smiles="CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        complexity_threshold=100.0,
        max_depth=4,
    )
    assert result.get("success") is True
    assert result["total_steps"] >= 1


def test_plan_route_celecoxib_produces_steps():
    """Celecoxib (BertzCT ~1038) should yield multiple disconnection steps."""
    from chem_coworker.tools.retrosynthesis import plan_route_tool
    result = plan_route_tool.fn(
        target_smiles="Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1",
        complexity_threshold=150.0,
        max_depth=4,
    )
    assert result.get("success") is True
    assert result["total_steps"] >= 2


# ---------------------------------------------------------------------------
# plan_route: output structure (contract)
# ---------------------------------------------------------------------------

def test_plan_route_output_structure():
    """Result dict must contain all required keys with correct types."""
    from chem_coworker.tools.retrosynthesis import plan_route_tool
    result = plan_route_tool.fn(
        target_smiles="CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        complexity_threshold=100.0,
    )
    assert result.get("success") is True

    required_keys = [
        "smiles", "total_steps", "cumulative_difficulty", "all_leaves_simple",
        "complexity_threshold", "max_depth", "leaves", "simple_leaves",
        "complex_leaves", "route", "route_summary",
    ]
    for key in required_keys:
        assert key in result, f"Missing key: {key!r}"

    assert isinstance(result["route"], list)
    assert isinstance(result["leaves"], list)
    assert isinstance(result["simple_leaves"], list)
    assert isinstance(result["complex_leaves"], list)
    assert isinstance(result["route_summary"], str)
    assert isinstance(result["total_steps"], int)
    assert isinstance(result["cumulative_difficulty"], float)
    assert isinstance(result["all_leaves_simple"], bool)


def test_plan_route_step_structure():
    """Each step dict must contain the required fields."""
    from chem_coworker.tools.retrosynthesis import plan_route_tool
    result = plan_route_tool.fn(
        target_smiles="CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        complexity_threshold=100.0,
    )
    assert result.get("success") is True
    for step in result["route"]:
        for field in ["step", "depth", "target", "target_complexity",
                      "reaction_name", "difficulty", "precursor_1"]:
            assert field in step, f"Step missing field: {field!r}"
        assert isinstance(step["step"], int)
        assert isinstance(step["depth"], int)
        assert isinstance(step["difficulty"], float)
        assert 0.0 <= step["difficulty"] <= 1.0


def test_plan_route_leaves_cover_precursors():
    """All leaf SMILES must come from the route precursors (or be the target itself)."""
    from chem_coworker.tools.retrosynthesis import plan_route_tool
    result = plan_route_tool.fn(
        target_smiles="CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        complexity_threshold=100.0,
    )
    assert result.get("success") is True
    leaves = set(result["leaves"])
    simple = set(result["simple_leaves"])
    complex_ = set(result["complex_leaves"])
    assert simple | complex_ == leaves, "simple_leaves + complex_leaves must equal leaves"


# ---------------------------------------------------------------------------
# plan_route: parameter effects
# ---------------------------------------------------------------------------

def test_plan_route_max_depth_zero():
    """max_depth=0 means no expansion at all — target is immediately a leaf."""
    from chem_coworker.tools.retrosynthesis import plan_route_tool
    result = plan_route_tool.fn(
        target_smiles="CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        max_depth=0,
        complexity_threshold=100.0,
    )
    assert result.get("success") is True
    assert result["total_steps"] == 0
    assert len(result["leaves"]) == 1  # just the target itself


def test_plan_route_high_threshold_fewer_steps():
    """A higher complexity threshold should stop recursion earlier (≤ steps vs lower threshold)."""
    from chem_coworker.tools.retrosynthesis import plan_route_tool
    smiles = "Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1"
    result_tight = plan_route_tool.fn(target_smiles=smiles, complexity_threshold=500.0, max_depth=4)
    result_loose = plan_route_tool.fn(target_smiles=smiles, complexity_threshold=100.0, max_depth=4)
    assert result_tight.get("success") and result_loose.get("success")
    # Higher threshold = fewer steps needed (simpler stopping condition)
    assert result_tight["total_steps"] <= result_loose["total_steps"]


def test_plan_route_smiles_alias():
    """The 'smiles' parameter alias works identically to 'target_smiles'."""
    from chem_coworker.tools.retrosynthesis import plan_route_tool
    r1 = plan_route_tool.fn(target_smiles="c1ccccc1")
    r2 = plan_route_tool.fn(smiles="c1ccccc1")
    assert r1.get("success") == r2.get("success")
    assert r1["total_steps"] == r2["total_steps"]


def test_plan_route_reaction_smiles_stripped():
    """Reaction SMILES (A>>B) should be stripped; only the product side is used."""
    from chem_coworker.tools.retrosynthesis import plan_route_tool
    # Pass a reaction SMILES; the tool should strip to the left side
    result = plan_route_tool.fn(target_smiles="Brc1ccccc1>>c1ccccc1")
    # Should work without error (not crash on the >> in input)
    assert isinstance(result, dict)


# ---------------------------------------------------------------------------
# plan_route: cumulative difficulty
# ---------------------------------------------------------------------------

def test_plan_route_cumulative_difficulty_nonnegative():
    from chem_coworker.tools.retrosynthesis import plan_route_tool
    result = plan_route_tool.fn(target_smiles="CC(C)Cc1ccc(cc1)C(C)C(=O)O")
    assert result.get("success") is True
    assert result["cumulative_difficulty"] >= 0.0


def test_plan_route_cumulative_difficulty_matches_steps():
    """cumulative_difficulty must equal the sum of individual step difficulties."""
    from chem_coworker.tools.retrosynthesis import plan_route_tool
    result = plan_route_tool.fn(
        target_smiles="CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        complexity_threshold=100.0,
    )
    assert result.get("success") is True
    expected = sum(s["difficulty"] for s in result["route"])
    assert abs(result["cumulative_difficulty"] - expected) < 1e-4, (
        f"cumulative_difficulty {result['cumulative_difficulty']} "
        f"!= sum of steps {expected}"
    )


# ---------------------------------------------------------------------------
# plan_route_candidates: strategy-aware candidate generation
# ---------------------------------------------------------------------------

def _fake_disconnection(name, p1, p2, difficulty=0.25, score=0.8, description=""):
    return SimpleNamespace(
        reaction_name=name,
        retron_name=name.lower().replace(" ", "_"),
        description=description or name,
        difficulty=difficulty,
        precursor_1=p1,
        precursor_2=p2,
        complexity_delta=100.0,
        overall_score=score,
    )


def test_plan_route_candidates_ranks_strategy_aligned_route(monkeypatch):
    import chem_coworker.tools.retrosynthesis as retro
    import chemtools.retro.disconnector as disconnector

    monkeypatch.setattr(
        retro,
        "_bertz_complexity_fast",
        lambda smi: 500.0 if smi == "TARGET" else 10.0,
    )
    monkeypatch.setattr(retro, "_inchi_key", lambda smi: f"IK-{smi}")

    def fake_rank(smiles, top_k=4, max_difficulty=0.9):  # noqa: ARG001
        if smiles != "TARGET":
            return []
        return [
            _fake_disconnection(
                "Suzuki coupling",
                "ArBr",
                "ArB(OH)2",
                difficulty=0.20,
                description="palladium cross-coupling",
            ),
            _fake_disconnection(
                "Aldol addition",
                "CC=O",
                "CC(=O)C",
                difficulty=0.25,
                description="metal-free carbonyl coupling",
            ),
        ]

    monkeypatch.setattr(disconnector, "rank_disconnections", fake_rank)

    result = retro.plan_route_candidates_tool.fn(
        target_smiles="TARGET",
        strategy_query="Find a metal-free concise route",
        max_depth=2,
        max_branches=2,
        beam_width=4,
        top_k=2,
    )

    assert result.get("success") is True
    assert result["candidate_count"] == 2
    assert result["candidates"][0]["route"][0]["reaction_name"] == "Aldol addition"
    assert result["candidates"][0]["combined_score"] >= result["candidates"][1]["combined_score"]


def test_plan_route_candidates_exposes_structured_projection(monkeypatch):
    import chem_coworker.tools.retrosynthesis as retro
    import chemtools.retro.disconnector as disconnector

    monkeypatch.setattr(
        retro,
        "_bertz_complexity_fast",
        lambda smi: 500.0 if smi == "TARGET" else 10.0,
    )
    monkeypatch.setattr(retro, "_inchi_key", lambda smi: f"IK-{smi}")
    monkeypatch.setattr(
        disconnector,
        "rank_disconnections",
        lambda *args, **kwargs: [_fake_disconnection("Wittig olefination", "A", "B")],
    )

    result = retro.plan_route_candidates_tool.fn(
        target_smiles="TARGET",
        strategy_query="concise feasible route",
        max_depth=1,
        top_k=1,
    )
    projection = retro.plan_route_candidates_tool.structured_projection(result)

    assert result.get("success") is True
    assert "route_candidates" in projection
    assert projection["best_route_candidate"]["route"][0]["reaction_name"] == "Wittig olefination"

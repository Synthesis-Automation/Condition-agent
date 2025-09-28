from __future__ import annotations

import pathlib
import sys

# Ensure repository root is on sys.path for imports when tests run in isolation.
ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from chemtools.condition_rules import recommend_rule_based, get_family_defaults  # noqa: E402


def _base_features() -> dict:
    return {
        "electrophile": {"class": "aryl chloride"},
        "nucleophile": {"class": "primary aniline"},
    }


def test_rule_based_returns_expected_playbook():
    features = _base_features()
    job_ctx = {"base": "K2CO3", "mode": "batch"}
    recs = recommend_rule_based("Buchwald_CN", features, job_ctx=job_ctx, top_n=1)
    assert recs, "Expected at least one playbook recommendation"
    top = recs[0]
    assert top["playbook_id"] == "BH-ARCL-PRIM-ANILINE-GENERAL"
    assert top["suggested_recipe"]["base"] in {"NaOtBu", "K2CO3"}


def test_guard_penalizes_kotbu_for_electron_poor():
    features = {
        "electrophile": {"class": "aryl chloride"},
        "electrophile_electronics": "electron-poor",
        "nucleophile": {"class": "primary aniline"},
    }
    job_ctx = {"base": "KOtBu", "mode": "batch"}
    recs = recommend_rule_based("Buchwald_CN", features, job_ctx=job_ctx, top_n=1)
    guard_messages = recs[0]["guard_messages"]
    assert any(msg["id"] == "GUARD-SNAR-KOTBU" for msg in guard_messages)
    assert recs[0]["score_breakdown"]["guard_adjust"] < 0


def test_priors_provide_boost():
    features = {
        "electrophile": {"class": "aryl chloride"},
        "electrophile_electronics": "electron-poor",
        "ortho_sub_count": 1,
        "nucleophile_class": "aniline",
        "solvent": "DMF",
    }
    job_ctx = {"base": "DBU", "mode": "batch"}
    recs = recommend_rule_based("Buchwald_CN", features, job_ctx=job_ctx, top_n=1)
    assert recs[0]["score"] > 1.0


def test_family_defaults_available():
    defaults = get_family_defaults("Buchwald_CN")
    assert defaults["temperature_C"] == 100
    assert "solvent_priority" in defaults

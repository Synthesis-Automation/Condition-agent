from __future__ import annotations

from pathlib import Path

import pytest

from chemtools.taxonomy.new import template_eval
from chemtools.util.rdkit_helpers import rdkit_available


pytestmark = pytest.mark.skipif(not rdkit_available(), reason="RDKit not available or disabled")


POC_DIR = Path(__file__).resolve().parents[1] / "chemtools" / "taxonomy" / "new"


def test_expand_templates_matches_expected_smarts_in_spec() -> None:
    spec = template_eval.load_template_spec(POC_DIR / "smarts_templates.suzuki_buchwald.poc.json")
    generated = template_eval.expand_templates(spec.fragments, spec.templates, spec.expansions)

    by_token = {g.feature_token: g.smarts for g in generated}
    for exp in spec.expansions:
        if exp.expected_smarts is None:
            continue
        assert by_token[exp.feature_token] == exp.expected_smarts


def test_expand_templates_is_deterministic() -> None:
    spec = template_eval.load_template_spec(POC_DIR / "smarts_templates.suzuki_buchwald.poc.json")
    g1 = template_eval.expand_templates(spec.fragments, spec.templates, spec.expansions)
    g2 = template_eval.expand_templates(spec.fragments, spec.templates, spec.expansions)
    assert g1 == g2


def test_generated_smarts_compile() -> None:
    spec = template_eval.load_template_spec(POC_DIR / "smarts_templates.suzuki_buchwald.poc.json")
    generated = template_eval.expand_templates(spec.fragments, spec.templates, spec.expansions)
    successes, failures = template_eval.compile_generated_smarts(generated)
    assert successes == len(generated)
    assert failures == ()


def test_cross_check_against_atomic_features() -> None:
    spec = template_eval.load_template_spec(POC_DIR / "smarts_templates.suzuki_buchwald.poc.json")
    generated = template_eval.expand_templates(spec.fragments, spec.templates, spec.expansions)
    atomic = template_eval.load_atomic_feature_smarts(POC_DIR / "calculable_features.atomic.suzuki_buchwald.poc.json")

    missing_tokens, missing_smarts, dup_tokens, dup_smarts = template_eval.cross_check_against_atomic_features(
        generated, atomic
    )
    assert missing_tokens == ()
    assert missing_smarts == ()
    assert dup_tokens == ()
    assert dup_smarts == ()


def test_functional_match_tests_all_pass() -> None:
    spec = template_eval.load_template_spec(POC_DIR / "smarts_templates.suzuki_buchwald.poc.json")
    generated = template_eval.expand_templates(spec.fragments, spec.templates, spec.expansions)
    results = template_eval.run_functional_match_tests(generated)
    failures = [r for r in results if not r.passed]
    assert failures == []


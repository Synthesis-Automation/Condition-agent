"""
Tests for the forward synthesis module.

Covers:
  1. chemtools.forward.reactor.ReactantAnalyzer — profile() and find_compatible_templates()
  2. chemtools.forward.reactor.ForwardReactor — generate()
  3. chemtools.forward.scoring.score_products()
  4. chem_coworker.tools.forward_synthesis — all 8 agent tools
  5. chem_coworker.workflow — forward_synthesis workflow registration
  6. retro hte_templates.json — forward_smarts field presence and validity
"""
import pytest
from types import SimpleNamespace

from chemtools.forward import ForwardReactor, ReactantAnalyzer, score_products
from chemtools.forward.contracts import ProductPrediction
from chem_coworker.tools.forward_synthesis import (
    _generate_products,
    _identify_reactions,
    _inspect_reactants,
    _rank_products,
    _find_forward_precedent,
    _search_reactant_precedent,
    _recommend_forward_conditions,
    _plan_forward_route,
    FORWARD_SYNTHESIS_TOOLS,
)
from chem_coworker.tools._helpers import _validate_reaction_smiles


# ---------------------------------------------------------------------------
# 1. retro hte_templates.json — forward_smarts field
# ---------------------------------------------------------------------------

class TestHTETemplatesForwardSmarts:
    def setup_method(self):
        import json, pathlib
        p = pathlib.Path("chemtools/retro/data/hte_templates.json")
        self.templates = json.loads(p.read_text(encoding="utf-8"))["templates"]

    def test_all_templates_have_forward_smarts(self):
        missing = [t["name"] for t in self.templates if "forward_smarts" not in t]
        assert missing == [], f"Templates missing forward_smarts: {missing}"

    def test_all_templates_have_n_reactants(self):
        missing = [t["name"] for t in self.templates if "n_reactants" not in t]
        assert missing == [], f"Templates missing n_reactants: {missing}"

    def test_forward_smarts_contain_arrow(self):
        for t in self.templates:
            fwd = t.get("forward_smarts", "")
            assert ">>" in fwd, f"Template {t['name']} forward_smarts missing '>>': {fwd!r}"

    def test_n_reactants_valid_values(self):
        for t in self.templates:
            assert t.get("n_reactants") in (1, 2), (
                f"Template {t['name']} has invalid n_reactants={t.get('n_reactants')}"
            )

    def test_key_templates_are_bimolecular(self):
        names_2 = {"buchwald_hartwig", "suzuki_miyaura", "amide_coupling",
                   "cuaac_triazole", "snar_amination"}
        by_name = {t["name"]: t for t in self.templates}
        for name in names_2:
            assert by_name[name]["n_reactants"] == 2, f"{name} should be bimolecular"

    def test_key_templates_are_unimolecular(self):
        names_1 = {"ketone_from_nabh4", "sandmeyer_bromide", "wacker_oxidation",
                   "aldehyde_from_oxidation", "deoxy_fluorination"}
        by_name = {t["name"]: t for t in self.templates}
        for name in names_1:
            assert by_name[name]["n_reactants"] == 1, f"{name} should be unimolecular"


# ---------------------------------------------------------------------------
# 2. ReactantAnalyzer.profile()
# ---------------------------------------------------------------------------

class TestReactantAnalyzerProfile:
    def setup_method(self):
        self.analyzer = ReactantAnalyzer()

    def test_aryl_bromide_detected(self):
        profile = self.analyzer.profile("Brc1ccccc1")
        assert "aryl_bromide" in profile.functional_groups

    def test_aryl_bromide_is_electrophile(self):
        profile = self.analyzer.profile("Brc1ccccc1")
        assert "aryl_bromide" in profile.electrophilic_sites

    def test_aniline_amine_detected(self):
        profile = self.analyzer.profile("Nc1ccccc1")
        assert "aniline" in profile.functional_groups or "primary_amine" in profile.functional_groups

    def test_boronic_acid_detected(self):
        profile = self.analyzer.profile("OB(O)c1ccccc1")
        assert "arylboronic_acid" in profile.functional_groups

    def test_leaving_group_aryl_bromide(self):
        profile = self.analyzer.profile("Brc1ccccc1")
        assert profile.leaving_group == "Br"

    def test_leaving_group_aryl_fluoride(self):
        profile = self.analyzer.profile("Fc1ccccc1")
        assert profile.leaving_group == "F"

    def test_canonical_smiles_roundtrip(self):
        profile = self.analyzer.profile("c1ccc(Br)cc1")   # non-canonical
        assert "c" in profile.canonical_smiles              # aromatic

    def test_mw_reasonable(self):
        profile = self.analyzer.profile("Brc1ccccc1")
        assert 155 < profile.molecular_weight < 160

    def test_heavy_atoms_count(self):
        profile = self.analyzer.profile("Brc1ccccc1")
        assert profile.heavy_atoms == 7   # 6 C + 1 Br

    def test_invalid_smiles_raises(self):
        with pytest.raises(ValueError):
            self.analyzer.profile("not_a_smiles!!")


# ---------------------------------------------------------------------------
# 3. ReactantAnalyzer.find_compatible_templates()
# ---------------------------------------------------------------------------

class TestFindCompatibleTemplates:
    def setup_method(self):
        self.analyzer = ReactantAnalyzer()

    def test_suzuki_detected(self):
        templates = self.analyzer.find_compatible_templates(
            "Brc1ccccc1", "OB(O)c1cccnc1"
        )
        names = [t.template_name for t in templates]
        assert "suzuki_miyaura" in names

    def test_buchwald_detected(self):
        templates = self.analyzer.find_compatible_templates(
            "Brc1ccccc1", "Nc1ccccc1"
        )
        names = [t.template_name for t in templates]
        assert "buchwald_hartwig" in names

    def test_amide_coupling_detected(self):
        templates = self.analyzer.find_compatible_templates(
            "CC(=O)O", "Nc1ccccc1"
        )
        names = [t.template_name for t in templates]
        assert "amide_coupling" in names

    def test_cuaac_detected(self):
        templates = self.analyzer.find_compatible_templates(
            "[N-]=[N+]=NCc1ccccc1", "C#Cc1ccccc1"
        )
        names = [t.template_name for t in templates]
        assert "cuaac_triazole" in names

    def test_unimolecular_templates_excluded_when_smiles_b_given(self):
        """When two reactants are supplied, n_reactants==1 templates must not fire."""
        templates = self.analyzer.find_compatible_templates(
            "CC(=O)c1ccccc1", "Nc1ccccc1"    # ketone + amine
        )
        names = [t.template_name for t in templates]
        # NaBH4 reduction is unimolecular — should not appear when smiles_b is given
        assert "ketone_from_nabh4" not in names

    def test_unimolecular_templates_included_when_single_reactant(self):
        """When only smiles_a is given, n_reactants==1 templates should fire."""
        templates = self.analyzer.find_compatible_templates("CC(=O)c1ccccc1")
        names = [t.template_name for t in templates]
        assert "ketone_from_nabh4" in names

    def test_sorted_by_difficulty(self):
        templates = self.analyzer.find_compatible_templates(
            "CC(=O)O", "Nc1ccccc1"
        )
        difficulties = [t.difficulty for t in templates]
        assert difficulties == sorted(difficulties)

    def test_empty_result_for_inert_pair(self):
        """Two molecules with no reactive FGs should return no matches."""
        templates = self.analyzer.find_compatible_templates("CCCC", "CCCC")
        assert templates == []


# ---------------------------------------------------------------------------
# 4. ForwardReactor.generate()
# ---------------------------------------------------------------------------

class TestForwardReactorGenerate:
    def setup_method(self):
        self.reactor = ForwardReactor()

    def test_suzuki_generates_biaryl(self):
        preds = self.reactor.generate("Brc1ccccc1", "OB(O)c1cccnc1")
        assert len(preds) >= 1
        products = [p.product_smiles for p in preds]
        # biphenyl-type product should be produced
        assert any("c1ccc" in smi and "n" in smi for smi in products)

    def test_buchwald_generates_aryl_amine(self):
        preds = self.reactor.generate("Brc1ccccc1", "Nc1ccccc1")
        assert len(preds) >= 1
        # N-aryl amine: nitrogen bridge between two aryl groups
        products = [p.product_smiles for p in preds]
        assert any("Nc" in smi or "N(c" in smi for smi in products)

    def test_amide_coupling_produces_amide(self):
        preds = self.reactor.generate("CC(=O)O", "Nc1ccccc1")
        products = [p.product_smiles for p in preds]
        # should contain C(=O)N motif
        assert any("C(=O)N" in smi or "NC(=O)" in smi for smi in products)

    def test_cuaac_generates_triazole(self):
        preds = self.reactor.generate("[N-]=[N+]=NCc1ccccc1", "C#Cc1ccccc1")
        assert len(preds) >= 1

    def test_nabh4_reduction_unimolecular(self):
        """Unimolecular reduction should work when no smiles_b given."""
        preds = self.reactor.generate("CC(=O)c1ccccc1")
        names = [p.template_name for p in preds]
        assert "ketone_from_nabh4" in names

    def test_reaction_name_filter(self):
        preds = self.reactor.generate(
            "Brc1ccccc1", "OB(O)c1cccnc1", reaction_name="suzuki_miyaura"
        )
        assert all(p.template_name == "suzuki_miyaura" for p in preds)

    def test_top_k_respected(self):
        preds = self.reactor.generate("CC(=O)O", "Nc1ccccc1", top_k=1)
        assert len(preds) <= 3   # may over-generate then trim in scoring step

    def test_prediction_has_reaction_smiles(self):
        preds = self.reactor.generate("Brc1ccccc1", "OB(O)c1cccnc1")
        for p in preds:
            assert ">>" in p.reaction_smiles

    def test_prediction_taxonomy_id_set(self):
        preds = self.reactor.generate("Brc1ccccc1", "OB(O)c1cccnc1")
        suzuki_preds = [p for p in preds if p.template_name == "suzuki_miyaura"]
        assert suzuki_preds
        assert suzuki_preds[0].taxonomy_id != ""


# ---------------------------------------------------------------------------
# 5. score_products()
# ---------------------------------------------------------------------------

class TestScoreProducts:
    def setup_method(self):
        self.reactor = ForwardReactor()

    def _make_prediction(self, template_name: str = "suzuki_miyaura") -> ProductPrediction:
        return ProductPrediction(
            product_smiles="c1ccc(-c2cccnc2)cc1",
            reactant_a="Brc1ccccc1",
            reactant_b="OB(O)c1cccnc1",
            template_name=template_name,
            taxonomy_id="Suzuki_miyaura",
            reaction_smiles="Brc1ccccc1.OB(O)c1cccnc1>>c1ccc(-c2cccnc2)cc1",
            hte_families=["suzuki_miyaura"],
        )

    def test_score_products_returns_sorted(self):
        preds = self.reactor.generate("Brc1ccccc1", "OB(O)c1cccnc1")
        ranked = score_products(preds)
        scores = [p.overall_score for p in ranked]
        assert scores == sorted(scores, reverse=True)

    def test_score_sets_confidence_label(self):
        pred = self._make_prediction()
        scored = score_products([pred])
        assert scored[0].confidence_label in ("high", "medium", "low")

    def test_score_sets_hte_yield_proxy(self):
        pred = self._make_prediction()
        scored = score_products([pred])
        assert scored[0].hte_yield_proxy > 0

    def test_chemoselectivity_penalty_zero_for_single(self):
        pred = self._make_prediction()
        scored = score_products([pred])
        assert scored[0].chemoselectivity_penalty == 0.0

    def test_chemoselectivity_penalty_nonzero_for_competing(self):
        preds = [
            self._make_prediction("suzuki_miyaura"),
            self._make_prediction("stille_coupling"),
        ]
        scored = score_products(preds)
        assert all(p.chemoselectivity_penalty > 0 for p in scored)

    def test_empty_list_returns_empty(self):
        assert score_products([]) == []


# ---------------------------------------------------------------------------
# 6. Agent tool: _inspect_reactants
# ---------------------------------------------------------------------------

class TestInspectReactants:
    def test_basic_pair(self):
        r = _inspect_reactants(smiles_a="Brc1ccccc1", smiles_b="Nc1ccccc1")
        assert r["success"]
        assert "reactant_a" in r
        assert "reactant_b" in r

    def test_fgs_populated(self):
        r = _inspect_reactants(smiles_a="Brc1ccccc1", smiles_b="Nc1ccccc1")
        assert len(r["reactant_a"]["functional_groups"]) > 0

    def test_compatibility_flag_present(self):
        r = _inspect_reactants(smiles_a="Brc1ccccc1", smiles_b="Nc1ccccc1")
        assert len(r["compatibility_flags"]) > 0

    def test_missing_smiles_a_returns_error(self):
        r = _inspect_reactants()
        assert not r["success"]

    def test_unimolecular_mode(self):
        r = _inspect_reactants(smiles_a="CC(=O)c1ccccc1")
        assert r["success"]
        assert "reactant_b" not in r   # single reactant

    def test_alias_reactant_a(self):
        r = _inspect_reactants(reactant_a="Brc1ccccc1", reactant_b="Nc1ccccc1")
        assert r["success"]


# ---------------------------------------------------------------------------
# 7. Agent tool: _identify_reactions
# ---------------------------------------------------------------------------

class TestIdentifyReactions:
    def test_suzuki_identified(self):
        r = _identify_reactions(smiles_a="Brc1ccccc1", smiles_b="OB(O)c1ccccc1")
        assert r["success"]
        names = [x["name"] for x in r.get("compatible_reactions", [])]
        assert "suzuki_miyaura" in names

    def test_no_fgs_returns_zero(self):
        r = _identify_reactions(smiles_a="CCCC", smiles_b="CCCC")
        assert r["success"]
        assert r.get("total_matched", 0) == 0

    def test_missing_smiles_a_returns_error(self):
        r = _identify_reactions()
        assert not r["success"]

    def test_sorted_by_difficulty(self):
        r = _identify_reactions(smiles_a="CC(=O)O", smiles_b="Nc1ccccc1")
        diffs = [x["difficulty"] for x in r.get("compatible_reactions", [])]
        assert diffs == sorted(diffs)


# ---------------------------------------------------------------------------
# 8. Agent tool: _generate_products
# ---------------------------------------------------------------------------

class TestGenerateProducts:
    def test_suzuki_product(self):
        r = _generate_products(smiles_a="Brc1ccccc1", smiles_b="OB(O)c1cccnc1")
        assert r["success"]
        assert r.get("total_products", 0) >= 1
        products = r.get("products", [])
        assert products[0]["template_name"] == "suzuki_miyaura"

    def test_buchwald_product(self):
        r = _generate_products(smiles_a="Brc1ccccc1", smiles_b="Nc1ccccc1")
        assert r["success"]
        products = r.get("products", [])
        assert len(products) >= 1
        names = [p["template_name"] for p in products]
        assert "buchwald_hartwig" in names

    def test_product_has_required_fields(self):
        r = _generate_products(smiles_a="Brc1ccccc1", smiles_b="OB(O)c1cccnc1")
        prod = r["products"][0]
        for field in ["rank", "product_smiles", "reaction_smiles", "template_name",
                      "overall_score", "confidence_label", "hte_yield_proxy"]:
            assert field in prod, f"Missing field: {field}"

    def test_reaction_smiles_contains_arrow(self):
        r = _generate_products(smiles_a="Brc1ccccc1", smiles_b="OB(O)c1cccnc1")
        for p in r.get("products", []):
            assert ">>" in p["reaction_smiles"]
            _, err = _validate_reaction_smiles(p["reaction_smiles"], require_product=True)
            assert err is None

    def test_missing_smiles_a_returns_error(self):
        r = _generate_products()
        assert not r["success"]

    def test_top_k_limits_output(self):
        r = _generate_products(smiles_a="Brc1ccccc1", smiles_b="Nc1ccccc1", top_k=1)
        assert r.get("total_products", 0) <= 1

    def test_unimolecular_nabh4(self):
        r = _generate_products(smiles_a="CC(=O)c1ccccc1")
        assert r["success"]
        names = [p["template_name"] for p in r.get("products", [])]
        assert "ketone_from_nabh4" in names

    def test_reaction_name_filter(self):
        r = _generate_products(
            smiles_a="Brc1ccccc1", smiles_b="OB(O)c1cccnc1",
            reaction_name="suzuki_miyaura"
        )
        for p in r.get("products", []):
            assert p["template_name"] == "suzuki_miyaura"

    def test_filters_invalid_generated_reaction_smiles(self, monkeypatch):
        bad = SimpleNamespace(
            product_smiles="c1ccccc1",
            reaction_smiles="Brc1ccccc1.Nc1ccccc1>>",
            template_name="buchwald_hartwig",
            taxonomy_id="buchwald_hartwig",
            description="bad",
            difficulty=0.2,
            hte_yield_proxy=55.0,
            overall_score=90.0,
            confidence_label="high",
            new_stereocenters=0,
            all_product_smiles=["c1ccccc1"],
            competing_templates=[],
            notes="",
            hte_families=["buchwald_hartwig"],
        )
        good = SimpleNamespace(
            product_smiles="Nc1ccccc1c1ccccc1",
            reaction_smiles="Brc1ccccc1.Nc1ccccc1>>Nc1ccccc1c1ccccc1",
            template_name="buchwald_hartwig",
            taxonomy_id="buchwald_hartwig",
            description="good",
            difficulty=0.3,
            hte_yield_proxy=50.0,
            overall_score=80.0,
            confidence_label="medium",
            new_stereocenters=0,
            all_product_smiles=["Nc1ccccc1c1ccccc1"],
            competing_templates=[],
            notes="",
            hte_families=["buchwald_hartwig"],
        )

        class _FakeReactor:
            def generate(self, *args, **kwargs):  # noqa: ANN002, ANN003
                return [bad, good]

        monkeypatch.setattr("chemtools.forward.reactor.ForwardReactor", _FakeReactor)
        monkeypatch.setattr(
            "chemtools.forward.scoring.score_products",
            lambda preds, smiles_a="", smiles_b="": preds,  # noqa: ARG005
        )

        r = _generate_products(smiles_a="Brc1ccccc1", smiles_b="Nc1ccccc1", top_k=5)
        assert r["success"]
        assert r.get("invalid_reaction_smiles_filtered", 0) == 1
        assert r.get("total_products", 0) == 1
        _, err = _validate_reaction_smiles(r["products"][0]["reaction_smiles"], require_product=True)
        assert err is None


class TestRecommendForwardConditions:
    def test_requires_product_smiles(self):
        r = _recommend_forward_conditions(
            smiles_a="Brc1ccccc1",
            smiles_b="Nc1ccccc1",
            product_smiles="",
        )
        assert not r["success"]
        assert "product_smiles is required" in r["error"]


# ---------------------------------------------------------------------------
# 9. Agent tool: _rank_products
# ---------------------------------------------------------------------------

class TestRankProducts:
    def test_ranks_from_generate(self):
        gen = _generate_products(smiles_a="CC(=O)O", smiles_b="Nc1ccccc1", top_k=5)
        products = gen.get("products", [])
        r = _rank_products(products=products, top_k=2)
        assert r["success"]
        assert len(r.get("ranked_products", [])) <= 2

    def test_selectivity_note_present(self):
        gen = _generate_products(smiles_a="CC(=O)O", smiles_b="Nc1ccccc1", top_k=5)
        r = _rank_products(products=gen.get("products", []))
        assert "selectivity_recommendation" in r

    def test_empty_products_returns_success(self):
        r = _rank_products(products=[])
        assert r["success"]


# ---------------------------------------------------------------------------
# 10. Agent tool: _plan_forward_route
# ---------------------------------------------------------------------------

class TestPlanForwardRoute:
    def test_single_step_route(self):
        r = _plan_forward_route(
            scaffold_smiles="Brc1ccccc1",
            reagent_smiles_list=["OB(O)c1cccnc1"],
            reaction_sequence=["suzuki_miyaura"],
        )
        assert r["success"]
        steps = r.get("route_steps", [])
        assert len(steps) >= 1
        assert steps[0]["step_number"] == 1

    def test_final_product_set(self):
        r = _plan_forward_route(
            scaffold_smiles="Brc1ccccc1",
            reagent_smiles_list=["OB(O)c1cccnc1"],
        )
        assert r["success"]
        assert r.get("final_product") is not None

    def test_missing_scaffold_returns_error(self):
        r = _plan_forward_route(scaffold_smiles="")
        assert not r["success"]

    def test_complexity_increases(self):
        r = _plan_forward_route(
            scaffold_smiles="Brc1ccccc1",
            reagent_smiles_list=["OB(O)c1cccnc1"],
        )
        assert r.get("complexity_end", 0) >= r.get("complexity_start", 0)


# ---------------------------------------------------------------------------
# 11. Tool registration
# ---------------------------------------------------------------------------

class TestToolRegistration:
    def test_all_8_tools_registered(self):
        from chem_coworker.tools import REGISTRY
        fwd_tools = REGISTRY.categories().get("forward_synthesis", [])
        expected = {
            "inspect_reactants", "identify_reactions", "generate_products",
            "rank_products", "find_forward_precedent", "search_reactant_precedent",
            "recommend_forward_conditions", "plan_forward_route",
        }
        assert expected.issubset(set(fwd_tools)), (
            f"Missing tools: {expected - set(fwd_tools)}"
        )

    def test_forward_synthesis_tools_list_length(self):
        assert len(FORWARD_SYNTHESIS_TOOLS) == 8

    def test_all_tools_have_fn(self):
        for tool in FORWARD_SYNTHESIS_TOOLS:
            assert callable(tool.fn), f"{tool.name} missing fn"

    def test_all_tools_have_category(self):
        for tool in FORWARD_SYNTHESIS_TOOLS:
            assert tool.category == "forward_synthesis"


# ---------------------------------------------------------------------------
# 12. Workflow registration
# ---------------------------------------------------------------------------

class TestWorkflowRegistration:
    def test_forward_synthesis_workflow_registered(self):
        from chem_coworker.workflow import WORKFLOW_REGISTRY
        names = WORKFLOW_REGISTRY.names()
        assert "forward_synthesis" in names

    def test_forward_synthesis_workflow_predicate(self):
        from chem_coworker.workflow import WORKFLOW_REGISTRY
        wf = WORKFLOW_REGISTRY.get_for_task("forward_synthesis")
        assert wf.name == "forward_synthesis"

    def test_retrosynthesis_still_works(self):
        from chem_coworker.workflow import WORKFLOW_REGISTRY
        wf = WORKFLOW_REGISTRY.get_for_task("retrosynthesis")
        assert wf.name == "retrosynthesis"

    def test_fallback_for_unknown_task(self):
        from chem_coworker.workflow import WORKFLOW_REGISTRY
        wf = WORKFLOW_REGISTRY.get_for_task("some_other_task")
        assert wf.name == "forward_chemistry"

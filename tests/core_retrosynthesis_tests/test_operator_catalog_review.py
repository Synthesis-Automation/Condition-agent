"""Graphical generic-operator catalog regressions."""

from __future__ import annotations

from pathlib import Path

from core_retrosynthesis.generic_models import (
    GenericCoreTemplate,
    GenericGraphOperator,
    GenericTemplateLibrary,
    GenericTemplatePrecedent,
)
from core_retrosynthesis.models import TemplateContext
from core_retrosynthesis.operator_catalog_review import (
    render_generic_operator_catalog_html,
    write_generic_operator_catalog_html,
)


def _operator(operator_id: str) -> GenericGraphOperator:
    return GenericGraphOperator(
        operator_id=operator_id,
        operator_signature=f"signature:{operator_id}",
        edit_tokens=("formed:C-N:NONE>SINGLE",),
        realization_ids=(f"realization:{operator_id}",),
        abstraction_levels=("L1", "L2"),
        observation_support=3,
        independent_reference_support=2,
    )


def _template(operator_id: str) -> GenericCoreTemplate:
    precedent = GenericTemplatePrecedent(
        reaction_id=f"reaction:{operator_id}",
        reference_id=f"reference:{operator_id}",
        product_smiles="CN",
        precursor_smiles="C.N",
        mapped_reaction_smiles="[CH4:1].[NH3:2]>>[CH3:1][NH2:2]",
        context=TemplateContext((), (), ()),
    )
    return GenericCoreTemplate(
        template_id=f"template:{operator_id}",
        operator_id=operator_id,
        transformation_kind="bond_formation",
        abstraction_level="L2",
        compiler_engine="rdchiral",
        reaction_smarts="[C:1].[N:2]>>[C:1]-[N:2]",
        product_smarts="[C]-[N]",
        precursor_smarts="[C].[N]",
        edit_tokens=("formed:C-N:NONE>SINGLE",),
        handle_signature="handle",
        stereo_policy="exact",
        observation_support=3,
        independent_reference_support=2,
        precedents=(precedent,),
        realization_id=f"realization:{operator_id}",
    )


def test_operator_catalog_renders_first_stable_id_with_graphic() -> None:
    library = GenericTemplateLibrary(
        templates=(_template("OP1:b"), _template("OP1:a")),
        source_row_count=2,
        accepted_observation_count=2,
        rejection_counts={},
        definition={},
        operators=(_operator("OP1:b"), _operator("OP1:a")),
    )

    page = render_generic_operator_catalog_html(
        library,
        limit=1,
        library_source="operators.json.gz",
    )

    assert page.count('<article class="operator-card"') == 1
    assert "OP1:a" in page
    assert "OP1:b" not in page
    assert "formed:C-N:NONE&gt;SINGLE" in page
    assert "<svg" in page
    assert "shortest_observed_precedent_per_operator" in page


def test_operator_catalog_writer_reports_counts(tmp_path: Path) -> None:
    library = GenericTemplateLibrary(
        templates=(_template("OP1:a"),),
        source_row_count=1,
        accepted_observation_count=1,
        rejection_counts={},
        definition={},
        operators=(_operator("OP1:a"),),
    )
    output = tmp_path / "catalog.html"

    summary = write_generic_operator_catalog_html(
        library,
        output,
        limit=100,
    )

    assert output.is_file()
    assert summary["rendered_operator_count"] == 1
    assert summary["library_operator_count"] == 1
    assert summary["html_bytes"] == output.stat().st_size

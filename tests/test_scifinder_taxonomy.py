import importlib.util
import json
import sys
from pathlib import Path

import pytest


pytest.importorskip("PyQt6")


DATA_PROCESSOR_DIR = Path(__file__).resolve().parents[1] / "data-processor"
if str(DATA_PROCESSOR_DIR) not in sys.path:
    sys.path.insert(0, str(DATA_PROCESSOR_DIR))


MODULE_PATH = DATA_PROCESSOR_DIR / "Scifinder_rdf_processer.py"
spec = importlib.util.spec_from_file_location("scifinder_rdf_processer", MODULE_PATH)
scifinder = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(scifinder)


def _make_taxonomy() -> scifinder._TaxonomyIndex:  # type: ignore[attr-defined]
    tax_dir = Path(__file__).resolve().parents[1] / "data" / "compound_taxonomy"
    return scifinder._TaxonomyIndex(str(tax_dir))


def test_taxonomy_roles_cover_new_reagent_classes():
    taxonomy = _make_taxonomy()

    assert taxonomy.role_for("Sulfuric acid", "7664-93-9") == "ACID"
    assert taxonomy.role_for("Tetrabutylammonium bromide", "1643-19-2") == "ADDITIVE"
    assert taxonomy.role_for("Oxygen (g)", "7782-44-7") == "OXIDANT"
    assert taxonomy.role_for("Zinc dust", "7440-66-6") == "REDUCTANT"


def test_condition_core_falls_back_to_acid_token_when_no_catalyst():
    taxonomy = _make_taxonomy()
    generator = scifinder.ReactionMarkdownGenerator(taxonomy=taxonomy)

    row = {
        "CatalystCoreGeneric": json.dumps([]),
        "CatalystCoreDetail": json.dumps([]),
        "Ligand": json.dumps([]),
        "Reagent": json.dumps(["Sulfuric acid|7664-93-9"]),
        "ReagentRole": json.dumps(["ACID"]),
    }

    assert generator._compute_condition_core(row) == "Acid: H2SO4"

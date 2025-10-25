from chemtools.taxonomy import load_registry, reset_registry
from chemtools.taxonomy.registry import TaxonomyRegistry
from scripts.taxonomy.convert_legacy_taxonomies import main as convert_legacy_main


def test_load_taxonomy_registry_placeholder():
    reset_registry()
    registry = load_registry()
    assert isinstance(registry, TaxonomyRegistry)
    assert registry.manifest.schema_version
    assert registry.manifest.taxonomy_version


def test_convert_legacy_taxonomies_dry_run():
    assert convert_legacy_main(["--dry-run"]) == 0

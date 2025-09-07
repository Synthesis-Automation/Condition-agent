import os
import importlib


def setup_module(module):
    # Point resolver to tiny fixture dataset and reload module to rebuild index
    os.environ["CHEMTOOLS_REGISTRY_PATH"] = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "data", "registry_tiny.jsonl")
    )
    global registry
    import chemtools.registry as registry  # type: ignore
    importlib.reload(registry)  # ensure index uses fixture path
    module.registry = registry


def test_resolve_by_cas_and_token():
    # By CAS
    r1 = registry.resolve("1310-58-3")
    assert r1.get("uid") == "1310-58-3"
    assert r1.get("role") == "BASE"
    assert "Potassium hydroxide" in r1.get("aliases", [])
    # By token
    r2 = registry.resolve("KOH")
    assert r2.get("uid") == "1310-58-3"
    assert r2.get("role") == "BASE"


def test_resolve_canonicalize_digits_only_cas():
    r = registry.resolve("1310583")
    assert r.get("uid") == "1310-58-3"
    assert r.get("role") == "BASE"


def test_resolve_catalyst_alias_and_uid():
    r1 = registry.resolve("CuI")
    assert r1.get("uid") == "7681-65-4"
    assert r1.get("role") == "CATALYST"
    r2 = registry.resolve("7681-65-4")
    assert r2.get("uid") == "7681-65-4"
    assert r2.get("role") == "CATALYST"


def test_resolve_ligand_punctuation_insensitive():
    # Name '1,10-Phenanthroline' should be found via punctuation-insensitive match
    r = registry.resolve("1,10-phenanthroline")
    assert r.get("uid") == "72-52-8"
    assert r.get("role") == "LIGAND"


def test_not_found_returns_error():
    r = registry.resolve("this_is_not_in_registry")
    assert r.get("error") == "NOT_FOUND"


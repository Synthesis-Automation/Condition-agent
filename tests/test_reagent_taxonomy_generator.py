import json
import shutil
import subprocess
import sys
from pathlib import Path
import importlib.util

import pytest

MODULE_PATH = Path(__file__).resolve().parents[1] / "data-processor" / "reagent_taxonomy_generator.py"
spec = importlib.util.spec_from_file_location("reagent_taxonomy_generator", MODULE_PATH)
rrg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rrg)

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "compound_taxonomy"


def _copy_taxonomy(tmp_path: Path) -> Path:
    dest = tmp_path / "taxonomy"
    shutil.copytree(DATA_DIR, dest)
    return dest


def _make_test_cas(digits: str) -> str:
    digits = ''.join(ch for ch in digits if ch.isdigit())
    if len(digits) < 3:
        raise ValueError("need at least 3 digits for CAS generation")
    total = sum(int(d) * (idx + 1) for idx, d in enumerate(reversed(digits)))
    check = total % 10
    prefix = str(int(digits[:-2]))
    mid = digits[-2:]
    return f"{prefix}-{mid}-{check}"


@pytest.mark.parametrize(
    "name,synonyms,expected_family",
    [
        ("XPhos", ["XPhos"], "dialkylbiaryl_phosphines"),
        ("Potassium tert-butoxide", ["KOtBu"], "alkoxides_hindered"),
        ("Sodium borohydride", ["NaBH4"], "complex_hydrides"),
    ],
)
def test_infer_family_from_tokens(tmp_path, name, synonyms, expected_family):
    taxonomy_dir = _copy_taxonomy(Path(tmp_path))
    store = rrg.TaxonomyStore(taxonomy_dir)
    heuristics = rrg.RoleHeuristics(store)
    inferred = heuristics.infer_family(name, synonyms)
    assert inferred is not None
    role, family, matches = inferred
    assert role in rrg.ROLE_FILES
    assert family == expected_family
    assert matches, "expected token overlap to justify inference"


def test_cli_adds_entry(tmp_path):
    taxonomy_dir = _copy_taxonomy(Path(tmp_path))
    cas = _make_test_cas("654321")
    result = subprocess.run(
        [
            "python",
            str(MODULE_PATH),
            "--cas",
            cas,
            "--name",
            "XPhos",
            "--abbr",
            "XPhos",
            "--taxonomy-dir",
            str(taxonomy_dir),
        ],
        check=True,
        capture_output=True,
        text=True,
        cwd=str(Path(__file__).resolve().parents[1]),
    )
    payload = json.loads(result.stdout)
    assert payload["status"] == "written"
    assert payload["family_id"] == "dialkylbiaryl_phosphines"
    target_file = Path(payload["written_to"])
    data = json.loads(target_file.read_text(encoding="utf-8"))
    family_map = {fam["family_id"]: fam for fam in data.get("families", [])}
    family = family_map[payload["family_id"]]
    members = family.get("example_members", [])
    assert any(mem.get("cas") == cas for mem in members)
    assert any(mem.get("name") == "XPhos" for mem in members)






def test_auto_resolve_populates_name(monkeypatch, tmp_path, capsys):
    taxonomy_dir = _copy_taxonomy(Path(tmp_path))
    cas = _make_test_cas("765432")

    def fake_resolver(_cas: str, *, timeout: float = 0.0, session=None):
        assert _cas == cas
        return {
            "name": "StubLigand",
            "synonyms": ["StubLigand", "SLig"],
            "smiles": "C=C",
            "source": "stub",
        }

    monkeypatch.setattr(rrg, "resolve_identity_from_cas", fake_resolver)
    argv = [
        "reagent_taxonomy_generator.py",
        "--cas",
        cas,
        "--taxonomy-dir",
        str(taxonomy_dir),
        "--family",
        "dialkylbiaryl_phosphines",
        "--dry-run",
    ]
    monkeypatch.setattr(sys, "argv", argv)
    rrg.main()
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert payload["status"] == "dry_run"
    assert payload["name"] == "StubLigand"
    assert payload.get("auto_resolve_source") == "stub"
    assert payload["family_id"] == "dialkylbiaryl_phosphines"
    entry_preview = payload.get("entry_preview")
    assert entry_preview
    assert entry_preview["synonyms"][0] == "StubLigand"


def test_default_family_for_reductant(tmp_path):
    taxonomy_dir = _copy_taxonomy(Path(tmp_path))
    store = rrg.TaxonomyStore(taxonomy_dir)
    heuristics = rrg.RoleHeuristics(store)
    assert heuristics.default_family_for_role("reductant") == "metal_powders"

def test_list_families_includes_reductant(tmp_path):
    taxonomy_dir = _copy_taxonomy(Path(tmp_path))
    store = rrg.TaxonomyStore(taxonomy_dir)
    families = [fid for role, fid, _ in store.list_families("reductant")]
    assert "metal_powders" in families
    assert len(families) >= 1



def test_resolver_prefers_pubchem_cas_match():
    cas = "2592-95-2"

    class FakeResponse:
        def __init__(self, payload):
            self.status_code = 200
            self._payload = payload
            self.text = ""

        def json(self):
            return self._payload

    property_payload = {
        "PropertyTable": {
            "Properties": [
                {
                    "CID": 807,
                    "Title": "Iodine",
                    "CanonicalSMILES": "II",
                },
                {
                    "CID": 75771,
                    "Title": "1-Hydroxybenzotriazole",
                    "IsomericSMILES": "C1=CC=C2C(=C1)N=NN2O",
                },
            ]
        }
    }
    synonyms_payload = {
        "InformationList": {
            "Information": [
                {"CID": 807, "Synonym": ["iodine", "7553-56-2"]},
                {"CID": 75771, "Synonym": ["1-Hydroxybenzotriazole", cas, "HOBt"]},
            ]
        }
    }

    class FakeSession:
        def __init__(self):
            self.headers = {}

        def get(self, url, timeout):
            if url.endswith('/property/Title,IUPACName,IsomericSMILES,CanonicalSMILES/JSON'):
                return FakeResponse(property_payload)
            if url.endswith('/synonyms/JSON'):
                return FakeResponse(synonyms_payload)
            raise AssertionError(f"Unexpected URL requested: {url}")

    identity = rrg.resolve_identity_from_cas(cas, timeout=1.0, session=FakeSession())
    assert identity is not None
    assert identity.get("source") == "pubchem"
    assert identity["name"] == "1-Hydroxybenzotriazole"
    assert cas in identity["synonyms"]
    assert "iodine" not in [syn.lower() for syn in identity["synonyms"]]
    assert identity["smiles"] == "C1=CC=C2C(=C1)N=NN2O"



def test_resolver_falls_back_when_pubchem_ambiguous(monkeypatch):
    cas = "2592-95-2"

    class FakeResponse:
        def __init__(self, payload):
            self.status_code = 200
            self._payload = payload
            self.text = ""

        def json(self):
            return self._payload

    property_payload = {
        "PropertyTable": {
            "Properties": [
                {"CID": 807, "Title": "Iodine", "CanonicalSMILES": "II"},
                {"CID": 75771, "Title": "1-Hydroxybenzotriazole", "IsomericSMILES": "C1=CC=C2C(=C1)N=NN2O"},
            ]
        }
    }

    class FakeSession:
        def __init__(self):
            self.headers = {}

        def get(self, url, timeout):
            if url.endswith('/property/Title,IUPACName,IsomericSMILES,CanonicalSMILES/JSON'):
                return FakeResponse(property_payload)
            if url.endswith('/synonyms/JSON'):
                raise RuntimeError('timeout')
            raise AssertionError(f"Unexpected URL requested: {url}")

    def fake_cactus(_session, cas_value, timeout):
        assert cas_value == cas
        return {"name": "CactusName", "synonyms": ["CactusName"], "smiles": "C"}

    monkeypatch.setattr(rrg, '_resolve_via_cactus', fake_cactus)

    identity = rrg.resolve_identity_from_cas(cas, timeout=1.0, session=FakeSession())
    assert identity is not None
    assert identity.get("source") == "cactus"
    assert identity["name"] == "CactusName"
    assert identity["synonyms"] == ["CactusName"]




def test_additive_family_inference_for_hobt(tmp_path):
    taxonomy_dir = _copy_taxonomy(Path(tmp_path))
    store = rrg.TaxonomyStore(taxonomy_dir)
    heuristics = rrg.RoleHeuristics(store)
    name = "1-Hydroxybenzotriazole"
    synonyms = ["HOBt", "1-Hydroxybenzotriazole"]
    inferred = heuristics.infer_family(name, synonyms)
    assert inferred is not None
    role, family_id, matches = inferred
    assert role == "additive"
    assert family_id == "amide_coupling_additives"
    assert any(token in matches for token in {"hobt", "1hydroxybenzotriazole", "hydroxybenzotriazole"})



def test_additive_default_requires_overlap(tmp_path):
    taxonomy_dir = _copy_taxonomy(Path(tmp_path))
    cas = _make_test_cas("987654")
    cmd = [
        "python",
        str(MODULE_PATH),
        "--cas",
        cas,
        "--name",
        "MysteryXYZ1234",
        "--role",
        "additive",
        "--taxonomy-dir",
        str(taxonomy_dir),
        "--allow-default-family",
        "--dry-run",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode != 0
    message = (result.stderr or result.stdout).lower()
    assert "default family" in message
    assert "rejected" in message


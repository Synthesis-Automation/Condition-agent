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

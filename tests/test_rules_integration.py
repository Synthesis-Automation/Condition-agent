from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import pytest

from clients.mcp_rules_client import McpRulesClient

ROOT = Path(__file__).resolve().parents[1]
SERVER = [str(sys.executable), str(ROOT / 'mcp_rules_server' / 'server.py')]
RULES_SRC = ROOT / "data" / "rules" / "buchwald_cn.json"


@pytest.fixture()
def temp_rules_path(tmp_path: Path) -> Path:
    copy_path = tmp_path / "rules.json"
    shutil.copyfile(RULES_SRC, copy_path)
    return copy_path


@pytest.fixture()
def client(temp_rules_path: Path) -> McpRulesClient:
    cli = McpRulesClient(SERVER, temp_rules_path, startup_timeout_s=5)
    cli.start()
    yield cli
    cli.stop()


def test_apply_merge_audit_roundtrip(client: McpRulesClient, temp_rules_path: Path, tmp_path: Path) -> None:
    features = {
        "electrophile": {"class": "aryl chloride"},
        "nucleophile": {"class": "primary aniline"},
        "mode": "batch",
    }
    result = client.apply("Clc1ccc(Cl)cc1.Nc1ccccc1>>", features)
    assert result["suggestions"], "Expected at least one suggestion"

    candidate = {
        "playbooks": [
            {
                "id": "BH-CAND-INTEG-1",
                "name": "Integration candidate",
                "status": "candidate",
                "when": {"electrophile": {"class": "aryl chloride"}},
                "recipe": {"ligands": ["XPhos"], "base": ["K3PO4"], "solvent": ["toluene"]},
            }
        ]
    }
    merge_res = client.merge(candidate)
    assert merge_res["message"] == "merge_ok"

    audit_res = client.audit("BH-CAND-INTEG-1")
    assert audit_res["rule"]["id"] == "BH-CAND-INTEG-1"

    written = json.loads(temp_rules_path.read_text(encoding="utf-8"))
    playbook_ids = {pb.get("id") for pb in written.get("playbooks", [])}
    assert "BH-CAND-INTEG-1" in playbook_ids

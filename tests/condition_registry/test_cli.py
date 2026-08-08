import json

from condition_registry.cli import main


def test_resolve_by_cas_and_alias(capsys) -> None:
    assert main(["resolve", "--cas", "14221-01-3", "--format", "json"]) == 0
    cas_result = json.loads(capsys.readouterr().out)
    assert cas_result["status"] == "resolved"
    assert cas_result["substance"]["canonical_name"] == "Tetrakis(triphenylphosphine)palladium"

    assert main(["resolve", "--name", "HCl"]) == 0
    assert "substance: Hydrochloric Acid" in capsys.readouterr().out

    assert main([
        "resolve",
        "--identifier",
        "Dipotassium carbonate",
        "--identifier-type",
        "systematic_name",
        "--format",
        "json",
    ]) == 0
    identifier_result = json.loads(capsys.readouterr().out)
    assert identifier_result["matched_identifier"]["identifier_type"] == (
        "systematic_name"
    )
    assert identifier_result["substance"]["substance_id"] == "cas:584-08-7"


def test_resolve_returns_failure_for_invalid_or_unknown_query(capsys) -> None:
    assert main(["resolve", "--cas", "7732-18-4"]) == 1
    assert "status: invalid_identifier" in capsys.readouterr().out

    assert main(["resolve", "--name", "not a registered substance"]) == 1
    assert "status: unresolved" in capsys.readouterr().out


def test_self_test_text_and_json(capsys) -> None:
    assert main(["self-test"]) == 0
    assert "overall: PASS" in capsys.readouterr().out

    assert main(["self-test", "--format", "json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["overall"] == "PASS"
    assert payload["failed"] == 0


def test_validate_reports_current_registry_state(capsys) -> None:
    assert main(["validate"]) == 0
    assert "duplicate CAS values: 0" in capsys.readouterr().out

    exit_code = main(["validate", "--format", "json"])
    report = json.loads(capsys.readouterr().out)
    assert exit_code == (1 if report["issue_rows"] else 0)
    assert report["accepted_rows"] + report["issue_rows"] == report["total_rows"]


def test_audit_writes_outputs(tmp_path, capsys) -> None:
    assert main(["audit", str(tmp_path)]) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["total_rows"] > 0
    assert {path.name for path in tmp_path.iterdir()} == {
        "accepted.jsonl",
        "issues.csv",
        "migration_report.json",
        "quarantine.jsonl",
    }

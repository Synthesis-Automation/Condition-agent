from pathlib import Path

import pytest

from chem_coworker.extractor import NotesExtractor


def _make_extractor() -> NotesExtractor:
    ex = object.__new__(NotesExtractor)
    ex.llm = object()
    ex.model_name = "fake"
    ex.verbose = False
    return ex


def test_plan_reaction_type_filing_canonicalizes_and_routes_unknown_to_general() -> None:
    ex = _make_extractor()
    plan = ex._plan_reaction_type_filing(
        reaction_type_hint="",
        detected_types=["not_a_real_reaction_label"],
        mismatch_policy="warn",
        unknown_reaction_policy="general",
        confirm_callback=None,
    )
    assert plan["reaction_types"] == ["general"]
    assert plan["reaction_types_unknown"] == ["not_a_real_reaction_label"]
    assert any("Unknown reaction label" in w for w in plan["warnings"])
    assert plan["warning_details"]
    assert plan["warning_details"][0]["severity"] == "warn"
    assert plan["warning_details"][0]["code"] == "unknown_reaction_type"
    assert "not_a_real_reaction_label" in plan["reaction_type_suggestions"]


def test_plan_reaction_type_filing_rejects_hint_detected_mismatch() -> None:
    ex = _make_extractor()
    plan = ex._plan_reaction_type_filing(
        reaction_type_hint="suzuki",
        detected_types=["heck"],
        mismatch_policy="reject",
        unknown_reaction_policy="general",
        confirm_callback=None,
    )
    assert "error" in plan
    assert "hint=suzuki_miyaura" in plan["error"]
    assert plan["mismatch"] is True


def test_plan_reaction_type_filing_force_uses_hint_only() -> None:
    ex = _make_extractor()
    plan = ex._plan_reaction_type_filing(
        reaction_type_hint="suzuki",
        detected_types=["heck"],
        mismatch_policy="force",
        unknown_reaction_policy="general",
        confirm_callback=None,
    )
    assert plan["reaction_types"] == ["suzuki_miyaura"]
    assert plan["mismatch"] is True


def test_intake_dry_run_does_not_write_files(monkeypatch, tmp_path: Path) -> None:
    ex = _make_extractor()

    monkeypatch.setattr(
        ex,
        "_load_source",
        lambda source, save: ("doc text", "demo source", [], ""),
    )
    monkeypatch.setattr(
        ex,
        "_extract_with_llm",
        lambda *args, **kwargs: ("## Source: demo\n\nnotes", ["suzuki"]),
    )

    append_calls = []

    def _fake_append_notes(*args, **kwargs):
        append_calls.append((args, kwargs))
        return tmp_path / "should_not_exist.md"

    def _fake_get_notes_path(reaction_type: str, note_type: str = "reactions"):
        return tmp_path / note_type / f"{reaction_type}.md"

    monkeypatch.setattr("chem_coworker.tools.notes.append_notes", _fake_append_notes)
    monkeypatch.setattr("chem_coworker.tools.notes.get_notes_path", _fake_get_notes_path)

    result = ex.intake("raw text", dry_run=True)
    assert result["success"] is True
    assert result["dry_run"] is True
    assert result["write_performed"] is False
    assert append_calls == []
    assert result["reaction_types"] == ["suzuki_miyaura"]
    assert str(tmp_path / "reactions" / "suzuki_miyaura.md") in result["notes_files"]


def test_intake_rejects_unknown_reaction_label(monkeypatch) -> None:
    ex = _make_extractor()
    monkeypatch.setattr(
        ex,
        "_load_source",
        lambda source, save: ("doc text", "demo source", [], ""),
    )
    monkeypatch.setattr(
        ex,
        "_extract_with_llm",
        lambda *args, **kwargs: ("## Source: demo\n\nnotes", ["definitely_not_real"]),
    )

    result = ex.intake("raw text", dry_run=True, unknown_reaction_policy="reject")
    assert result["success"] is False
    assert "Unknown reaction label" in result["error"]


def test_intake_quarantine_policy_routes_to_quarantine_store(monkeypatch, tmp_path: Path) -> None:
    ex = _make_extractor()
    monkeypatch.setattr(
        ex,
        "_load_source",
        lambda source, save: ("doc text", "demo source", [], ""),
    )
    monkeypatch.setattr(
        ex,
        "_extract_with_llm",
        lambda *args, **kwargs: ("## Source: demo\n\nnotes", ["definitely_not_real"]),
    )

    def _fake_quarantine_path(note_type: str = "reactions", bucket: str = "unknown_reaction_labels"):
        return tmp_path / "_quarantine" / note_type / f"{bucket}.md"

    monkeypatch.setattr("chem_coworker.tools.notes.get_quarantine_notes_path", _fake_quarantine_path)

    result = ex.intake("raw text", dry_run=True, unknown_reaction_policy="quarantine")
    assert result["success"] is True
    assert result["quarantined"] is True
    assert result["reaction_types"] == []
    assert "_quarantine" in result["quarantine_file"]
    assert any(d["code"] == "quarantine_routing" for d in result["warning_details"])

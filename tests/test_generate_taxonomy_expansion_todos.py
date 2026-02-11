from __future__ import annotations

from pathlib import Path

from scripts import generate_taxonomy_expansion_todos as todo_mod


def test_build_todo_payload_generates_tracks_and_stubs() -> None:
    discovery_payload = {
        "clusters": [
            {
                "cluster_id": "Ar-X|R-NH2 -> Ar-NR2 || events:c_n_bond_formation",
                "count": 12,
                "motif_signature": "Ar-X|R-NH2 -> Ar-NR2",
                "event_signature": "c_n_bond_formation",
                "reasons": [["unknown_reaction_type", 9], ["motif_outside_reaction_taxonomy", 9]],
                "recommended_action": "expand_existing_slots_or_aliases",
                "top_taxonomy_candidates": [
                    {"reaction_id": "c_n_cross_coupling", "score": 0.86, "reactant_overlap": 0.5, "product_overlap": 1.0}
                ],
                "samples": [{"reaction_smiles": "A>>B"}],
            },
            {
                "cluster_id": "Foo-Bar -> Baz-Qux || events:none",
                "count": 3,
                "motif_signature": "Foo-Bar -> Baz-Qux",
                "event_signature": "none",
                "reasons": [["unknown_reaction_type", 3], ["unclassified_motif", 1]],
                "recommended_action": "consider_new_reaction_family_or_multistep_bucket",
                "top_taxonomy_candidates": [
                    {"reaction_id": "c_n_cross_coupling", "score": 0.21, "reactant_overlap": 0.0, "product_overlap": 0.0}
                ],
                "samples": [{"reaction_smiles": "C>>D"}],
            },
        ]
    }
    taxonomy_context = {
        "reacted_slot_motifs": {"Ar-Br", "R-NH2"},
        "formed_slot_motifs": {"Ar-NR2"},
        "compound_ids": {"Ar-X", "R-NH2", "Ar-NR2"},
        "group_ids": {"Ar", "R", "-X", "-NH2", "-NR2"},
    }

    payload = todo_mod.build_todo_payload(
        discovery_payload,
        taxonomy_context=taxonomy_context,
        top_clusters=10,
        top_motifs_per_cluster=5,
    )

    assert payload["summary"]["todo_entries"] == 2
    todos = payload["todos"]
    assert todos

    expand_todo = next(item for item in todos if item["cluster_id"].startswith("Ar-X|R-NH2"))
    assert expand_todo["action_track"] == "expand_existing_reaction_type"
    patch = expand_todo["patch_stubs"]["reaction_type_update"]
    assert patch is not None
    assert patch["reaction_id"] == "c_n_cross_coupling"
    assert "Ar-X" in (patch["reactant_slot_additions_todo"]["electrophile_or_substrate"] or [])
    assert expand_todo["patch_stubs"]["new_reaction_family"] is None

    new_family_todo = next(item for item in todos if item["cluster_id"].startswith("Foo-Bar"))
    assert new_family_todo["action_track"] == "propose_new_reaction_family"
    assert new_family_todo["patch_stubs"]["new_reaction_family"] is not None
    compounds = new_family_todo["patch_stubs"]["new_compound_entries"]
    assert compounds and isinstance(compounds, list)


def test_write_markdown_renders_sections(tmp_path: Path) -> None:
    payload = {
        "summary": {"input_clusters": 1, "todo_entries": 1, "tracks": {}},
        "todos": [
            {
                "priority_rank": 1,
                "cluster_id": "cluster-1",
                "cluster_count": 5,
                "priority_score": 9,
                "action_track": "expand_existing_reaction_type",
                "motif_signature": "Ar-X -> Ar-NR2",
                "event_signature": "c_n_bond_formation",
                "reason_counts": {"unknown_reaction_type": 5},
                "top_candidate": {"reaction_id": "c_n_cross_coupling", "score": 0.8},
                "gap_analysis": {
                    "missing_reacted_slot_motifs": ["Ar-X"],
                    "missing_formed_slot_motifs": [],
                    "existing_compounds_missing_slots": ["Ar-X"],
                    "unknown_compound_ids": [],
                },
                "patch_stubs": {
                    "reaction_type_update": {"reaction_id": "c_n_cross_coupling"},
                    "new_reaction_family": None,
                    "new_compound_entries": [],
                },
            }
        ],
    }
    out = tmp_path / "todos.md"
    todo_mod.write_markdown(out, payload)
    text = out.read_text(encoding="utf-8")
    assert "# Taxonomy Expansion TODOs" in text
    assert "Priority 1" in text
    assert "Patch Stubs" in text

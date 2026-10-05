"""Released route qualification preserves uncertainty and ignores array order."""

from copy import deepcopy
import gzip
import json
import sqlite3
import zlib

import pytest

from core_retrosynthesis.processed_routes import normalize_released_route, build_processed_route_catalog, ProcessedRouteCatalog
from core_retrosynthesis.route_conversion import build_observed_route_tree, iter_route_trees
from core_retrosynthesis.route_curation import RouteQualityError


def _wrapper():
    reactions = [
        ("root", "[CH3:1][OH:2].[CH2:3]([CH3:4])[Br:5]>>[CH3:1][O:2][CH2:3][CH3:4]"),
        ("methanol", "[CH3:1][Br:2].[OH2:3]>>[CH3:1][OH:3]"),
        ("ethyl_bromide", "[CH2:1]([CH3:2])[OH:3].[BrH:4]>>[CH2:1]([CH3:2])[Br:4]"),
    ]
    return {"raw_route": {"route_id": "US_TEST_0", "original_tree": {
        "num_reactions": 3, "reaction_ids": [key for key, _ in reactions]},
        "subtrees": [{"subtree_id": "sub", "reactions": [
            {"_id": f"sub_{key}", "reaction_smiles": smiles,
             "abstracted_reaction_smiles": "invalid algorithmic label"}
            for key, smiles in reversed(reactions)]}]}}


def test_branched_route_reconstructs_topology_and_archives_weak_labels():
    wrapper = _wrapper()
    tree = build_observed_route_tree(normalize_released_route(wrapper))
    assert tree.reaction_count == 3 and tree.maximum_depth == 2
    assert tree.root.reaction.evidence.source_reaction_id == "root"
    assert tree.root.reaction.evidence.abstraction_status == "archived_algorithmic_annotation"
    assert wrapper["raw_route"]["subtrees"][0]["reactions"][0]["abstracted_reaction_smiles"] == "invalid algorithmic label"


def test_duplicate_subtree_memberships_deduplicate_only_agreeing_steps():
    wrapper = _wrapper()
    wrapper["raw_route"]["subtrees"].append(deepcopy(wrapper["raw_route"]["subtrees"][0]))
    assert len(normalize_released_route(wrapper)["steps"]) == 3
    wrapper["raw_route"]["subtrees"][1]["reactions"][0]["reaction_smiles"] = "[CH3:1][OH:2]>>[CH2:1]=[O:2]"
    with pytest.raises(RouteQualityError, match="conflicting_original_reaction"):
        normalize_released_route(wrapper)


def test_invalid_supplied_maps_do_not_become_validated_trees():
    wrapper = _wrapper()
    wrapper["raw_route"]["subtrees"][0]["reactions"][0]["reaction_smiles"] = "CO>>C=O"
    with pytest.raises(RouteQualityError, match="mapping"):
        normalize_released_route(wrapper)


def _write(path, rows):
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row) + "\n")


def test_catalog_keeps_rejected_sources_and_exact_step_joins(tmp_path):
    wrapper = _wrapper()
    invalid = deepcopy(wrapper)
    invalid["raw_route"]["route_id"] = "US_OTHER_0"
    invalid["raw_route"]["subtrees"][0]["reactions"][0]["reaction_smiles"] = "CO>>C=O"
    routes, steps, output = (tmp_path / name for name in ("routes.jsonl.gz", "steps.jsonl.gz", "routes.sqlite"))
    _write(routes, [wrapper, invalid])
    _write(steps, [{"observation_id": "obs-" + key, "source": {"source_record_id": key}}
                   for key in wrapper["raw_route"]["original_tree"]["reaction_ids"]])
    report = build_processed_route_catalog(routes, steps, output)
    assert report["route_count"] == 2
    assert report["counts"] == {"unresolved": 1, "validated_tree": 1}
    with sqlite3.connect(output) as db:
        raw = db.execute("SELECT source FROM routes WHERE id='US_OTHER_0'").fetchone()[0]
        assert json.loads(zlib.decompress(raw)) == invalid
        assert db.execute("SELECT count(*) FROM route_steps").fetchone()[0] == 6
    assert build_processed_route_catalog(routes, steps, output) == report
    catalog = ProcessedRouteCatalog(output)
    summary = catalog.route("US_TEST_0")
    assert all(m["join_status"] == "resolved" for m in summary["memberships"])
    assert "source" not in summary and "tree" not in summary
    assert catalog.route("US_OTHER_0", include_source=True)["source"] == invalid
    assert len(list(iter_route_trees(output))) == 1


def test_missing_step_observations_remain_review_evidence(tmp_path):
    routes, steps, output = (tmp_path / name for name in ("routes.jsonl.gz", "steps.jsonl.gz", "routes.sqlite"))
    _write(routes, [_wrapper()])
    _write(steps, [])
    report = build_processed_route_catalog(routes, steps, output)
    assert report["missing_memberships"] == 3
    row = ProcessedRouteCatalog(output).route("US_TEST_0", include_tree=True)
    assert row["status"] == "unresolved" and row["tree_available"]
    assert not list(iter_route_trees(output))

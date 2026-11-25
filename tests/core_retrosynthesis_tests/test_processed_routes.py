"""Released route qualification preserves uncertainty and ignores array order."""

from copy import deepcopy
import gzip
import hashlib
import json
import sqlite3
import zlib

import pytest

from core_retrosynthesis.processed_routes import (
    build_processed_route_catalog, build_route_bundle_catalog, normalize_released_route, ProcessedRouteCatalog,
)
from core_retrosynthesis.route_conversion import build_observed_route_tree, iter_route_trees
from core_retrosynthesis.route_curation import RouteQualityError
from condition_recommender.corpus_io import file_sha256


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


def _literature_wrapper():
    route_id = "cbyd:vTest:route"
    rows = [{"reaction_id": "duplicated author id", "reaction_smiles_mapped": smiles} for smiles in (
        "[CH3:1][OH:2]>>[CH2:1]=[O:2]", "[CH3:1][Br:2].[OH2:3]>>[CH3:1][OH:3]")]
    rows[0]["reaction_smiles"] = "CO>>C=O"
    rows[1]["reaction_smiles"] = "CBr.O>>CO"
    nodes = [{"id": name, "t": "m", "s": smiles, "r": role} for name, smiles, role in (
        ("target", "C=O", "target"), ("methanol", "CO", "intermediate"),
        ("start1", "CBr", "start"), ("start2", "O", "start"))]
    graph = {"rx": [{"rid": r["reaction_id"], "m": r["reaction_smiles_mapped"]} for r in rows],
             "nodes": nodes + [{"id": "rx1", "t": "r", "i": 0}, {"id": "rx2", "t": "r", "i": 1}],
             "edges": [["rx1", "target"], ["methanol", "rx1"], ["rx2", "methanol"],
                       ["start1", "rx2"], ["start2", "rx2"]], "warn": []}
    return {"schema_version": "cbyd_route_source.v1", "route_id": route_id,
            "source": {"corpus_id": "cbyd", "release_id": "vTest", "reference": "Test article"},
            "raw_route": {"synthesis": rows}, "publisher_graph": graph,
            "step_memberships": [{"source_row_number": i, "source_record_id": f"{route_id}:step:{i}",
                                  "observation_id": f"obs-{i}", "source_reaction_id": "duplicated author id"}
                                 for i in (1, 2)]}


def test_mixed_catalog_reconstructs_literature_tree_and_scopes_memberships(tmp_path):
    patent, literature = _wrapper(), _literature_wrapper()
    pairs = []
    for name, route in [("patent", patent), ("literature", literature)]:
        routes, steps = tmp_path / f"{name}.gz", tmp_path / f"{name}_steps.gz"
        pairs.append((routes, steps))
        _write(routes, [route])
    members = literature["step_memberships"]
    legacy_ids = patent["raw_route"]["original_tree"]["reaction_ids"]
    # The same source-record alias is present in another bundle. Its observation cannot join here.
    _write(pairs[0][1], [{"source": {"source_record_id": i}, "observation_id": "patent-" + i}
                         for i in legacy_ids + [members[0]["source_record_id"]]])
    _write(pairs[1][1], [{"source": m, "observation_id": m["observation_id"]} for m in members])
    destination = tmp_path / "catalog.sqlite"

    def interrupt(event):
        if event["bundle"] == 1:
            raise InterruptedError()

    with pytest.raises(InterruptedError):
        build_route_bundle_catalog(pairs, destination, progress=interrupt)
    assert not destination.exists()
    report = build_route_bundle_catalog(pairs, destination)
    assert report["counts"] == {"validated_tree": 2}
    catalog = ProcessedRouteCatalog(destination)
    summary = catalog.route(literature["route_id"], include_tree=True, include_source=True)
    assert summary["source"]["publisher_graph"] == literature["publisher_graph"]
    assert summary["publisher_graph_available"]
    assert [m["observation_ids"] for m in summary["memberships"]] == [["obs-1"], ["obs-2"]]
    assert catalog.routes_for_observation("obs-2") == [literature["route_id"]]
    assert not catalog.routes_for_observation("patent-" + members[0]["source_record_id"])
    tree = summary["tree"]
    assert tree["patent_id"] is None and tree["source_reference"] == "Test article"
    assert tree["source_dataset_id"] == "cbyd:vTest" and tree["maximum_depth"] == 2
    assert "patent" not in json.dumps(tree["warnings"]).lower()
    assert tree["root"]["reaction"]["evidence"]["abstraction_status"] == "none"
    assert len(list(catalog.iter_trees())) == 2
    with pytest.raises(ValueError, match="different sources"):
        build_route_bundle_catalog(list(reversed(pairs)), destination)


def test_literature_source_graph_conflicts_and_observation_mismatch_are_preserved(tmp_path):
    route = _literature_wrapper()
    route["publisher_graph"]["nodes"][0]["s"] = "C"
    routes, steps, output = [tmp_path / name for name in ("routes.gz", "steps.gz", "routes.sqlite")]
    _write(routes, [route])
    _write(steps, [{"source": m, "observation_id": "inconsistent-" + m["observation_id"]}
                   for m in route["step_memberships"]])
    report = build_processed_route_catalog(routes, steps, output)
    assert report["conflicting_memberships"] == 2
    catalog = ProcessedRouteCatalog(output)
    result = catalog.route(route["route_id"], include_source=True)
    assert result["status"] == "unresolved" and result["unresolved_reason"] == "source_graph_product_mismatch"
    assert result["source"] == route and result["publisher_graph_available"]
    assert all(m["join_status"] == "conflicting_expected_observation" for m in result["memberships"])
    assert not catalog.routes_for_observation("inconsistent-obs-1")
    assert not list(catalog.iter_trees())


def test_literature_raw_chemistry_and_repeated_edges_cannot_silently_qualify():
    route = _literature_wrapper()
    route["raw_route"]["synthesis"][0]["reaction_smiles"] = "CO>>C"
    with pytest.raises(RouteQualityError, match="source_mapped_products_disagree"):
        normalize_released_route(route)
    route = _literature_wrapper()
    route["publisher_graph"]["edges"].append(route["publisher_graph"]["edges"][0])
    with pytest.raises(RouteQualityError, match="publisher_graph_repeated_edge"):
        normalize_released_route(route)


@pytest.mark.parametrize(("issues", "signature", "disposition"), [
    ([], True, "structural_precedent"),
    ([], False, "source_lookup_only"),
    (["source_mapped_products_disagree"], True, "quarantined"),
    (None, True, "quarantined"),
])
def test_screen_preserves_route_positions_and_blocks_trees_across_gaps(
    tmp_path, issues, signature, disposition,
):
    from condition_recommender.processed_catalog import build_processed_catalog
    from condition_recommender.record_storage import write_object_records

    wrapper = _literature_wrapper()
    members = wrapper["step_memberships"]
    routes, steps, evidence, records, output = [tmp_path / p for p in
        ("route.gz", "steps.gz", "evidence.sqlite", "records.jsonl", "routes.sqlite")]
    _write(routes, [wrapper])
    _write(steps, [{"observation_id": m["observation_id"], "source": m} for m in members])
    rows = []
    for i, member in enumerate(members):
        mapping = issues if i == 0 else []
        assessment = {"mapping": {"issues": mapping}} if mapping is not None else None
        rows.append({"observation_id": member["observation_id"], "reaction_id": member["source_record_id"],
                     "source_dataset": "cbyd", "source": {"source_transfer_assessment": assessment},
                     "reaction_signature": {"id": "signature"} if signature or i else None,
                     "admission_reasons": ["partial_atom_mapping", "source_model_mapping_requires_review"]})
    with records.open("w", encoding="utf-8") as stream:
        write_object_records(stream, rows)
    build_processed_catalog(records, evidence)
    build_route_bundle_catalog([(routes, steps)], output, evidence_catalog=evidence)
    catalog = ProcessedRouteCatalog(output)
    route = catalog.route(wrapper["route_id"], include_tree=True, include_source=True)
    assert route["source"] == wrapper
    assert [m["source_row_number"] for m in route["memberships"]] == [1, 2]
    assert [m["structural_disposition"] for m in route["memberships"]] == [disposition, "structural_precedent"]
    eligible = disposition == "structural_precedent"
    assert ("tree" in route) == eligible
    assert len(list(catalog.iter_trees())) == int(eligible)


def test_unchanged_verified_legacy_catalog_can_seed_new_build(tmp_path, monkeypatch):
    from core_retrosynthesis import processed_routes

    wrapper = _wrapper()
    routes, steps = tmp_path / "routes.gz", tmp_path / "steps.gz"
    aliases = wrapper["raw_route"]["original_tree"]["reaction_ids"]
    _write(routes, [wrapper])
    _write(steps, [{"observation_id": a, "source": {"source_record_id": a}} for a in aliases])
    old_version = processed_routes.LEGACY_ROUTE_CATALOG_VERSION
    binding = hashlib.sha256((file_sha256(routes) + file_sha256(steps) + old_version).encode()).hexdigest()
    metadata = {"schema_version": old_version, "build_complete": True, "source_binding": binding}
    release = tmp_path / "old"
    release.mkdir()
    legacy = release / "routes.sqlite"
    tree = build_observed_route_tree(normalize_released_route(wrapper)).to_dict()
    with sqlite3.connect(legacy) as db:
        db.executescript("""
            CREATE TABLE metadata(payload);
            CREATE TABLE routes(id,status,reason,source,tree);
            CREATE TABLE route_steps(route_id,reaction_id);
            CREATE TABLE step_observations(reaction_id,observation_id);
        """)
        db.execute("INSERT INTO metadata VALUES (?)", (json.dumps(metadata),))
        db.execute("INSERT INTO routes VALUES ('US_TEST_0','validated_tree',NULL,?,?)",
                   (processed_routes._encode(wrapper), processed_routes._encode(tree)))
        db.executemany("INSERT INTO route_steps VALUES ('US_TEST_0',?)", [(a,) for a in aliases])
        db.executemany("INSERT INTO step_observations VALUES (?,?)", [(a, a) for a in aliases])
    manifest = {"schema_version": "processed_reaction_release.v1", "build_complete": True,
                "artifacts": {"route_catalog": {"relative_path": legacy.name, "size_bytes": legacy.stat().st_size,
                                                "sha256": file_sha256(legacy)}}}
    (release / "manifest.json").write_text(json.dumps(manifest))

    def unexpected_conversion(*args):
        pytest.fail("Checksum-bound patent bundle should have been reused")

    monkeypatch.setattr(processed_routes, "_convert_wrapper", unexpected_conversion)
    dest = tmp_path / "new.sqlite"
    report = build_route_bundle_catalog([(routes, steps)], dest, reuse_release=release)
    assert report["counts"] == {"validated_tree": 1}
    assert ProcessedRouteCatalog(legacy).routes_for_observation("root") == ["US_TEST_0"]
    assert ProcessedRouteCatalog(dest).route("US_TEST_0", include_tree=True)["tree"] == tree

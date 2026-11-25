"""Checks for occurrence identity, source scope and release-level integrity."""

import base64
import gzip
import hashlib
import io
import json
import tarfile
from pathlib import Path
from typing import Any

import pytest

from condition_recommender.ingestion.artifacts import PreprocessingCancelled
from condition_recommender.ingestion.artifacts import preprocess_files
from condition_recommender.ingestion.migration import (
    regenerate_intermediate_datasets,
    validate_intermediate_manifest,
)
from condition_recommender.conversion.input_schema import discover_conversion_datasets
from condition_recommender.ingestion.cbyd import (
    CBYD_ROUTE_SCHEMA_VERSION,
    _outcomes,
    prepare_cbyd_release,
)


def _dataset(root: Path, *, bad_graph: bool = False) -> tuple[dict, dict]:
    root.mkdir(exist_ok=True)
    salt = "[Na+].[Cl-]"
    first = {
        "reaction_id": "duplicate-author-id",
        "reaction_smiles": "CN->[Cu].Cl>>CN",
        "reaction_smiles_mapped": "[CH3:1][NH2:2]->[Cu]>Cl>[CH3:1][NH2:2]",
        "mapping_confidence": 0.78,
        "conditions": {
            "raw_description": "salt; 25 °C; 72%",
            "TEMPERATURE": ["25"],
            "SOLVENTS": ["salt"],
            "SOLVENTS_SMILES": [salt],
            "YIELD": ["72%"],
        },
    }
    second = {
        "reaction_id": first["reaction_id"],
        "reaction_smiles": "CO>>CO.CN",
        "reaction_smiles_mapped": "[CH3:1][OH:2]>>[CH3:1][OH:2].CN",
        "mapping_confidence": None,
        "conditions": {
            "raw_description": "warm from 0 °C to rt; 86% over two steps",
            "GENERAL_REAGENTS": ["A", "B"],
            "GENERAL_REAGENTS_SMILES": [salt],
            "MODIFIERS": ["air"],
            "TEMPERATURE": ["0", "rt"],
            "YIELD": ["86%"],
        },
    }
    route = {
        "background": {
            "name": "Reported target",
            "doi": "https://doi.org/10.1234/ABC",
            "citation": "Original reference",
            "extra": None,
        },
        "synthesis": [first, second],
        "preserved_field": ["extra"],
    }
    graph = {
        "rx": [
            {"rid": s["reaction_id"], "m": s["reaction_smiles_mapped"]}
            for s in route["synthesis"]
        ],
        "nodes": [{"id": "r0", "t": "r", "i": 0, "inv": True}],
        "edges": [["r0", "r0"], ["r0", "r0"]],
        "warn": ["disconnected step"],
    }
    if bad_graph:
        graph["rx"][1]["rid"] = "unrelated-row"
    encoded = base64.b64encode(
        gzip.compress(json.dumps(graph).encode(), mtime=0)
    ).decode()
    index = [{"name": "A", "n_steps": 2}, {"name": "B", "n_steps": 2}]
    embedded = {"route-index": index, "route-data": {"A": encoded, "B": encoded}}
    dashboard = "".join(
        f'<script id="{k}" type="application/json">{json.dumps(v)}</script>'
        for k, v in embedded.items()
    ).encode()
    member = json.dumps(route, indent=2).encode()
    tarpath = root / "cbyd-v1.tar.gz"
    with tarfile.open(tarpath, "w:gz") as tar:
        for name, content in (
            ("dashboard.html", dashboard),
            ("data/routes/A.json", member),
            ("data/routes/B.json", member),
        ):
            item = tarfile.TarInfo("cbyd-1/" + name)
            item.size = len(content)
            tar.addfile(item, io.BytesIO(content))
    (root / "provenance.json").write_text(
        json.dumps(
            {
                "version": "v1",
                "archive_sha256": hashlib.sha256(tarpath.read_bytes()).hexdigest(),
                "archive_bytes": tarpath.stat().st_size,
                "verified_route_blobs": 2,
                "authors": "Authors",
                "license": "CC BY 4.0",
            }
        )
    )
    return route, graph


def _rows(path: Path) -> list[dict[str, Any]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle]


def test_lossless_occurrences_scope_and_linkage(tmp_path):
    root, out = tmp_path / "release", tmp_path / "bundle"
    route, graph = _dataset(root)
    report = prepare_cbyd_release(root, out)
    assert report["coverage_complete"] and not report["reused"]
    assert report["counts"]["route_count"] == 2
    assert report["counts"]["output_step_observations"] == 4
    assert report["counts"]["distinct_source_reaction_ids"] == 1
    routes = _rows(out / "routes.source.jsonl.gz")
    steps = _rows(out / "route_steps.observations.jsonl.gz")
    assert len({r["observation_id"] for r in steps}) == 4
    for wrapper in routes:
        assert wrapper["schema_version"] == CBYD_ROUTE_SCHEMA_VERSION
        assert wrapper["publisher_graph"] == graph and wrapper["raw_route"] == route
        assert wrapper["source"]["reference"] == "10.1234/abc"
        for member in wrapper["step_memberships"]:
            pos = member["source_row_number"]
            row = next(
                s for s in steps if s["observation_id"] == member["observation_id"]
            )
            assert member["source_reaction_id"] == "duplicate-author-id"
            assert row["source"]["source_record_id"] == member["source_record_id"]
            assert member["source_record_id"] == f"{wrapper['route_id']}:step:{pos}"
            assert row["raw_fields"]["source_step"] == route["synthesis"][pos - 1]
            assert (
                row["source"]["source_file_sha256"]
                == wrapper["source"]["source_file_sha256"]
            )

    first, pooled = steps[:2]
    reaction = first["reaction"]
    assert reaction["reaction_smiles"] == "[CH3:1][NH2:2]->[Cu]>>[CH3:1][NH2:2]"
    assert reaction["source_labels"]["source_agent_smiles"] == "Cl"
    assert reaction["source_labels"]["mapping_origin"] == "RXNMapper"
    assert reaction["supplied_mapping_status"] == "supplied_unvalidated"
    (claim,) = first["conditions"]["components"]
    assert [i["value"] for i in claim["identifiers"]] == ["salt", "[Na+].[Cl-]"]
    assert (
        claim["introduced_in_stage"] is None
        and claim["provenance"]["requires_participation_check"]
    )
    assert first["conditions"]["stages"][0]["temperature_c"] == 25
    assert first["conditions"]["procedure_text"] == ""
    assert first["outcomes"][0]["outcome_type"] == "reported_yield_pct"
    assert first["outcomes"][0]["metadata"]["primary_eligible"]
    assert (
        pooled["conditions"]["components"] == []
        and pooled["conditions"]["stages"] == []
    )
    (group,) = pooled["conditions"]["component_groups"]
    assert group["provenance"]["names"] == ["A", "B"]
    assert (
        "unpaired_condition_evidence:GENERAL_REAGENTS"
        in pooled["conditions"]["warnings"]
    )
    assert "unresolved_source_modifiers" in pooled["conditions"]["warnings"]
    (outcome,) = pooled["outcomes"]
    assert (outcome["outcome_type"], outcome["value"]) == (
        "reported_sequence_yield",
        86,
    )
    assert (
        outcome["metadata"]["scope"] == "sequence"
        and not outcome["metadata"]["primary_eligible"]
    )


def test_reproducible_bundle_retry_and_cancel(tmp_path):
    root, out = tmp_path / "release", tmp_path / "bundle"
    _dataset(root)
    first = prepare_cbyd_release(root, out)
    assert prepare_cbyd_release(root, out)["reused"]
    rebuilt = prepare_cbyd_release(root, tmp_path / "separate", force=True)
    assert first["output_sha256"] == rebuilt["output_sha256"]
    (out / "route_steps.observations.jsonl.gz").write_bytes(b"incomplete")
    repaired = prepare_cbyd_release(root, out)
    assert (
        not repaired["reused"] and repaired["output_sha256"] == first["output_sha256"]
    )
    originals = {path.name: path.read_bytes() for path in out.iterdir()}
    checks = iter([False, False, True])
    with pytest.raises(PreprocessingCancelled):
        prepare_cbyd_release(root, out, force=True, cancel_check=lambda: next(checks))
    assert originals == {path.name: path.read_bytes() for path in out.iterdir()}


def test_release_verifies_archive_and_dashboard_before_publish(tmp_path):
    root, out = tmp_path / "release", tmp_path / "bundle"
    _dataset(root)
    prepare_cbyd_release(root, out)
    original = (out / "route_release.manifest.json").read_bytes()
    _dataset(root, bad_graph=True)
    with pytest.raises(ValueError, match="step correspondence"):
        prepare_cbyd_release(root, out)
    assert original == (out / "route_release.manifest.json").read_bytes()
    tarpath = root / "cbyd-v1.tar.gz"
    with tarpath.open("ab") as stream:
        stream.write(b"unverified")
    with pytest.raises(ValueError, match="pinned"):
        prepare_cbyd_release(root, out)


def test_migration_and_selected_file_publish_the_same_release(tmp_path):
    raw = tmp_path / "raw"
    relative = Path("routes/cbyd/v1")
    release = raw / relative
    release.mkdir(parents=True)
    _dataset(release)
    migrated, selected = tmp_path / "migrated", tmp_path / "selected"
    report = regenerate_intermediate_datasets(raw, migrated, cbyd_release="v1")
    assert report["route_count"] == 2 and report["single_step_observation_count"] == 4
    (entry,) = report["sources"]
    assert entry["bundle_relative_path"] == relative.as_posix()
    assert entry["source_relative_path"] == "routes/cbyd/v1/cbyd-v1.tar.gz"
    assert validate_intermediate_manifest(migrated) == {
        "coverage_complete": True,
        "verified_artifact_count": 2,
    }
    batch = preprocess_files([release / "cbyd-v1.tar.gz"], selected, source_root=raw)
    assert batch["file_count"] == 1 and batch["row_count"] == 4
    (file_report,) = batch["files"]
    assert file_report["output_sha256"] == entry["report"]["output_sha256"]
    assert file_report["definition_version"] == "cbyd_release.v1"
    assert (
        Path(file_report["output_path"])
        == selected / relative / "route_steps.observations.jsonl.gz"
    )
    for output in (selected, migrated):
        (obs,) = discover_conversion_datasets(output)
        assert obs == output / relative / "route_steps.observations.jsonl.gz"
        row = _rows(obs)[0]
        assert row["source"]["source_record_id"] == "cbyd:v1:A:step:1"


@pytest.mark.parametrize(
    ("field", "arrow", "product", "kind", "value", "basis"),
    [
        (
            "86%",
            "salt; 86% overall",
            "CC",
            "reported_aggregate_yield",
            86,
            "source_unspecified",
        ),
        ("80%", "solvent; 80% (NMR)", "CC", "reported_yield_pct", 80, "nmr"),
        (
            "isolated yield 80%",
            "isolated yield 80%",
            "CC",
            "reported_yield_pct",
            80,
            "isolated",
        ),
        ("80%", "80%", "CC.O", "reported_outcome_text", 80, "source_unspecified"),
        (
            "80%",
            "80% based on recovered substrate",
            "CC",
            "reported_recovered_material_yield",
            80,
            "source_unspecified",
        ),
        (
            "75% (E/Z=80:20)",
            "75% (E/Z=80:20)",
            "CC",
            "reported_selectivity_or_yield",
            None,
            "source_unspecified",
        ),
        ("92% ee", "92% ee", "CC", "reported_selectivity", 92, "source_unspecified"),
        (">90%", ">90%", "CC", "reported_outcome_text", None, "source_unspecified"),
        (
            "96%",
            "50%; 96% ee",
            "CC",
            "reported_selectivity_or_yield",
            96,
            "source_unspecified",
        ),
        ("150%", "150%", "CC", "reported_outcome_text", None, "source_unspecified"),
    ],
)
def test_outcome_scope(field, arrow, product, kind, value, basis):
    (result,) = _outcomes({"YIELD": [field], "raw_description": arrow}, product)
    assert (result.outcome_type, result.value, result.raw_value) == (kind, value, field)
    assert result.metadata["measurement_basis"] == basis
    assert result.metadata["primary_eligible"] == (kind == "reported_yield_pct")

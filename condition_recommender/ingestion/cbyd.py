"""Prepare occurrence-linked Chemistry by Design routes and step observations.

The tagged archive is the authority for both the route members and the publisher's
display graphs. Reported order, roles, process conditions and atom maps are source
evidence. Chemistry qualification belongs to the downstream converter.
"""

from __future__ import annotations

import base64
import gzip
import hashlib
import json
import re
import tarfile
from collections import Counter
from pathlib import Path, PurePosixPath
from tempfile import TemporaryDirectory
from typing import Any, Callable

from reactive_taxonomy.chemistry.reaction_smiles import split_reaction_fields

from .adapters.base import observation_id, supplied_mapping_status
from .artifacts import PreprocessingCancelled, PreprocessingProgress, _sha256
from .models import (
    INTERMEDIATE_OBSERVATION_SCHEMA_VERSION,
    CanonicalSourceObservation,
    ConditionComponentClaim,
    ConditionComponentGroup,
    ConditionInput,
    ConditionStageInput,
    OutcomeInput,
    ReactionEvidenceInput,
    SourceIdentifier,
    SourceProvenance,
)
from .route_release import _compressed_writer, _json

CBYD_DEFINITION_VERSION = "cbyd_release.v1"
CBYD_ADAPTER_VERSION = "1.1"
CBYD_ROUTE_SCHEMA_VERSION = "cbyd_route_source.v1"
_ROUTES = "routes.source.jsonl.gz"
_STEPS = "route_steps.observations.jsonl.gz"
_MANIFEST = "route_release.manifest.json"
_SLOTS = {
    "GENERAL_REAGENTS": "reagent",
    "CATALYSTS": "catalyst",
    "LIGANDS": "ligand",
    "SOLVENTS": "solvent",
}
_NUMBER = r"[+-]?(?:\d+(?:\.\d+)?|\.\d+)"
_PERCENT = re.compile(rf"({_NUMBER})\s*%")
_CELSIUS = re.compile(rf"({_NUMBER})\s*°\s*C\b", re.I)
_SEQUENCE = re.compile(
    r"\b(?:[2-9]|\d{2,}|two|three|four|five|six|seven|eight|nine|ten|multiple|several)"
    r"\s*[- ]?steps?\b|\b(?:multi[- ]?step|(?:over|overall)\b.{0,12}\bsteps?)\b",
    re.I,
)
_AGGREGATE = re.compile(r"\b(?:overall|pooled|combined|total yield)\b", re.I)
_RECOVERED = re.compile(
    r"\b(?:brsm|borsm|recovered|recovery|(?:consumed|unreacted) (?:starting material|substrate))\b",
    re.I,
)
_BASIS = r"(?:isolated|q?NMR|HPLC|GC)"
_YIELD_SEGMENT = re.compile(
    rf"(?:(?P<prefix>{_BASIS})(?:\s+yield)?\s*[:=]?\s*)?"
    rf"(?P<value>{_NUMBER})\s*%(?:\s*\(?(?P<suffix>{_BASIS})(?:\s+yield)?\)?)?",
    re.I,
)
_SELECTIVITY = re.compile(
    r"\b(?:e\.?e\.?|d\.?e\.?|d\.?r\.?|e\.?r\.?|E/Z|enantioselectivity|"
    r"diastereoselectivity|enantiomeric excess|diastereomeric ratio)\b",
    re.I,
)
_ORDERED = re.compile(r"(?:^|[;\s])(?:[1-9]|i{1,3})\)|\b(?:then|followed by)\b", re.I)
_TEMP_CHANGE = re.compile(
    r"\b(?:to|rt|r\.t|reflux|ambient|room temp\w*|warming|cooling|ramp\w*)\b", re.I
)


def _values(fields: dict[str, Any], name: str) -> list[str]:
    value = fields.get(name, [])
    if not isinstance(value, list) or any(
        item is not None and not isinstance(item, str) for item in value
    ):
        raise ValueError(f"Invalid CbyD source array {name}")
    return [item or "" for item in value]


def _reference(background: dict[str, Any]) -> str:
    # Import here so inventory discovery does not load the conversion engine.
    from ..conversion.references import normalize_reference

    doi, citation = (str(background.get(k) or "") for k in ("doi", "citation"))
    identities = [normalize_reference(value) for value in (doi, citation)]
    for key in ("doi", "patent_number"):
        for identity in identities:
            if value := getattr(identity, key):
                return value
    return citation


def _reaction_parts(value: str) -> tuple[str, str, str] | None:
    parts = split_reaction_fields(value)
    if len(parts) == 3 and parts[0] and parts[2]:
        return parts[0], parts[1], parts[2]
    return None


def _temperature(
    fields: dict[str, Any], description: str
) -> tuple[ConditionStageInput, ...]:
    values = _values(fields, "TEMPERATURE")
    if len(values) != 1 or not re.fullmatch(_NUMBER, values[0].strip()):
        return ()
    text = description.replace("−", "-")
    matches = _CELSIUS.findall(text)
    # Require a standalone source segment, rather than the endpoint of a range.
    if (
        len(matches) != 1
        or float(matches[0]) != float(values[0])
        or _TEMP_CHANGE.search(text)
        or _ORDERED.search(text)
        or _SEQUENCE.search(text)
    ):
        return ()
    segment = next(
        (s.strip() for s in text.split(";") if _CELSIUS.fullmatch(s.strip())), ""
    )
    if not segment:
        return ()
    return (
        ConditionStageInput(
            stage_index=1,
            temperature_c=float(matches[0]),
            source_text=segment,
            provenance={
                "source": "TEMPERATURE and raw_description",
                "component_order": "unknown",
            },
        ),
    )


def _conditions(fields: dict[str, Any]) -> ConditionInput:
    components, groups, warnings = [], [], []
    represented = []
    for slot, role in _SLOTS.items():
        names, structures = _values(fields, slot), _values(fields, slot + "_SMILES")
        represented.extend(names)
        if len(names) != len(structures) or any(
            not name.strip() or not smiles.strip()
            for name, smiles in zip(names, structures)
        ):
            warnings.append(f"unpaired_condition_evidence:{slot}")
            groups.append(
                ConditionComponentGroup(
                    group_key=slot,
                    source_slot=slot,
                    provenance={
                        "requires_review": True,
                        "names": names,
                        "structures": structures,
                        "alignment": "independent_source_arrays",
                    },
                )
            )
            continue
        for position, (name, smiles) in enumerate(zip(names, structures), 1):
            components.append(
                ConditionComponentClaim(
                    component_key=f"{slot}:{position}",
                    source_slot=slot,
                    source_role_hint=role,
                    identifiers=(
                        SourceIdentifier("name", name, slot),
                        SourceIdentifier("smiles", smiles, slot + "_SMILES"),
                    ),
                    introduced_in_stage=None,
                    provenance={
                        "source_position": position,
                        "requires_participation_check": True,
                        "alignment": "source_position_only",
                    },
                )
            )
    description = str(fields.get("raw_description") or "")
    modifiers = _values(fields, "MODIFIERS")
    if any(modifiers):
        warnings.append("unresolved_source_modifiers")
    if description.strip() and not components:
        warnings.append("raw_only_condition_evidence")
    # Only exact source labels and clearly typed operating/outcome segments are
    # accounted for. Additional arrow text may identify an unlisted material.
    represented.extend(modifiers + _values(fields, "YIELD"))
    represented_text = {s.strip().casefold() for s in represented}
    for segment in (s.strip() for s in description.split(";")):
        if (
            segment
            and segment.casefold() not in represented_text
            and not _CELSIUS.fullmatch(segment.replace("−", "-"))
            and not re.fullmatch(r"(?:r\.?t\.?|reflux|room temperature)", segment, re.I)
            and not re.fullmatch(
                rf"{_NUMBER}\s*(?:h|hr|hours?|min|minutes?|days?)", segment, re.I
            )
        ):
            warnings.append("unresolved_raw_condition_evidence")
            break
    stages = _temperature(fields, description)
    if any(_values(fields, "TEMPERATURE")) and not stages:
        warnings.append("unresolved_temperature_scope")
    return ConditionInput(
        components=tuple(components),
        component_groups=tuple(groups),
        stages=stages,
        warnings=tuple(warnings),
    )


def _outcomes(fields: dict[str, Any], product: str) -> tuple[OutcomeInput, ...]:
    values = [v for v in _values(fields, "YIELD") if v.strip()]
    description = str(fields.get("raw_description") or "")
    results = []
    for raw in values:
        token = raw.strip()
        numeric_source = _YIELD_SEGMENT.fullmatch(token)
        matches = _PERCENT.findall(token)
        number = (
            float(matches[0])
            if len(matches) == 1 and (numeric_source or _PERCENT.match(token))
            else None
        )
        if number is not None and not 0 <= number <= 100:
            number = None
        combined = token + " " + description
        basis = "source_unspecified"
        # An isolated percentage in an arrow summary can be an assay, loading
        # or selectivity. Match its whole source segment before assigning it.
        corroborated = [
            m
            for segment in description.split(";")
            if (m := _YIELD_SEGMENT.fullmatch(segment.strip()))
            and number is not None
            and float(m["value"]) == number
        ]
        bases = {
            label.casefold()
            for match in [numeric_source, *corroborated]
            if match
            for label in (match["prefix"], match["suffix"])
            if label
        }
        if len(bases) == 1:
            basis = next(iter(bases))
        scope, kind = "unassigned", "reported_outcome_text"
        if _SELECTIVITY.search(token):
            # The leading percent in "75% (E/Z=80:20)" is not the E/Z ratio.
            if not matches or re.fullmatch(
                rf"{_NUMBER}\s*%\s*(?:e\.?e\.?|d\.?e\.?)", token, re.I
            ):
                kind = "reported_selectivity"
            else:
                kind, number = "reported_selectivity_or_yield", None
        elif _RECOVERED.search(combined):
            kind = "reported_recovered_material_yield"
        elif _SEQUENCE.search(combined):
            scope, kind = "sequence", "reported_sequence_yield"
        elif _AGGREGATE.search(combined):
            scope, kind = "aggregate", "reported_aggregate_yield"
        elif _SELECTIVITY.search(combined):
            kind = "reported_selectivity_or_yield"
        elif (
            numeric_source
            and corroborated
            and len(bases) <= 1
            and number is not None
            and len(values) == 1
            and product
            and "." not in product
            and not _ORDERED.search(description)
            and _PERCENT.findall(description) == matches
        ):
            scope, kind = "step", "reported_yield_pct"
        metadata = {
            "scope": scope,
            "primary_eligible": scope == "step",
            "measurement_basis": basis,
        }
        if scope == "step":
            metadata["product_smiles"] = product
        results.append(
            OutcomeInput(
                kind,
                number,
                "percent" if number is not None else "",
                raw,
                "YIELD",
                metadata,
            )
        )
    return tuple(results)


def _step_observation(
    step: dict[str, Any], position: int, route_id: str, source: dict[str, Any]
) -> CanonicalSourceObservation:
    raw = str(step.get("reaction_smiles") or "")
    mapped = str(step.get("reaction_smiles_mapped") or "")
    mapped_parts, raw_parts = _reaction_parts(mapped), _reaction_parts(raw)
    left, middle, right = mapped_parts or raw_parts or ("", "", "")
    transformation = left + ">>" + right if left and right else ""
    mapping_status = supplied_mapping_status(transformation)
    warnings = []
    if not mapped_parts:
        warnings.append("missing_or_incomplete_mapped_reaction")
    if not raw_parts:
        warnings.append("missing_or_incomplete_unmapped_reaction")
    occurrence = f"{route_id}:step:{position}"
    provenance = SourceProvenance(
        corpus_id="cbyd",
        release_id=source["release_id"],
        adapter_id=CBYD_DEFINITION_VERSION,
        adapter_version=CBYD_ADAPTER_VERSION,
        source_file=source["source_file"],
        source_file_sha256=source["source_file_sha256"],
        source_row_number=position,
        source_record_id=occurrence,
        reference=source["reference"],
        source_groups={
            "route_id": route_id,
            "target_name": source["target_name"],
            "source_doi": source["source_doi"],
            "source_citation": source["source_citation"],
            "source_reaction_id": str(step.get("reaction_id") or ""),
            "route_membership_artifact": _ROUTES,
        },
    )
    fields = step.get("conditions") or {}
    return CanonicalSourceObservation(
        observation_id=observation_id(
            adapter_id=provenance.adapter_id,
            source_sha256=provenance.source_file_sha256,
            row_number=position,
            record_id=occurrence,
        ),
        observation_kind="structure_backed" if transformation else "source_label_only",
        source=provenance,
        reaction=ReactionEvidenceInput(
            evidence_kind="source_structure",
            reaction_smiles=transformation or None,
            structure_available=bool(transformation),
            supplied_mapping_status=mapping_status,
            source_labels={
                "source_reaction_smiles": mapped,
                "source_unmapped_reaction_smiles": raw,
                "source_agent_smiles": middle,
                "mapping_confidence": step.get("mapping_confidence"),
                "mapping_origin": "RXNMapper" if mapped_parts else "",
            },
        ),
        conditions=_conditions(fields),
        outcomes=_outcomes(fields, right),
        ingestion_status="review" if warnings else "accepted",
        warnings=tuple(warnings),
        raw_fields={"source_step": step},
    )


def _embedded(dashboard: bytes, element: str) -> Any:
    start = f'<script id="{element}" type="application/json">'.encode()
    begin = dashboard.find(start)
    end = dashboard.find(b"</script>", begin + len(start))
    if begin == -1 or end == -1:
        raise ValueError(f"Missing CbyD dashboard element {element}")
    return json.loads(dashboard[begin + len(start) : end])


def _member_bytes(archive: tarfile.TarFile, member: tarfile.TarInfo) -> bytes:
    if not member.isfile():
        raise ValueError(f"Expected regular archive member {member.name}")
    with archive.extractfile(member) as handle:
        return handle.read()


def _write_bundle(
    archive: tarfile.TarFile,
    staging: Path,
    provenance: dict[str, Any],
    check_cancelled: Callable[[], None],
    emit: Callable[[str, int, str], None],
) -> dict[str, Any]:
    members = archive.getmembers()
    dashboards = [m for m in members if re.fullmatch(r"[^/]+/dashboard\.html", m.name)]
    if len(dashboards) != 1:
        raise ValueError("Expected exactly one deposited CbyD dashboard")
    dashboard = _member_bytes(archive, dashboards[0])
    dashboard_hash = hashlib.sha256(dashboard).hexdigest()
    payloads = _embedded(dashboard, "route-data")
    index_rows = _embedded(dashboard, "route-index")
    index = {row["name"]: row for row in index_rows}
    del dashboard
    prefix = dashboards[0].name.split("/")[0] + "/data/routes/"
    routes = [
        m for m in members if m.name.startswith(prefix) and m.name.endswith(".json")
    ]
    keys = {PurePosixPath(m.name).stem for m in routes}
    if (
        not routes
        or len(keys) != len(routes)
        or len(index) != len(index_rows)
        or keys != set(index)
        or keys != set(payloads)
        or provenance.get("verified_route_blobs", len(routes)) != len(routes)
    ):
        raise ValueError("CbyD archive/dashboard route inventory mismatch")

    counts: Counter[str] = Counter()
    seen_observations, author_ids = set(), set()
    release = provenance["version"]
    emit("routes", 0, f"Preparing {len(routes):,} Chemistry by Design routes.")
    with (
        _compressed_writer(staging / (_ROUTES + ".pending")) as route_out,
        _compressed_writer(staging / (_STEPS + ".pending")) as step_out,
    ):
        for member in routes:
            check_cancelled()
            key = PurePosixPath(member.name).stem
            route_id = f"cbyd:{release}:{key}"
            member_bytes = _member_bytes(archive, member)
            route = json.loads(member_bytes)
            graph = json.loads(
                gzip.decompress(base64.b64decode(payloads[key], validate=True))
            )
            steps = route["synthesis"]
            if len(steps) != index[key]["n_steps"] or [
                (s.get("reaction_id"), s.get("reaction_smiles_mapped") or "")
                for s in steps
            ] != [(r["rid"], r["m"]) for r in graph["rx"]]:
                raise ValueError(f"CbyD dashboard step correspondence failed for {key}")
            source = {
                "corpus_id": "cbyd",
                "release_id": release,
                "source_file": member.name.split("/", 1)[1],
                "source_file_sha256": hashlib.sha256(member_bytes).hexdigest(),
                "reference": _reference(route["background"]),
                "target_name": str(route["background"].get("name") or ""),
                "source_doi": str(route["background"].get("doi") or ""),
                "source_citation": str(route["background"].get("citation") or ""),
            }
            memberships = []
            for position, step in enumerate(steps, 1):
                obs = _step_observation(step, position, route_id, source)
                if obs.observation_id in seen_observations:
                    raise ValueError(
                        f"Duplicate CbyD occurrence in {key} at {position}"
                    )
                seen_observations.add(obs.observation_id)
                author_ids.add(step.get("reaction_id"))
                step_out.write(_json(obs.to_dict()) + "\n")
                memberships.append(
                    {
                        "source_row_number": position,
                        "source_record_id": obs.source.source_record_id,
                        "observation_id": obs.observation_id,
                        "source_reaction_id": step.get("reaction_id"),
                    }
                )
                counts["output_step_observations"] += 1
                counts["step_status_" + obs.ingestion_status] += 1
                counts.update(set(obs.conditions.warnings))
                counts.update({o.outcome_type for o in obs.outcomes})
            route_out.write(
                _json(
                    {
                        "schema_version": CBYD_ROUTE_SCHEMA_VERSION,
                        "route_id": route_id,
                        "source": source,
                        "raw_route": route,
                        "publisher_graph": graph,
                        "publisher_index": index[key],
                        "graph_derivation": "CbyD deposited build_route_tree.py; inferred display connectivity",
                        "step_memberships": memberships,
                    }
                )
                + "\n"
            )
            counts["route_count"] += 1
            counts["step_memberships"] += len(memberships)
            if counts["route_count"] % 100 == 0:
                emit(
                    "routes",
                    counts["output_step_observations"],
                    f"Prepared {counts['route_count']:,} of {len(routes):,} routes.",
                )
    counts["distinct_source_reaction_ids"] = len(author_ids - {None, ""})
    return {"counts": dict(sorted(counts.items())), "dashboard_sha256": dashboard_hash}


def prepare_cbyd_release(
    release_root: str | Path,
    output_dir: str | Path,
    *,
    force: bool = False,
    progress_callback: Callable[[PreprocessingProgress], None] | None = None,
    cancel_check: Callable[[], bool] | None = None,
) -> dict[str, Any]:
    """Verify a local tagged archive and publish both linked intermediate streams."""
    root, destination = Path(release_root), Path(output_dir)
    provenance_path = root / "provenance.json"
    provenance_bytes = provenance_path.read_bytes()
    provenance = json.loads(provenance_bytes)
    archives = list(root.glob("cbyd-*.tar.gz"))
    if len(archives) != 1 or archives[0].name != f"cbyd-{provenance['version']}.tar.gz":
        raise ValueError(
            "Expected a single CbyD archive matching its release provenance"
        )
    archive = archives[0]
    archive_sha = _sha256(archive)
    if archive_sha != provenance[
        "archive_sha256"
    ] or archive.stat().st_size != provenance.get(
        "archive_bytes", archive.stat().st_size
    ):
        raise ValueError("CbyD archive does not match pinned release provenance")
    source_identity = {
        "definition_version": CBYD_DEFINITION_VERSION,
        "adapter_version": CBYD_ADAPTER_VERSION,
        "corpus_id": "cbyd",
        "source_sha256": archive_sha,
        "source_path": str(archive.resolve()),
        "provenance_sha256": hashlib.sha256(provenance_bytes).hexdigest(),
        "release_id": provenance["version"],
        "provenance": provenance,
    }

    def check_cancelled() -> None:
        if cancel_check and cancel_check():
            raise PreprocessingCancelled("Chemistry by Design preparation cancelled")

    def emit(phase: str, count: int, message: str) -> None:
        if progress_callback:
            progress_callback(
                PreprocessingProgress(phase, str(archive), 1, 1, count, message)
            )

    check_cancelled()
    manifest = destination / _MANIFEST
    if not force and manifest.is_file():
        try:
            previous = json.loads(manifest.read_bytes())
        except (UnicodeError, ValueError):
            previous = {}
        if (
            isinstance(previous, dict)
            and previous.get("coverage_complete") is True
            and all(previous.get(k) == v for k, v in source_identity.items())
            and set(previous.get("output_sha256", {})) == {_ROUTES, _STEPS}
            and all(
                (destination / name).is_file() and _sha256(destination / name) == digest
                for name, digest in previous["output_sha256"].items()
            )
        ):
            return {**previous, "reused": True}
    destination.mkdir(parents=True, exist_ok=True)
    # Even an interrupted scratch directory is outside the intermediate file globs.
    with TemporaryDirectory(prefix=".cbyd-", dir=destination.parent) as temporary:
        staging = Path(temporary)
        with tarfile.open(archive, "r:gz") as source:
            report = _write_bundle(source, staging, provenance, check_cancelled, emit)
        report.update(source_identity)
        report.update(
            {
                "coverage_complete": True,
                "reused": False,
                "schemas": {
                    _ROUTES: CBYD_ROUTE_SCHEMA_VERSION,
                    _STEPS: INTERMEDIATE_OBSERVATION_SCHEMA_VERSION,
                },
                "output_sha256": {
                    name: _sha256(staging / (name + ".pending"))
                    for name in (_ROUTES, _STEPS)
                },
            }
        )
        (staging / (_MANIFEST + ".pending")).write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        check_cancelled()
        for name in (_ROUTES, _STEPS, _MANIFEST):
            (staging / (name + ".pending")).replace(destination / name)
    emit(
        "complete",
        report["counts"]["output_step_observations"],
        "Prepared Chemistry by Design.",
    )
    return report

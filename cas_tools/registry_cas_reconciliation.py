"""Reconcile curator-supplied CAS candidates with the condition registry.

This application-layer utility composes the network-backed CAS lookup client with
the deterministic condition registry.  The registry package itself remains free
of network access.
"""

from __future__ import annotations

import csv
import json
import os
import re
import shutil
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass, replace
from difflib import SequenceMatcher
from pathlib import Path
from typing import Callable, Dict, Iterable, Mapping, Optional, Sequence, Tuple

from condition_registry import (
    CompoundAdditionRequest,
    CompoundAliasInput,
    ConditionRegistry,
    SubstanceAliasAdditionRequest,
    add_compound,
    add_substance_aliases,
)
from condition_registry.loader import (
    ADDITIONS_PATH,
    IDENTIFIERS_PATH,
    SUBSTANCES_PATH,
)
from condition_registry.normalization import normalize_cas, normalize_chemical_name

from .compound_lookup import CompoundLookupResult, lookup_compound_by_cas


CAS_PATTERN = re.compile(r"(?<!\d)(\d{2,7}-\d{2}-\d)(?!\d)")
REQUIRED_COLUMNS = {
    "name",
    "aliases",
    "roles",
    "normalized_identity",
    "possible_cas_no",
    "cas_match_status",
    "cas_confidence",
    "source_basis",
    "source_url",
    "match_notes",
}
UNSAFE_MATCH_STATUSES = {
    "ambiguous_multiple_candidates",
    "component_or_parent_only",
    "no_defensible_cas_found",
    "no_unique_cas",
    "possible_formula_match",
    "unresolved",
}
METAL_NAME_PATTERNS = {
    "Pd": re.compile(r"(?<![A-Za-z])Pd(?![a-z])|palladium", re.I),
    "Ni": re.compile(r"(?<![A-Za-z])Ni(?![a-z])|nickel", re.I),
    "Ru": re.compile(r"(?<![A-Za-z])Ru(?![a-z])|ruthenium", re.I),
    "Rh": re.compile(r"(?<![A-Za-z])Rh(?![a-z])|rhodium", re.I),
    "Ir": re.compile(r"(?<![A-Za-z])Ir(?![a-z])|iridium", re.I),
    "Pt": re.compile(r"(?<![A-Za-z])Pt(?![a-z])|platinum", re.I),
    "Au": re.compile(r"(?<![A-Za-z])Au(?![a-z])|gold", re.I),
    "Ag": re.compile(r"(?<![A-Za-z])Ag(?![a-z])|silver", re.I),
    "Cu": re.compile(r"(?<![A-Za-z])Cu(?![a-z])|copper", re.I),
    "Co": re.compile(r"(?<![A-Za-z])Co(?![a-z])|cobalt", re.I),
    "Fe": re.compile(r"(?<![A-Za-z])Fe(?![a-z])|iron", re.I),
}
ALIAS_SOURCE = "weak_label_v2.1:user_possible_cas_reconciliation"
AUDIT_FIELDNAMES = (
    "name",
    "aliases",
    "roles",
    "mention_count",
    "normalized_identity",
    "possible_cas_no",
    "cas_match_status",
    "cas_confidence",
    "cas_checksum_status",
    "registry_status_before",
    "registry_substance_id",
    "registry_canonical_name",
    "lookup_status",
    "lookup_canonical_name",
    "identity_match_method",
    "identity_match_score",
    "decision",
    "decision_reason",
    "aliases_added_or_planned",
    "source_basis",
    "source_url",
    "match_notes",
    "lookup_sources",
    "lookup_warnings",
)


LookupFunction = Callable[[str], CompoundLookupResult]


@dataclass(frozen=True)
class CasReconciliationSummary:
    """Counts and output paths from one CAS reconciliation run."""

    input_rows: int
    rows_with_one_valid_cas: int
    existing_substances_updated: int
    aliases_added: int
    new_compounds_added: int
    review_rows: int
    audit_path: str
    review_path: str
    lookup_cache_path: str


@dataclass(frozen=True)
class _InputRow:
    raw: Mapping[str, str]
    cas_values: Tuple[str, ...]

    @property
    def name(self) -> str:
        return str(self.raw.get("name") or "").strip()

    @property
    def normalized_identity(self) -> str:
        return str(self.raw.get("normalized_identity") or "").strip()

    @property
    def confidence(self) -> str:
        return str(self.raw.get("cas_confidence") or "").strip().casefold()

    @property
    def match_status(self) -> str:
        return str(self.raw.get("cas_match_status") or "").strip().casefold()


@dataclass(frozen=True)
class _GroupPlan:
    cas: str
    rows: Tuple[_InputRow, ...]
    decision: str
    reason: str
    registry_status: str
    registry_substance_id: Optional[str]
    registry_canonical_name: Optional[str]
    lookup: CompoundLookupResult
    match_method: str
    match_score: float
    aliases: Tuple[str, ...] = ()
    addition_request: Optional[CompoundAdditionRequest] = None


def _split_aliases(value: str) -> Tuple[str, ...]:
    return tuple(
        part.strip()
        for part in str(value or "").split("|")
        if part.strip()
    )


def _name_key(value: str) -> str:
    normalized = normalize_chemical_name(str(value or ""))
    return "".join(character for character in normalized if character.isalnum())


def _registry_name_key(value: str) -> str:
    """Return the exact normalization key used by name resolution."""
    return normalize_chemical_name(str(value or ""))


def _distinct_names(values: Iterable[str]) -> Tuple[str, ...]:
    selected = []
    seen = set()
    for raw in values:
        value = " ".join(str(raw or "").split())
        key = _registry_name_key(value)
        if not value or not key or key in seen:
            continue
        seen.add(key)
        selected.append(value)
    return tuple(selected)


def _parse_rows(path: Path) -> Tuple[_InputRow, ...]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or ())
        if missing:
            raise ValueError(
                "CAS reconciliation input is missing columns: "
                + ", ".join(sorted(missing))
            )
        rows = []
        for raw in reader:
            extracted = CAS_PATTERN.findall(str(raw.get("possible_cas_no") or ""))
            rows.append(
                _InputRow(
                    raw=dict(raw),
                    cas_values=tuple(
                        cas for item in extracted if (cas := normalize_cas(item))
                    ),
                )
            )
    return tuple(rows)


def _lookup_to_dict(result: CompoundLookupResult) -> Dict[str, object]:
    return asdict(result)


def _lookup_from_dict(value: Mapping[str, object]) -> CompoundLookupResult:
    tuple_fields = {"synonyms", "source_ids", "source_urls", "warnings"}
    prepared = {
        key: tuple(item) if key in tuple_fields and isinstance(item, list) else item
        for key, item in value.items()
    }
    return CompoundLookupResult(**prepared)


def _load_lookup_cache(path: Path) -> Dict[str, CompoundLookupResult]:
    if not path.exists():
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, Mapping):
        return {}
    return {
        str(cas): _lookup_from_dict(value)
        for cas, value in payload.items()
        if isinstance(value, Mapping)
    }


def _write_lookup_cache(
    path: Path,
    values: Mapping[str, CompoundLookupResult],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        cas: _lookup_to_dict(result)
        for cas, result in sorted(values.items())
    }
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _resolve_lookups(
    cas_values: Iterable[str],
    *,
    cache_path: Path,
    lookup: LookupFunction,
    workers: int,
) -> Dict[str, CompoundLookupResult]:
    cache = _load_lookup_cache(cache_path)
    requested = set(cas_values)
    retryable = {
        cas
        for cas, result in cache.items()
        if cas in requested
        and not result.smiles
        and "CORE_LOOKUP_RETRIED" not in result.warnings
        and (
            result.status == "lookup_error"
            or any(
                warning.startswith("PUBCHEM_PROPERTY_LOOKUP_FAILED")
                for warning in result.warnings
            )
        )
    }
    missing = sorted((requested - set(cache)) | retryable)
    with ThreadPoolExecutor(max_workers=max(1, int(workers))) as executor:
        futures = {executor.submit(lookup, cas): cas for cas in missing}
        for future in as_completed(futures):
            cas = futures[future]
            try:
                result = future.result()
                if cas in retryable and not result.smiles:
                    result = replace(
                        result,
                        warnings=tuple(
                            dict.fromkeys((*result.warnings, "CORE_LOOKUP_RETRIED"))
                        ),
                    )
                cache[cas] = result
            except Exception as error:  # pragma: no cover - defensive network edge
                cache[cas] = CompoundLookupResult(
                    cas=cas,
                    status="lookup_error",
                    warnings=(f"LOOKUP_ERROR:{type(error).__name__}:{error}",),
                )
    _write_lookup_cache(cache_path, cache)
    return cache


def _identity_match(
    proposed: Sequence[str],
    known: Sequence[str],
) -> Tuple[str, float]:
    best_method = "no_name_evidence"
    best_score = 0.0
    for proposed_name in proposed:
        proposed_key = _name_key(proposed_name)
        if not proposed_key:
            continue
        for known_name in known:
            known_key = _name_key(known_name)
            if not known_key:
                continue
            if proposed_key == known_key:
                return "exact_normalized_name", 1.0
            shorter, longer = sorted((proposed_key, known_key), key=len)
            if len(shorter) >= 8 and shorter in longer:
                score = len(shorter) / len(longer)
                if score > best_score:
                    best_method, best_score = "normalized_name_containment", score
            score = SequenceMatcher(None, proposed_key, known_key).ratio()
            if score > best_score:
                best_method, best_score = "normalized_name_similarity", score
    return best_method, round(best_score, 4)


def _structure_element_contradiction(
    proposed_names: Sequence[str],
    formula: Optional[str],
) -> Optional[str]:
    """Detect an explicit metal in the proposed identity missing from formula."""
    if not formula:
        return None
    combined = " | ".join(proposed_names)
    formula_elements = set(re.findall(r"[A-Z][a-z]?", formula))
    for element, pattern in METAL_NAME_PATTERNS.items():
        if pattern.search(combined) and element not in formula_elements:
            return element
    return None


def _group_names(rows: Sequence[_InputRow]) -> Tuple[str, ...]:
    values = []
    for row in rows:
        values.extend((row.normalized_identity, row.name))
        values.extend(_split_aliases(str(row.raw.get("aliases") or "")))
    return _distinct_names(values)


def _normalized_identity_keys(rows: Sequence[_InputRow]) -> set[str]:
    return {
        key
        for row in rows
        if (key := _name_key(row.normalized_identity))
    }


def _addition_aliases(
    *,
    canonical_name: str,
    rows: Sequence[_InputRow],
    lookup: CompoundLookupResult,
    registry: ConditionRegistry,
) -> Tuple[CompoundAliasInput, ...]:
    candidates = []
    for row in rows:
        candidates.append(("common_name", row.normalized_identity))
        candidates.append(("legacy_name", row.name))
        candidates.extend(
            ("legacy_name", alias)
            for alias in _split_aliases(str(row.raw.get("aliases") or ""))
        )
    candidates.append(("systematic_name", lookup.iupac_name or ""))
    candidates.extend(("common_name", value) for value in lookup.synonyms)
    aliases = []
    seen = {_registry_name_key(canonical_name)}
    for identifier_type, value in candidates:
        value = " ".join(str(value or "").split())
        key = _registry_name_key(value)
        if not value or not key or key in seen:
            continue
        seen.add(key)
        if registry.resolve(name=value).status != "unresolved":
            continue
        aliases.append(
            CompoundAliasInput(
                identifier_type=identifier_type,
                value=value,
                language="en",
            )
        )
    return tuple(aliases)


def _plan_group(
    cas: str,
    rows: Tuple[_InputRow, ...],
    *,
    registry: ConditionRegistry,
    lookup: CompoundLookupResult,
) -> _GroupPlan:
    cas_resolution = registry.resolve(cas=cas)
    substance = cas_resolution.substance
    registry_names = (
        tuple(identifier.value for identifier in substance.identifiers)
        if substance is not None
        else ()
    )
    web_names = _distinct_names(
        (
            lookup.canonical_name or "",
            lookup.iupac_name or "",
            *lookup.synonyms,
        )
    )
    proposed_names = _group_names(rows)
    method, score = _identity_match(
        proposed_names,
        (*registry_names, *web_names),
    )
    shared_identity_keys = _normalized_identity_keys(rows)
    if len(rows) > 1 and len(shared_identity_keys) > 1:
        return _GroupPlan(
            cas, rows, "review", "SHARED_CAS_DIFFERENT_IDENTITIES",
            cas_resolution.status, substance.substance_id if substance else None,
            substance.canonical_name if substance else None, lookup, method, score,
        )
    unsafe = sorted({row.match_status for row in rows} & UNSAFE_MATCH_STATUSES)
    if unsafe:
        return _GroupPlan(
            cas, rows, "review", "UNSAFE_MATCH_STATUS:" + "|".join(unsafe),
            cas_resolution.status, substance.substance_id if substance else None,
            substance.canonical_name if substance else None, lookup, method, score,
        )
    confidence = "high" if any(row.confidence == "high" for row in rows) else "medium"
    if cas_resolution.status == "ambiguous":
        return _GroupPlan(
            cas, rows, "review", "REGISTRY_CAS_AMBIGUOUS",
            cas_resolution.status, None, None, lookup, method, score,
        )
    threshold = 0.72 if confidence == "high" else 0.82
    has_curator_source = any(
        str(row.raw.get("source_basis") or "").strip()
        and str(row.raw.get("source_url") or "").strip()
        for row in rows
    )
    if substance is not None and (
        score >= threshold or (confidence == "high" and has_curator_source)
    ):
        return _GroupPlan(
            cas=cas,
            rows=rows,
            decision="add_aliases",
            reason="CAS_RESOLVES_EXISTING_IDENTITY_AND_NAME_CORROBORATED",
            registry_status=cas_resolution.status,
            registry_substance_id=substance.substance_id,
            registry_canonical_name=substance.canonical_name,
            lookup=lookup,
            match_method=method,
            match_score=score,
            aliases=proposed_names,
        )
    if score < threshold:
        return _GroupPlan(
            cas, rows, "review", f"IDENTITY_NOT_CORROBORATED:{score:.4f}<{threshold:.2f}",
            cas_resolution.status, substance.substance_id if substance else None,
            substance.canonical_name if substance else None, lookup, method, score,
        )
    if not lookup.found or not lookup.smiles:
        return _GroupPlan(
            cas, rows, "review", "WEB_LOOKUP_LACKS_NAME_OR_STRUCTURE",
            cas_resolution.status, None, None, lookup, method, score,
        )
    contradicted_element = _structure_element_contradiction(
        proposed_names,
        lookup.formula,
    )
    if contradicted_element:
        return _GroupPlan(
            cas,
            rows,
            "review",
            f"STRUCTURE_ELEMENT_CONTRADICTION:EXPECTED_{contradicted_element}",
            cas_resolution.status,
            None,
            None,
            lookup,
            method,
            score,
        )
    canonical_name = (
        lookup.canonical_name
        or next(
            (
                row.normalized_identity
                for row in rows
                if row.normalized_identity
            ),
            "",
        )
        or lookup.iupac_name
        or rows[0].name
    )
    aliases = _addition_aliases(
        canonical_name=canonical_name,
        rows=rows,
        lookup=lookup,
        registry=registry,
    )
    source_ids = "+".join(lookup.source_ids) or "web_lookup"
    roles = sorted(
        {
            role.strip()
            for row in rows
            for role in str(row.raw.get("roles") or "").split("|")
            if role.strip()
        }
    )
    evidence = "; ".join(
        dict.fromkeys(
            str(row.raw.get("source_basis") or "").strip()
            for row in rows
            if str(row.raw.get("source_basis") or "").strip()
        )
    )
    notes = (
        "Imported from curator-updated weak-label unresolved-compounds CSV. "
        f"Observed source roles (contextual only): {' | '.join(roles) or 'none'}. "
        f"CAS evidence: {evidence or 'not supplied'}. "
        "Intrinsic registry roles were intentionally not inferred from dataset columns."
    )
    request = CompoundAdditionRequest(
        canonical_name=canonical_name,
        cas=cas,
        source=f"{ALIAS_SOURCE}+{source_ids}",
        smiles=lookup.smiles,
        substance_kind=lookup.substance_kind,
        aliases=aliases,
        curator_notes=notes,
    )
    return _GroupPlan(
        cas=cas,
        rows=rows,
        decision="add_compound",
        reason="UNREGISTERED_CAS_WITH_CORROBORATED_WEB_IDENTITY_AND_STRUCTURE",
        registry_status=cas_resolution.status,
        registry_substance_id=f"cas:{cas}",
        registry_canonical_name=canonical_name,
        lookup=lookup,
        match_method=method,
        match_score=score,
        aliases=tuple(alias.value for alias in aliases),
        addition_request=request,
    )


def _alias_requests(
    plans: Sequence[_GroupPlan],
    registry: ConditionRegistry,
) -> Tuple[SubstanceAliasAdditionRequest, ...]:
    requests = []
    seen = set()
    for plan in plans:
        if plan.decision != "add_aliases" or not plan.registry_substance_id:
            continue
        normalized_identities = {
            _registry_name_key(row.normalized_identity): row.normalized_identity
            for row in plan.rows
            if _registry_name_key(row.normalized_identity)
        }
        for value in plan.aliases:
            identifier_type = (
                "common_name"
                if _registry_name_key(value) in normalized_identities
                else "legacy_name"
            )
            key = (plan.registry_substance_id, _registry_name_key(value))
            if not key[1] or key in seen:
                continue
            seen.add(key)
            existing = registry.resolve(name=value)
            if (
                existing.status == "resolved"
                and existing.substance is not None
                and existing.substance.substance_id == plan.registry_substance_id
            ):
                continue
            if existing.status != "unresolved":
                continue
            requests.append(
                SubstanceAliasAdditionRequest(
                    substance_id=plan.registry_substance_id,
                    identifier_type=identifier_type,
                    value=value,
                    source=f"{ALIAS_SOURCE}:{plan.cas}",
                    language="en",
                )
            )
    return tuple(requests)


def _copy_definitions(
    directory: Path,
    paths: Sequence[Path],
) -> Tuple[Path, Path, Path]:
    directory.mkdir(parents=True, exist_ok=True)
    copied = []
    for path in paths:
        target = directory / path.name
        shutil.copy2(path, target)
        copied.append(target)
    return tuple(copied)  # type: ignore[return-value]


def _apply_requests(
    additions: Sequence[CompoundAdditionRequest],
    aliases: Sequence[SubstanceAliasAdditionRequest],
    *,
    substances_path: Path,
    additions_path: Path,
    identifiers_path: Path,
) -> Tuple[int, int]:
    compound_count = 0
    alias_count = 0
    for request in additions:
        add_compound(
            request,
            substances_path=substances_path,
            additions_path=additions_path,
            identifiers_path=identifiers_path,
        )
        compound_count += 1
    if aliases:
        result = add_substance_aliases(
            aliases,
            substances_path=substances_path,
            additions_path=additions_path,
            identifiers_path=identifiers_path,
        )
        alias_count = len(result.added)
    return compound_count, alias_count


def _write_csv(
    path: Path,
    fieldnames: Sequence[str],
    rows: Iterable[Mapping[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _audit_row(plan: _GroupPlan, row: _InputRow) -> Dict[str, object]:
    return {
        "name": row.name,
        "aliases": row.raw.get("aliases") or "",
        "roles": row.raw.get("roles") or "",
        "mention_count": row.raw.get("mention_count") or "",
        "normalized_identity": row.normalized_identity,
        "possible_cas_no": row.raw.get("possible_cas_no") or "",
        "cas_match_status": row.match_status,
        "cas_confidence": row.confidence,
        "cas_checksum_status": "valid",
        "registry_status_before": plan.registry_status,
        "registry_substance_id": plan.registry_substance_id or "",
        "registry_canonical_name": plan.registry_canonical_name or "",
        "lookup_status": plan.lookup.status,
        "lookup_canonical_name": plan.lookup.canonical_name or "",
        "identity_match_method": plan.match_method,
        "identity_match_score": f"{plan.match_score:.4f}",
        "decision": plan.decision,
        "decision_reason": plan.reason,
        "aliases_added_or_planned": " | ".join(plan.aliases),
        "source_basis": row.raw.get("source_basis") or "",
        "source_url": row.raw.get("source_url") or "",
        "match_notes": row.raw.get("match_notes") or "",
        "lookup_sources": " | ".join(plan.lookup.source_urls),
        "lookup_warnings": " | ".join(plan.lookup.warnings),
    }


def reconcile_registry_from_cas_csv(
    input_path: str | Path,
    output_dir: str | Path,
    *,
    apply_changes: bool = False,
    lookup: LookupFunction = lookup_compound_by_cas,
    lookup_workers: int = 4,
    substances_path: str | Path = SUBSTANCES_PATH,
    additions_path: str | Path = ADDITIONS_PATH,
    identifiers_path: str | Path = IDENTIFIERS_PATH,
) -> CasReconciliationSummary:
    """Plan, verify, optionally apply, and audit curator-supplied CAS mappings."""
    input_path = Path(input_path)
    output_dir = Path(output_dir)
    substances_path = Path(substances_path)
    additions_path = Path(additions_path)
    identifiers_path = Path(identifiers_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = _parse_rows(input_path)
    grouped: Dict[str, list[_InputRow]] = {}
    for row in rows:
        if len(row.cas_values) == 1:
            grouped.setdefault(row.cas_values[0], []).append(row)
    cache_path = output_dir / "cas_lookup_cache.json"
    lookups = _resolve_lookups(
        grouped,
        cache_path=cache_path,
        lookup=lookup,
        workers=lookup_workers,
    )
    registry = ConditionRegistry(
        substances_path=substances_path,
        additions_path=additions_path,
        identifiers_path=identifiers_path,
    )
    plans = tuple(
        _plan_group(
            cas,
            tuple(group_rows),
            registry=registry,
            lookup=lookups[cas],
        )
        for cas, group_rows in sorted(grouped.items())
    )
    plan_by_cas = {plan.cas: plan for plan in plans}
    audit_rows = []
    for row in rows:
        if len(row.cas_values) == 1:
            audit_rows.append(_audit_row(plan_by_cas[row.cas_values[0]], row))
            continue
        raw_cas = str(row.raw.get("possible_cas_no") or "").strip()
        if not raw_cas:
            reason = "NO_CAS_SUPPLIED"
            checksum = "not_applicable"
        elif len(row.cas_values) > 1:
            reason = "MULTIPLE_VALID_CAS_CANDIDATES"
            checksum = "valid_multiple"
        else:
            reason = "NO_VALID_CAS_CHECKSUM"
            checksum = "invalid"
        empty_lookup = CompoundLookupResult(cas=raw_cas, status="not_run")
        placeholder = _GroupPlan(
            raw_cas, (row,), "review", reason, "not_checked", None, None,
            empty_lookup, "not_checked", 0.0,
        )
        audit = _audit_row(placeholder, row)
        audit["cas_checksum_status"] = checksum
        audit_rows.append(audit)

    alias_requests = _alias_requests(plans, registry)
    addition_requests = tuple(
        plan.addition_request
        for plan in plans
        if plan.addition_request is not None
    )
    compounds_added = 0
    aliases_added = 0
    if apply_changes and (addition_requests or alias_requests):
        paths = (substances_path, additions_path, identifiers_path)
        with tempfile.TemporaryDirectory(prefix="condition-registry-cas-dry-run-") as temp:
            dry_paths = _copy_definitions(Path(temp), paths)
            _apply_requests(
                addition_requests,
                alias_requests,
                substances_path=dry_paths[0],
                additions_path=dry_paths[1],
                identifiers_path=dry_paths[2],
            )
        originals = {path: path.read_bytes() for path in paths}
        try:
            compounds_added, aliases_added = _apply_requests(
                addition_requests,
                alias_requests,
                substances_path=substances_path,
                additions_path=additions_path,
                identifiers_path=identifiers_path,
            )
        except Exception:
            for path, content in originals.items():
                with tempfile.NamedTemporaryFile(
                    mode="wb", dir=path.parent, delete=False
                ) as handle:
                    handle.write(content)
                    temporary = handle.name
                os.replace(temporary, path)
            raise

    audit_path = output_dir / "weak_label_cas_registry_audit.csv"
    review_path = output_dir / "weak_label_cas_still_unresolved.csv"
    _write_csv(audit_path, AUDIT_FIELDNAMES, audit_rows)
    _write_csv(
        review_path,
        AUDIT_FIELDNAMES,
        (row for row in audit_rows if row["decision"] == "review"),
    )
    updated_existing = len(
        {
            request.substance_id
            for request in alias_requests
        }
    )
    return CasReconciliationSummary(
        input_rows=len(rows),
        rows_with_one_valid_cas=sum(len(value) for value in grouped.values()),
        existing_substances_updated=updated_existing if apply_changes else 0,
        aliases_added=aliases_added,
        new_compounds_added=compounds_added,
        review_rows=sum(row["decision"] == "review" for row in audit_rows),
        audit_path=str(audit_path),
        review_path=str(review_path),
        lookup_cache_path=str(cache_path),
    )


__all__ = [
    "AUDIT_FIELDNAMES",
    "CasReconciliationSummary",
    "reconcile_registry_from_cas_csv",
]

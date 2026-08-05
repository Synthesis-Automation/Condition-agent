"""Reconcile weak-label condition names with registry substance identities."""

from __future__ import annotations

import argparse
import csv
import re
import unicodedata
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional, Sequence, Tuple

from .curation import SubstanceAliasAdditionRequest, add_substance_aliases
from .loader import ADDITIONS_PATH, IDENTIFIERS_PATH, SUBSTANCES_PATH, load_substances
from .models import CONDITION_NAME_IDENTIFIER_TYPES, Substance
from .normalization import normalize_chemical_name
from .resolver import ConditionRegistry

WEAK_LABEL_ROLE_COLUMNS = {
    "Base": "base",
    "Catalyst": "catalyst",
    "Solvent": "solvent",
    "Ligand": "ligand",
    "Additive": "additive",
    "Coupling Reagent": "coupling_reagent",
    "Secondary Solvent": "solvent",
    "Tertiary Solvent": "solvent",
}
MISSING_COMPONENT_LABELS = {"", "missing", "none", "nan", "n/a"}
ALIAS_SOURCE = "weak_label_v2.1:condition_name_reconciliation"


def _curated_alias_targets() -> Dict[str, Tuple[str, str]]:
    mappings = {
        "tAmOH": ("cas:75-85-4", "established abbreviation: tert-amyl alcohol"),
        "NEt3": ("cas:121-44-8", "established abbreviation: triethylamine"),
        "NMM": ("cas:109-02-4", "established abbreviation: N-methylmorpholine"),
        "TBME": ("cas:1634-04-4", "established abbreviation: MTBE"),
        "TBAC": ("cas:1112-67-0", "established abbreviation: tetrabutylammonium chloride"),
        "TMG": ("cas:80-70-6", "established abbreviation: tetramethylguanidine"),
        "LiOtBu": ("cas:1907-33-1", "established abbreviation: lithium tert-butoxide"),
        "nBuOAc": ("cas:123-86-4", "established abbreviation: n-butyl acetate"),
        "2-BuOH": ("cas:78-92-2", "established abbreviation: 2-butanol"),
        "trans-NN-Dimethylcyclohexane-12-diamine": (
            "cas:67579-81-1",
            "punctuation-normalized chemical name",
        ),
        "AgF": ("cas:7775-41-9", "molecular formula abbreviation: silver fluoride"),
        "CuF2": ("cas:7789-19-7", "molecular formula abbreviation: copper(II) fluoride"),
        "AdCOOH": ("cas:828-51-3", "established abbreviation: 1-adamantanecarboxylic acid"),
        "NiBr2-glyme": (
            "cas:28923-39-9",
            "established abbreviation: nickel(II) bromide DME complex",
        ),
        "dtbpy": ("cas:72914-19-3", "established ligand abbreviation: dtbbpy"),
        "Me4phen": ("cas:1660-93-1", "established ligand abbreviation: tmphen"),
        "8-Quinolol": ("cas:148-24-3", "established name: 8-hydroxyquinoline"),
        "LDA": ("cas:4111-54-0", "established abbreviation: lithium diisopropylamide"),
        "TrichloroACN": ("cas:545-06-2", "established abbreviation: trichloroacetonitrile"),
        "NaHMDS": ("cas:1070-89-9", "established abbreviation: sodium HMDS"),
        "NaOtPent": ("cas:14593-46-5", "established abbreviation: sodium tert-amylate"),
        "K2S2O5": ("cas:16731-55-8", "molecular formula: potassium metabisulfite"),
    }
    return {
        normalize_chemical_name(alias): target
        for alias, target in mappings.items()
    }


CURATED_ALIAS_TARGETS = _curated_alias_targets()


@dataclass(frozen=True)
class ExtractedWeakLabelReagent:
    """One normalized reagent with observed spellings, roles, and counts."""

    normalized_name: str
    preferred_name: str
    aliases: Tuple[str, ...]
    role_counts: Mapping[str, int]
    source_columns: Tuple[str, ...]

    @property
    def mention_count(self) -> int:
        return sum(self.role_counts.values())


@dataclass(frozen=True)
class ReconciliationMatch:
    """Registry reconciliation decision for one weak-label reagent."""

    reagent: ExtractedWeakLabelReagent
    status: str
    match_method: str
    substance_id: Optional[str] = None
    canonical_name: Optional[str] = None
    evidence: Optional[str] = None
    candidates: Tuple[str, ...] = ()
    alias_added: bool = False


def _display_text(value: str) -> str:
    return " ".join(unicodedata.normalize("NFKC", str(value or "")).split())


def _split_component_value(value: str) -> Tuple[str, ...]:
    """Split source list delimiters without splitting systematic-name commas."""
    return tuple(
        component
        for raw_component in re.split(r",\s+", _display_text(value))
        if (component := _display_text(raw_component))
        and component.casefold() not in MISSING_COMPONENT_LABELS
        and not component.casefold().startswith("no ")
    )


def extract_weak_label_reagents(
    source_path: str | Path,
) -> Tuple[ExtractedWeakLabelReagent, ...]:
    """Extract and deduplicate all individual condition reagents and roles."""
    variants: Dict[str, Counter[str]] = defaultdict(Counter)
    role_counts: Dict[str, Counter[str]] = defaultdict(Counter)
    source_columns: Dict[str, set[str]] = defaultdict(set)
    with Path(source_path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        missing_columns = set(WEAK_LABEL_ROLE_COLUMNS) - set(reader.fieldnames or ())
        if missing_columns:
            raise ValueError(
                "Weak-label source is missing role columns: "
                + ", ".join(sorted(missing_columns))
            )
        for row in reader:
            row_mentions: set[tuple[str, str]] = set()
            for column, role in WEAK_LABEL_ROLE_COLUMNS.items():
                for component in _split_component_value(row.get(column) or ""):
                    normalized = normalize_chemical_name(component)
                    if not normalized:
                        continue
                    variants[normalized][component] += 1
                    source_columns[normalized].add(column)
                    row_mentions.add((normalized, role))
            for normalized, role in row_mentions:
                role_counts[normalized][role] += 1
    reagents = []
    for normalized, observed_variants in variants.items():
        ordered_variants = tuple(
            value
            for value, _count in sorted(
                observed_variants.items(),
                key=lambda item: (-item[1], len(item[0]), item[0].casefold(), item[0]),
            )
        )
        reagents.append(
            ExtractedWeakLabelReagent(
                normalized_name=normalized,
                preferred_name=ordered_variants[0],
                aliases=ordered_variants,
                role_counts=dict(sorted(role_counts[normalized].items())),
                source_columns=tuple(sorted(source_columns[normalized])),
            )
        )
    return tuple(sorted(reagents, key=lambda item: item.normalized_name))


def _compact_name(value: str) -> str:
    return normalize_chemical_name(value).replace(" ", "")


def _compact_registry_index(
    substances: Iterable[Substance],
) -> Dict[str, Dict[str, set[str]]]:
    index: Dict[str, Dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for substance in substances:
        for identifier in substance.identifiers:
            if (
                identifier.status == "active"
                and identifier.identifier_type in CONDITION_NAME_IDENTIFIER_TYPES
            ):
                compact = _compact_name(identifier.value)
                if compact:
                    index[compact][substance.substance_id].add(identifier.value)
    return index


def reconcile_weak_label_reagents(
    reagents: Iterable[ExtractedWeakLabelReagent],
    *,
    registry: ConditionRegistry,
    substances: Iterable[Substance],
) -> Tuple[ReconciliationMatch, ...]:
    """Reconcile exact, formatting-only, and curated high-confidence aliases."""
    compact_index = _compact_registry_index(substances)
    matches = []
    for reagent in reagents:
        resolved_ids = set()
        ambiguous_candidates = set()
        for alias in reagent.aliases:
            result = registry.resolve(name=alias)
            if result.status == "resolved" and result.substance is not None:
                resolved_ids.add(result.substance.substance_id)
            elif result.status == "ambiguous":
                ambiguous_candidates.update(result.candidates)
        if len(resolved_ids) == 1 and not ambiguous_candidates:
            substance_id = next(iter(resolved_ids))
            substance = registry.resolve_id(substance_id).substance
            matches.append(
                ReconciliationMatch(
                    reagent=reagent,
                    status="matched",
                    match_method="exact_registry_name_or_alias",
                    substance_id=substance_id,
                    canonical_name=(substance.canonical_name if substance else None),
                    evidence="registry exact normalized identifier",
                )
            )
            continue
        if len(resolved_ids) > 1 or ambiguous_candidates:
            matches.append(
                ReconciliationMatch(
                    reagent=reagent,
                    status="ambiguous",
                    match_method="ambiguous_registry_identifier",
                    candidates=tuple(sorted(resolved_ids | ambiguous_candidates)),
                )
            )
            continue
        compact_owners = compact_index.get(_compact_name(reagent.preferred_name), {})
        if len(compact_owners) == 1:
            substance_id = next(iter(compact_owners))
            substance = registry.resolve_id(substance_id).substance
            matched_values = sorted(compact_owners[substance_id])
            matches.append(
                ReconciliationMatch(
                    reagent=reagent,
                    status="alias_candidate",
                    match_method="unique_compact_identifier",
                    substance_id=substance_id,
                    canonical_name=(substance.canonical_name if substance else None),
                    evidence="formatting-only match to: " + " | ".join(matched_values),
                )
            )
            continue
        curated = CURATED_ALIAS_TARGETS.get(reagent.normalized_name)
        if curated is not None:
            substance_id, evidence = curated
            substance = registry.resolve_id(substance_id).substance
            if substance is not None:
                matches.append(
                    ReconciliationMatch(
                        reagent=reagent,
                        status="alias_candidate",
                        match_method="curated_chemical_alias",
                        substance_id=substance_id,
                        canonical_name=substance.canonical_name,
                        evidence=evidence,
                    )
                )
                continue
        matches.append(
            ReconciliationMatch(
                reagent=reagent,
                status="unresolved",
                match_method="no_defensible_registry_match",
            )
        )
    return tuple(matches)


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


def _role_rows(
    reagents: Iterable[ExtractedWeakLabelReagent],
) -> Iterable[Mapping[str, object]]:
    for reagent in reagents:
        for role, count in sorted(reagent.role_counts.items()):
            yield {
                "reagent_name": reagent.preferred_name,
                "observed_aliases": " | ".join(reagent.aliases),
                "role": role,
                "mention_count": count,
                "source_columns": " | ".join(reagent.source_columns),
                "normalized_name": reagent.normalized_name,
            }


def _reconciliation_rows(
    matches: Iterable[ReconciliationMatch],
) -> Iterable[Mapping[str, object]]:
    for match in matches:
        yield {
            "reagent_name": match.reagent.preferred_name,
            "observed_aliases": " | ".join(match.reagent.aliases),
            "observed_roles": " | ".join(match.reagent.role_counts),
            "mention_count": match.reagent.mention_count,
            "registry_status": match.status,
            "match_method": match.match_method,
            "substance_id": match.substance_id or "",
            "registry_canonical_name": match.canonical_name or "",
            "alias_added": str(match.alias_added).lower(),
            "evidence": match.evidence or "",
            "candidate_substance_ids": " | ".join(match.candidates),
        }


def _unresolved_rows(
    matches: Iterable[ReconciliationMatch],
) -> Iterable[Mapping[str, object]]:
    for match in matches:
        if match.status != "unresolved":
            continue
        yield {
            "name": match.reagent.preferred_name,
            "aliases": " | ".join(match.reagent.aliases),
            "roles": " | ".join(match.reagent.role_counts),
            "mention_count": match.reagent.mention_count,
            "source_columns": " | ".join(match.reagent.source_columns),
            "reconciliation_status": match.status,
        }


def run_weak_label_reconciliation(
    source_path: str | Path,
    output_dir: str | Path,
    *,
    apply_aliases: bool = False,
    substances_path: str | Path = SUBSTANCES_PATH,
    additions_path: str | Path = ADDITIONS_PATH,
    identifiers_path: str | Path = IDENTIFIERS_PATH,
) -> Mapping[str, int]:
    """Extract, reconcile, optionally apply aliases, and write review CSVs."""
    source_path = Path(source_path)
    output_dir = Path(output_dir)
    substances_path = Path(substances_path)
    additions_path = Path(additions_path)
    identifiers_path = Path(identifiers_path)
    reagents = extract_weak_label_reagents(source_path)
    substances = load_substances(
        substances_path=substances_path,
        additions_path=additions_path,
        identifiers_path=identifiers_path,
    )
    registry = ConditionRegistry(substances=substances)
    matches = reconcile_weak_label_reagents(
        reagents,
        registry=registry,
        substances=substances,
    )
    alias_candidates = tuple(
        match for match in matches if match.status == "alias_candidate"
    )
    added_values: set[tuple[str, str]] = set()
    if apply_aliases and alias_candidates:
        result = add_substance_aliases(
            (
                SubstanceAliasAdditionRequest(
                    substance_id=match.substance_id or "",
                    identifier_type="legacy_name",
                    value=match.reagent.preferred_name,
                    source=ALIAS_SOURCE,
                )
                for match in alias_candidates
            ),
            substances_path=substances_path,
            additions_path=additions_path,
            identifiers_path=identifiers_path,
        )
        added_values = {
            (identifier.substance_id, identifier.value)
            for identifier in result.added
        }
        matches = tuple(
            ReconciliationMatch(
                reagent=match.reagent,
                status=("matched" if match.status == "alias_candidate" else match.status),
                match_method=match.match_method,
                substance_id=match.substance_id,
                canonical_name=match.canonical_name,
                evidence=match.evidence,
                candidates=match.candidates,
                alias_added=(
                    match.substance_id or "",
                    match.reagent.preferred_name,
                )
                in added_values,
            )
            for match in matches
        )

    role_path = output_dir / "weak_label_reagent_roles_deduplicated.csv"
    reconciliation_path = output_dir / "weak_label_registry_reconciliation.csv"
    unresolved_path = output_dir / "weak_label_unresolved_reagents.csv"
    summary_path = output_dir / "weak_label_reconciliation_summary.csv"
    _write_csv(
        role_path,
        (
            "reagent_name",
            "observed_aliases",
            "role",
            "mention_count",
            "source_columns",
            "normalized_name",
        ),
        _role_rows(reagents),
    )
    _write_csv(
        reconciliation_path,
        (
            "reagent_name",
            "observed_aliases",
            "observed_roles",
            "mention_count",
            "registry_status",
            "match_method",
            "substance_id",
            "registry_canonical_name",
            "alias_added",
            "evidence",
            "candidate_substance_ids",
        ),
        _reconciliation_rows(matches),
    )
    _write_csv(
        unresolved_path,
        (
            "name",
            "aliases",
            "roles",
            "mention_count",
            "source_columns",
            "reconciliation_status",
        ),
        _unresolved_rows(matches),
    )
    counts = Counter(match.status for match in matches)
    summary = {
        "source_rows": sum(1 for _ in source_path.open(encoding="utf-8-sig")) - 1,
        "unique_reagents": len(reagents),
        "unique_reagent_role_rows": sum(len(item.role_counts) for item in reagents),
        "matched_reagents": counts["matched"],
        "ambiguous_reagents": counts["ambiguous"],
        "unresolved_reagents": counts["unresolved"],
        "aliases_added": len(added_values),
    }
    _write_csv(
        summary_path,
        ("metric", "value"),
        ({"metric": key, "value": value} for key, value in summary.items()),
    )
    return summary


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the weak-label registry reconciliation command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_path", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--apply-aliases", action="store_true")
    arguments = parser.parse_args(argv)
    summary = run_weak_label_reconciliation(
        arguments.source_path,
        arguments.output_dir,
        apply_aliases=arguments.apply_aliases,
    )
    for key, value in summary.items():
        print(f"{key}: {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "ALIAS_SOURCE",
    "CURATED_ALIAS_TARGETS",
    "ExtractedWeakLabelReagent",
    "ReconciliationMatch",
    "extract_weak_label_reagents",
    "reconcile_weak_label_reagents",
    "run_weak_label_reconciliation",
]

"""Conservative exact-evidence substance and identifier resolver."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Optional, Tuple

from .loader import (
    ADDITIONS_PATH,
    IDENTIFIERS_PATH,
    SUBSTANCES_PATH,
    load_substances,
)
from .models import (
    CONDITION_IDENTIFIER_TYPES,
    CONDITION_NAME_IDENTIFIER_TYPES,
    ResolutionResult,
    Substance,
    SubstanceIdentifier,
)
from .normalization import normalize_cas, normalize_identifier

_CAS_SHAPE_RE = re.compile(r"\d{2,7}-\d{2}-\d")


class ConditionRegistry:
    """Immutable in-memory index over versioned condition definitions."""

    def __init__(
        self,
        *,
        substances_path: str | Path = SUBSTANCES_PATH,
        additions_path: str | Path = ADDITIONS_PATH,
        identifiers_path: str | Path = IDENTIFIERS_PATH,
        substances: Optional[Iterable[Substance]] = None,
    ) -> None:
        loaded = tuple(substances) if substances is not None else load_substances(
            substances_path=substances_path,
            additions_path=additions_path,
            identifiers_path=identifiers_path,
        )
        self._by_id: Dict[str, Substance] = {}
        self._by_cas: DefaultDict[str, List[Substance]] = defaultdict(list)
        self._by_name: DefaultDict[str, List[SubstanceIdentifier]] = defaultdict(list)
        self._by_identifier_type: Dict[
            str, DefaultDict[str, List[SubstanceIdentifier]]
        ] = {
            identifier_type: defaultdict(list)
            for identifier_type in CONDITION_IDENTIFIER_TYPES
        }
        self._identifier_substances: Dict[str, Substance] = {}
        identifier_ids: set[str] = set()
        for substance in loaded:
            if substance.substance_id in self._by_id:
                raise ValueError(f"Duplicate substance ID: {substance.substance_id}")
            self._by_id[substance.substance_id] = substance
            if substance.cas:
                self._by_cas[substance.cas].append(substance)
            for identifier in substance.identifiers:
                if identifier.identifier_id in identifier_ids:
                    raise ValueError(
                        f"Duplicate substance identifier ID: {identifier.identifier_id}"
                    )
                identifier_ids.add(identifier.identifier_id)
                if identifier.status != "active":
                    continue
                normalized = normalize_identifier(
                    identifier.value, identifier.identifier_type
                )
                if not normalized:
                    continue
                self._by_identifier_type[identifier.identifier_type][normalized].append(
                    identifier
                )
                self._identifier_substances[identifier.identifier_id] = substance
                if identifier.identifier_type in CONDITION_NAME_IDENTIFIER_TYPES:
                    self._by_name[normalized].append(identifier)

    @staticmethod
    def _identifier_sort_key(
        identifier: SubstanceIdentifier,
    ) -> Tuple[int, int, str]:
        type_priority = CONDITION_IDENTIFIER_TYPES.index(identifier.identifier_type)
        return (-int(identifier.is_preferred), type_priority, identifier.identifier_id)

    def _result_from_identifiers(
        self,
        *,
        query: str,
        identifiers: Iterable[SubstanceIdentifier],
        match_kind: Optional[str] = None,
    ) -> ResolutionResult:
        matched = tuple(sorted(identifiers, key=self._identifier_sort_key))
        substances = {
            self._identifier_substances[item.identifier_id].substance_id:
            self._identifier_substances[item.identifier_id]
            for item in matched
        }
        if len(substances) == 1:
            identifier = matched[0]
            return ResolutionResult(
                query=query,
                status="resolved",
                substance=next(iter(substances.values())),
                match_kind=match_kind or f"exact_{identifier.identifier_type}",
                matched_identifier=identifier,
            )
        if len(substances) > 1:
            return ResolutionResult(
                query=query,
                status="ambiguous",
                match_kind=match_kind,
                candidates=tuple(sorted(substances)),
                warnings=("AMBIGUOUS_IDENTIFIER",),
            )
        return ResolutionResult(query=query, status="unresolved")

    def resolve_identifier(
        self,
        value: str,
        *,
        identifier_type: str = "auto",
    ) -> ResolutionResult:
        """Resolve an exact typed identifier while preserving ambiguity."""
        query = str(value or "").strip()
        if not query:
            return ResolutionResult(query=query, status="empty_query")
        if identifier_type == "substance_id":
            return self.resolve_id(query)
        if identifier_type == "name":
            normalized = normalize_identifier(query, "common_name")
            return self._result_from_identifiers(
                query=query,
                identifiers=self._by_name.get(normalized or "", ()),
                match_kind="exact_alias",
            )
        if identifier_type == "auto":
            if _CAS_SHAPE_RE.fullmatch(query):
                return self.resolve(cas=query)
            normalized = normalize_identifier(query, "common_name")
            return self._result_from_identifiers(
                query=query,
                identifiers=self._by_name.get(normalized or "", ()),
                match_kind="exact_alias",
            )
        if identifier_type not in CONDITION_IDENTIFIER_TYPES:
            raise ValueError(f"Unsupported condition identifier type: {identifier_type}")
        normalized = normalize_identifier(query, identifier_type)
        if identifier_type == "cas" and normalized is None:
            return ResolutionResult(
                query=query,
                status="invalid_identifier",
                warnings=("INVALID_CAS",),
            )
        return self._result_from_identifiers(
            query=query,
            identifiers=self._by_identifier_type[identifier_type].get(
                normalized or "", ()
            ),
        )

    def resolve(
        self,
        *,
        cas: Optional[str] = None,
        name: Optional[str] = None,
    ) -> ResolutionResult:
        if cas:
            normalized_cas = normalize_cas(cas)
            if normalized_cas is None:
                return ResolutionResult(
                    query=cas,
                    status="invalid_identifier",
                    warnings=("INVALID_CAS",),
                )
            substances = {
                substance.substance_id: substance
                for substance in self._by_cas.get(normalized_cas, ())
            }
            if len(substances) > 1:
                return ResolutionResult(
                    query=cas,
                    status="ambiguous",
                    match_kind="exact_cas",
                    candidates=tuple(sorted(substances)),
                    warnings=("AMBIGUOUS_IDENTIFIER",),
                )
            identifiers = self._by_identifier_type["cas"].get(
                normalized_cas, ()
            )
            return self._result_from_identifiers(
                query=cas,
                identifiers=identifiers,
                match_kind="exact_cas",
            )
        return self.resolve_identifier(name or "", identifier_type="name")

    def resolve_id(self, substance_id: str) -> ResolutionResult:
        """Resolve an exact stable registry identity without alias inference."""
        query = str(substance_id or "").strip()
        if not query:
            return ResolutionResult(query=query, status="empty_query")
        substance = self._by_id.get(query)
        if substance is None:
            return ResolutionResult(query=query, status="unresolved")
        return ResolutionResult(
            query=query,
            status="resolved",
            substance=substance,
            match_kind="exact_substance_id",
        )

    @property
    def size(self) -> int:
        return len(self._by_id)


__all__ = ["ConditionRegistry"]

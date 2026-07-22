"""Conservative exact-evidence substance resolver."""

from __future__ import annotations

from collections import defaultdict
from typing import DefaultDict, Dict, List, Optional

from .loader import iter_substance_rows, row_to_substance
from .models import ResolutionResult, Substance
from .normalization import normalize_cas, normalize_name


class ConditionRegistry:
    def __init__(self) -> None:
        self._by_id: Dict[str, Substance] = {}
        self._by_cas: Dict[str, Substance] = {}
        self._by_name: DefaultDict[str, List[Substance]] = defaultdict(list)
        for row_number, row in enumerate(iter_substance_rows(), start=2):
            substance = row_to_substance(row, row_number)
            self._by_id[substance.substance_id] = substance
            if substance.cas:
                self._by_cas[substance.cas] = substance
            for value in (substance.canonical_name, *substance.aliases):
                normalized = normalize_name(value)
                if normalized:
                    self._by_name[normalized].append(substance)

    def resolve(self, *, cas: Optional[str] = None, name: Optional[str] = None) -> ResolutionResult:
        if cas:
            normalized_cas = normalize_cas(cas)
            if normalized_cas is None:
                return ResolutionResult(cas, "invalid_identifier", warnings=("INVALID_CAS",))
            substance = self._by_cas.get(normalized_cas)
            return ResolutionResult(cas, "resolved", substance, "exact_cas") if substance else ResolutionResult(cas, "unresolved")
        normalized_name = normalize_name(name or "")
        if not normalized_name:
            return ResolutionResult(name or "", "empty_query")
        candidates = self._by_name.get(normalized_name, [])
        unique = {item.substance_id: item for item in candidates}
        if len(unique) == 1:
            substance = next(iter(unique.values()))
            return ResolutionResult(name or "", "resolved", substance, "exact_alias")
        if len(unique) > 1:
            return ResolutionResult(name or "", "ambiguous", candidates=tuple(sorted(unique)))
        return ResolutionResult(name or "", "unresolved")

    def resolve_id(self, substance_id: str) -> ResolutionResult:
        """Resolve an exact stable registry identity without alias inference."""
        query = str(substance_id or "").strip()
        if not query:
            return ResolutionResult(query, "empty_query")
        substance = self._by_id.get(query)
        if substance is None:
            return ResolutionResult(query, "unresolved")
        return ResolutionResult(query, "resolved", substance, "exact_substance_id")

    @property
    def size(self) -> int:
        return len(self._by_id)


__all__ = ["ConditionRegistry"]

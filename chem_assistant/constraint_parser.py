"""
Constraint parsing utilities for LangChain ChemTools integration.

This module converts loose natural language preferences into structured
constraint specifications that can be applied to ChemTools recommendation
results. The parsing is intentionally deterministic and keyword-driven so it
can run safely inside tool wrappers and the CLI without any LLM involvement.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set

# Canonical metal tokens that commonly appear in catalyst core labels.
_METAL_ALIASES: Dict[str, Set[str]] = {
    "PD": {"pd", "palladium"},
    "CU": {"cu", "copper"},
    "NI": {"ni", "nickel"},
    "FE": {"fe", "iron", "ferric", "ferrous"},
    "CO": {"co", "cobalt"},
    "ZN": {"zn", "zinc"},
    "AG": {"ag", "silver"},
    "AU": {"au", "gold"},
    "RU": {"ru", "ruthenium"},
}

# Environmental / operational keywords for constraint_rules passthrough.
_RULE_KEYWORDS = {
    "no_chlorinated": {
        "no chlorinated",
        "avoid chlorinated",
        "ban chlorinated",
    },
    "no_HMPA": {
        "no hmpa",
        "ban hmpa",
        "avoid hmpa",
    },
    "aqueous_only": {
        "aqueous only",
        "water only",
        "must be aqueous",
        "use water solvent only",
        "solvent must be water",
    },
}

_ALLOW_PATTERNS = [
    r"\bonly\s+(?:use|allow|permit|keep|stick to)\s+{token}",
    r"\b(use|require|keep|stick to)\s+{token}\b",
    r"{token}\s+only\b",
]

_EXCLUDE_PATTERNS = [
    r"\bno\s+{token}\b",
    r"\bavoid\s+{token}\b",
    r"\bwithout\s+{token}\b",
    r"{token}[- ]?free\b",
    r"\bexclude\s+{token}\b",
]

_PREFER_PATTERNS = [
    r"\bprefer\s+{token}\b",
    r"\bfavor\s+{token}\b",
    r"\blean towards\s+{token}\b",
    r"\bprioritize\s+{token}\b",
]

_CROSS_FAMILY_PATTERNS = [
    "cross-family",
    "cross family",
    "search all families",
    "all reaction families",
    "ignore reaction family",
    "search across families",
]


def _normalize_metals(values: Optional[Sequence[str]]) -> Set[str]:
    """Normalize user-provided metal strings into canonical tokens."""
    out: Set[str] = set()
    if not values:
        return out
    for value in values:
        token = (value or "").strip().lower()
        if not token:
            continue
        for canonical, aliases in _METAL_ALIASES.items():
            if token == canonical.lower() or token in aliases:
                out.add(canonical)
                break
        else:
            # If not matched, take uppercase abbreviation.
            out.add(token.upper())
    return out


def _scan_text_patterns(text: str, patterns: Iterable[str], token: str) -> bool:
    """Return True when any regex pattern matches the text for a given token."""
    token_re = re.escape(token)
    for pattern in patterns:
        regex = pattern.format(token=token_re)
        if re.search(regex, text, flags=re.IGNORECASE):
            return True
    return False


@dataclass
class ConstraintSpec:
    """Structured constraints extracted from user input."""

    allow_metals: Set[str] = field(default_factory=set)
    exclude_metals: Set[str] = field(default_factory=set)
    prefer_metals: Set[str] = field(default_factory=set)
    search_all_families: bool = False
    constraint_rules: Dict[str, object] = field(default_factory=dict)

    def merge(self, other: "ConstraintSpec") -> "ConstraintSpec":
        """Merge another spec into this one (in-place) and return self."""
        self.allow_metals |= other.allow_metals
        self.exclude_metals |= other.exclude_metals
        self.prefer_metals |= other.prefer_metals
        self.search_all_families = self.search_all_families or other.search_all_families
        if other.constraint_rules:
            merged = dict(self.constraint_rules)
            merged.update(other.constraint_rules)
            self.constraint_rules = merged
        # Ensure allow/exclude sets are consistent.
        self.exclude_metals -= self.allow_metals
        return self

    def formatted_summary(self) -> str:
        """Return a human-readable summary for logging or CLI display."""
        pieces: List[str] = []
        if self.allow_metals:
            pieces.append(f"allow: {', '.join(sorted(self.allow_metals))}")
        if self.exclude_metals:
            pieces.append(f"avoid: {', '.join(sorted(self.exclude_metals))}")
        if self.prefer_metals:
            pieces.append(f"prefer: {', '.join(sorted(self.prefer_metals))}")
        if self.search_all_families:
            pieces.append("search all families")
        if not pieces and not self.constraint_rules:
            return "none"
        if self.constraint_rules:
            rule_bits = ", ".join(sorted(self.constraint_rules.keys()))
            pieces.append(f"rules: {rule_bits}")
        return "; ".join(pieces)


def parse_constraint_text(text: Optional[str]) -> ConstraintSpec:
    """Parse free-form constraint text into a structured specification."""
    spec = ConstraintSpec()
    if not text:
        return spec

    lowered = " ".join(text.strip().split()).lower()
    if not lowered:
        return spec

    # Search for cross-family requests.
    if any(phrase in lowered for phrase in _CROSS_FAMILY_PATTERNS):
        spec.search_all_families = True

    # Scan for metal preferences.
    for canonical, aliases in _METAL_ALIASES.items():
        for alias in aliases | {canonical.lower()}:
            if _scan_text_patterns(lowered, _EXCLUDE_PATTERNS, alias):
                spec.exclude_metals.add(canonical)
            if _scan_text_patterns(lowered, _ALLOW_PATTERNS, alias):
                spec.allow_metals.add(canonical)
            if _scan_text_patterns(lowered, _PREFER_PATTERNS, alias):
                spec.prefer_metals.add(canonical)

        if re.search(rf"\b{canonical.lower()}[- ]?free\b", lowered):
            spec.exclude_metals.add(canonical)

    # Environmental/operational rules.
    for rule_key, triggers in _RULE_KEYWORDS.items():
        if any(trigger in lowered for trigger in triggers):
            spec.constraint_rules[rule_key] = True

    return spec


def build_constraint_spec(
    *,
    text: Optional[str] = None,
    allow_metals: Optional[Sequence[str]] = None,
    exclude_metals: Optional[Sequence[str]] = None,
    prefer_metals: Optional[Sequence[str]] = None,
    search_all_families: Optional[bool] = None,
    base_constraint_rules: Optional[Dict[str, object]] = None,
) -> ConstraintSpec:
    """Combine explicit parameters and text into a single constraint spec."""
    spec = parse_constraint_text(text)
    spec.allow_metals |= _normalize_metals(allow_metals)
    spec.exclude_metals |= _normalize_metals(exclude_metals)
    spec.prefer_metals |= _normalize_metals(prefer_metals)

    if search_all_families is not None:
        spec.search_all_families = bool(search_all_families) or spec.search_all_families

    if base_constraint_rules:
        merged = dict(spec.constraint_rules)
        merged.update(base_constraint_rules)
        spec.constraint_rules = merged

    spec.exclude_metals -= spec.allow_metals
    return spec


def merge_specs(*specs: ConstraintSpec) -> ConstraintSpec:
    """Merge multiple specs into a single combined spec."""
    combined = ConstraintSpec()
    for spec in specs:
        combined.merge(spec)
    return combined


def format_constraints_for_prompt(spec: ConstraintSpec) -> Optional[str]:
    """Return a prompt-friendly string describing the constraint spec."""
    if not spec:
        return None
    summary = spec.formatted_summary()
    if summary == "none":
        return None
    return (
        "Apply the following user constraints to any recommendations: "
        f"{summary}."
    )


def filter_cores_by_constraints(
    candidates: Sequence[str],
    *,
    allow_metals: Optional[Set[str]] = None,
    exclude_metals: Optional[Set[str]] = None,
    prefer_metals: Optional[Set[str]] = None,
) -> List[str]:
    """
    Reorder catalyst core candidates according to metal constraints.

    Args:
        candidates: Ordered list of catalyst core strings.
        allow_metals: Allowed metal tokens (whitelist). If present, candidates
            must contain at least one allowed metal token.
        exclude_metals: Metals that should be removed from candidates.
        prefer_metals: Metals to prioritize (moved to the front if present).

    Returns:
        Filtered list of candidate cores. If filtering removes every option,
        the original list is returned to avoid empty recommendations.
    """
    if not candidates:
        return []

    def matches_any(core: str, metals: Set[str]) -> bool:
        upper = core.upper()
        return any(metal in upper for metal in metals)

    allow_metals = set(allow_metals or set())
    exclude_metals = set(exclude_metals or set())
    prefer_metals = set(prefer_metals or set())

    filtered = []
    for core in candidates:
        if exclude_metals and matches_any(core, exclude_metals):
            continue
        if allow_metals and not matches_any(core, allow_metals):
            continue
        filtered.append(core)

    if not filtered:
        # No matches; fall back to original candidates to avoid empty output.
        return list(candidates)

    if not prefer_metals:
        return filtered

    preferred = [core for core in filtered if matches_any(core, prefer_metals)]
    others = [core for core in filtered if core not in preferred]
    return preferred + others


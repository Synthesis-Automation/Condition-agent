"""
LLM Output Validator
=====================

MOSAIC-inspired structured output validation for chemistry LLM responses.

Provides ``pass_check()`` — a schema-driven validator that examines the
parsed JSON dict returned by any chemistry LLM call and:

1. Asserts all **required fields** are present and non-empty.
2. Detects **filler / placeholder text** that indicates the model refused to
   commit (e.g. "N/A", "not specified", "TBD", "null", …).
3. Validates **field types** and plausible **value ranges** where known
   (e.g. confidence ∈ [0, 1], temperature in reasonable °C range).
4. Performs lightweight **nested-schema** validation for complex outputs like
   MULTI_SOURCE_SYNTHESIS consensus blocks.

Usage::

    from chemtools.util.llm_output_validator import pass_check

    ok, failures, warns = pass_check(synthesis_dict, "multi_source_synthesis")
    if not ok:
        raise ValueError(f"LLM output failed validation: {failures}")

Design notes:
    - Schemas are declarative dicts — easy to extend.
    - ``pass_check`` never raises; it returns ``(bool, list[str], list[str])``.
    - All string comparisons are case-insensitive and strip-whitespace-safe.
"""

from __future__ import annotations

import re
from typing import Any

# ---------------------------------------------------------------------------
# Filler / placeholder patterns
# ---------------------------------------------------------------------------

_FILLER_PATTERN = re.compile(
    r"""^\s*(
        n/?a |
        not\s+specified |
        not\s+applicable |
        not\s+provided |
        unknown |
        tbd |
        to\s+be\s+determined |
        none\s+applicable |
        none\s+specified |
        none\s+provided |
        placeholder |
        string |                # literal "string" from schema examples
        \.\.\. |
        <[^>]+> |               # <value>, <catalyst>
        \[.*?\]                 # lone bracket from template artefacts
    )\s*$""",
    re.VERBOSE | re.IGNORECASE,
)

# Text usually indicates a hallucinated/generic response that will fail downstream
_WEAK_FILLER = re.compile(
    r"\b(as\s+needed|if\s+applicable|optional|varies|standard|typical|usual)\b",
    re.IGNORECASE,
)

# ---------------------------------------------------------------------------
# Schema definitions
# ---------------------------------------------------------------------------

# Each schema is a dict with keys:
#   required       – fields that MUST be present and non-empty (use dot-notation for nested)
#   filler_check   – fields where filler text is flagged as a FAILURE (not just warning)
#   typed          – {field: type_or_tuple}  (only top-level keys)
#   ranges         – {field: (min, max)} for numeric fields
#   nested         – {field: nested_schema_name}  for recursive validation

_SCHEMAS: dict[str, dict[str, Any]] = {
    # ------------------------------------------------------------------
    # Direct condition recommendation (CONDITION_RECOMMENDATION prompt)
    # ------------------------------------------------------------------
    "condition_recommendation": {
        "required": [
            "catalyst",
            "solvent",
            "base",
            "temperature_C",
            "time_h",
            "rationale",
            "confidence",
        ],
        "filler_check": ["catalyst", "solvent", "base", "rationale"],
        "typed": {
            "temperature_C": (int, float),
            "time_h": (int, float),
            "confidence": (int, float),
        },
        "ranges": {
            "temperature_C": (-80.0, 300.0),
            "time_h": (0.0, 168.0),
            "confidence": (0.0, 1.0),
        },
    },

    # ------------------------------------------------------------------
    # Multi-source synthesis output (MULTI_SOURCE_SYNTHESIS / V2)
    # ------------------------------------------------------------------
    "multi_source_synthesis": {
        "required": [
            "recommended_condition",
            "confidence_level",
            "consensus_analysis",
            "warnings",
            "recommended_condition.catalyst",
            "recommended_condition.solvent",
            "recommended_condition.base",
            "recommended_condition.rationale",
        ],
        "filler_check": [
            "recommended_condition.catalyst",
            "recommended_condition.solvent",
            "confidence_reasoning",
        ],
        "typed": {
            "warnings": list,
            "backup_conditions": list,
            "confidence_level": str,
        },
        "ranges": {},
        "enum_fields": {
            "confidence_level": {"high", "medium", "low"},
        },
    },

    # ------------------------------------------------------------------
    # Reaction featurization review (REACTION_FEATURIZATION_REVIEW)
    # ------------------------------------------------------------------
    "reaction_featurization_review": {
        "required": [
            "suggested_reaction_type",
            "confidence",
            "rationale",
            "requires_human_review",
        ],
        "filler_check": ["suggested_reaction_type", "rationale"],
        "typed": {
            "confidence": (int, float),
            "requires_human_review": bool,
            "tautomer_or_representation_issue": bool,
            "taxonomy_gap_suspected": bool,
        },
        "ranges": {
            "confidence": (0.0, 1.0),
        },
    },

    # ------------------------------------------------------------------
    # Reagent role classification (REAGENT_ROLE_CLASSIFICATION)
    # ------------------------------------------------------------------
    "reagent_role_classification": {
        "required": ["role", "confidence", "reasoning"],
        "filler_check": ["role", "reasoning"],
        "typed": {
            "confidence": (int, float),
            "role": str,
        },
        "ranges": {
            "confidence": (0.0, 1.0),
        },
        "enum_fields": {
            "role": {
                "ligand", "metal_catalyst", "base", "acid",
                "condensation_agent", "oxidant", "reductant",
                "additive", "solvent", "organo_catalyst",
                "enzyme", "other_reagent",
            }
        },
    },

    # ------------------------------------------------------------------
    # Reagent registry review (REAGENT_REGISTRY_REVIEW)
    # ------------------------------------------------------------------
    "reagent_registry_review": {
        "required": ["status", "proposed_role", "proposed_family", "confidence", "justification"],
        "filler_check": ["justification"],
        "typed": {
            "confidence": (int, float),
            "alerts": list,
        },
        "ranges": {"confidence": (0.0, 1.0)},
        "enum_fields": {"status": {"confirm", "needs_review", "reject"}},
    },

    # ------------------------------------------------------------------
    # Rule builder extraction (RULE_BUILDER_EXTRACTION)
    # ------------------------------------------------------------------
    "rule_builder_extraction": {
        "required": [
            "notes",
            "scope",
            "applies_if",
            "default_rule",
            "base_rules",
            "default_rule.conditions",
        ],
        "filler_check": ["notes"],
        "typed": {
            "base_rules": list,
            "modifiers": list,
        },
        "ranges": {},
    },
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _get_nested(obj: dict, dotted_key: str) -> tuple[bool, Any]:
    """Resolve ``"a.b.c"`` style key into ``obj``.

    Returns (found: bool, value).
    """
    parts = dotted_key.split(".", 1)
    if not isinstance(obj, dict):
        return False, None
    if parts[0] not in obj:
        return False, None
    if len(parts) == 1:
        return True, obj[parts[0]]
    return _get_nested(obj[parts[0]], parts[1])


def _is_empty(value: Any) -> bool:
    """True when a value is semantically absent or meaningless.

    Notes
    -----
    * ``None`` and empty strings are always empty.
    * Empty **dicts** are treated as empty (e.g. a ``consensus_analysis: {}``
      block that the model skipped entirely is a real problem).
    * Empty **lists** are NOT empty — ``"warnings": []`` is a valid "no
      warnings" response, not a missing field.
    """
    if value is None:
        return True
    if isinstance(value, str) and value.strip() == "":
        return True
    if isinstance(value, dict) and len(value) == 0:
        return True
    return False


def _is_filler(value: Any) -> bool:
    """True when a string value matches a placeholder / filler pattern."""
    if not isinstance(value, str):
        return False
    return bool(_FILLER_PATTERN.match(value.strip()))


def _check_type(value: Any, expected: type | tuple) -> bool:
    """Return True if *value* matches *expected* type(s)."""
    if isinstance(expected, tuple):
        return isinstance(value, expected)
    return isinstance(value, expected)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def pass_check(
    response_dict: dict,
    schema_name: str,
    *,
    strict: bool = False,
) -> tuple[bool, list[str], list[str]]:
    """
    Validate a parsed LLM output dict against a named schema.

    Parameters
    ----------
    response_dict:
        The parsed JSON dict returned by the LLM.
    schema_name:
        One of the keys in ``_SCHEMAS`` (e.g. ``"multi_source_synthesis"``).
    strict:
        If *True*, weak-filler text in filler-checked fields causes a
        *failure* rather than a *warning*.

    Returns
    -------
    (passed, failures, warnings):
        ``passed`` is ``True`` only when ``failures`` is empty.
        ``failures`` blocks the result; ``warnings`` are advisory.
    """
    if not isinstance(response_dict, dict):
        return False, ["response_dict is not a dict"], []

    schema = _SCHEMAS.get(schema_name)
    if schema is None:
        available = ", ".join(sorted(_SCHEMAS))
        return False, [f"Unknown schema '{schema_name}'. Available: {available}"], []

    failures: list[str] = []
    warnings: list[str] = []

    # 1. Required-field presence
    for field in schema.get("required", []):
        found, value = _get_nested(response_dict, field)
        if not found:
            failures.append(f"Missing required field: '{field}'")
        elif _is_empty(value):
            failures.append(f"Required field '{field}' is empty")

    # 2. Filler-text checks on designated fields
    for field in schema.get("filler_check", []):
        found, value = _get_nested(response_dict, field)
        if not found or value is None:
            continue  # already caught in step 1
        if _is_filler(str(value)):
            failures.append(f"Field '{field}' contains filler/placeholder text: {value!r}")
        elif strict and isinstance(value, str) and _WEAK_FILLER.search(value):
            failures.append(f"Field '{field}' contains weak/generic text: {value!r}")
        elif not strict and isinstance(value, str) and _WEAK_FILLER.search(value):
            warnings.append(f"Field '{field}' may contain generic text: {value!r}")

    # 3. Type checks
    for field, expected_type in schema.get("typed", {}).items():
        found, value = _get_nested(response_dict, field)
        if not found or value is None:
            continue
        if not _check_type(value, expected_type):
            failures.append(
                f"Field '{field}' has wrong type: expected {expected_type}, "
                f"got {type(value).__name__} ({value!r})"
            )

    # 4. Range checks
    for field, (lo, hi) in schema.get("ranges", {}).items():
        found, value = _get_nested(response_dict, field)
        if not found or value is None:
            continue
        try:
            fval = float(value)
        except (TypeError, ValueError):
            continue
        if not (lo <= fval <= hi):
            failures.append(
                f"Field '{field}' value {fval} is outside expected range [{lo}, {hi}]"
            )

    # 5. Enum checks
    for field, allowed in schema.get("enum_fields", {}).items():
        found, value = _get_nested(response_dict, field)
        if not found or value is None:
            continue
        if isinstance(value, str) and value.lower() not in {v.lower() for v in allowed}:
            failures.append(
                f"Field '{field}' value {value!r} not in allowed set {sorted(allowed)}"
            )

    passed = len(failures) == 0
    return passed, failures, warnings


def validate_batch(
    responses: list[dict],
    schema_name: str,
    *,
    strict: bool = False,
) -> list[dict[str, Any]]:
    """
    Run ``pass_check`` over a list of response dicts.

    Returns a parallel list of ``{"passed", "failures", "warnings", "index"}``
    dicts, one per input.
    """
    results = []
    for i, resp in enumerate(responses):
        passed, failures, warns = pass_check(resp, schema_name, strict=strict)
        results.append({
            "index": i,
            "passed": passed,
            "failures": failures,
            "warnings": warns,
        })
    return results


def list_schemas() -> list[str]:
    """Return the names of all registered validation schemas."""
    return sorted(_SCHEMAS.keys())


def register_schema(name: str, schema: dict[str, Any]) -> None:
    """
    Register a custom validation schema at runtime.

    Example::

        register_schema("my_tool", {
            "required": ["answer", "confidence"],
            "typed": {"confidence": float},
            "ranges": {"confidence": (0.0, 1.0)},
            "filler_check": ["answer"],
        })
    """
    _SCHEMAS[name] = schema


__all__ = [
    "pass_check",
    "validate_batch",
    "list_schemas",
    "register_schema",
]

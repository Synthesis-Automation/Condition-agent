"""
Deterministic boolean-expression evaluation for taxonomy rules.

This module is shared by:
- calculable feature derived_shortcuts (calculable_features.json)
- higher-level chemistry term rules (chem_terms.json)

Supported syntax (space-delimited operators):
  - AND / OR / NOT
  - Parentheses for grouping (must be preceded by whitespace or start)
  - Integer comparisons: token >= 2, token == 0, etc.

Notes:
  - The parser is intentionally conservative and deterministic.
  - It avoids treating parentheses embedded in token names as grouping
    (e.g., "ArB(OH)2_reactant") by requiring grouping parens to be preceded by
    whitespace or start-of-string.
"""

from __future__ import annotations

import re
from typing import Any, Mapping

_COMPARISON_PATTERN = re.compile(r"(\w+)\s*(>=|<=|>|<|==|!=)\s*(\d+)")
def _reduce_grouping_parentheses(expression: str, values: Mapping[str, Any]) -> str:
    """
    Evaluate grouping parentheses iteratively while ignoring parentheses inside tokens.

    Heuristic:
      - Treat '(' as grouping only if it is at start or preceded by whitespace.
      - Treat ')' as grouping only if it is at end or followed by whitespace.

    This prevents token-internal parentheses such as "ArB(OH)2_present" from
    being interpreted as expression grouping.
    """
    max_iterations = 100
    current = expression

    for _ in range(max_iterations):
        stack: list[int] = []
        pairs: list[tuple[int, int]] = []
        length = len(current)

        for idx, ch in enumerate(current):
            if ch == "(":
                if idx == 0 or current[idx - 1].isspace() or current[idx - 1] == "(":
                    stack.append(idx)
            elif ch == ")":
                if stack and (idx == length - 1 or current[idx + 1].isspace() or current[idx + 1] == ")"):
                    start = stack.pop()
                    pairs.append((start, idx))

        if not pairs:
            return current

        # Pick the last discovered pair (innermost in most common cases).
        start, end = pairs[-1]
        inner = current[start + 1 : end]
        inner_value = evaluate(inner, values)
        current = current[:start] + str(inner_value) + current[end + 1 :]

    return current


def evaluate(expr: str, values: Mapping[str, Any]) -> bool:
    """
    Evaluate a boolean expression against a mapping of token values.

    Args:
        expr: Expression like "a_present AND NOT b_present" or "count >= 2".
        values: Mapping from token to bool/int/etc. Missing tokens are falsey.

    Returns:
        Boolean result.
    """
    if not expr:
        return False

    expression = str(expr).strip()
    if not expression:
        return False

    def _coerce_number(value: Any) -> float:
        if value is None:
            return 0.0
        if isinstance(value, bool):
            return 1.0 if value else 0.0
        try:
            return float(value)
        except Exception:
            return 0.0

    def _eval_comparison(match: re.Match[str]) -> str:
        token = match.group(1)
        operator = match.group(2)
        target = float(match.group(3))
        current = _coerce_number(values.get(token))

        if operator == ">=":
            ok = current >= target
        elif operator == "<=":
            ok = current <= target
        elif operator == ">":
            ok = current > target
        elif operator == "<":
            ok = current < target
        elif operator == "==":
            ok = current == target
        elif operator == "!=":
            ok = current != target
        else:
            ok = False

        return "True" if ok else "False"

    # Replace integer comparisons with boolean literals.
    expression = _COMPARISON_PATTERN.sub(_eval_comparison, expression)

    # Reduce grouping parentheses iteratively while ignoring token-internal parens.
    expression = _reduce_grouping_parentheses(expression, values)

    # Protect boolean literals from being interpreted as token names in recursion.
    expression = expression.replace("True", "true_token").replace("False", "false_token")

    # OR has lower precedence than AND.
    if " OR " in expression:
        for part in (p.strip() for p in expression.split(" OR ")):
            if evaluate(part, values):
                return True
        return False

    # Evaluate AND chain.
    parts = [p.strip() for p in expression.split(" AND ") if p.strip()]
    for part in parts:
        if part == "true_token":
            continue
        if part == "false_token":
            return False

        if part.startswith("NOT "):
            token = part[4:].strip()
            if bool(values.get(token, False)):
                return False
            continue

        if not bool(values.get(part, False)):
            return False

    return True

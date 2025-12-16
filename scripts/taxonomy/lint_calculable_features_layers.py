"""
Lint the layered calculable feature specification.

Checks:
  - Layer invariants (foundation vs properties vs derived)
  - No smartsrx-style reactant_types expansion (reactant_types must have reactant_metadata)
  - Derived rules only reference existing tokens (features + derived_shortcuts)
  - SMARTS strings stay small (guardrail against PFAS/super-pattern reintroduction)
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple


_OPERATORS = {"AND", "OR", "NOT", ">=", "<=", ">", "<", "==", "!="}
_COMPARISON_RE = re.compile(r"(\w+)\s*(>=|<=|>|<|==|!=)\s*(\d+)")


def _load(path: Path) -> Dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Expected JSON object at {path}")
    return data


def _iter_features(spec: Dict[str, Any]) -> Iterable[Dict[str, Any]]:
    features = spec.get("features") or []
    if not isinstance(features, list):
        return []
    return (f for f in features if isinstance(f, dict))


def _extract_rule_tokens(expr: str) -> List[str]:
    tokens: List[str] = []
    for part in str(expr).split():
        if part in _OPERATORS:
            continue
        cleaned = part.strip().strip("()")
        if not cleaned or cleaned in {"True", "False"}:
            continue
        if re.fullmatch(r"\d+", cleaned):
            continue
        # For comparisons like "halogen_count >= 2", keep just the token.
        match = _COMPARISON_RE.fullmatch(cleaned)
        if match:
            tokens.append(match.group(1))
            continue
        tokens.append(cleaned)
    return tokens


def _collect_tokens(features: Iterable[Dict[str, Any]]) -> List[str]:
    tokens: List[str] = []
    for feat in features:
        token = feat.get("token")
        if isinstance(token, str) and token:
            tokens.append(token)
    return tokens


def _collect_smarts_strings(features: Iterable[Dict[str, Any]]) -> List[Tuple[str, str]]:
    smarts: List[Tuple[str, str]] = []
    for feat in features:
        token = feat.get("token")
        if not isinstance(token, str) or not token:
            continue
        detect = feat.get("detect")
        if not isinstance(detect, dict):
            continue
        for key in ("smarts_any", "smarts_count"):
            value = detect.get(key)
            if isinstance(value, str) and value:
                smarts.append((token, value))
            elif isinstance(value, list):
                for s in value:
                    if isinstance(s, str) and s:
                        smarts.append((token, s))
    return smarts


def lint_layers(*, root: Path, max_smarts_len: int) -> List[str]:
    errors: List[str] = []

    base_path = root / "calculable_features.json"
    properties_path = root / "calculable_features_properties.json"
    derived_path = root / "calculable_features_derived.json"

    base = _load(base_path)
    props = _load(properties_path)
    derived = _load(derived_path)

    base_features = list(_iter_features(base))
    props_features = list(_iter_features(props))
    derived_features = list(_iter_features(derived))

    # ------------------------------------------------------------------ layer invariants
    if "derived_shortcuts" in base and (base.get("derived_shortcuts") or []):
        errors.append("Base layer should not define derived_shortcuts (keep in derived layer).")

    for feat in base_features:
        if isinstance(feat.get("derive") or feat.get("derived"), str):
            errors.append(f"Base layer feature should not be derived: {feat.get('token')}")
        detect = feat.get("detect")
        if isinstance(detect, dict) and "heuristic" in detect:
            errors.append(f"Base layer feature should not be heuristic: {feat.get('token')}")

    for feat in props_features:
        detect = feat.get("detect")
        if not (isinstance(detect, dict) and "heuristic" in detect):
            errors.append(f"Properties layer must use detect.heuristic: {feat.get('token')}")
        if isinstance(feat.get("derive") or feat.get("derived"), str):
            errors.append(f"Properties layer feature should not be derived: {feat.get('token')}")

    for feat in derived_features:
        derive_expr = feat.get("derive") or feat.get("derived")
        if not (isinstance(derive_expr, str) and derive_expr.strip()):
            errors.append(f"Derived layer must include derive/derived: {feat.get('token')}")

    # ------------------------------------------------------------------ merged-token validation
    all_features = base_features + props_features + derived_features
    feature_tokens = _collect_tokens(all_features)
    duplicates = sorted({t for t in feature_tokens if feature_tokens.count(t) > 1})
    if duplicates:
        errors.append(f"Duplicate feature tokens: {duplicates[:20]}")

    derived_shortcuts = derived.get("derived_shortcuts") or []
    if not isinstance(derived_shortcuts, list):
        errors.append("derived_shortcuts must be a list in derived layer.")
        derived_shortcuts = []

    shortcut_tokens: List[str] = []
    shortcut_rules: List[Tuple[str, str]] = []
    for entry in derived_shortcuts:
        if not isinstance(entry, dict):
            continue
        token = entry.get("token")
        expr = entry.get("derive") or ""
        if isinstance(token, str) and token:
            shortcut_tokens.append(token)
        if isinstance(token, str) and token and isinstance(expr, str) and expr.strip():
            shortcut_rules.append((token, expr.strip()))

    shortcut_dupes = sorted({t for t in shortcut_tokens if shortcut_tokens.count(t) > 1})
    if shortcut_dupes:
        errors.append(f"Duplicate derived_shortcuts tokens: {shortcut_dupes[:20]}")

    overlap = set(feature_tokens) & set(shortcut_tokens)
    # Overlap is allowed only as a "category-level broadening" pattern, where the
    # derived shortcut is self-inclusive (e.g., `RSH_present OR ArSH_present`).
    if overlap:
        by_token: Dict[str, Dict[str, Any]] = {}
        for entry in derived_shortcuts:
            if not isinstance(entry, dict):
                continue
            token = entry.get("token")
            if isinstance(token, str) and token and token not in by_token:
                by_token[token] = entry
        for token in sorted(overlap):
            entry = by_token.get(token) or {}
            expr = entry.get("derive") or ""
            if not isinstance(expr, str) or not expr.strip():
                errors.append(f"Overlapping derived_shortcuts token missing derive: {token}")
                continue
            refs = _extract_rule_tokens(expr)
            if token not in refs:
                errors.append(
                    f"Overlapping derived_shortcuts token must be self-inclusive ({token} not in refs): {token} => {expr}"
                )

    known_tokens = set(feature_tokens) | set(shortcut_tokens)

    # ------------------------------------------------------------------ derived rule references
    for feat in derived_features:
        token = feat.get("token")
        expr = feat.get("derive") or feat.get("derived") or ""
        if not isinstance(expr, str) or not expr.strip():
            continue
        if "+" in expr:
            errors.append(f"Unsupported arithmetic in derived rule for {token!r}: {expr}")
        refs = _extract_rule_tokens(expr)
        missing = sorted({r for r in refs if r not in known_tokens})
        if missing:
            errors.append(f"Derived rule for {token!r} references missing tokens: {missing[:20]}")

    for token, expr in shortcut_rules:
        if "+" in expr:
            errors.append(f"Unsupported arithmetic in derived_shortcuts for {token!r}: {expr}")
        refs = _extract_rule_tokens(expr)
        missing = sorted({r for r in refs if r not in known_tokens})
        if missing:
            errors.append(f"derived_shortcuts rule for {token!r} references missing tokens: {missing[:20]}")

    # ------------------------------------------------------------------ guardrails against heavy SMARTS
    for token, smarts in _collect_smarts_strings(all_features):
        if len(smarts) > max_smarts_len:
            errors.append(f"SMARTS too long ({len(smarts)}>{max_smarts_len}) for {token}: {smarts[:60]}...")
        if "**" in smarts:
            errors.append(f"SMARTS contains disallowed '**' for {token}: {smarts}")

    # Prevent reintroducing the smartsrx expansion layer.
    for feat in all_features:
        if feat.get("category") == "reactant_types" and not isinstance(feat.get("reactant_metadata"), dict):
            errors.append(f"reactant_types feature missing reactant_metadata: {feat.get('token')}")

    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("chemtools/taxonomy/data"),
        help="Taxonomy data directory containing the layered spec files.",
    )
    parser.add_argument(
        "--max-smarts-len",
        type=int,
        default=160,
        help="Maximum allowed SMARTS string length (guardrail against super-patterns).",
    )
    args = parser.parse_args()

    root = args.root.resolve()
    errors = lint_layers(root=root, max_smarts_len=int(args.max_smarts_len))

    if errors:
        print("FAIL: calculable feature spec lint errors:")
        for err in errors:
            print(f"- {err}")
        return 1

    print("OK: calculable feature spec layers are consistent.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

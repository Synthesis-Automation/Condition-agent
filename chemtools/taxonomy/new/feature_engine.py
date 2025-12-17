"""
Suzuki + Buchwald POC feature engine.

Implements the runtime evaluation algorithm described in:
`chemtools/taxonomy/new/CODEX_TEST_BRIEF_Suzuki_Buchwald_POC.md`.

This module intentionally uses the repo-wide SMARTS cache
(`chemtools.util.smarts_cache.compile_smarts`) and RDKit helpers
(`chemtools.util.rdkit_helpers`) for consistency and performance.
"""

from __future__ import annotations

import json
from collections import defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

from chemtools.util.rdkit_helpers import parse_smiles, rdkit_available
from chemtools.util.smarts_cache import compile_smarts


FeatureValues = Dict[str, bool]


class FeatureValidationError(ValueError):
    """Raised when feature definitions fail validation."""


@dataclass(frozen=True)
class ReactionTypeMatch:
    """A matched reaction type with minimal explainability."""

    id: str
    name: str
    description: str
    why_all_of: Tuple[str, ...]
    why_any_of: Tuple[str, ...]
    why_none_of: Tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "name": self.name,
            "description": self.description,
            "why_all_of": list(self.why_all_of),
            "why_any_of": list(self.why_any_of),
            "why_none_of": list(self.why_none_of),
        }


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def load_feature_definitions(
    *,
    atomic_path: Optional[Path] = None,
    derived_path: Optional[Path] = None,
    v2_path: Optional[Path] = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """
    Load atomic/derived feature definitions from JSON.

    Args:
        atomic_path: Path to `calculable_features.atomic.*.json`
        derived_path: Path to `calculable_features.derived.*.json`
        v2_path: Path to `calculable_features.v2.*.json` (combined dict)

    Returns:
        (atomic_defs, derived_defs)
    """
    paths = [p is not None for p in (atomic_path, derived_path, v2_path)]
    if sum(paths) not in {1, 2}:
        raise ValueError("Provide either v2_path, or atomic_path+derived_path.")

    if v2_path is not None:
        payload = _read_json(v2_path)
        if not isinstance(payload, dict):
            raise ValueError(f"Expected object in {v2_path}, got {type(payload).__name__}")
        atomic = payload.get("atomic", [])
        derived = payload.get("derived", [])
        if not isinstance(atomic, list) or not isinstance(derived, list):
            raise ValueError(f"Invalid v2 schema in {v2_path}")
        return atomic, derived

    if atomic_path is None or derived_path is None:
        raise ValueError("Both atomic_path and derived_path are required when v2_path is not provided.")

    atomic = _read_json(atomic_path)
    derived = _read_json(derived_path)
    if not isinstance(atomic, list) or not isinstance(derived, list):
        raise ValueError("Atomic and derived JSON must be lists.")
    return atomic, derived


def compute_atomic_features(mol: Any, atomic_feature_defs: Sequence[Mapping[str, Any]]) -> FeatureValues:
    """
    Compute atomic SMARTS features for a single molecule.

    Args:
        mol: RDKit Mol (or None)
        atomic_feature_defs: Feature definitions with `detect.smarts_any`.

    Returns:
        Mapping token -> bool.
    """
    results: FeatureValues = {}
    if mol is None:
        for feature in atomic_feature_defs:
            token = str(feature.get("token", ""))
            if token:
                results[token] = False
        return results

    for feature in atomic_feature_defs:
        token = str(feature.get("token", ""))
        if not token:
            continue
        detect = feature.get("detect") or {}
        smarts_any = detect.get("smarts_any") or []
        if not isinstance(smarts_any, list):
            raise ValueError(f"{token}: detect.smarts_any must be a list")

        matched = False
        for smarts in smarts_any:
            if not isinstance(smarts, str) or not smarts.strip():
                continue
            pattern = compile_smarts(smarts, validate=False)
            if pattern is None:
                continue
            try:
                if mol.HasSubstructMatch(pattern):
                    matched = True
                    break
            except Exception:
                continue
        results[token] = matched
    return results


def _derive_dependencies(derive_expr: Mapping[str, Any]) -> list[str]:
    deps: list[str] = []
    for key in ("all_of", "any_of", "none_of"):
        values = derive_expr.get(key, [])
        if not values:
            continue
        if not isinstance(values, list):
            raise ValueError(f"derive.{key} must be a list")
        deps.extend(str(t) for t in values)
    return deps


def _toposort_derived(derived_defs: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    token_to_def: dict[str, dict[str, Any]] = {}
    for entry in derived_defs:
        token = str(entry.get("token", ""))
        if not token:
            continue
        token_to_def[token] = dict(entry)

    derived_tokens = set(token_to_def.keys())
    indegree: dict[str, int] = {t: 0 for t in derived_tokens}
    graph: dict[str, set[str]] = {t: set() for t in derived_tokens}

    for token, entry in token_to_def.items():
        derive_expr = entry.get("derive") or {}
        if not isinstance(derive_expr, dict):
            raise ValueError(f"{token}: derive must be an object")
        deps = _derive_dependencies(derive_expr)
        for dep in deps:
            if dep in derived_tokens:
                if token not in graph[dep]:
                    graph[dep].add(token)
                    indegree[token] += 1

    queue = deque(sorted([t for t, deg in indegree.items() if deg == 0]))
    ordered: list[str] = []
    while queue:
        current = queue.popleft()
        ordered.append(current)
        for nxt in sorted(graph[current]):
            indegree[nxt] -= 1
            if indegree[nxt] == 0:
                queue.append(nxt)

    if len(ordered) != len(derived_tokens):
        remaining = sorted([t for t in derived_tokens if t not in set(ordered)])
        raise FeatureValidationError(f"Derived feature graph has a cycle involving: {remaining}")

    return [token_to_def[t] for t in ordered]


def compute_derived_features(
    base_values: Mapping[str, bool], derived_feature_defs: Sequence[Mapping[str, Any]]
) -> FeatureValues:
    """
    Compute derived features using a dependency-aware evaluation order.

    Args:
        base_values: Existing token -> bool mapping (usually atomic results).
        derived_feature_defs: Derived definitions with `derive.any_of/all_of/none_of`.

    Returns:
        New dict containing base_values plus derived results.
    """
    values: FeatureValues = dict(base_values)
    ordered_defs = _toposort_derived(derived_feature_defs)

    for feature in ordered_defs:
        token = str(feature.get("token", ""))
        if not token:
            continue
        derive_expr = feature.get("derive") or {}
        if not isinstance(derive_expr, dict):
            raise ValueError(f"{token}: derive must be an object")

        result = True

        if "all_of" in derive_expr:
            deps = derive_expr.get("all_of") or []
            if not isinstance(deps, list):
                raise ValueError(f"{token}: derive.all_of must be a list")
            result = result and all(bool(values.get(str(t), False)) for t in deps)

        if "any_of" in derive_expr:
            deps = derive_expr.get("any_of") or []
            if not isinstance(deps, list):
                raise ValueError(f"{token}: derive.any_of must be a list")
            result = result and any(bool(values.get(str(t), False)) for t in deps)

        if "none_of" in derive_expr:
            deps = derive_expr.get("none_of") or []
            if not isinstance(deps, list):
                raise ValueError(f"{token}: derive.none_of must be a list")
            result = result and (not any(bool(values.get(str(t), False)) for t in deps))

        values[token] = bool(result)

    return values


def compute_all_features(
    mol: Any,
    *,
    atomic_defs: Sequence[Mapping[str, Any]],
    derived_defs: Sequence[Mapping[str, Any]],
) -> FeatureValues:
    """Compute atomic + derived features for a single molecule."""
    atomic = compute_atomic_features(mol, atomic_defs)
    return compute_derived_features(atomic, derived_defs)


def validate_features(
    atomic_defs: Sequence[Mapping[str, Any]],
    derived_defs: Sequence[Mapping[str, Any]],
    *,
    reactant_types: Optional[Sequence[Mapping[str, Any]]] = None,
    reaction_types: Optional[Sequence[Mapping[str, Any]]] = None,
) -> list[str]:
    """
    Validate feature definitions for common errors.

    Returns:
        List of warning strings. Raises FeatureValidationError on hard failures.
    """
    warnings: list[str] = []

    atomic_tokens: list[str] = []
    for entry in atomic_defs:
        token = str(entry.get("token", "")).strip()
        if token:
            atomic_tokens.append(token)

    derived_tokens: list[str] = []
    for entry in derived_defs:
        token = str(entry.get("token", "")).strip()
        if token:
            derived_tokens.append(token)

    all_tokens = atomic_tokens + derived_tokens
    if not all_tokens:
        raise FeatureValidationError("No feature tokens found.")

    seen: set[str] = set()
    duplicates_set: set[str] = set()
    for token in all_tokens:
        if token in seen:
            duplicates_set.add(token)
        seen.add(token)
    duplicates = sorted(duplicates_set)
    if duplicates:
        raise FeatureValidationError(f"Duplicate feature tokens: {duplicates}")

    token_set = set(all_tokens)

    if not rdkit_available():
        raise FeatureValidationError("RDKit is not available; cannot validate SMARTS patterns.")

    seen_smarts: dict[str, list[str]] = defaultdict(list)
    for entry in atomic_defs:
        token = str(entry.get("token", "")).strip()
        detect = entry.get("detect") or {}
        if not isinstance(detect, dict):
            raise FeatureValidationError(f"{token}: detect must be an object")
        smarts_any = detect.get("smarts_any") or []
        if not isinstance(smarts_any, list):
            raise FeatureValidationError(f"{token}: detect.smarts_any must be a list")
        for smarts in smarts_any:
            if not isinstance(smarts, str) or not smarts.strip():
                continue
            compile_smarts(smarts, validate=True)
            seen_smarts[smarts].append(token)

    redundant = sorted([(s, toks) for s, toks in seen_smarts.items() if len(toks) > 1])
    if redundant:
        for smarts, toks in redundant:
            warnings.append(f"SMARTS pattern reused across tokens {sorted(toks)}: {smarts}")

    for entry in derived_defs:
        token = str(entry.get("token", "")).strip()
        derive_expr = entry.get("derive") or {}
        if not isinstance(derive_expr, dict):
            raise FeatureValidationError(f"{token}: derive must be an object")
        for dep in _derive_dependencies(derive_expr):
            if dep not in token_set:
                raise FeatureValidationError(f"{token}: references unknown token '{dep}'")

    _toposort_derived(derived_defs)

    if reaction_types is not None:
        for rt in reaction_types:
            rt_id = str(rt.get("id", "")).strip()
            when = rt.get("when") or {}
            if not isinstance(when, dict):
                raise FeatureValidationError(f"{rt_id}: when must be an object")
            for key in ("all_of", "any_of", "none_of"):
                for token in when.get(key, []) or []:
                    token_str = str(token)
                    if token_str not in token_set:
                        raise FeatureValidationError(f"{rt_id}: when.{key} references unknown token '{token_str}'")

    if reactant_types is not None:
        missing: list[tuple[str, str]] = []

        def _walk(node: Mapping[str, Any]) -> None:
            node_id = str(node.get("id", "")).strip()
            token = str(node.get("feature_token", "")).strip()
            if token and token not in token_set:
                missing.append((node_id, token))
            for child in node.get("members") or []:
                if isinstance(child, dict):
                    _walk(child)

        for root in reactant_types:
            _walk(root)

        if missing:
            missing_sorted = sorted(missing)
            raise FeatureValidationError(f"Reactant type nodes reference unknown tokens: {missing_sorted}")

    return warnings


def match_reactant_types(
    feature_values: Mapping[str, bool],
    reactant_type_tree: Sequence[Mapping[str, Any]],
    *,
    leaf_only: bool = True,
) -> list[dict[str, Any]]:
    """
    Match reactant taxonomy nodes against computed features.

    Args:
        feature_values: token -> bool mapping.
        reactant_type_tree: Root nodes from `reactant_types.*.json`.
        leaf_only: If True, return only most-specific matches.

    Returns:
        List of matched nodes (each includes id/name/feature_token/path_ids).
    """
    matches: list[dict[str, Any]] = []

    def _walk(node: Mapping[str, Any], path: tuple[str, ...]) -> list[dict[str, Any]]:
        token = str(node.get("feature_token", "")).strip()
        if not token or not bool(feature_values.get(token, False)):
            return []

        node_id = str(node.get("id", "")).strip()
        current_path = path + (node_id,) if node_id else path

        children = node.get("members") or []
        child_matches: list[dict[str, Any]] = []
        for child in children:
            if isinstance(child, dict):
                child_matches.extend(_walk(child, current_path))

        node_match = {
            "id": node_id,
            "name": str(node.get("name", "")),
            "feature_token": token,
            "path_ids": current_path,
        }

        if leaf_only:
            return child_matches if child_matches else [node_match]
        return [node_match] + child_matches

    for root in reactant_type_tree:
        if isinstance(root, dict):
            matches.extend(_walk(root, ()))

    return matches


def aggregate_reaction_features(per_reactant: Sequence[Mapping[str, bool]]) -> tuple[FeatureValues, dict[str, int]]:
    """
    Aggregate reactant-level booleans into reaction-level booleans and counts.
    """
    reaction_has: FeatureValues = {}
    reaction_count: dict[str, int] = defaultdict(int)

    for feats in per_reactant:
        for token, present in feats.items():
            if bool(present):
                reaction_has[token] = True
                reaction_count[token] += 1
            else:
                reaction_has.setdefault(token, False)
                reaction_count.setdefault(token, 0)

    return reaction_has, dict(reaction_count)


def match_reaction_types(
    reaction_has: Mapping[str, bool], reaction_types: Sequence[Mapping[str, Any]]
) -> list[ReactionTypeMatch]:
    """
    Match reaction types based on a reaction-level feature mapping.
    """
    hits: list[ReactionTypeMatch] = []
    for rt in reaction_types:
        rt_id = str(rt.get("id", "")).strip()
        name = str(rt.get("name", "")).strip()
        desc = str(rt.get("description", "")).strip()
        when = rt.get("when") or {}
        if not isinstance(when, dict):
            raise ValueError(f"{rt_id}: when must be an object")

        all_of = [str(t) for t in (when.get("all_of") or [])]
        any_of = [str(t) for t in (when.get("any_of") or [])]
        none_of = [str(t) for t in (when.get("none_of") or [])]

        ok = True
        if all_of:
            ok = ok and all(bool(reaction_has.get(t, False)) for t in all_of)
        if "any_of" in when:
            ok = ok and any(bool(reaction_has.get(t, False)) for t in any_of)
        if none_of:
            ok = ok and (not any(bool(reaction_has.get(t, False)) for t in none_of))

        if not ok:
            continue

        why_all = tuple(sorted([t for t in all_of if bool(reaction_has.get(t, False))]))
        why_any = tuple(sorted([t for t in any_of if bool(reaction_has.get(t, False))]))
        why_none = tuple(sorted([t for t in none_of if not bool(reaction_has.get(t, False))]))

        hits.append(
            ReactionTypeMatch(
                id=rt_id,
                name=name,
                description=desc,
                why_all_of=why_all,
                why_any_of=why_any,
                why_none_of=why_none,
            )
        )
    return hits


def classify_reaction_smiles(
    reactant_smiles: Sequence[str],
    *,
    atomic_defs: Sequence[Mapping[str, Any]],
    derived_defs: Sequence[Mapping[str, Any]],
    reactant_types: Optional[Sequence[Mapping[str, Any]]] = None,
    reaction_types: Optional[Sequence[Mapping[str, Any]]] = None,
) -> dict[str, Any]:
    """
    Classify a reaction (list of reactant SMILES) into reaction types.

    Returns a dict suitable for CLI printing or tests.
    """
    if not rdkit_available():
        raise RuntimeError("RDKit is not available; cannot classify SMILES.")

    per_reactant: list[dict[str, Any]] = []
    per_features: list[FeatureValues] = []

    for smiles in reactant_smiles:
        mol = parse_smiles(smiles)
        if mol is None:
            raise ValueError(f"Bad SMILES: {smiles}")
        feats = compute_all_features(mol, atomic_defs=atomic_defs, derived_defs=derived_defs)
        per_features.append(feats)

        reactant_matches = None
        if reactant_types is not None:
            reactant_matches = match_reactant_types(feats, reactant_types, leaf_only=True)

        per_reactant.append(
            {
                "smiles": smiles,
                "features": feats,
                "reactant_type_matches": reactant_matches,
            }
        )

    reaction_has, reaction_count = aggregate_reaction_features(per_features)

    rxn_hits: Optional[list[ReactionTypeMatch]] = None
    if reaction_types is not None:
        rxn_hits = match_reaction_types(reaction_has, reaction_types)

    return {
        "reactants": per_reactant,
        "reaction_has": reaction_has,
        "reaction_count": reaction_count,
        "reaction_type_matches": rxn_hits,
    }

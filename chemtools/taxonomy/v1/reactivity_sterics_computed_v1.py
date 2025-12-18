#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computed steric reactivity features (POC v1).

Implements the `reactivity_features.computed.v1.json` spec by:
- extracting Ar-* anchor sites from templated atomic SMARTS (atom-map `:1`)
- computing per-site ortho substitution counts
- aggregating into max/sum summary features
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    from .ar_context_sterics_v1 import (
        DEFAULT_DATA_DIR,
        find_ar_anchor_sites,
        load_compiled_features,
        load_groups,
        ortho_substitution_count,
        parse_smiles,
        rdkit_available,
        resolve_data_path,
    )
except Exception:
    from ar_context_sterics_v1 import (
        DEFAULT_DATA_DIR,
        find_ar_anchor_sites,
        load_compiled_features,
        load_groups,
        ortho_substitution_count,
        parse_smiles,
        rdkit_available,
        resolve_data_path,
    )


def _read_json(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object in {path}, got {type(payload).__name__}")
    return payload


def load_computed_feature_spec(path: Path) -> Dict[str, Any]:
    return _read_json(path)


def _aggregate_max(values: List[int]) -> int:
    return max(values) if values else 0


def _aggregate_sum(values: List[int]) -> int:
    return int(sum(values)) if values else 0


@dataclass(frozen=True)
class _AnchorContext:
    context_group_id: str
    context_atom_map_num: int


def compute_reactivity_sterics_poc_v1(
    smiles: str,
    *,
    compiled_features_path: Optional[Path] = None,
    groups_path: Optional[Path] = None,
    computed_spec_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """
    Compute steric features defined in `reactivity_features.computed.v1.json`.

    Returns a dict that includes:
    - computed feature tokens
    - `sites`: per-anchor site debug/explainability payload
    """
    if compiled_features_path is None:
        compiled_features_path = DEFAULT_DATA_DIR / "calculable_features.compiled.v1.json"
    if groups_path is None:
        groups_path = DEFAULT_DATA_DIR / "organic_groups.v1.json"
    if computed_spec_path is None:
        computed_spec_path = DEFAULT_DATA_DIR / "reactivity_features.computed.v1.json"

    compiled_path = resolve_data_path(compiled_features_path)
    groups_path_resolved = resolve_data_path(groups_path)
    computed_path = resolve_data_path(computed_spec_path)

    if not rdkit_available():
        return {"smiles": smiles, "error": "RDKit is not available"}

    compiled = load_compiled_features(compiled_path)
    groups_doc = load_groups(groups_path_resolved)
    spec = load_computed_feature_spec(computed_path)

    mol = parse_smiles(smiles)
    if mol is None:
        return {"smiles": smiles, "error": "Invalid SMILES"}

    computed_defs = spec.get("computed_features") or []
    if not isinstance(computed_defs, list):
        raise ValueError(f"{computed_path}: computed_features must be a list")

    anchor_specs = [
        (d.get("compute") or {}).get("anchor")
        for d in computed_defs
        if isinstance(d, dict) and isinstance((d.get("compute") or {}).get("anchor"), dict)
    ]
    anchor_spec = anchor_specs[0] if anchor_specs else {}

    ctx = _AnchorContext(
        context_group_id=str(anchor_spec.get("context_group_id", "Ar")),
        context_atom_map_num=int(anchor_spec.get("context_atom_map_num", 1)),
    )

    sites = find_ar_anchor_sites(
        mol,
        compiled,
        groups_doc,
        context_group_id=ctx.context_group_id,
        context_atom_map_num=ctx.context_atom_map_num,
    )

    site_payloads: List[Dict[str, Any]] = []
    ortho_counts: List[int] = []
    for s in sites:
        c = ortho_substitution_count(mol, s.ipso_idx)
        # Defensive fallback: for a valid Ar-* anchor, we expect ring/ortho atoms to exist.
        ortho_counts.append(int(c) if c is not None else 0)
        site_payloads.append(
            {
                "label": s.label,
                "feature_token": s.token,
                "ipso_atom_idx": s.ipso_idx,
                "ortho_sub_count": int(c) if c is not None else 0,
            }
        )

    values: Dict[str, Any] = {"smiles": smiles, "sites": site_payloads}

    for feature_def in computed_defs:
        if not isinstance(feature_def, dict):
            continue
        token = str(feature_def.get("token", "")).strip()
        compute = feature_def.get("compute") or {}
        if not token or not isinstance(compute, dict):
            continue

        method = str(compute.get("method", "")).strip()
        if method == "count_anchor_sites":
            values[token] = int(len(site_payloads))
        elif method == "aryl_ortho_substitution_count":
            values[token] = list(ortho_counts)
        elif method == "aggregate_max":
            from_token = str(compute.get("from", "")).strip()
            src = values.get(from_token, [])
            values[token] = _aggregate_max(list(src) if isinstance(src, list) else [])
        elif method == "aggregate_sum":
            from_token = str(compute.get("from", "")).strip()
            src = values.get(from_token, [])
            values[token] = _aggregate_sum(list(src) if isinstance(src, list) else [])
        else:
            raise ValueError(f"Unknown computed feature method: {method} (token={token})")

    return values

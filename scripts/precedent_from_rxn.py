from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

# Ensure project root on path when running from source
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from chemtools.smiles import normalize_reaction
from chemtools.router import detect_family
from chemtools.featurizers import molecular as feat_molecular
from chemtools.precedent import knn

try:
    from chemtools.integrations import molpipeline as molpipe  # type: ignore

    _MOLPIPE_ENV = molpipe.environment_snapshot()
    _MOLPIPE_AVAILABLE = bool(_MOLPIPE_ENV.available)
except Exception:
    molpipe = None  # type: ignore
    _MOLPIPE_ENV = None
    _MOLPIPE_AVAILABLE = False


def _pick_electrophile_nucleophile(reactants: List[str]) -> Tuple[str, str]:
    def is_electrophile(s: str) -> bool:
        t = (s or "").lower()
        return (
            ("br" in t) or ("cl" in t) or (" i" in t)
            or ("os(=o)(=o)c(f)(f)f" in t) or ("otf" in t)
        )

    if not reactants:
        return "", ""
    if len(reactants) == 1:
        return reactants[0], ""
    r0, r1 = reactants[0], reactants[1]
    if is_electrophile(r0):
        return r0, r1
    if is_electrophile(r1):
        return r1, r0
    return r0, r1


def _load_json_arg(value: str | None) -> Dict[str, Any]:
    if not value:
        return {}
    candidate = value.strip()
    if not candidate:
        return {}
    data = candidate
    path_candidate: Path | None = None
    if candidate.startswith("@"):
        path_candidate = Path(candidate[1:])
    else:
        maybe_path = Path(candidate)
        if maybe_path.exists() and maybe_path.is_file():
            path_candidate = maybe_path
    if path_candidate is not None and path_candidate.is_file():
        try:
            data = path_candidate.read_text(encoding="utf-8")
        except Exception:
            return {}
    try:
        obj = json.loads(data)
    except Exception:
        return {}
    return obj if isinstance(obj, dict) else {}


def _sanitize_molpipeline_config(raw: Dict[str, Any]) -> Dict[str, Any]:
    cfg: Dict[str, Any] = {}
    roles = raw.get("roles")
    if isinstance(roles, (list, tuple, set)):
        clean_roles = []
        for role in roles:
            if isinstance(role, str) and role.strip():
                clean_roles.append(role.strip().upper())
        if clean_roles:
            cfg["roles"] = clean_roles
    include_role = raw.get("include_role_features")
    cfg["include_role_features"] = bool(True if include_role is None else include_role)
    include_concat = raw.get("include_concat")
    cfg["include_concat"] = bool(True if include_concat is None else include_concat)
    cfg["suppress_errors"] = bool(raw.get("suppress_errors", True))
    aggregate = raw.get("aggregate")
    if isinstance(aggregate, str) and aggregate.strip():
        cfg["aggregate"] = aggregate.strip()
    missing_strategy = raw.get("missing_strategy")
    if isinstance(missing_strategy, str) and missing_strategy.strip():
        cfg["missing_strategy"] = missing_strategy.strip()
    for key in ("n_jobs", "ligand_n_bits", "ligand_radius"):
        if key in raw:
            try:
                cfg[key] = int(raw[key])
            except Exception:
                continue
    return cfg


def _sanitize_query_map(raw: Dict[str, Any]) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    if not isinstance(raw, dict):
        return out
    for role, value in raw.items():
        if not isinstance(role, str):
            continue
        role_norm = role.strip().upper()
        if not role_norm:
            continue
        smiles: List[str] = []
        if isinstance(value, str):
            trimmed = value.strip()
            if trimmed:
                smiles.append(trimmed)
        elif isinstance(value, (list, tuple, set)):
            for item in value:
                if isinstance(item, str) and item.strip():
                    smiles.append(item.strip())
        if smiles:
            out[role_norm] = smiles
    return out


def main(argv: List[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Run precedent (kNN) search from a reaction SMILES")
    p.add_argument("reaction", help="Reaction SMILES, e.g. 'Brc1ccccc1.NCCOC>>COCCNc1ccccc1'")
    p.add_argument("--k", type=int, default=10, help="Number of neighbors (default: 10)")
    p.add_argument("--strict-bin", dest="strict_bin", action="store_true", help="Require exact bin matches before fallbacks")
    p.add_argument("--no-strict-bin", dest="strict_bin", action="store_false", help="Allow fallbacks when exact bin scarce")
    p.set_defaults(strict_bin=False)
    p.add_argument("--min-candidates", type=int, default=5, help="Minimum candidates to gather with fallbacks")
    p.add_argument("--molpipeline", action="store_true", help="Attach MolPipeline role features when the optional extras are installed")
    p.add_argument("--molpipeline-config", default="", help="MolPipeline config JSON string or @/path/to/config.json")
    p.add_argument("--molpipeline-query", default="", help="MolPipeline query role SMILES JSON string or @/path/to/query.json")
    p.add_argument("--jsonl", action="store_true", help="Emit a compact JSON object")
    p.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    args = p.parse_args(argv)

    # 1) Normalize reaction and extract reactant SMILES
    norm = normalize_reaction(args.reaction)
    reactants = [
        (r.get("smiles_norm") or r.get("largest_smiles") or r.get("input") or "")
        for r in (norm.get("reactants") or [])
    ]

    # 2) Detect family
    fam = detect_family(reactants).get("family") or "Unknown"

    # 3) Pick electrophile/nucleophile and featurize (Ullmann C-N supported)
    elec, nuc = _pick_electrophile_nucleophile(reactants)
    features: Dict[str, Any] = feat_molecular.featurize(elec, nuc)

    # 4) Run kNN
    # Pass normalized reaction SMILES and enable DRFP re-ranking when available
    relax: Dict[str, Any] = {
        "strict_bin": bool(args.strict_bin),
        "min_candidates": int(args.min_candidates),
        "reaction_smiles": norm.get("normalized") or args.reaction,
        "use_drfp": True,  # best-effort; falls back gracefully if DRFP not installed
        "precompute_drfp": True,
        "precompute_scope": "candidates",  # or "dataset" to warm all rows (slower)
    }

    molpipeline_summary: Dict[str, Any] | None = None
    if args.molpipeline:
        molpipeline_summary = {
            "available": bool(_MOLPIPE_AVAILABLE),
        }
        if _MOLPIPE_ENV is not None:
            molpipeline_summary.update(
                {
                    "version": getattr(_MOLPIPE_ENV, "version", None),
                    "rdkit": getattr(_MOLPIPE_ENV, "rdkit_version", None),
                    "sklearn": getattr(_MOLPIPE_ENV, "sklearn_version", None),
                    "shap": getattr(_MOLPIPE_ENV, "shap_version", None),
                }
            )
        if not _MOLPIPE_AVAILABLE:
            molpipeline_summary["enabled"] = False
            molpipeline_summary["error"] = (
                "MolPipeline integration not available. Install chemtools with the molpipeline extras."
            )
        else:
            cfg_raw = _load_json_arg(args.molpipeline_config)
            cfg = _sanitize_molpipeline_config(cfg_raw)
            query_raw = _load_json_arg(args.molpipeline_query)
            query_map = _sanitize_query_map(query_raw)
            if query_map:
                cfg["query_role_smiles"] = query_map
            relax["molpipeline"] = cfg
            molpipeline_summary["enabled"] = True
            summary_cfg: Dict[str, Any] = {
                key: cfg[key]
                for key in (
                    "roles",
                    "include_role_features",
                    "include_concat",
                    "aggregate",
                    "missing_strategy",
                    "n_jobs",
                    "ligand_n_bits",
                    "ligand_radius",
                )
                if key in cfg
            }
            if "query_role_smiles" in cfg:
                summary_cfg["query_roles"] = sorted(cfg["query_role_smiles"].keys())
            if summary_cfg:
                molpipeline_summary["config"] = summary_cfg

    result = knn(family=fam, features=features, k=int(args.k), relax=relax)

    if molpipeline_summary is not None and isinstance(result, dict):
        warnings = result.get("molpipeline_warnings")
        if warnings:
            molpipeline_summary["warnings"] = warnings
        query_vec = result.get("molpipeline_query_vector")
        if isinstance(query_vec, list):
            molpipeline_summary["query_vector_length"] = len(query_vec)
        first_vec = None
        for row in result.get("precedents") or []:
            if not isinstance(row, dict):
                continue
            vec = row.get("molpipeline_feature_vector")
            if isinstance(vec, list):
                first_vec = vec
                break
        if first_vec is not None:
            molpipeline_summary["feature_length"] = len(first_vec)

    payload = {
        "input_reaction": args.reaction,
        "family": fam,
        "electrophile": elec,
        "nucleophile": nuc,
        "features": features,
        "result": result,
    }
    if molpipeline_summary is not None:
        payload["molpipeline"] = molpipeline_summary

    if args.jsonl and not args.pretty:
        print(json.dumps(payload, ensure_ascii=False))
    else:
        print(json.dumps(payload, ensure_ascii=False, indent=(2 if args.pretty else None)))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional
import json

try:  # Optional dependency
    from jsonschema import validate  # type: ignore
except Exception:  # pragma: no cover - jsonschema may be absent
    validate = None  # type: ignore


PACKAGE_DIR = Path(__file__).resolve().parent
DEFAULT_CRL_PATH = PACKAGE_DIR / "crl.json"
DEFAULT_SCHEMA_PATH = PACKAGE_DIR / "schema.json"

_STATUS_RANK = {"stable": 0, "candidate": 1, "experimental": 2, "deprecated": 3}

_SECTION_KEY_MAPPINGS: Dict[str, Dict[str, str]] = {
    "electrophile": {
        "electrophile_class": "class",
        "electrophile_subclass": "subclass",
        "electrophile_electronics": "electronics",
        "electrophile_aromaticity": "aromaticity",
        "electrophile_hybridization": "hybridization",
        "electrophile_form": "form",
        "electrophile_state": "state",
        "ortho_sub_count": "ortho_sub_count",
        "para_sub_count": "para_sub_count",
    },
    "nucleophile": {
        "nucleophile_class": "class",
        "nucleophile_subclass": "subclass",
        "nucleophile_type": "type",
        "nucleophile_identity": "identity",
    },
    "ligand": {"ligand_class": "class"},
}


@dataclass
class GuardResult:
    ok: bool
    messages: List[Dict[str, Any]]
    patch: Dict[str, Any]
    score_adjust: float


def _merge_nested(base: Dict[str, Any], extra: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    result: Dict[str, Any] = {**(base or {})}
    for key, value in (extra or {}).items():
        if isinstance(result.get(key), dict) and isinstance(value, dict):
            result[key] = _merge_nested(result[key], value)
        else:
            result[key] = value
    return result


def _augment_features(features: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    ctx: Dict[str, Any] = _merge_nested({}, features or {})

    for section, key_map in _SECTION_KEY_MAPPINGS.items():
        section_value = ctx.get(section)
        if isinstance(section_value, dict):
            section_dict = dict(section_value)
        else:
            section_dict = {} if section_value is None else {"value": section_value}

        updated = False
        for source_key, target_key in key_map.items():
            if source_key in ctx and ctx[source_key] is not None:
                if not isinstance(section_dict, dict):
                    section_dict = {}
                if target_key not in section_dict:
                    section_dict[target_key] = ctx[source_key]
                    updated = True

        if updated:
            ctx[section] = section_dict

    return ctx


def _flatten_context(data: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    flat: Dict[str, Any] = {}
    if not isinstance(data, dict):
        return flat

    def visit(value: Any, path: List[str]) -> None:
        if isinstance(value, dict):
            for k, v in value.items():
                visit(v, path + [str(k)])
        else:
            if not path:
                return
            key_dot = ".".join(path)
            key_us = "_".join(path)
            flat[key_dot] = value
            flat[key_us] = value

    for k, v in data.items():
        flat[k] = v
        visit(v, [str(k)])
    return flat


def _to_float(value: Any) -> Optional[float]:
    try:
        if isinstance(value, (int, float)):
            return float(value)
        if isinstance(value, str):
            return float(value.strip())
    except Exception:
        return None
    return None


def _match_operator_dict(value: Any, condition: Dict[str, Any]) -> bool:
    for op, target in condition.items():
        op_key = str(op).lower()
        if op_key == "in":
            if isinstance(value, (list, tuple, set)):
                value_set = {str(v).lower() for v in value}
            else:
                value_set = {str(value).lower()}
            allowed = {str(v).lower() for v in (target or [])}
            if not value_set & allowed:
                return False
        elif op_key == "not_in":
            if isinstance(value, (list, tuple, set)):
                value_set = {str(v).lower() for v in value}
            else:
                value_set = {str(value).lower()}
            forbidden = {str(v).lower() for v in (target or [])}
            if value_set & forbidden:
                return False
        elif op_key in {">=", "<=", ">", "<"}:
            value_num = _to_float(value)
            target_num = _to_float(target)
            if value_num is None or target_num is None:
                return False
            if op_key == ">=" and not (value_num >= target_num):
                return False
            if op_key == "<=" and not (value_num <= target_num):
                return False
            if op_key == ">" and not (value_num > target_num):
                return False
            if op_key == "<" and not (value_num < target_num):
                return False
        elif op_key in {"eq", "equals"}:
            if str(value).lower() != str(target).lower():
                return False
        else:
            # Unknown operator; fallback to equality on nested field
            if isinstance(value, dict):
                if not _match_condition(value.get(op), target):
                    return False
            else:
                return False
    return True


def _match_condition(value: Any, condition: Any) -> bool:
    if isinstance(condition, dict):
        operator_keys = {"in", "not_in", ">", "<", ">=", "<=", "eq", "equals"}
        if any(str(k).lower() in operator_keys for k in condition.keys()):
            return _match_operator_dict(value, condition)
        if not isinstance(value, dict):
            return False
        return all(_match_condition(value.get(k), v) for k, v in condition.items())

    if isinstance(condition, list):
        if isinstance(value, (list, tuple, set)):
            value_norm = {str(v).lower() for v in value}
            cond_norm = {str(v).lower() for v in condition}
            return bool(value_norm & cond_norm)
        return str(value).lower() in {str(v).lower() for v in condition}

    if isinstance(value, (list, tuple, set)):
        return str(condition).lower() in {str(v).lower() for v in value}

    if value is None:
        return False
    return str(value).lower() == str(condition).lower()


def _match_when(when: Dict[str, Any], context: Dict[str, Any]) -> bool:
    if not when:
        return True
    for key, condition in when.items():
        value = context.get(key)
        if not _match_condition(value, condition):
            return False
    return True


def load_crl(
    path: Optional[Path | str] = None,
    schema_path: Optional[Path | str] = None,
    *,
    validate_schema: bool = True,
) -> Dict[str, Any]:
    crl_path = Path(path or DEFAULT_CRL_PATH)
    schema = Path(schema_path or DEFAULT_SCHEMA_PATH)

    data = json.loads(crl_path.read_text(encoding="utf-8"))
    if validate_schema and validate is not None and schema.exists():
        schema_data = json.loads(schema.read_text(encoding="utf-8"))
        validate(instance=data, schema=schema_data)  # type: ignore[arg-type]
    return data


def load_default_crl() -> Dict[str, Any]:
    return load_crl()


def _playbook_sort_key(playbook: Dict[str, Any]) -> tuple[int, float]:
    status = playbook.get("status", "candidate")
    rank = _STATUS_RANK.get(str(status).lower(), 1)
    metrics = playbook.get("metrics") or {}
    runs = float(metrics.get("n_runs", 0)) if isinstance(metrics, dict) else 0.0
    return (rank, -runs)


def select_playbooks(
    crl: Dict[str, Any],
    family: str,
    features: Dict[str, Any],
    *,
    fallback_to_all: bool = True,
) -> List[Dict[str, Any]]:
    fam = (crl.get("families") or {}).get(family, {})
    playbooks: List[Dict[str, Any]] = list(fam.get("playbooks") or [])
    if not playbooks:
        return []

    context = _augment_features(features)
    matched = [pb for pb in playbooks if _match_when(pb.get("when") or {}, context)]
    if not matched and fallback_to_all:
        matched = playbooks

    return sorted(matched, key=_playbook_sort_key)


def _combine_contexts(*contexts: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    combined: Dict[str, Any] = {}
    for ctx in contexts:
        if isinstance(ctx, dict):
            combined = _merge_nested(combined, ctx)
    return combined


def apply_guards(
    crl: Dict[str, Any],
    family: str,
    candidate: Dict[str, Any],
    guard_context: Optional[Dict[str, Any]] = None,
) -> GuardResult:
    fam = (crl.get("families") or {}).get(family, {})
    guards: Iterable[Dict[str, Any]] = fam.get("guards") or []
    ctx = _merge_nested({}, guard_context or {})

    messages: List[Dict[str, Any]] = []
    patch: Dict[str, Any] = {}
    score_adjust = 0.0
    ok = True

    for guard in guards:
        condition = guard.get("if") or {}
        if not _match_when(condition, ctx):
            continue

        action = str(guard.get("action", "warn")).lower()
        message = {
            "id": guard.get("id"),
            "action": action,
            "rationale": guard.get("rationale"),
            "params": guard.get("params"),
        }
        messages.append(message)

        if action == "enforce":
            ok = False
            patch = _merge_nested(patch, guard.get("params") or {})
        elif action in {"avoid_or_penalize", "penalize"}:
            score_adjust -= 0.4 if action == "avoid_or_penalize" else 0.3
        elif action == "boost":
            score_adjust += 0.3

    return GuardResult(ok=ok, messages=messages, patch=patch, score_adjust=score_adjust)


def _lookup(flat_ctx: Dict[str, Any], key: str) -> Any:
    if key in flat_ctx:
        return flat_ctx[key]
    key_dot = key.replace("_", ".")
    if key_dot in flat_ctx:
        return flat_ctx[key_dot]
    key_us = key.replace(".", "_")
    return flat_ctx.get(key_us)


def _evaluate_single_feature(expr: str, flat_ctx: Dict[str, Any]) -> bool:
    expr = expr.strip()
    if not expr:
        return False

    if " in [" in expr:
        key, rest = expr.split(" in [", 1)
        options_block = rest.strip()
        if options_block.endswith("]"):
            options_block = options_block[:-1]
        options = [opt.strip() for opt in options_block.split(",") if opt.strip()]
        value = _lookup(flat_ctx, key.strip())
        if value is None:
            return False
        if isinstance(value, (list, tuple, set)):
            value_set = {str(v).lower() for v in value}
        else:
            value_set = {str(value).lower()}
        options_set = {opt.lower() for opt in options}
        return not options_set.isdisjoint(value_set)

    for op in (">=", "<=", ">", "<"):
        if op in expr:
            key, val = expr.split(op, 1)
            value = _lookup(flat_ctx, key.strip())
            if value is None:
                return False
            value_num = _to_float(value)
            target_num = _to_float(val.strip())
            if value_num is None or target_num is None:
                return False
            if op == ">=" and not (value_num >= target_num):
                return False
            if op == "<=" and not (value_num <= target_num):
                return False
            if op == ">" and not (value_num > target_num):
                return False
            if op == "<" and not (value_num < target_num):
                return False
            return True

    if "=" in expr:
        key, val = expr.split("=", 1)
        value = _lookup(flat_ctx, key.strip())
        if value is None:
            return False
        if isinstance(value, (list, tuple, set)):
            return val.strip().lower() in {str(v).lower() for v in value}
        return str(value).lower() == val.strip().lower()

    return False


def _evaluate_prior(feature_expr: str, flat_ctx: Dict[str, Any]) -> bool:
    parts = [part.strip() for part in feature_expr.split("&") if part.strip()]
    if not parts:
        return False
    return all(_evaluate_single_feature(part, flat_ctx) for part in parts)


def score_with_priors(
    crl: Dict[str, Any],
    family: str,
    candidate: Dict[str, Any],
    features: Dict[str, Any],
) -> float:
    fam = (crl.get("families") or {}).get(family, {})
    priors: Iterable[Dict[str, Any]] = fam.get("priors") or []
    flat_ctx = _flatten_context(_augment_features(features))

    bump = 0.0
    for prior in priors:
        feature_expr = prior.get("feature")
        if not feature_expr:
            continue
        if _evaluate_prior(str(feature_expr), flat_ctx):
            effect = str(prior.get("effect", "boost")).lower()
            if effect == "boost":
                bump += 0.4
            elif effect == "slight_boost":
                bump += 0.2
            elif effect == "penalize":
                bump -= 0.3
    return bump


def _materialize_recipe(recipe: Dict[str, Any]) -> Dict[str, Any]:
    concrete: Dict[str, Any] = {}
    for key, value in (recipe or {}).items():
        if isinstance(value, list) and value:
            concrete[key] = value[0]
        else:
            concrete[key] = value
    return concrete


def recommend(
    crl: Dict[str, Any],
    family: str,
    features: Dict[str, Any],
    *,
    job_ctx: Optional[Dict[str, Any]] = None,
    top_n: int = 3,
    include_all: bool = False,
) -> List[Dict[str, Any]]:
    augmented_features = _augment_features(features)
    combined_context = _combine_contexts(augmented_features, job_ctx)
    playbooks = select_playbooks(crl, family, combined_context, fallback_to_all=include_all)
    if not playbooks:
        return []

    defaults = crl.get("defaults") or {}
    fam_defaults = (crl.get("families") or {}).get(family, {}).get("defaults") or {}
    base_defaults = _merge_nested(defaults, fam_defaults)

    guard_context = _merge_nested({}, combined_context)

    results: List[Dict[str, Any]] = []
    for pb in playbooks:
        guard_res = apply_guards(crl, family, pb, guard_context)
        prior_bump = score_with_priors(crl, family, pb, augmented_features)
        score = 1.0 + guard_res.score_adjust + prior_bump
        concrete = _materialize_recipe(pb.get("recipe") or {})
        result = {
            "playbook_id": pb.get("id"),
            "name": pb.get("name"),
            "status": pb.get("status"),
            "score": round(score, 3),
            "score_breakdown": {
                "base": 1.0,
                "guard_adjust": round(guard_res.score_adjust, 3),
                "prior": round(prior_bump, 3),
            },
            "ok": guard_res.ok,
            "guard_messages": guard_res.messages,
            "patch": guard_res.patch,
            "recipe_options": pb.get("recipe"),
            "suggested_recipe": _merge_nested(base_defaults, concrete),
        }
        results.append(result)

    results.sort(key=lambda r: (not r["ok"], -r["score"]))
    return results[: max(1, top_n)] if not include_all else results

#!/usr/bin/env python3
"""Lightweight MCP-style JSON-RPC server for rule suggestions.

This server reads and writes JSON messages from stdin/stdout so it can be
invoked as a subprocess by the orchestrator. It intentionally avoids network
sockets to keep deployments hermetic.

Supported methods
=================
- ``ping`` → ``{"status": "ok"}``
- ``rules.apply``
- ``rules.merge``
- ``rules.audit``

The payload format aligns with the integration plan documented in
``plan/MCP_Server_Integration_Plan.md``.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import threading
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

LOCK = threading.Lock()


def _now_iso() -> str:
    return datetime.utcnow().replace(microsecond=0).isoformat() + "Z"


def _is_sequence(value: Any) -> bool:
    return isinstance(value, (list, tuple, set))


def _normalize_sequence(value: Any) -> List[Any]:
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, set):
        return list(value)
    return [value]


def _to_number(value: Any) -> Optional[float]:
    try:
        if isinstance(value, (int, float)):
            return float(value)
        if isinstance(value, str):
            return float(value.strip())
    except Exception:
        return None
    return None


def _match_operator(value: Any, operator: str, target: Any) -> bool:
    operator = operator.lower()
    if operator == "in":
        value_set = {_normalize_case(v) for v in _normalize_sequence(value)}
        target_set = {_normalize_case(v) for v in _normalize_sequence(target)}
        return not value_set.isdisjoint(target_set)
    if operator == "not_in":
        value_set = {_normalize_case(v) for v in _normalize_sequence(value)}
        target_set = {_normalize_case(v) for v in _normalize_sequence(target)}
        return value_set.isdisjoint(target_set)
    if operator in {">", ">=", "<", "<="}:
        value_num = _to_number(value)
        target_num = _to_number(target)
        if value_num is None or target_num is None:
            return False
        if operator == ">":
            return value_num > target_num
        if operator == ">=":
            return value_num >= target_num
        if operator == "<":
            return value_num < target_num
        if operator == "<=":
            return value_num <= target_num
    if operator in {"eq", "equals"}:
        return _normalize_case(value) == _normalize_case(target)
    return False


def _normalize_case(value: Any) -> str:
    return str(value).strip().lower()


def _match_condition(candidate_value: Any, rule_condition: Any) -> bool:
    if isinstance(rule_condition, dict):
        operator_keys = {"in", "not_in", ">", "<", ">=", "<=", "eq", "equals"}
        if any(_normalize_case(k) in operator_keys for k in rule_condition.keys()):
            return all(_match_operator(candidate_value, k, rule_condition[k]) for k in rule_condition)
        if not isinstance(candidate_value, dict):
            return False
        return all(_match_condition(candidate_value.get(k), v) for k, v in rule_condition.items())

    if isinstance(rule_condition, list):
        if isinstance(candidate_value, (list, tuple, set)):
            rule_norm = {_normalize_case(v) for v in rule_condition}
            cand_norm = {_normalize_case(v) for v in candidate_value}
            return not rule_norm.isdisjoint(cand_norm)
        return _normalize_case(candidate_value) in {_normalize_case(v) for v in rule_condition}

    if isinstance(candidate_value, (list, tuple, set)):
        return _normalize_case(rule_condition) in {_normalize_case(v) for v in candidate_value}

    if candidate_value is None:
        return False
    return _normalize_case(candidate_value) == _normalize_case(rule_condition)


def _match_when(when: Dict[str, Any], features: Dict[str, Any]) -> bool:
    if not when:
        return True
    for key, cond in when.items():
        if key not in features:
            return False
        if not _match_condition(features.get(key), cond):
            return False
    return True


def _materialize_recipe(recipe: Dict[str, Any]) -> Dict[str, Any]:
    concrete: Dict[str, Any] = {}
    for key, value in (recipe or {}).items():
        if _is_sequence(value):
            seq = _normalize_sequence(value)
            concrete[key] = seq[0] if seq else None
        else:
            concrete[key] = value
    return concrete


@dataclass
class RulesRepository:
    rules_path: Path
    data: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.reload()

    def reload(self) -> None:
        with LOCK:
            self.data = json.loads(self.rules_path.read_text(encoding="utf-8"))
            if "playbooks" not in self.data:
                raise ValueError(f"Rules file '{self.rules_path}' is missing 'playbooks'.")

    @property
    def family_name(self) -> str:
        return str(self.data.get("reaction") or self.rules_path.stem)

    def playbooks(self) -> Iterable[Dict[str, Any]]:
        return self.data.get("playbooks", [])

    def defaults(self) -> Dict[str, Any]:
        defaults = dict(self.data.get("defaults") or {})
        defaults.setdefault("reaction", self.family_name)
        return defaults

    def find_playbook(self, rule_id: str) -> Optional[Dict[str, Any]]:
        for pb in self.playbooks():
            if str(pb.get("id")) == str(rule_id):
                return pb
        return None

    def apply(self, features: Dict[str, Any], max_suggestions: int) -> Dict[str, Any]:
        matched: List[Dict[str, Any]] = []
        for pb in self.playbooks():
            if _match_when(pb.get("when") or {}, features):
                matched.append(pb)

        suggestions: List[Dict[str, Any]] = []
        for pb in matched[: max_suggestions or 5]:
            suggestion = {
                "playbook_id": pb.get("id"),
                "name": pb.get("name"),
                "status": pb.get("status", "candidate"),
                "score": 1.0,
                "conditions": _materialize_recipe(pb.get("recipe") or {}),
                "notes": pb.get("notes", []),
                "diagnostics": pb.get("diagnostics", []),
            }
            suggestions.append(suggestion)

        return {
            "input": {
                "features": features,
                "max_suggestions": max_suggestions,
            },
            "defaults": self.defaults(),
            "matched_playbooks": [pb.get("id") for pb in matched],
            "suggestions": suggestions,
        }

    def merge(self, candidates: Dict[str, Any]) -> Dict[str, Any]:
        if not isinstance(candidates, dict):
            raise ValueError("'candidates' must be a dict")
        updated = False
        with LOCK:
            playbooks = list(self.data.get("playbooks") or [])
            candidate_playbooks = candidates.get("playbooks") or []
            if not candidate_playbooks:
                raise ValueError("Merge payload must contain playbooks")

            existing_index = {str(pb.get("id")): idx for idx, pb in enumerate(playbooks) if pb.get("id")}
            for pb in candidate_playbooks:
                rule_id = str(pb.get("id")) if pb.get("id") is not None else None
                if not rule_id:
                    raise ValueError("Each playbook must define an 'id'")
                if rule_id in existing_index:
                    playbooks[existing_index[rule_id]] = pb
                else:
                    playbooks.append(pb)
                updated = True

            if candidates.get("guards"):
                guards = list(self.data.get("guards") or [])
                guards.extend(candidates["guards"])
                self.data["guards"] = guards
                updated = True

            if candidates.get("priors"):
                priors = list(self.data.get("priors") or [])
                priors.extend(candidates["priors"])
                self.data["priors"] = priors
                updated = True

            if not updated:
                return {"message": "no_changes"}

            self.data["playbooks"] = playbooks
            self.data["version"] = _now_iso()
            self.rules_path.write_text(json.dumps(self.data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        return {
            "message": "merge_ok",
            "rules_path": str(self.rules_path),
            "version": self.data.get("version"),
        }

    def audit(self, rule_id: str) -> Dict[str, Any]:
        pb = self.find_playbook(rule_id)
        if pb is None:
            raise KeyError(f"Rule '{rule_id}' not found")
        return {
            "section": "playbooks",
            "rule": pb,
            "defaults": self.defaults(),
        }


class JsonRpcServer:
    def __init__(self, repo: RulesRepository, log_level: str = "info") -> None:
        self.repo = repo
        self.log_level = log_level.lower()

    def log(self, message: str) -> None:
        if self.log_level in {"debug", "info"}:
            ts = _now_iso()
            print(json.dumps({"ts": ts, "level": self.log_level, "message": message}), file=sys.stderr)

    def handle_request(self, request: Dict[str, Any]) -> Dict[str, Any]:
        method = request.get("method")
        params = request.get("params") or {}
        request_id = request.get("id")

        start = time.time()
        try:
            if method == "ping":
                result = {"status": "ok", "ts": _now_iso()}
            elif method == "rules.apply":
                features = params.get("features") or {}
                max_suggestions = int(params.get("max_suggestions", 5))
                result = self.repo.apply(features, max_suggestions)
                result["latency_ms"] = int((time.time() - start) * 1000)
            elif method == "rules.merge":
                candidates = params.get("candidates") or {}
                result = self.repo.merge(candidates)
            elif method == "rules.audit":
                rule_id = params.get("rule_id")
                if not rule_id:
                    raise ValueError("'rule_id' is required")
                result = self.repo.audit(rule_id)
            else:
                raise ValueError(f"Unknown method '{method}'")
            response = {"jsonrpc": "2.0", "id": request_id, "result": result}
        except Exception as exc:  # pylint: disable=broad-except
            response = {
                "jsonrpc": "2.0",
                "id": request_id,
                "error": {
                    "message": str(exc),
                    "type": exc.__class__.__name__,
                },
            }
        finally:
            if method and method != "ping":
                duration = int((time.time() - start) * 1000)
                self.log(f"method={method} duration_ms={duration}")
        return response

    def serve_forever(self) -> None:
        for line in sys.stdin:
            line = line.strip()
            if not line:
                continue
            try:
                request = json.loads(line)
            except json.JSONDecodeError as exc:
                error = {"jsonrpc": "2.0", "id": None, "error": {"message": f"Invalid JSON: {exc}"}}
                print(json.dumps(error), flush=True)
                continue
            response = self.handle_request(request)
            print(json.dumps(response, ensure_ascii=False), flush=True)


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the MCP rules JSON-RPC server")
    parser.add_argument("--rules", required=True, help="Path to the rules JSON file")
    parser.add_argument("--log-level", default=os.environ.get("MCP_LOG_LEVEL", "info"), help="Logging level")
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> None:
    args = parse_args(argv)
    rules_path = Path(args.rules).resolve()
    if not rules_path.exists():
        raise SystemExit(f"Rules file not found: {rules_path}")
    repo = RulesRepository(rules_path)
    server = JsonRpcServer(repo, log_level=args.log_level)
    server.log(f"server_start rules={rules_path}")
    server.serve_forever()


if __name__ == "__main__":
    main()

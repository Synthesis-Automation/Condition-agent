#!/usr/bin/env python3
"""
Reliability benchmark gate for reaction-key/type quality reports.

The gate can validate an absolute threshold and optional regression delta against
another baseline report.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


def _rate(report: Dict[str, Any], key: str) -> Optional[float]:
    payload = report.get(key)
    if not isinstance(payload, dict):
        return None
    try:
        return float(payload.get("rate"))
    except Exception:
        return None


def evaluate_gate(
    current: Dict[str, Any],
    *,
    max_unknown_rate: float,
    max_low_quality_rate: float,
    max_empty_key_rate: float,
    baseline: Optional[Dict[str, Any]] = None,
    max_unknown_delta: Optional[float] = None,
    max_low_quality_delta: Optional[float] = None,
    max_empty_key_delta: Optional[float] = None,
) -> Tuple[bool, List[str], Dict[str, Any]]:
    errors: List[str] = []
    details: Dict[str, Any] = {"current": {}, "baseline": {}}

    thresholds = {
        "unknown_reaction_type": max_unknown_rate,
        "low_reaction_key_quality": max_low_quality_rate,
        "empty_reaction_key": max_empty_key_rate,
    }
    for key, threshold in thresholds.items():
        value = _rate(current, key)
        details["current"][key] = value
        if value is None:
            errors.append(f"missing_metric:{key}")
            continue
        if value > float(threshold):
            errors.append(f"threshold_failed:{key}:{value:.4f}>{float(threshold):.4f}")

    if baseline is not None:
        delta_limits = {
            "unknown_reaction_type": max_unknown_delta,
            "low_reaction_key_quality": max_low_quality_delta,
            "empty_reaction_key": max_empty_key_delta,
        }
        deltas: Dict[str, float] = {}
        for key, limit in delta_limits.items():
            current_rate = _rate(current, key)
            baseline_rate = _rate(baseline, key)
            details["baseline"][key] = baseline_rate
            if current_rate is None or baseline_rate is None:
                errors.append(f"missing_baseline_metric:{key}")
                continue
            delta = float(current_rate) - float(baseline_rate)
            deltas[key] = round(delta, 4)
            if limit is not None and delta > float(limit):
                errors.append(
                    f"delta_failed:{key}:{delta:.4f}>{float(limit):.4f}"
                )
        details["deltas"] = deltas

    return len(errors) == 0, errors, details


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Gate reaction reliability report against thresholds")
    parser.add_argument("--report", required=True, help="Current report JSON path")
    parser.add_argument("--baseline", default="", help="Optional baseline report JSON path")

    parser.add_argument("--max-unknown-rate", type=float, default=0.30)
    parser.add_argument("--max-low-quality-rate", type=float, default=0.02)
    parser.add_argument("--max-empty-key-rate", type=float, default=0.01)

    parser.add_argument("--max-unknown-delta", type=float, default=None)
    parser.add_argument("--max-low-quality-delta", type=float, default=None)
    parser.add_argument("--max-empty-key-delta", type=float, default=None)

    parser.add_argument("--output-json", default="", help="Optional gate result JSON path")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    report_path = Path(args.report).expanduser().resolve()
    if not report_path.exists():
        raise FileNotFoundError(f"Report not found: {report_path}")
    current = json.loads(report_path.read_text(encoding="utf-8"))

    baseline = None
    baseline_path = str(args.baseline or "").strip()
    if baseline_path:
        baseline_file = Path(baseline_path).expanduser().resolve()
        if not baseline_file.exists():
            raise FileNotFoundError(f"Baseline report not found: {baseline_file}")
        baseline = json.loads(baseline_file.read_text(encoding="utf-8"))

    passed, errors, details = evaluate_gate(
        current,
        max_unknown_rate=float(args.max_unknown_rate),
        max_low_quality_rate=float(args.max_low_quality_rate),
        max_empty_key_rate=float(args.max_empty_key_rate),
        baseline=baseline,
        max_unknown_delta=args.max_unknown_delta,
        max_low_quality_delta=args.max_low_quality_delta,
        max_empty_key_delta=args.max_empty_key_delta,
    )

    payload = {
        "passed": passed,
        "errors": errors,
        "details": details,
        "thresholds": {
            "max_unknown_rate": args.max_unknown_rate,
            "max_low_quality_rate": args.max_low_quality_rate,
            "max_empty_key_rate": args.max_empty_key_rate,
            "max_unknown_delta": args.max_unknown_delta,
            "max_low_quality_delta": args.max_low_quality_delta,
            "max_empty_key_delta": args.max_empty_key_delta,
        },
    }
    print(json.dumps(payload, indent=2, sort_keys=True))

    output = str(args.output_json or "").strip()
    if output:
        out_path = Path(output).expanduser().resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")

    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())

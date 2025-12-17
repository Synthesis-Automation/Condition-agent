"""
Attachment-point SMARTS templating evaluator for the Suzuki/Buchwald POC.

Implements the evaluation workflow described in:
`chemtools/taxonomy/new/CODEX_EVAL_AttachmentPoint_Templating.md`.

Key responsibilities:
  - Parse the template spec (`smarts_templates...json`)
  - Expand templates deterministically into SMARTS strings
  - Compile generated SMARTS (via the centralized cache)
  - Cross-check generated patterns against the runtime atomic feature file
  - Run small functional match checks on representative SMILES
  - Render a short Markdown evaluation report
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence

from chemtools.util.rdkit_helpers import parse_smiles, rdkit_available
from chemtools.util.smarts_cache import compile_smarts


class TemplateEvalError(ValueError):
    """Raised when the template spec or evaluation inputs are invalid."""


_PLACEHOLDER_RE = re.compile(r"{([a-zA-Z_][a-zA-Z0-9_]*)}")


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _infer_placeholders(pattern: str) -> list[str]:
    placeholders: list[str] = []
    seen: set[str] = set()
    for match in _PLACEHOLDER_RE.finditer(pattern):
        name = match.group(1)
        if name in seen:
            continue
        placeholders.append(name)
        seen.add(name)
    return placeholders


@dataclass(frozen=True)
class ExpansionSpec:
    feature_token: str
    template: str
    fragments: tuple[str, ...]
    expected_smarts: Optional[str] = None


@dataclass(frozen=True)
class GeneratedPattern:
    feature_token: str
    smarts: str
    template: str
    fragments: tuple[str, ...]
    placeholder_to_fragment: Dict[str, str]


@dataclass(frozen=True)
class TemplateSpec:
    fragments: Dict[str, str]
    templates: Dict[str, str]
    expansions: tuple[ExpansionSpec, ...]


@dataclass(frozen=True)
class CompileFailure:
    feature_token: str
    smarts: str
    reason: str


@dataclass(frozen=True)
class FunctionalTestCase:
    name: str
    smiles: str
    expected_tokens: tuple[str, ...]
    xfail: bool = False
    note: str = ""


@dataclass(frozen=True)
class FunctionalTestResult:
    case: FunctionalTestCase
    passed: bool
    missing_tokens: tuple[str, ...]
    unexpected_tokens: tuple[str, ...]


@dataclass(frozen=True)
class TemplateEvalResult:
    fragments_count: int
    templates_count: int
    expansions_count: int
    generated_count: int
    compile_successes: int
    compile_failures: tuple[CompileFailure, ...]
    generated_token_duplicates: tuple[str, ...]
    generated_smarts_duplicates: tuple[str, ...]
    tokens_missing_in_atomic: tuple[str, ...]
    smarts_missing_in_atomic: tuple[str, ...]
    functional_tests: tuple[FunctionalTestResult, ...]

    @property
    def compile_failure_count(self) -> int:
        return len(self.compile_failures)


def load_template_spec(path: Path) -> TemplateSpec:
    """
    Parse `smarts_templates.suzuki_buchwald.poc.json` into a normalized TemplateSpec.
    """
    payload = _read_json(path)
    if not isinstance(payload, dict):
        raise TemplateEvalError(f"Expected object in {path}, got {type(payload).__name__}")

    raw_frags = payload.get("fragments", {})
    raw_templates = payload.get("templates", {})
    raw_generated = payload.get("generated_atomic_features", [])

    if not isinstance(raw_frags, dict) or not isinstance(raw_templates, dict) or not isinstance(raw_generated, list):
        raise TemplateEvalError("Invalid template spec schema (fragments/templates/generated_atomic_features).")

    fragments: dict[str, str] = {}
    for name, meta in raw_frags.items():
        if not isinstance(name, str):
            continue
        if isinstance(meta, dict) and isinstance(meta.get("smarts"), str):
            fragments[name] = meta["smarts"]

    templates: dict[str, str] = {}
    for name, meta in raw_templates.items():
        if not isinstance(name, str):
            continue
        if isinstance(meta, dict) and isinstance(meta.get("pattern"), str):
            templates[name] = meta["pattern"]

    expansions: list[ExpansionSpec] = []
    for entry in raw_generated:
        if not isinstance(entry, dict):
            continue
        token = str(entry.get("token", "")).strip()
        template = str(entry.get("template", "")).strip()
        frags = entry.get("fragments", [])
        expected_smarts = entry.get("smarts")
        if not token or not template:
            continue
        if not isinstance(frags, list) or not all(isinstance(f, str) for f in frags):
            raise TemplateEvalError(f"{token}: fragments must be a list[str]")
        expansions.append(
            ExpansionSpec(
                feature_token=token,
                template=template,
                fragments=tuple(frags),
                expected_smarts=str(expected_smarts).strip() if isinstance(expected_smarts, str) else None,
            )
        )

    # Deterministic order: sort by token, then template, then fragment tuple.
    expansions.sort(key=lambda e: (e.feature_token, e.template, e.fragments))

    return TemplateSpec(fragments=fragments, templates=templates, expansions=tuple(expansions))


def expand_templates(
    fragments: Mapping[str, str],
    templates: Mapping[str, str],
    expansions: Sequence[ExpansionSpec],
) -> list[GeneratedPattern]:
    """
    Deterministically expand templates into SMARTS patterns.

    Args:
        fragments: fragment id -> SMARTS snippet (may include attachment labels like ":1")
        templates: template id -> format string (e.g. "{core}{lg}")
        expansions: expansion plan entries (token + template + fragment ids)
    """
    generated: list[GeneratedPattern] = []
    for exp in expansions:
        pattern = templates.get(exp.template)
        if not pattern:
            raise TemplateEvalError(f"{exp.feature_token}: unknown template '{exp.template}'")

        placeholders = _infer_placeholders(pattern)
        if len(placeholders) != len(exp.fragments):
            raise TemplateEvalError(
                f"{exp.feature_token}: placeholder/fragment count mismatch "
                f"(placeholders={placeholders}, fragments={list(exp.fragments)})"
            )

        placeholder_to_fragment: dict[str, str] = {}
        placeholder_to_smarts: dict[str, str] = {}
        for placeholder, frag_id in zip(placeholders, exp.fragments, strict=True):
            frag_smarts = fragments.get(frag_id)
            if frag_smarts is None:
                raise TemplateEvalError(f"{exp.feature_token}: unknown fragment '{frag_id}'")
            placeholder_to_fragment[placeholder] = frag_id
            placeholder_to_smarts[placeholder] = frag_smarts

        try:
            smarts = pattern.format(**placeholder_to_smarts)
        except Exception as e:
            raise TemplateEvalError(f"{exp.feature_token}: failed to expand template '{exp.template}': {e}") from e

        generated.append(
            GeneratedPattern(
                feature_token=exp.feature_token,
                smarts=smarts,
                template=exp.template,
                fragments=exp.fragments,
                placeholder_to_fragment=placeholder_to_fragment,
            )
        )

    # Deterministic output ordering, independent of input order.
    generated.sort(key=lambda g: (g.feature_token, g.template, g.fragments, g.smarts))
    return generated


def compile_generated_smarts(generated: Sequence[GeneratedPattern]) -> tuple[int, tuple[CompileFailure, ...]]:
    """
    Compile generated SMARTS (requires RDKit).
    """
    if not rdkit_available():
        raise TemplateEvalError("RDKit is not available; cannot compile SMARTS patterns.")

    failures: list[CompileFailure] = []
    successes = 0
    for entry in generated:
        try:
            patt = compile_smarts(entry.smarts, validate=True)
            if patt is None:
                failures.append(
                    CompileFailure(feature_token=entry.feature_token, smarts=entry.smarts, reason="compiled to None")
                )
            else:
                successes += 1
        except Exception as e:
            failures.append(CompileFailure(feature_token=entry.feature_token, smarts=entry.smarts, reason=str(e)))
    return successes, tuple(failures)


def load_atomic_feature_smarts(path: Path) -> Dict[str, tuple[str, ...]]:
    """
    Load the runtime atomic feature file as token -> tuple(smarts_any...).
    """
    payload = _read_json(path)
    if not isinstance(payload, list):
        raise TemplateEvalError(f"Expected list in {path}, got {type(payload).__name__}")

    out: dict[str, tuple[str, ...]] = {}
    for entry in payload:
        if not isinstance(entry, dict):
            continue
        token = str(entry.get("token", "")).strip()
        detect = entry.get("detect") or {}
        if not token or not isinstance(detect, dict):
            continue
        smarts_any = detect.get("smarts_any") or []
        if not isinstance(smarts_any, list):
            continue
        out[token] = tuple(str(s) for s in smarts_any if isinstance(s, str))

    return out


def _find_duplicates(items: Iterable[str]) -> tuple[str, ...]:
    seen: set[str] = set()
    dupes: set[str] = set()
    for item in items:
        if item in seen:
            dupes.add(item)
        seen.add(item)
    return tuple(sorted(dupes))


def cross_check_against_atomic_features(
    generated: Sequence[GeneratedPattern], atomic_smarts: Mapping[str, Sequence[str]]
) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    """
    Cross-check generated patterns vs runtime atomic features.

    Returns:
        (missing_tokens, missing_smarts, duplicate_tokens, duplicate_smarts)
    """
    tokens = [g.feature_token for g in generated]
    duplicate_tokens = _find_duplicates(tokens)

    smarts_to_tokens: dict[str, set[str]] = {}
    for g in generated:
        smarts_to_tokens.setdefault(g.smarts, set()).add(g.feature_token)
    duplicate_smarts = tuple(sorted([s for s, toks in smarts_to_tokens.items() if len(toks) > 1]))

    missing_tokens: list[str] = []
    missing_smarts: list[str] = []
    for g in generated:
        if g.feature_token not in atomic_smarts:
            missing_tokens.append(g.feature_token)
            continue
        if g.smarts not in set(atomic_smarts[g.feature_token]):
            missing_smarts.append(g.feature_token)

    return (
        tuple(sorted(set(missing_tokens))),
        tuple(sorted(set(missing_smarts))),
        duplicate_tokens,
        duplicate_smarts,
    )


def functional_test_cases() -> tuple[FunctionalTestCase, ...]:
    """
    Default functional test cases for this POC.

    Note: The evaluation brief uses "aryl_*" token names in examples; this POC uses
    "aromatic_*" tokens.
    """
    return (
        FunctionalTestCase(
            name="bromobenzene electrophile",
            smiles="Brc1ccccc1",
            expected_tokens=("aromatic_bromide_present",),
        ),
        FunctionalTestCase(
            name="chlorobenzene electrophile",
            smiles="Clc1ccccc1",
            expected_tokens=("aromatic_chloride_present",),
        ),
        FunctionalTestCase(
            name="iodobenzene electrophile",
            smiles="Ic1ccccc1",
            expected_tokens=("aromatic_iodide_present",),
        ),
        FunctionalTestCase(
            name="phenylboronic acid",
            smiles="OB(O)c1ccccc1",
            expected_tokens=("aromatic_boronic_acid_present",),
        ),
        FunctionalTestCase(
            name="phenyl Bpin (one common SMILES)",
            smiles="c1ccccc1B1OC(C)(C)C(C)(C)O1",
            expected_tokens=("aromatic_bpin_present",),
        ),
        FunctionalTestCase(
            name="vinyl bromide (optional)",
            smiles="C=CBr",
            expected_tokens=("vinyl_bromide_present",),
        ),
    )


def run_functional_match_tests(
    generated: Sequence[GeneratedPattern], *, strict: bool = False
) -> tuple[FunctionalTestResult, ...]:
    """
    Run small molecule match checks using RDKit HasSubstructMatch.

    Args:
        generated: Generated atomic SMARTS patterns.
        strict: If True, treat any unexpected matches as failures.
                If False (default), only missing expected tokens fail; unexpected
                matches are still reported in the result for review.
    """
    if not rdkit_available():
        raise TemplateEvalError("RDKit is not available; cannot run functional match tests.")

    compiled: dict[str, Any] = {}
    for g in generated:
        patt = compile_smarts(g.smarts, validate=True)
        if patt is None:
            raise TemplateEvalError(f"Generated SMARTS compiled to None for {g.feature_token}: {g.smarts}")
        compiled[g.feature_token] = patt

    results: list[FunctionalTestResult] = []
    for case in functional_test_cases():
        mol = parse_smiles(case.smiles)
        if mol is None:
            # Treat bad SMILES as failure, unless xfail.
            results.append(
                FunctionalTestResult(
                    case=case,
                    passed=bool(case.xfail),
                    missing_tokens=case.expected_tokens,
                    unexpected_tokens=tuple(),
                )
            )
            continue

        hit_tokens = {token for token, patt in compiled.items() if mol.HasSubstructMatch(patt)}
        expected = set(case.expected_tokens)
        missing = tuple(sorted(expected - hit_tokens))
        unexpected = tuple(sorted(hit_tokens - expected))
        passed = (not missing) and ((not unexpected) if strict else True)
        if case.xfail:
            passed = True
        results.append(
            FunctionalTestResult(
                case=case,
                passed=passed,
                missing_tokens=missing,
                unexpected_tokens=unexpected,
            )
        )

    return tuple(results)


def evaluate_attachment_point_templating(
    *,
    templates_path: Path,
    atomic_features_path: Path,
) -> TemplateEvalResult:
    """
    End-to-end evaluation: expand -> compile -> cross-check -> functional tests.
    """
    spec = load_template_spec(templates_path)
    generated = expand_templates(spec.fragments, spec.templates, spec.expansions)

    compile_successes, compile_failures = compile_generated_smarts(generated)

    atomic_smarts = load_atomic_feature_smarts(atomic_features_path)
    missing_tokens, missing_smarts, dup_tokens, dup_smarts = cross_check_against_atomic_features(
        generated, atomic_smarts
    )

    functional = run_functional_match_tests(generated)

    return TemplateEvalResult(
        fragments_count=len(spec.fragments),
        templates_count=len(spec.templates),
        expansions_count=len(spec.expansions),
        generated_count=len(generated),
        compile_successes=compile_successes,
        compile_failures=compile_failures,
        generated_token_duplicates=dup_tokens,
        generated_smarts_duplicates=dup_smarts,
        tokens_missing_in_atomic=missing_tokens,
        smarts_missing_in_atomic=missing_smarts,
        functional_tests=functional,
    )


def render_markdown_report(result: TemplateEvalResult) -> str:
    """
    Render a short Markdown report with requested metrics.
    """
    lines: list[str] = []
    lines.append("# TEMPLATE_EVAL_REPORT")
    lines.append("")
    lines.append("## Build integrity")
    lines.append("")
    lines.append(f"- fragments: {result.fragments_count}")
    lines.append(f"- templates: {result.templates_count}")
    lines.append(f"- expansions: {result.expansions_count}")
    lines.append(f"- generated patterns: {result.generated_count}")
    lines.append(f"- SMARTS compile successes: {result.compile_successes}")
    lines.append(f"- SMARTS compile failures: {result.compile_failure_count}")
    if result.compile_failures:
        lines.append("")
        lines.append("### Compile failures")
        lines.append("")
        for f in result.compile_failures:
            lines.append(f"- {f.feature_token}: `{f.smarts}` ({f.reason})")

    lines.append("")
    lines.append("## Consistency with atomic features")
    lines.append("")
    token_ok = result.generated_count - len(result.tokens_missing_in_atomic)
    smarts_ok = result.generated_count - len(result.smarts_missing_in_atomic)
    token_pct = (100.0 * token_ok / result.generated_count) if result.generated_count else 0.0
    smarts_pct = (100.0 * smarts_ok / result.generated_count) if result.generated_count else 0.0
    lines.append(f"- generated tokens found in atomic features: {token_ok}/{result.generated_count} ({token_pct:.1f}%)")
    lines.append(f"- generated SMARTS found in atomic features: {smarts_ok}/{result.generated_count} ({smarts_pct:.1f}%)")
    lines.append(f"- duplicate generated tokens: {len(result.generated_token_duplicates)}")
    lines.append(f"- duplicate generated SMARTS (under multiple tokens): {len(result.generated_smarts_duplicates)}")

    if result.tokens_missing_in_atomic:
        lines.append("")
        lines.append("### Missing tokens in atomic features")
        lines.append("")
        for t in result.tokens_missing_in_atomic:
            lines.append(f"- {t}")

    if result.smarts_missing_in_atomic:
        lines.append("")
        lines.append("### Missing SMARTS in atomic features")
        lines.append("")
        for t in result.smarts_missing_in_atomic:
            lines.append(f"- {t}")

    if result.generated_token_duplicates:
        lines.append("")
        lines.append("### Duplicate generated tokens")
        lines.append("")
        for t in result.generated_token_duplicates:
            lines.append(f"- {t}")

    if result.generated_smarts_duplicates:
        lines.append("")
        lines.append("### Duplicate SMARTS strings")
        lines.append("")
        for s in result.generated_smarts_duplicates:
            lines.append(f"- `{s}`")

    lines.append("")
    lines.append("## Functional match tests")
    lines.append("")
    passed = sum(1 for r in result.functional_tests if r.passed)
    total = len(result.functional_tests)
    lines.append(f"- passed: {passed}/{total}")
    lines.append("")
    lines.append("| test | smiles | expected tokens | status | notes |")
    lines.append("|---|---|---|---|---|")
    for r in result.functional_tests:
        case = r.case
        status = "PASS" if r.passed else "FAIL"
        notes: list[str] = []
        if r.missing_tokens:
            notes.append(f"missing={list(r.missing_tokens)}")
        if r.unexpected_tokens:
            notes.append(f"unexpected={list(r.unexpected_tokens)}")
        if case.xfail:
            notes.append("xfail")
        if case.note:
            notes.append(case.note)
        notes_str = "; ".join(notes)
        lines.append(
            f"| {case.name} | `{case.smiles}` | {', '.join(case.expected_tokens)} | {status} | {notes_str} |"
        )

    return "\n".join(lines) + "\n"


def write_markdown_report(result: TemplateEvalResult, path: Path) -> None:
    path.write_text(render_markdown_report(result), encoding="utf-8")


__all__ = [
    "TemplateEvalError",
    "ExpansionSpec",
    "GeneratedPattern",
    "TemplateSpec",
    "CompileFailure",
    "FunctionalTestCase",
    "FunctionalTestResult",
    "TemplateEvalResult",
    "load_template_spec",
    "expand_templates",
    "compile_generated_smarts",
    "load_atomic_feature_smarts",
    "cross_check_against_atomic_features",
    "run_functional_match_tests",
    "evaluate_attachment_point_templating",
    "render_markdown_report",
    "write_markdown_report",
]

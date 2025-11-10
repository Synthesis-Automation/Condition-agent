"""Utilities for drafting protocol JSON files from freeform inputs."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
import json
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from chemtools.smiles import normalize_reaction

try:
    from llmtools.clients import LLMClient
except Exception:  # pragma: no cover - optional dependency
    LLMClient = None

LLM_AVAILABLE = LLMClient is not None


@dataclass
class ProtocolDraftInput:
    """Input payload for `generate_protocol_draft`."""

    reaction_smiles: str
    procedure_text: str
    title: Optional[str] = None
    metadata_id: Optional[str] = None
    reaction_family: Optional[str] = None
    tags: Optional[List[str]] = None
    solvent_hint: Optional[str] = None
    base_hint: Optional[str] = None
    catalyst_hint: Optional[str] = None
    additive_hint: Optional[str] = None
    temperature_hint: Optional[str] = None
    time_hint: Optional[str] = None
    scale_mmol: float = 2.0
    atmosphere_hint: Optional[str] = None
    use_llm: bool = False
    llm_provider: Optional[str] = None
    llm_model: Optional[str] = None


@dataclass
class ProtocolDraftResult:
    """Result payload returned by `generate_protocol_draft`."""

    draft: Dict[str, Any]
    addition_sequence: List[Dict[str, Any]] = field(default_factory=list)
    issues: List[str] = field(default_factory=list)
    llm_used: bool = False
    llm_metadata: Dict[str, Any] = field(default_factory=dict)


class ProtocolDraftError(ValueError):
    """Raised when draft generation fails."""


_ACTION_KEYWORDS: List[Tuple[str, str]] = [
    (r"\bcharge\b", "Charge"),
    (r"\badd\b", "Add"),
    (r"\bintroduce\b", "Introduce"),
    (r"\bcombine\b", "Combine"),
    (r"\bstir\b", "Stir"),
    (r"\bheated?\b", "Heat"),
    (r"\bcool\b", "Cool"),
    (r"\bquench\b", "Quench"),
    (r"\bwash\b", "Wash"),
    (r"\bdry\b", "Dry"),
    (r"\bfilter\b", "Filter"),
    (r"\bpurif(y|ication)\b", "Purify"),
]

_ROLE_HINTS = {
    "solvent": {"solvent", "dmf", "thf", "dcm", "acetonitrile", "methanol", "ethanol"},
    "base": {"base", "carbonate", "amine", "amine-base", "hydroxide"},
    "ligand": {"ligand", "xphos", "binap", "dppf"},
    "metal_catalyst": {"palladium", "nickel", "copper", "catalyst", "pd(", "ni(", "cu("},
    "additive": {"additive", "salt", "halide"},
}

LLM_PROTOCOL_SYSTEM = (
    "You are an expert synthetic chemist. Extract a structured protocol from the user's "
    "procedure text. Respond with compact JSON only, no prose. The JSON MUST follow:\n"
    "{\n"
    '  "steps": [\n'
    '     {"order":1,"action":"Add","material":"DMF","amount":"40 mL","role":"solvent",'
    '"summary":"Add 40 mL DMF to the flask."}, ...\n'
    "  ],\n"
    '  "conditions": {"temperature_C":60,"time_h":16,"atmosphere":"N2"},\n'
    '  "workup": [{"step":"Quench","details":"Quench with sat. NH4Cl."}],\n'
    '  "issues": ["Missing reagent amounts"]\n'
    "}\n"
    "Always include all four top-level keys even if empty (use [] or {})."
)

LLM_PROTOCOL_PROMPT = """Reaction context:
- Reaction SMILES: {reaction_smiles}
- Scale (mmol): {scale_mmol}
- Optional hints: {hints}

Procedure text:
\"\"\"{procedure_text}\"\"\"

Return structured JSON as specified. Avoid commentary."""


def generate_protocol_draft(params: ProtocolDraftInput) -> ProtocolDraftResult:
    """Convert user text into a protocol draft structure."""

    reaction_smiles = (params.reaction_smiles or "").strip()
    procedure_text = (params.procedure_text or "").strip()
    if not reaction_smiles:
        raise ProtocolDraftError("Reaction SMILES is required.")
    if not procedure_text:
        raise ProtocolDraftError("Operation procedure text is required.")

    try:
        normalized = normalize_reaction(reaction_smiles)
    except Exception as exc:  # pragma: no cover - defensive
        raise ProtocolDraftError(f"Failed to normalize reaction SMILES: {exc}") from exc

    metadata_id = params.metadata_id or _build_metadata_id(params.title, normalized)
    title = params.title or f"Protocol draft for {metadata_id}"
    tags = params.tags or []

    issues: List[str] = []
    addition_sequence: List[Dict[str, Any]] = []
    workup: List[Dict[str, Any]] = []
    llm_used = False
    llm_metadata: Dict[str, Any] = {}
    llm_conditions: Dict[str, Any] = {}

    if params.use_llm:
        if not LLM_AVAILABLE:
            issues.append("LLM parser requested but llmtools is not installed/configured.")
        else:
            try:
                llm_payload, llm_metadata = _run_llm_parser(
                    params,
                    procedure_text,
                    normalized.get("normalized") or reaction_smiles,
                )
                addition_sequence = _convert_llm_steps(llm_payload.get("steps", []))
                llm_conditions = llm_payload.get("conditions") or {}
                workup = llm_payload.get("workup") or []
                llm_reported_issues = llm_payload.get("issues") or []
                if llm_reported_issues:
                    issues.extend([f"LLM: {msg}" for msg in llm_reported_issues])
                llm_used = True
            except ProtocolDraftError as exc:
                issues.append(str(exc))
            except Exception as exc:  # pragma: no cover - defensive
                issues.append(f"LLM parsing failed: {exc}")

    if not addition_sequence:
        addition_sequence = _derive_addition_sequence(procedure_text)
    if not workup:
        workup = _derive_workup(procedure_text)

    chemicals = _build_chemicals(params, addition_sequence)
    conditions = _build_conditions(params, procedure_text)
    if llm_conditions:
        conditions = _merge_conditions(conditions, llm_conditions)

    if normalized.get("errors"):
        issues.append(
            f"Unrecognized fragments in reaction SMILES: {', '.join(normalized['errors'])}"
        )
    if not addition_sequence:
        issues.append("Unable to detect addition steps. Provide explicit 'Add ...' statements.")
    if not chemicals:
        issues.append("No reagents or solvents inferred. Provide hints or structured inputs.")
    if not procedure_text:
        issues.append("Procedure text missing.")

    draft = {
        "schema_version": "2.0",
        "source_type": "protocol_draft",
        "metadata": {
            "id": metadata_id,
            "name": title,
            "version": datetime.utcnow().strftime("%Y.%m.%d.%H%M"),
            "tags": tags,
            "created_date": datetime.utcnow().isoformat(),
            "last_modified": datetime.utcnow().isoformat(),
            "generation": {
                "llm_used": llm_used,
                "llm_provider": llm_metadata.get("provider"),
                "llm_model": llm_metadata.get("model"),
                "llm_tokens": llm_metadata.get("tokens"),
            },
        },
        "reaction": {
            "reaction_smiles": normalized.get("normalized") or reaction_smiles,
            "family": params.reaction_family or "Unknown",
            "notes": _summarize_text(procedure_text),
            "addition_sequence": addition_sequence,
        },
        "reaction_setup": [
            {
                "chemicals": chemicals,
                "conditions": [conditions] if conditions else [],
            }
        ],
        "workup_and_purification": workup,
        "original_procedure": procedure_text,
    }

    return ProtocolDraftResult(
        draft=draft,
        addition_sequence=addition_sequence,
        issues=issues,
        llm_used=llm_used,
        llm_metadata=llm_metadata,
    )


def _build_metadata_id(title: Optional[str], normalized: Dict[str, Any]) -> str:
    base = title or normalized.get("normalized") or "protocol_draft"
    slug = _slugify(base)[:48] or "protocol_draft"
    return f"{slug}_{datetime.utcnow().strftime('%Y%m%d%H%M%S')}"


def _slugify(value: str) -> str:
    return re.sub(r"[^a-zA-Z0-9]+", "_", value or "").strip("_").lower()


def _summarize_text(text: str, limit: int = 240) -> str:
    clean = " ".join(text.split())
    if len(clean) <= limit:
        return clean
    return clean[: limit - 3] + "..."


def _derive_addition_sequence(text: str) -> List[Dict[str, Any]]:
    segments = _split_segments(text)
    sequence: List[Dict[str, Any]] = []
    for idx, segment in enumerate(segments, 1):
        action = _detect_action(segment)
        material = _extract_material(segment)
        amount = _extract_amount(segment)
        role = _guess_role(segment, material)
        if not action and not material:
            continue
        sequence.append(
            {
                "order": idx,
                "action": action or "Step",
                "summary": segment.strip(),
                "material": material,
                "amount": amount,
                "role": role,
            }
        )
    return sequence


def _split_segments(text: str) -> List[str]:
    lines = []
    for raw in text.splitlines():
        stripped = raw.strip(" -*\t")
        if not stripped:
            continue
        parts = re.split(r"(?<=[.!?])\s+", stripped)
        lines.extend([p.strip() for p in parts if p.strip()])
    return lines


def _detect_action(segment: str) -> Optional[str]:
    lowered = segment.lower()
    for pattern, label in _ACTION_KEYWORDS:
        if re.search(pattern, lowered):
            return label
    return None


def _extract_material(segment: str) -> Optional[str]:
    match = re.search(
        r"(?:add|charge|introduce|combine|treat(?:ed)? with)\s+([^\d,.;]+?)"
        r"(?:\s+\d|\s+at|\s+and|,|\.|;)",
        segment,
        re.IGNORECASE,
    )
    if match:
        return match.group(1).strip()
    match = re.search(r"with\s+([A-Za-z0-9\-\(\)\/,\s]+)$", segment, re.IGNORECASE)
    if match:
        return match.group(1).strip()
    return None


def _extract_amount(segment: str) -> Optional[str]:
    match = re.search(
        r"(\d+(?:\.\d+)?)\s*(mL|L|µL|g|mg|kg|mmol|mol|equiv|mol%)",
        segment,
        re.IGNORECASE,
    )
    if match:
        return f"{match.group(1)} {match.group(2)}"
    return None


def _guess_role(segment: str, material: Optional[str]) -> Optional[str]:
    target = (segment or "").lower()
    if material:
        target += " " + material.lower()
    for role, hints in _ROLE_HINTS.items():
        if any(hint in target for hint in hints):
            return role
    return None


def _build_chemicals(
    params: ProtocolDraftInput,
    addition_sequence: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    entries: List[Dict[str, Any]] = []

    def add_entry(name: Optional[str], role: str) -> None:
        label = (name or "").strip()
        if not label:
            return
        for existing in entries:
            if existing["name"].lower() == label.lower() and existing["role"] == role:
                return
        entries.append(
            {
                "name": label,
                "abbreviation": None,
                "cas": None,
                "smiles": None,
                "amount": {
                    "weight_g": None,
                    "mmol": None,
                    "volume_ml": None,
                    "equivalents": None,
                },
                "role": role,
            }
        )

    add_entry(params.solvent_hint, "solvent")
    add_entry(params.base_hint, "base")
    add_entry(params.catalyst_hint, "metal_catalyst")
    add_entry(params.additive_hint, "additive")

    for step in addition_sequence:
        if not step.get("material"):
            continue
        inferred_role = step.get("role") or "reagent"
        add_entry(step["material"], inferred_role)

    return entries


def _build_conditions(params: ProtocolDraftInput, procedure_text: str) -> Dict[str, Any]:
    conditions: Dict[str, Any] = {}

    temperature = _extract_temperature(params.temperature_hint, procedure_text)
    if temperature is not None:
        conditions["temperature_C"] = temperature
    time_h = _extract_time(params.time_hint, procedure_text)
    if time_h is not None:
        conditions["time_h"] = time_h
    atmosphere = params.atmosphere_hint or _extract_atmosphere(procedure_text)
    if atmosphere:
        conditions["atmosphere"] = atmosphere
    if params.scale_mmol:
        conditions["scale_mmol"] = round(float(params.scale_mmol), 3)
    return conditions


def _extract_temperature(hint: Optional[str], text: str) -> Optional[float]:
    if hint:
        value = _parse_numeric(hint)
        if value is not None:
            return value
    match = re.search(r"(\d+(?:\.\d+)?)\s*(?:°C|deg C)", text, re.IGNORECASE)
    if match:
        return float(match.group(1))
    return None


def _extract_time(hint: Optional[str], text: str) -> Optional[float]:
    if hint:
        value = _parse_numeric(hint)
        if value is not None:
            return value
    match = re.search(r"(\d+(?:\.\d+)?)\s*(h|hr|hrs|hours)", text, re.IGNORECASE)
    if match:
        return float(match.group(1))
    match = re.search(r"(\d+(?:\.\d+)?)\s*(min|minutes)", text, re.IGNORECASE)
    if match:
        return round(float(match.group(1)) / 60.0, 3)
    return None


def _extract_atmosphere(text: str) -> Optional[str]:
    if re.search(r"\bN2\b|\bnitrogen\b", text, re.IGNORECASE):
        return "N2"
    if re.search(r"\bAr\b|\bargon\b", text, re.IGNORECASE):
        return "Ar"
    if re.search(r"\bair\b", text, re.IGNORECASE):
        return "air"
    return None


def _parse_numeric(value: str) -> Optional[float]:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


def _derive_workup(procedure_text: str) -> List[Dict[str, Any]]:
    quench_steps: List[Dict[str, str]] = []
    workup_steps: List[Dict[str, str]] = []
    purification_steps: List[Dict[str, str]] = []
    notes: List[str] = []

    for segment in _split_segments(procedure_text):
        lowered = segment.lower()
        if "quench" in lowered:
            quench_steps.append({"step": "Quench", "details": segment})
        elif "wash" in lowered or "extract" in lowered:
            workup_steps.append({"step": "Workup", "details": segment})
        elif "dry" in lowered or "filter" in lowered:
            workup_steps.append({"step": "Dry/Filter", "details": segment})
        elif "column" in lowered or "purify" in lowered:
            purification_steps.append({"method": "Chromatography", "details": segment})
        else:
            notes.append(segment)

    if not (quench_steps or workup_steps or purification_steps):
        return []

    entry: Dict[str, Any] = {
        "quench": quench_steps,
        "workup": workup_steps,
        "purification": purification_steps,
    }
    if notes:
        entry["notes"] = notes[:5]
    return [entry]


def to_json_file(draft: Dict[str, Any], destination: Path) -> None:
    """Persist a draft to JSON."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8") as handle:
        json.dump(draft, handle, indent=2, ensure_ascii=False)


def _run_llm_parser(
    params: ProtocolDraftInput,
    procedure_text: str,
    normalized_smiles: str,
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    if not LLMClient:
        raise ProtocolDraftError("LLM client is unavailable.")

    provider = (params.llm_provider or os.getenv("LLM_PROVIDER") or "openai").lower()
    model = params.llm_model or os.getenv("LLM_MODEL")

    client = LLMClient(
        provider=provider,
        model=model,
        temperature=0.1,
        max_tokens=1500,
    )
    prompt = LLM_PROTOCOL_PROMPT.format(
        reaction_smiles=normalized_smiles,
        scale_mmol=params.scale_mmol,
        hints=_format_llm_hints(params),
        procedure_text=procedure_text,
    )
    response = client.chat(
        prompt=prompt,
        system=LLM_PROTOCOL_SYSTEM,
        temperature=0.1,
        max_tokens=1500,
    )
    payload = _parse_llm_json(response.content)
    if "steps" not in payload:
        raise ProtocolDraftError("LLM response missing 'steps'.")
    metadata = {
        "provider": provider,
        "model": client.model,
        "tokens": response.total_tokens,
    }
    return payload, metadata


def _convert_llm_steps(steps: Any) -> List[Dict[str, Any]]:
    sequence: List[Dict[str, Any]] = []
    if not isinstance(steps, list):
        return sequence
    for idx, entry in enumerate(steps, 1):
        if not isinstance(entry, dict):
            continue
        summary = str(entry.get("summary") or "").strip()
        material = str(entry.get("material") or "").strip()
        action = str(entry.get("action") or "").strip() or "Step"
        amount = str(entry.get("amount") or "").strip() or None
        role = str(entry.get("role") or "").strip() or None
        order = entry.get("order") or idx
        if not (summary or material):
            continue
        sequence.append(
            {
                "order": int(order),
                "action": action,
                "summary": summary or material,
                "material": material or None,
                "amount": amount,
                "role": role,
            }
        )
    return sequence


def _merge_conditions(base: Dict[str, Any], overrides: Dict[str, Any]) -> Dict[str, Any]:
    merged = dict(base)
    for key, value in overrides.items():
        if value in (None, "", []):
            continue
        if key not in merged or not merged[key]:
            merged[key] = value
        else:
            merged[key] = value
    return merged


def _format_llm_hints(params: ProtocolDraftInput) -> str:
    hints = []
    for label, value in (
        ("solvent", params.solvent_hint),
        ("base", params.base_hint),
        ("catalyst", params.catalyst_hint),
        ("additive", params.additive_hint),
        ("temperature_C", params.temperature_hint),
        ("time_h", params.time_hint),
        ("atmosphere", params.atmosphere_hint),
    ):
        if value:
            hints.append(f"{label}={value}")
    return ", ".join(hints) if hints else "None"


def _parse_llm_json(text: str) -> Dict[str, Any]:
    cleaned = _strip_code_fences(text.strip())
    try:
        return json.loads(cleaned)
    except json.JSONDecodeError:
        match = re.search(r"\{.*\}", cleaned, re.DOTALL)
        if match:
            try:
                return json.loads(match.group(0))
            except json.JSONDecodeError as exc:  # pragma: no cover - defensive
                raise ProtocolDraftError(f"Unable to parse LLM JSON: {exc}") from exc
    raise ProtocolDraftError("LLM response did not contain parseable JSON.")


def _strip_code_fences(text: str) -> str:
    if text.startswith("```"):
        text = re.sub(r"^```(?:json)?", "", text, flags=re.IGNORECASE).strip()
        if text.endswith("```"):
            text = text[: -3].strip()
    return text

"""Fragment-based substituent composition utilities.

This module composes additional substituent groups from linker/terminal
fragments while preserving the runtime A-B motif model.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from chemtools.util.rdkit_helpers import rdkit_available
from chemtools.util.smarts_cache import compile_smarts

SUBSTITUENT_FRAGMENTS_FILE = "substituent_fragments.v1.json"
_MAP_RE = re.compile(r":\d+(?=\])")
_LINKER_ALLOWED_KEYS = {
    "label",
    "priority",
    "smarts_template",
    "description",
    "allowed_terminal_groups",
    "blocked_terminal_groups",
}
_LINKER_REQUIRED_KEYS = ("label", "smarts_template")
_GROUP_ALLOWED_KEYS = {"id", "label", "terminal_group", "aliases", "description", "priority"}
_GROUP_REQUIRED_KEYS = ("id", "label", "terminal_group")


def _safe_load_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    return payload if isinstance(payload, dict) else {}


def _default_fragments_path(groups_path: Path) -> Path:
    return groups_path.resolve().parent / SUBSTITUENT_FRAGMENTS_FILE


def _strip_atom_maps(smarts: str) -> str:
    return _MAP_RE.sub("", smarts)


def validate_substituent_fragments_payload(fragments_payload: Dict[str, Any]) -> List[str]:
    """Validate substituent fragment payload schema."""
    errors: List[str] = []
    if not isinstance(fragments_payload, dict):
        return ["Payload must be a JSON object."]

    linkers = fragments_payload.get("linkers")
    groups = fragments_payload.get("groups")
    if not isinstance(linkers, list) or not linkers:
        errors.append("'linkers' must be a non-empty list.")
        linkers = []
    if not isinstance(groups, list) or not groups:
        errors.append("'groups' must be a non-empty list.")
        groups = []

    seen_linkers: set[str] = set()
    for idx, entry in enumerate(linkers):
        if not isinstance(entry, dict):
            errors.append(f"linkers[{idx}] must be an object.")
            continue
        unexpected = sorted(set(entry.keys()) - _LINKER_ALLOWED_KEYS)
        if unexpected:
            errors.append(f"linkers[{idx}] has unexpected keys: {', '.join(unexpected)}")
        for key in _LINKER_REQUIRED_KEYS:
            if not str(entry.get(key) or "").strip():
                errors.append(f"linkers[{idx}] missing required key '{key}'.")
        label = str(entry.get("label") or "").strip()
        if label and label != label.upper():
            errors.append(f"linkers[{idx}] label '{label}' must be uppercase (e.g., CO, SO2).")
        if label:
            upper = label.upper()
            if upper in seen_linkers:
                errors.append(f"Duplicate linker label '{upper}'.")
            seen_linkers.add(upper)
        smarts_template = str(entry.get("smarts_template") or "").strip()
        if smarts_template and "{TAIL}" not in smarts_template:
            errors.append(f"linkers[{idx}] smarts_template must include '{{TAIL}}'.")
        allowed_terminals = entry.get("allowed_terminal_groups")
        blocked_terminals = entry.get("blocked_terminal_groups")
        if allowed_terminals is not None and not isinstance(allowed_terminals, list):
            errors.append(f"linkers[{idx}] allowed_terminal_groups must be a list when provided.")
            allowed_terminals = []
        if blocked_terminals is not None and not isinstance(blocked_terminals, list):
            errors.append(f"linkers[{idx}] blocked_terminal_groups must be a list when provided.")
            blocked_terminals = []
        allowed_set = {str(v).strip() for v in (allowed_terminals or []) if str(v).strip()}
        blocked_set = {str(v).strip() for v in (blocked_terminals or []) if str(v).strip()}
        if allowed_terminals is not None and not allowed_set:
            errors.append(f"linkers[{idx}] allowed_terminal_groups must contain non-empty strings.")
        if blocked_terminals is not None and not blocked_set:
            errors.append(f"linkers[{idx}] blocked_terminal_groups must contain non-empty strings.")
        overlap = sorted(allowed_set & blocked_set)
        if overlap:
            errors.append(
                f"linkers[{idx}] terminal policy overlap in allowed/blocked: {', '.join(overlap)}"
            )

    seen_ids: set[str] = set()
    for idx, entry in enumerate(groups):
        if not isinstance(entry, dict):
            errors.append(f"groups[{idx}] must be an object.")
            continue
        unexpected = sorted(set(entry.keys()) - _GROUP_ALLOWED_KEYS)
        if unexpected:
            errors.append(f"groups[{idx}] has unexpected keys: {', '.join(unexpected)}")
        if "terminal" in entry and "terminal_group" not in entry:
            errors.append(f"groups[{idx}] uses deprecated key 'terminal'; use 'terminal_group'.")
        for key in _GROUP_REQUIRED_KEYS:
            if not str(entry.get(key) or "").strip():
                errors.append(f"groups[{idx}] missing required key '{key}'.")
        group_id = str(entry.get("id") or "").strip()
        if group_id and not group_id.startswith("-"):
            errors.append(f"groups[{idx}] id '{group_id}' must start with '-'.")
        if group_id:
            if group_id in seen_ids:
                errors.append(f"Duplicate group id '{group_id}'.")
            seen_ids.add(group_id)
        label = str(entry.get("label") or "").strip()
        if label and label != label.upper():
            errors.append(f"groups[{idx}] label '{label}' must be uppercase linker label.")
        if label and label.upper() not in seen_linkers:
            errors.append(f"groups[{idx}] references unknown linker label '{label}'.")
        aliases = entry.get("aliases")
        if aliases is not None and not isinstance(aliases, list):
            errors.append(f"groups[{idx}] aliases must be a list when provided.")
        elif isinstance(aliases, list):
            bad_alias = next((a for a in aliases if not str(a).strip()), None)
            if bad_alias is not None:
                errors.append(f"groups[{idx}] aliases must contain non-empty strings only.")

    return errors


def compose_groups_from_fragments(
    *,
    base_groups: List[Dict[str, Any]],
    fragments_payload: Dict[str, Any],
) -> Tuple[List[Dict[str, Any]], List[str]]:
    """Compose substituent groups from fragment definitions.

    Returns:
        Tuple of (generated_groups, errors).
    """
    errors: List[str] = []
    generated: List[Dict[str, Any]] = []

    existing_ids = {
        str(entry.get("id") or "").strip()
        for entry in base_groups
        if isinstance(entry, dict) and entry.get("id")
    }

    errors.extend(validate_substituent_fragments_payload(fragments_payload))

    raw_linkers = fragments_payload.get("linkers", []) or []

    linkers = {
        str(entry.get("label") or "").strip().upper(): entry
        for entry in raw_linkers
        if isinstance(entry, dict) and entry.get("label")
    }
    terminal_groups = {
        str(entry.get("id") or "").strip(): entry
        for entry in base_groups
        if isinstance(entry, dict) and entry.get("id")
    }

    compositions = fragments_payload.get("groups", []) or []
    for item in compositions:
        if not isinstance(item, dict):
            continue
        linker_label = str(item.get("label") or "").strip().upper()
        terminal_group_id = str(item.get("terminal_group") or "").strip()
        group_id = str(item.get("id") or "").strip()
        if not linker_label or not terminal_group_id or not group_id:
            errors.append("Invalid group entry: missing label/terminal_group/id.")
            continue
        if group_id in existing_ids:
            continue

        linker = linkers.get(linker_label)
        terminal = terminal_groups.get(terminal_group_id)
        if not linker or not terminal:
            errors.append(f"Composition {group_id}: unknown fragment reference.")
            continue
        allowed_terminals = {
            str(v).strip()
            for v in (linker.get("allowed_terminal_groups") or [])
            if str(v).strip()
        }
        blocked_terminals = {
            str(v).strip()
            for v in (linker.get("blocked_terminal_groups") or [])
            if str(v).strip()
        }
        if allowed_terminals and terminal_group_id not in allowed_terminals:
            errors.append(
                f"Composition {group_id}: terminal_group '{terminal_group_id}' is not allowed for linker '{linker_label}'."
            )
            continue
        if terminal_group_id in blocked_terminals:
            errors.append(
                f"Composition {group_id}: terminal_group '{terminal_group_id}' is blocked for linker '{linker_label}'."
            )
            continue

        template = str(linker.get("smarts_template") or "")
        tail_smarts = str(terminal.get("tail_smarts") or "")
        if not tail_smarts:
            tail_smarts = _strip_atom_maps(str(terminal.get("smarts") or ""))
        if not template or "{TAIL}" not in template or not tail_smarts:
            errors.append(f"Composition {group_id}: invalid SMARTS template/terminal.")
            continue

        smarts = template.replace("{TAIL}", tail_smarts)
        if rdkit_available() and compile_smarts(smarts, validate=False) is None:
            errors.append(f"Composition {group_id}: SMARTS failed to compile.")
            continue

        linker_priority = int(linker.get("priority", 1))
        terminal_priority = int(terminal.get("priority", 1))
        priority = int(item.get("priority", max(linker_priority, terminal_priority)))

        entry: Dict[str, Any] = {
            "id": group_id,
            "kind": "substituent",
            "priority": priority,
            "smarts": smarts,
            "description": str(item.get("description") or "").strip()
            or f"Composed substituent {group_id} from {linker_label}+{terminal_group_id}.",
            "generated": {
                "source": "substituent_fragments.v1",
                "linker": linker_label,
                "terminal_group": terminal_group_id,
            },
        }

        aliases = item.get("aliases")
        if isinstance(aliases, list):
            cleaned_aliases = [str(a).strip() for a in aliases if str(a).strip()]
            if cleaned_aliases:
                entry["aliases"] = cleaned_aliases

        generated.append(entry)
        existing_ids.add(group_id)

    return generated, errors


def load_organic_groups_with_compositions(
    groups_path: Path,
    *,
    fragments_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """Load organic groups and append composed groups when configured."""
    payload = _safe_load_json(groups_path)
    groups = payload.get("groups", []) or []
    if not isinstance(groups, list):
        payload["groups"] = []
        return payload

    base_groups = [g for g in groups if isinstance(g, dict)]
    merged_groups = list(base_groups)

    fragments_file = fragments_path or _default_fragments_path(groups_path)
    if not fragments_file.exists():
        payload["groups"] = merged_groups
        return payload

    fragments_payload = _safe_load_json(fragments_file)
    generated, errors = compose_groups_from_fragments(
        base_groups=merged_groups,
        fragments_payload=fragments_payload,
    )
    merged_groups.extend(generated)
    payload["groups"] = merged_groups
    payload["composed_groups"] = {
        "source_file": str(fragments_file),
        "generated_count": len(generated),
        "errors": errors,
    }
    return payload


__all__ = [
    "SUBSTITUENT_FRAGMENTS_FILE",
    "compose_groups_from_fragments",
    "load_organic_groups_with_compositions",
    "validate_substituent_fragments_payload",
]

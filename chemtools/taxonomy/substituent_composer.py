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


def _safe_load_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    return payload if isinstance(payload, dict) else {}


def _default_fragments_path(groups_path: Path) -> Path:
    return groups_path.resolve().parent / SUBSTITUENT_FRAGMENTS_FILE


def _strip_atom_maps(smarts: str) -> str:
    return _MAP_RE.sub("", smarts)


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
        terminal_group_id = str(item.get("terminal_group") or item.get("terminal") or "").strip()
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
]

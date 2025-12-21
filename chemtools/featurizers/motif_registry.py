"""
Build a motif registry for compound SMARTS derived from organic group templates.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import re
from typing import Any, Dict, Iterable, List, Mapping

from chemtools.util.smarts_cache import compile_smarts

_MAP_RE = re.compile(r":\d+")


@dataclass(frozen=True)
class CompoundPattern:
    compound_id: str
    smarts: str
    query: Any
    group_a: str
    group_b: str
    b_tags: List[str]


def build_compound_registry(registry_paths: Mapping[str, str | Path]) -> Dict[str, Any]:
    """
    Load group/template/compound files and compile compound SMARTS queries.
    """
    groups_path = Path(registry_paths["groups"])
    templates_path = Path(registry_paths["templates"])
    compounds_path = Path(registry_paths["compounds"])

    groups = _load_groups(groups_path)
    templates = _load_templates(templates_path)
    compounds = _load_compounds(compounds_path)

    compiled: List[CompoundPattern] = []
    compound_map: Dict[str, CompoundPattern] = {}

    for entry in compounds:
        compound_id = entry.get("id")
        if not compound_id:
            continue
        direct_smarts = entry.get("smarts")
        smarts_any = entry.get("smarts_any")
        smarts_list: List[str] = []
        if isinstance(direct_smarts, str):
            smarts_list = [direct_smarts]
        elif isinstance(direct_smarts, list):
            smarts_list = [s for s in direct_smarts if isinstance(s, str)]
        elif isinstance(smarts_any, list):
            smarts_list = [s for s in smarts_any if isinstance(s, str)]

        if smarts_list:
            for smarts in smarts_list:
                smarts = smarts.strip()
                if not smarts:
                    continue
                if not _has_atom_map(smarts):
                    smarts = _inject_map_on_first_atom(smarts, map_num=1)
                query = compile_smarts(smarts, validate=False)
                if query is None:
                    continue
                pattern = CompoundPattern(
                    compound_id=compound_id,
                    smarts=smarts,
                    query=query,
                    group_a=str(entry.get("A") or ""),
                    group_b=str(entry.get("B") or ""),
                    b_tags=list(entry.get("tags") or []),
                )
                compiled.append(pattern)
                compound_map[compound_id] = pattern
            continue

        template_id = entry.get("template", "")
        template_format = templates.get(template_id)
        if not template_format:
            continue

        group_a = entry.get("A")
        group_b = entry.get("B")
        if not group_a or not group_b:
            continue
        group_a_record = groups.get(group_a)
        group_b_record = groups.get(group_b)
        if not group_a_record or not group_b_record:
            continue

        a_smarts = group_a_record.get("smarts", "")
        b_smarts = group_b_record.get("smarts", "")
        if not a_smarts or not b_smarts:
            continue

        compound_smarts = _format_compound_smarts(
            template_id=template_id,
            template_format=template_format,
            a_smarts=a_smarts,
            b_smarts=b_smarts,
        )
        query = compile_smarts(compound_smarts, validate=False)
        if query is None:
            continue

        pattern = CompoundPattern(
            compound_id=compound_id,
            smarts=compound_smarts,
            query=query,
            group_a=group_a,
            group_b=group_b,
            b_tags=list(group_b_record.get("tags") or []),
        )
        compiled.append(pattern)
        compound_map[compound_id] = pattern

    return {
        "groups": groups,
        "templates": templates,
        "compiled_compounds": compiled,
        "compound_map": compound_map,
    }


def _load_groups(path: Path) -> Dict[str, Dict[str, Any]]:
    payload = _load_json(path)
    groups = payload.get("groups", [])
    return {group["id"]: group for group in groups if isinstance(group, dict) and group.get("id")}


def _load_templates(path: Path) -> Dict[str, str]:
    payload = _load_json(path)
    templates = payload.get("templates", {})
    formatted: Dict[str, str] = {}
    for template_id, entry in templates.items():
        if isinstance(entry, dict) and "format" in entry:
            formatted[template_id] = str(entry["format"])
    return formatted


def _load_compounds(path: Path) -> List[Dict[str, Any]]:
    payload = _load_json(path)
    return list(payload.get("compounds") or [])


def _load_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Motif registry file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _format_compound_smarts(
    *,
    template_id: str,
    template_format: str,
    a_smarts: str,
    b_smarts: str,
) -> str:
    if template_id == "via_oxygen" and not _has_atom_map(b_smarts):
        return f"{a_smarts}[O:2]{b_smarts}"

    if template_id == "single_bond" and not _has_atom_map(b_smarts):
        b_smarts = _inject_map_on_first_atom(b_smarts, map_num=2)

    return template_format.format(A=a_smarts, B=b_smarts)


def _has_atom_map(smarts: str) -> bool:
    return bool(_MAP_RE.search(smarts))


def _inject_map_on_first_atom(smarts: str, *, map_num: int) -> str:
    if not smarts:
        return smarts
    if smarts.startswith("["):
        end = smarts.find("]")
        if end == -1:
            return smarts
        return f"{smarts[:end]}:{map_num}{smarts[end:]}"

    match = re.match(r"([A-Za-z][a-z]?)", smarts)
    if not match:
        return smarts
    atom = match.group(1)
    return f"[{atom}:{map_num}]{smarts[len(atom):]}"

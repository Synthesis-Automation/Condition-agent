"""
Build a motif registry for compound SMARTS derived from organic group templates.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import re
from typing import Any, Dict, Iterable, List, Mapping, Optional, Set

from chemtools.util.smarts_cache import compile_smarts
from chemtools.util.rdkit_helpers import parse_smiles, rdkit_available

_MAP_RE = re.compile(r":\d+(?=\])")


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
    _validate_group_maps(groups)
    _validate_compound_templates(compounds, templates)

    compiled: List[CompoundPattern] = []
    compound_map: Dict[str, CompoundPattern] = {}

    for entry in compounds:
        compound_id = entry.get("id")
        if not compound_id:
            continue
        smarts_list = _extract_compound_smarts(entry)
        if smarts_list:
            for smarts in smarts_list:
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


def build_compound_detect_registry(registry_paths: Mapping[str, str | Path]) -> Dict[str, Any]:
    """
    Load group/template/compound files and compile detect SMARTS queries.

    Detect SMARTS strips atom-map labels for pure substructure matching.
    """
    groups_path = Path(registry_paths["groups"])
    templates_path = Path(registry_paths["templates"])
    compounds_path = Path(registry_paths["compounds"])

    groups = _load_groups(groups_path)
    templates = _load_templates(templates_path)
    compounds = _load_compounds(compounds_path)
    _validate_group_maps(groups)
    _validate_compound_templates(compounds, templates)

    compiled: Dict[str, List[Any]] = {}
    direct_compounds: Set[str] = set()

    for entry in compounds:
        compound_id = entry.get("id")
        if not compound_id:
            continue

        smarts_list = _extract_compound_smarts(entry)
        if smarts_list:
            direct_compounds.add(compound_id)
        else:
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
            smarts_list = [
                _format_compound_smarts(
                    template_format=template_format,
                    a_smarts=a_smarts,
                    b_smarts=b_smarts,
                )
            ]

        detect_smarts_list = [strip_atom_maps(s) for s in smarts_list]
        queries: List[Any] = []
        for smarts in detect_smarts_list:
            query = compile_smarts(smarts, validate=False)
            if query is not None:
                queries.append(query)
        if queries:
            compiled[compound_id] = queries

    return {
        "groups": groups,
        "templates": templates,
        "compiled_compounds": compiled,
        "direct_compounds": direct_compounds,
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
    compounds = list(payload.get("compounds") or [])
    for entry in compounds:
        if not isinstance(entry, dict):
            continue
        if "anchors" not in entry:
            entry["anchors"] = {"scaffold": "A", "substituent": "B"}
    return compounds


def _load_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Motif registry file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _format_compound_smarts(
    *,
    template_format: str,
    a_smarts: str,
    b_smarts: str,
) -> str:
    return template_format.format(A=a_smarts, B=b_smarts)


def _has_atom_map(smarts: str) -> bool:
    return bool(_MAP_RE.search(smarts))


def strip_atom_maps(smarts: str) -> str:
    return _MAP_RE.sub("", smarts)


def classify_compound_smiles(
    smiles: str,
    *,
    registry: Optional[Mapping[str, Any]] = None,
    registry_paths: Optional[Mapping[str, str | Path]] = None,
    include_best: bool = True,
    prefer_direct: bool = False,
) -> Dict[str, Any]:
    """
    Classify a SMILES string into compound motif IDs using detect SMARTS.
    """
    result: Dict[str, Any] = {"smiles": smiles, "ok": False, "hits": [], "best": None}
    if not rdkit_available():
        result["error"] = "rdkit_unavailable"
        return result

    mol = parse_smiles(smiles)
    if mol is None:
        result["error"] = "invalid_smiles"
        return result

    if registry is None:
        registry_paths = registry_paths or _default_registry_paths()
        registry = build_compound_detect_registry(registry_paths)

    compiled = registry.get("compiled_compounds") or {}
    hits: List[str] = []
    for compound_id, queries in compiled.items():
        if any(mol.HasSubstructMatch(query) for query in queries):
            hits.append(compound_id)

    result["ok"] = True
    result["hits"] = hits
    if include_best:
        best = choose_best_compound_hit(
            hits,
            direct_ids=set(registry.get("direct_compounds") or []),
            prefer_direct=prefer_direct,
        )
        result["best"] = best
    return result


def classify_compound_batch(
    smiles_list: Iterable[str],
    *,
    registry: Optional[Mapping[str, Any]] = None,
    registry_paths: Optional[Mapping[str, str | Path]] = None,
    include_best: bool = True,
    prefer_direct: bool = False,
) -> List[Dict[str, Any]]:
    """
    Classify a batch of SMILES strings into compound motif IDs.
    """
    if registry is None:
        registry_paths = registry_paths or _default_registry_paths()
        registry = build_compound_detect_registry(registry_paths)

    return [
        classify_compound_smiles(
            smiles,
            registry=registry,
            include_best=include_best,
            prefer_direct=prefer_direct,
        )
        for smiles in smiles_list
    ]


def choose_best_compound_hit(
    hits: Iterable[str],
    *,
    direct_ids: Optional[Set[str]] = None,
    prefer_direct: bool = False,
) -> Optional[str]:
    """
    Choose a single motif label using simple precedence heuristics.
    """
    hits_list = [h for h in hits if h]
    if not hits_list:
        return None
    direct_ids = direct_ids or set()

    def prefix_rank(hit: str) -> int:
        if hit.startswith("Ar-"):
            return 0
        if hit.startswith("Vinyl-"):
            return 1
        if hit.startswith(("R-", "Bn-", "Allyl-")):
            return 2
        return 3

    def rank(hit: str) -> tuple[int, int, int, str]:
        direct_rank = 0 if prefer_direct and hit in direct_ids else 1
        return (direct_rank, prefix_rank(hit), -len(hit), hit)

    return sorted(hits_list, key=rank)[0]


def _has_map(smarts: str, *, map_num: int) -> bool:
    if not smarts:
        return False
    return bool(re.search(rf":{map_num}(?=\])", smarts))


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


def _extract_compound_smarts(entry: Mapping[str, Any]) -> List[str]:
    direct_smarts = entry.get("smarts")
    smarts_any = entry.get("smarts_any")
    smarts_list: List[str] = []
    if isinstance(direct_smarts, str):
        smarts_list = [direct_smarts]
    elif isinstance(direct_smarts, list):
        smarts_list = [s for s in direct_smarts if isinstance(s, str)]
    elif isinstance(smarts_any, list):
        smarts_list = [s for s in smarts_any if isinstance(s, str)]
    return [s.strip() for s in smarts_list if isinstance(s, str) and s.strip()]


def _validate_group_maps(groups: Mapping[str, Mapping[str, Any]]) -> None:
    errors: List[str] = []
    for group_id, group in groups.items():
        if not isinstance(group, Mapping):
            continue
        kind = str(group.get("kind") or "")
        smarts = str(group.get("smarts") or "")
        if not kind or not smarts:
            continue
        expected_map = 1 if kind == "scaffold" else 2
        if not _has_map(smarts, map_num=expected_map):
            errors.append(f"{group_id} (kind={kind}, expected :{expected_map})")
    if errors:
        joined = ", ".join(sorted(errors))
        raise ValueError(f"Group SMARTS missing required attach maps: {joined}")


def _validate_compound_templates(
    compounds: Iterable[Dict[str, Any]],
    templates: Mapping[str, str],
) -> None:
    missing: set[str] = set()
    for entry in compounds:
        if not isinstance(entry, dict):
            continue
        if _extract_compound_smarts(entry):
            continue
        template_id = entry.get("template")
        if template_id and template_id not in templates:
            missing.add(str(template_id))
    if missing:
        joined = ", ".join(sorted(missing))
        raise ValueError(f"Compound templates missing from registry: {joined}")


def _default_registry_paths() -> Dict[str, Path]:
    base = Path(__file__).resolve().parent.parent / "taxonomy" / "v2_data"
    return {
        "groups": base / "organic_groups.v1.3.json",
        "compounds": base / "organic_compounds.v1.3.json",
        "templates": base / "smarts_templates.v1.json",
    }

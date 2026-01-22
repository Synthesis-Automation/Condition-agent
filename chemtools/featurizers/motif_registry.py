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

_DEFAULT_TEMPLATES = {
    "single_bond": "{A}-{B}",
    "via_oxygen": "{A}O{B}",
}


@dataclass(frozen=True)
class CompoundPattern:
    compound_id: str
    smarts: str
    query: Any
    group_a: str
    group_b: str
    b_tags: List[str]
    priority: int = 1
    complexity: int = 0
    reactivity_weight: float = 0.0


def build_compound_registry(registry_paths: Mapping[str, str | Path]) -> Dict[str, Any]:
    """
    Load group/template/compound files and compile compound SMARTS queries.
    """
    groups_path = Path(registry_paths["groups"])
    templates_path = Path(registry_paths["templates"]) if "templates" in registry_paths else None
    compounds_path = Path(registry_paths["compounds"])
    logic_path = Path(registry_paths["logic"]) if "logic" in registry_paths else None

    groups = _load_groups(groups_path)
    templates = _load_templates(templates_path)
    compounds = _load_compounds(compounds_path)
    
    if "scaffold_motifs" in registry_paths:
        scaffold_path = Path(registry_paths["scaffold_motifs"])
        compounds.extend(_load_compounds(scaffold_path))
    
    # Load priorities from groups and logic sets
    priorities = {g_id: g.get("priority", 1) for g_id, g in groups.items()}
    if logic_path and logic_path.exists():
        logic_data = _load_json(logic_path)
        group_sets = logic_data.get("group_sets", {})
        for set_id, set_data in group_sets.items():
            if "priority" in set_data:
                priorities[set_id] = set_data["priority"]

    _validate_group_maps(groups)
    _validate_compound_templates(compounds, templates)

    compiled: List[CompoundPattern] = []
    compound_map: Dict[str, CompoundPattern] = {}
    combination_map: Dict[Tuple[str, str], str] = {}

    for entry in compounds:
        compound_id = entry.get("id")
        group_a = str(entry.get("A") or "")
        group_b = str(entry.get("B") or "")
        
        # Auto-generate ID as A+B when missing (v1.3 simplified format)
        if not compound_id:
            if group_a and group_b:
                compound_id = f"{group_a}{group_b}"
            else:
                continue
        if group_a and group_b:
            combination_map[(group_a, group_b)] = compound_id

        # Calculate priority: Priority(A) + Priority(B), unless overridden
        priority = entry.get("priority")
        if priority is None:
            priority = priorities.get(group_a, 1) + priorities.get(group_b, 1)
        reactivity_weight = float(entry.get("reactivity_weight") or 0.0)

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
                    group_a=group_a,
                    group_b=group_b,
                    b_tags=list(entry.get("tags") or []),
                    priority=priority,
                    complexity=calculate_smarts_complexity(query),
                    reactivity_weight=reactivity_weight,
                )
                compiled.append(pattern)
                compound_map[compound_id] = pattern
            continue

        template_id = entry.get("template", "")
        template_format = templates.get(template_id)
        if not template_format:
            continue

        if not group_a or not group_b:
            continue
        group_a_record = groups.get(group_a)
        group_b_record = groups.get(group_b)
        
        # If not in groups, check if it's a logic set (we already have its priority)
        # But we need SMARTS to build the compound SMARTS.
        # For now, we assume A and B are either in groups or we skip if SMARTS can't be built.
        # (The existing logic already does this)

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
            priority=priority,
            complexity=calculate_smarts_complexity(query),
            reactivity_weight=reactivity_weight,
        )
        compiled.append(pattern)
        compound_map[compound_id] = pattern

    # Re-sort compiled patterns by priority and complexity
    # Higher priority first. If priority tied, more complex (narrower) SMARTS first.
    def sort_key(p: CompoundPattern) -> tuple[int, int, str]:
        return (-p.priority, -p.complexity, p.compound_id)

    compiled.sort(key=sort_key)

    compiled_groups: Dict[str, Dict[str, Any]] = {}
    for group_id, group in groups.items():
        smarts = group.get("smarts")
        if smarts:
            query = compile_smarts(smarts, validate=False)
            if query:
                compiled_groups[group_id] = {
                    "id": group_id,
                    "kind": group.get("kind"),
                    "query": query,
                    "priority": priorities.get(group_id, 1),
                    "complexity": calculate_smarts_complexity(query),
                    "tags": list(group.get("tags") or []),
                }

    return {
        "groups": groups,
        "templates": templates,
        "compiled_compounds": compiled,
        "compound_map": compound_map,
        "combination_map": combination_map,
        "priorities": priorities,
        "compiled_groups": compiled_groups,
    }


def _load_groups(path: Path) -> Dict[str, Dict[str, Any]]:
    payload = _load_json(path)
    groups = payload.get("groups", [])
    return {group["id"]: group for group in groups if isinstance(group, dict) and group.get("id")}


def _load_templates(path: Optional[Path]) -> Dict[str, str]:
    if path is None or not path.exists():
        return dict(_DEFAULT_TEMPLATES)

    payload = _load_json(path)
    templates = payload.get("templates", {})
    formatted: Dict[str, str] = dict(_DEFAULT_TEMPLATES)
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
        # Auto-generate compound ID from A+B
        if "id" not in entry:
            group_a = str(entry.get("A") or "")
            group_b = str(entry.get("B") or "")
            if group_a and group_b:
                # Add hyphen separator if B doesn't already start with "-"
                # Substituents have "-" prefix (e.g., "-Cl"), scaffolds don't (e.g., "Alkenyl")
                separator = "" if group_b.startswith("-") else "-"
                entry["id"] = f"{group_a}{separator}{group_b}"
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


def classify_compound_smiles(
    smiles: str,
    *,
    registry: Optional[Mapping[str, Any]] = None,
    registry_paths: Optional[Mapping[str, str | Path]] = None,
    include_best: bool = True,
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
        registry = build_compound_registry(registry_paths)

    compiled = registry.get("compiled_compounds") or []
    hits: List[str] = []
    
    # Use the compiled patterns which now include priority
    for pattern in compiled:
        if mol.HasSubstructMatch(pattern.query):
            hits.append(pattern.compound_id)

    result["ok"] = True
    result["hits"] = list(set(hits))
    if include_best:
        best = choose_best_compound_hit(
            result["hits"],
            compound_map=registry.get("compound_map"),
        )
        result["best"] = best
    return result


def classify_compound_batch(
    smiles_list: Iterable[str],
    *,
    registry: Optional[Mapping[str, Any]] = None,
    registry_paths: Optional[Mapping[str, str | Path]] = None,
    include_best: bool = True,
) -> List[Dict[str, Any]]:
    """
    Classify a batch of SMILES strings into compound motif IDs.
    """
    if registry is None:
        registry_paths = registry_paths or _default_registry_paths()
        registry = build_compound_registry(registry_paths)

    return [
        classify_compound_smiles(
            smiles,
            registry=registry,
            include_best=include_best,
        )
        for smiles in smiles_list
    ]


def choose_best_compound_hit(
    hits: Iterable[str],
    *,
    compound_map: Optional[Mapping[str, CompoundPattern]] = None,
) -> Optional[str]:
    """
    Choose the most specific motif label based on priority score.
    """
    hits_list = [h for h in hits if h]
    if not hits_list:
        return None
    
    if not compound_map:
        # Fallback to alphabetical if no map provided
        return sorted(hits_list)[0]

    def score_rank(hit: str) -> tuple[int, int, int, str]:
        pattern = compound_map.get(hit)
        priority = pattern.priority if pattern else 0
        complexity = calculate_smarts_complexity(pattern.query) if pattern else 0
        # Higher priority first, then higher complexity (narrower), 
        # then longer ID as tie-breaker, then alphabetical
        return (-priority, -complexity, -len(hit), hit)

    return sorted(hits_list, key=score_rank)[0]


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


def calculate_smarts_complexity(query_mol: Any) -> int:
    """
    Calculate a structural complexity score for a SMARTS query molecule.
    Higher scores indicate more specific patterns.
    """
    if not query_mol:
        return 0
    score = 0
    # Score atoms
    for atom in query_mol.GetAtoms():
        score += 10 # Base atom
        symbol = atom.GetSymbol()
        # Non-wildcard/Non-H atoms are more specific
        if symbol != "*" and symbol != "H":
            score += 10 # Heavy atom bonus
        if symbol != "*":
            score += 5 # Explicit element
        
        # Parse SMARTS constraints for specificity
        smarts = atom.GetSmarts()
        score += len(smarts)
        if "X" in smarts: score += 5  # Hybridization constraint
        if "R" in smarts or "r" in smarts: score += 5 # Ring constraint
        if "H" in smarts or "h" in smarts: score += 3 # Explicit H-count
        if "!" in smarts: score += 2 # Negations
        score += smarts.count(";") * 2 # AND connections
        score -= smarts.count(",") * 2 # OR connections (Genericness)

    # Score bonds
    for bond in query_mol.GetBonds():
        score += 5 # Base bond
        bond_smarts = bond.GetSmarts()
        score += len(bond_smarts)
        if bond.GetBondTypeAsDouble() > 1:
            score += 3 # Multiple bonds are more specific
        if bond_smarts and "@" in bond_smarts:
            score += 5 # Stereochemical constraints
            
    return score


def _default_registry_paths() -> Dict[str, Path]:
    base = Path(__file__).resolve().parent.parent / "taxonomy" / "data"
    return {
        "groups": base / "organic_groups.v1.3.json",
        "compounds": base / "organic_compounds.v1.3.json",
        "scaffold_motifs": base / "scaffold_motifs.v1.3.json",
        "logic": base / "group_logic.json",
    }

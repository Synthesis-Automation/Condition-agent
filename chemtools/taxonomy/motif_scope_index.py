"""
Build and persist taxonomy-driven motif scope index for recommendation matching.

The index is generated offline and cached as JSON. Runtime consumers can load it
without recomputing scaffold hierarchy relations.
"""

from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Set, Tuple

from chemtools.core.rdkit import parse_smiles, rdkit_available
from chemtools.core.smarts import compile_smarts

from . import loader


DEFAULT_PROBE_SMILES: Tuple[str, ...] = (
    "CC",
    "CCC",
    "CC(C)C",
    "CC(C)(C)C",
    "Cc1ccccc1",
    "CC=C",
    "C=CC",
    "CC=CC",
    "CC=CCO",
    "C=CC=C",
    "c1ccccc1",
    "c1ccncc1",
    "c1ncccc1",
    "C#C",
    "CC#C",
    "C#CC",
    "CC#CC",
    "CC#CCO",
    "CCCO",
    "COC",
    "CC(=O)O",
    "NCC",
)

_PREFERRED_PARENT_SCAFFOLDS: Set[str] = {
    "Ar",
    "Alkyl",
    "Alkenyl",
    "Alkynyl",
}
_BLOCKED_CHILDREN_BY_PARENT: Dict[str, Set[str]] = {
    "Alkyl": {"Acyl", "R_acidic"},
}


def _compound_entry_id(entry: Dict[str, Any]) -> str:
    motif_id = str(entry.get("id") or "").strip()
    if motif_id:
        return motif_id
    scaffold = str(entry.get("A") or "").strip()
    substituent = str(entry.get("B") or "").strip()
    if not scaffold or not substituent:
        return ""
    if substituent.startswith("-"):
        substituent = substituent[1:]
    if not substituent:
        return ""
    return f"{scaffold}-{substituent}"


def _probe_molecules(probe_smiles: Iterable[str]) -> List[Any]:
    mols: List[Any] = []
    for smiles in probe_smiles:
        mol = parse_smiles(str(smiles))
        if mol is not None:
            mols.append(mol)
    return mols


def _vector_subset(child_vec: Tuple[bool, ...], parent_vec: Tuple[bool, ...]) -> bool:
    return not any(c and not p for c, p in zip(child_vec, parent_vec))


def _infer_scaffold_parent_map(
    scaffold_payload: Dict[str, Any],
    probe_smiles: Iterable[str],
) -> Dict[str, Set[str]]:
    if not rdkit_available():
        return {}

    scaffold_defs: Dict[str, Dict[str, Any]] = {}
    for entry in scaffold_payload.get("groups", []) or []:
        if not isinstance(entry, dict):
            continue
        if str(entry.get("kind") or "").strip() != "scaffold":
            continue
        group_id = str(entry.get("id") or "").strip()
        smarts = str(entry.get("smarts") or "").strip()
        if not group_id or not smarts:
            continue
        patt = compile_smarts(smarts, validate=False)
        if patt is None:
            continue
        priority = int(entry.get("priority") or 0)
        scaffold_defs[group_id] = {"pattern": patt, "priority": priority}

    if not scaffold_defs:
        return {}

    probes = _probe_molecules(probe_smiles)
    if not probes:
        return {}

    match_vectors: Dict[str, Tuple[bool, ...]] = {}
    for scaffold_id, data in scaffold_defs.items():
        patt = data["pattern"]
        vector = tuple(bool(mol.HasSubstructMatch(patt)) for mol in probes)
        if any(vector):
            match_vectors[scaffold_id] = vector

    parent_map: Dict[str, Set[str]] = defaultdict(set)
    scaffold_ids = sorted(match_vectors.keys())
    for child_id in scaffold_ids:
        child_vec = match_vectors[child_id]
        child_priority = scaffold_defs[child_id]["priority"]
        for parent_id in scaffold_ids:
            if parent_id == child_id:
                continue
            if parent_id not in _PREFERRED_PARENT_SCAFFOLDS:
                continue
            parent_priority = scaffold_defs[parent_id]["priority"]
            if parent_priority > child_priority:
                continue
            parent_vec = match_vectors[parent_id]
            if not _vector_subset(child_vec, parent_vec):
                continue
            if not any(p and not c for c, p in zip(child_vec, parent_vec)):
                continue
            parent_map[child_id].add(parent_id)

    return dict(parent_map)


def build_motif_scope_index(
    *,
    probe_smiles: Iterable[str] = DEFAULT_PROBE_SMILES,
) -> Dict[str, Any]:
    compounds_payload = loader.load_organic_compounds()
    groups_payload = loader.load_organic_groups()

    scaffold_parents = _infer_scaffold_parent_map(groups_payload, probe_smiles)

    compound_records: List[Dict[str, str]] = []
    by_key: Dict[Tuple[str, str, str], List[str]] = defaultdict(list)
    for entry in compounds_payload.get("compounds", []) or []:
        if not isinstance(entry, dict):
            continue
        motif_id = _compound_entry_id(entry)
        scaffold = str(entry.get("A") or "").strip()
        substituent = str(entry.get("B") or "").strip()
        template = str(entry.get("template") or "").strip()
        if not motif_id or not scaffold or not substituent:
            continue
        record = {"id": motif_id, "A": scaffold, "B": substituent, "template": template}
        compound_records.append(record)
        by_key[(scaffold, substituent, template)].append(motif_id)

    scope_map: Dict[str, Set[str]] = defaultdict(set)
    for record in compound_records:
        child_id = record["id"]
        child_scaffold = record["A"]
        child_substituent = record["B"]
        child_template = record["template"]
        for parent_scaffold in sorted(scaffold_parents.get(child_scaffold, set())):
            if child_scaffold in _BLOCKED_CHILDREN_BY_PARENT.get(parent_scaffold, set()):
                continue
            parent_ids = by_key.get((parent_scaffold, child_substituent, child_template), [])
            for parent_id in parent_ids:
                if parent_id != child_id:
                    scope_map[parent_id].add(child_id)

    serial_scope_map = {k: sorted(v) for k, v in sorted(scope_map.items()) if v}
    serial_scaffold_parents = {
        child: sorted(parents) for child, parents in sorted(scaffold_parents.items()) if parents
    }

    return {
        "version": "v1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "chemtools.taxonomy.motif_scope_index",
        "rdkit_available": bool(rdkit_available()),
        "probe_smiles": [str(s) for s in probe_smiles],
        "scaffold_parent_map": serial_scaffold_parents,
        "scope_map": serial_scope_map,
    }


def save_motif_scope_index(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=True)
        handle.write("\n")


def build_and_save_motif_scope_index(path: Path) -> Dict[str, Any]:
    payload = build_motif_scope_index()
    save_motif_scope_index(path, payload)
    return payload

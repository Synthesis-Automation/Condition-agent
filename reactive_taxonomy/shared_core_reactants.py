"""Reactant-side retrieval views of validated, persisted reaction-core graphs.

These views nominate observed precedents. They neither predict products nor
establish reaction compatibility; the complete before/after core does that.
"""

from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from pathlib import Path
from typing import Any

from rdkit import Chem, rdBase

from .shared_core_graph import canonical_graph


REACTANT_PROJECTION_VERSION = "shared_core_reactant_projection.v1"


@lru_cache(maxsize=1)
def reactant_projection_definition_hash() -> str:
    """Validate the retrieval-only projection contract and return its identity."""
    path = (
        Path(__file__).with_name("definitions") / f"{REACTANT_PROJECTION_VERSION}.json"
    )
    rules = json.loads(path.read_text(encoding="utf-8"))
    if rules != {
        "definition_id": REACTANT_PROJECTION_VERSION,
        "schema_version": "1.0",
        "source": "shared_reaction_core.v2",
        "retained_side": "before",
        "retain_center_anchors": True,
        "retain_center_distances": True,
        "include_after_states": False,
        "include_edit_events": False,
        "usage": "candidate_retrieval_only",
    }:
        raise ValueError("invalid reactant-side projection definition")
    payload = json.dumps([rules, rdBase.rdkitVersion], sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()


@lru_cache(maxsize=8192)
def reactant_core_key(graph_payload: str) -> str:
    """Project the before-side of a canonical colored incidence graph.

    Center anchors, qualified hydrogen-port labels, before-state atom/bond labels,
    and within-component center distances survive. After states, product context
    and edit events are excluded. No molecular correspondence is inferred.
    """
    vocabulary, encoded = json.loads(graph_payload)
    # This is the serialized colored incidence graph, not a chemical SMARTS.
    graph = Chem.MolFromSmiles(encoded)
    if graph is None:
        raise ValueError("invalid canonical core graph")
    colors = [json.loads(value) for value in vocabulary]
    if any(
        atom.GetAtomicNum() != 0 or not 1 <= atom.GetIsotope() <= len(colors)
        for atom in graph.GetAtoms()
    ):
        raise ValueError("invalid core graph color")
    labels: list[Any] = []
    edges: list[tuple[int, int, Any]] = []
    retained: dict[int, int] = {}
    for atom in graph.GetAtoms():
        kind, label = colors[atom.GetIsotope() - 1]
        if kind != "atom" or not isinstance(label, list):
            continue
        if label[0] == "center" and label[1] is not None:
            before = ["center", label[1]]
        elif label[0] == "before":
            before = label
        else:
            continue
        retained[atom.GetIdx()] = len(labels)
        labels.append(before)
    for atom in graph.GetAtoms():
        kind, label = colors[atom.GetIsotope() - 1]
        if kind != "edge" or not isinstance(label, list) or label[0] != "before":
            continue
        endpoints = [neighbor.GetIdx() for neighbor in atom.GetNeighbors()]
        if len(endpoints) != 2:
            raise ValueError("invalid core incidence edge")
        if all(endpoint in retained for endpoint in endpoints):
            edges.append((retained[endpoints[0]], retained[endpoints[1]], label))
    if not labels:
        raise ValueError("reactant-side core has no observed atoms")
    payload = json.dumps(
        [reactant_projection_definition_hash(), canonical_graph(labels, edges)],
        separators=(",", ":"),
    )
    return "SCRS1:" + hashlib.sha256(payload.encode()).hexdigest()

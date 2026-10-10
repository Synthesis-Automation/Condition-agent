"""Bounded core-first discovery over the canonical fragment retrieval path."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from time import monotonic
from typing import Any, Callable

from rdkit import Chem

from reactive_taxonomy.fragment_search import compile_fragment_query, fragment_embeddings
from reactive_taxonomy.precedent_queries import discovery_policy, plan_precedent_queries
from .fragment_search import FragmentSearchSession, lookup_exact_product, search_fragment_precedents


def _explain_hit(hit: dict[str, Any], step: dict[str, Any], target: Any) -> dict[str, Any]:
    """Rank construction of the original core, not an added peripheral bond."""
    hit = deepcopy(hit)
    positions = [i for i, atom in enumerate(step["target_atom_ids"]) if atom in step["core_atom_ids"]]
    relationships = set()
    for match in hit["matches"]:
        selected = {match["query_to_original_atoms"][i] for i in positions}
        unresolved = "unresolved" in match["relationships"]
        if unresolved:
            relationships.add("unresolved")
            continue
        local = []
        for witness in match["witnesses"]:
            endpoints = witness["product_atoms"]
            if not any(i in selected for i in endpoints):
                continue
            relation = ("boundary_changed" if not all(i in selected for i in endpoints)
                        else "constructed" if witness["kind"] == "formed" else "modified")
            relationships.add(relation)
            local.append(witness)
        if not local:
            relationships.add("carried_through")
    precedence = discovery_policy()["relationship_order"]
    relationship = min(relationships or {"unresolved"}, key=precedence.index)
    # Independently inspect matched structures; do not infer source/target atom correspondence.
    product = Chem.MolFromSmiles(hit["matched_molecule_smiles"])
    query = compile_fragment_query(step["query"], "smarts", step["topology"])
    embeddings, truncated = fragment_embeddings(query, product, maximum=8)
    differences = []
    if embeddings:
        for i in positions:
            atom = target.GetAtomWithIdx(step["target_atom_ids"][i])
            if atom.GetAtomicNum() == 7:
                values = sorted({product.GetAtomWithIdx(e[i]).GetTotalNumHs() for e in embeddings})
                if values != [atom.GetTotalNumHs()]:
                    differences.append({"target_atom_id": atom.GetIdx(), "element": "N",
                                        "target_hydrogens": atom.GetTotalNumHs(), "source_hydrogen_counts": values})
    hit["discovery"] = {
        "core_id": step["core_id"], "level": step["level"], "query": step["query"],
        "query_format": "smarts", "topology": step["topology"],
        "core_relationship": relationship, "core_atom_count": len(positions),
        "context_atom_count": len(step["target_atom_ids"]),
        "target_core_atom_ids": step["core_atom_ids"], "nitrogen_hydrogen_differences": differences,
        "alignment_ambiguous": len(embeddings) > 1 or truncated,
        "explanation": f"Same selected core ({len(positions)} atoms); {relationship.replace('_', ' ')} evidence. "
                       f"Matches {len(step['target_atom_ids'])} target-context atoms. Peripheral substitution, nonstereo hydrogen counts and extra ring fusion were unconstrained; specified stereo hydrogen counts were retained.",
    }
    return hit


def find_synthesis_precedents(
    index_path: str | Path, target_smiles: str, *, limit: int = 10,
    timeout_seconds: int = 90, progress: Callable[[dict[str, Any]], None] | None = None,
) -> dict[str, Any]:
    """Search informative cores, refine dense results, and retain useful earlier hits.

    Deadlines are cooperative between bounded search batches, including library
    loading time. A complete result means the bounded policy finished, not that
    all possible cores, analogues or constructions have been exhausted.
    """
    policy = discovery_policy()
    if type(limit) is not int or not 1 <= limit <= 10:
        raise ValueError("limit must be an integer between 1 and 10")
    if type(timeout_seconds) is not int or not 1 <= timeout_seconds <= policy["max_timeout_seconds"]:
        raise ValueError("timeout_seconds must be an integer between 1 and 120")
    started = monotonic()
    plan = plan_precedent_queries(target_smiles)
    target = Chem.MolFromSmiles(plan["target_smiles"])
    exact = lookup_exact_product(index_path, plan["target_smiles"], limit=limit)
    attempts, candidates = [], {}
    timed_out = False
    with FragmentSearchSession(index_path) as session:
        if monotonic() - started < timeout_seconds:
            session.load_library()
        for ladder in plan["ladders"]:
            # Focused contexts are siblings. Only a fully enumerated query
            # whose selected target atoms are a subset may constrain another.
            enumerated = []
            focus = [step for step in plan["focused_queries"]
                     if step["core_id"] == ladder[0]["core_id"]]
            insertion = next((i for i, step in enumerate(ladder)
                              if step["level"] in {"neighbor_context", "target_context"}), len(ladder))
            schedule = ladder[:insertion] + focus + ladder[insertion:]
            for position, step in enumerate(schedule):
                selected_atoms = set(step["target_atom_ids"])
                reusable = [(atoms, ids) for atoms, ids in enumerated
                            if atoms <= selected_atoms]
                parent_ids = max(reusable, key=lambda item: len(item[0]))[1] if reusable else None
                remaining = timeout_seconds - (monotonic() - started)
                if remaining < 1:
                    timed_out = True
                    break
                search = search_fragment_precedents(
                    index_path, step["query"], "smarts", step["topology"], 10,
                    min(30, int(remaining)), target_smiles=plan["target_smiles"],
                    progress=progress, _session=session, _candidate_ids=parent_ids, _sample_broad_matches=True,
                )
                query_id = search["query"]["query_id"]
                current_ids = session.candidates.get(query_id)
                construction_refs = set()
                for hit in search["hits"]:
                    explained = _explain_hit(hit, step, target)
                    details = explained["discovery"]
                    if details["core_relationship"] == "constructed" and hit["reference_id"]:
                        construction_refs.add(hit["reference_id"])
                    key = hit["observation_id"]
                    existing = candidates.get(key)
                    precedence = policy["relationship_order"]
                    def rank(item: dict[str, Any]) -> tuple[Any, ...]:
                        d = item["discovery"]
                        return (precedence.index(d["core_relationship"]), -d["core_atom_count"],
                                -d["context_atom_count"], item["observation_id"])
                    if existing is None or rank(explained) < rank(existing):
                        candidates[key] = explained
                count = search["counts"]["products"]
                if search["stop_reason"] == "deadline":
                    decision = "deadline"
                    timed_out = True
                elif current_ids is not None and not current_ids:
                    has_alternative = any(not selected_atoms <= set(later["target_atom_ids"])
                                          for later in schedule[position + 1:])
                    decision = ("no_matches_try_next_context" if has_alternative
                                else "no_matches_try_next_core")
                elif search["search_status"] == "too_broad" or count["value"] > policy["refine_above_products"]:
                    decision = "add_context"
                elif len(construction_refs) >= policy["sufficient_construction_references"]:
                    decision = "enough_construction_references"
                elif construction_refs:
                    decision = "retain_construction_and_check_context"
                else:
                    decision = "inspect_more_specific_context"
                attempts.append({**step, "query_id": query_id, "search_status": search["search_status"],
                                 "stop_reason": search["stop_reason"], "counts": search["counts"],
                                 "construction_references_in_returned_hits": len(construction_refs),
                                 "decision": decision, "execution": search["execution"],
                                 "returned_observation_ids": [h["observation_id"] for h in search["hits"]]})
                if decision in {"deadline", "no_matches_try_next_core", "enough_construction_references"}:
                    break
                if current_ids is not None:
                    enumerated.append((selected_atoms, current_ids))
            if timed_out:
                break
        loads = session.library_loads
        manifest = session.manifest
    precedence = policy["relationship_order"]
    ordered = sorted(candidates.values(), key=lambda h: (
        precedence.index(h["discovery"]["core_relationship"]), -h["discovery"]["core_atom_count"],
        -h["discovery"]["context_atom_count"], h["observation_id"]))
    chosen, repeats, seen = [], [], set()
    for hit in ordered:
        # Within each evidence tier, prevent one reference from occupying the list.
        key = (hit["discovery"]["core_relationship"], hit["reference_id"] or hit["observation_id"])
        (repeats if key in seen else chosen).append(hit)
        seen.add(key)
    selected = sorted(chosen + repeats, key=lambda h: precedence.index(h["discovery"]["core_relationship"]))[:limit]
    for i, hit in enumerate(selected):
        hit["inspect_paths"] = {key: ["result", "hits", i, key] for key in ("matches", "record", "procedures")}
    partial = timed_out or any(a["search_status"] != "complete" for a in attempts)
    return {"schema_version": "synthesis_precedent_discovery.v1", "definition_version": policy["definition_version"],
            "target_smiles": plan["target_smiles"], "search_status": "partial" if partial else "complete",
            "stop_reason": "deadline" if timed_out else "bounded_policy_finished",
            "index_id": manifest["index_id"], "source_scope": manifest["source_scope"],
            "source_coverage_complete": manifest["source_coverage_complete"], "plan": plan,
            "exact_target": exact, "attempts": attempts, "hits": selected, "returned_count": len(selected),
            "ranking_scope": "retrieved_candidates", "ranking_policy": policy["relationship_order"],
            "execution": {"elapsed_seconds": round(monotonic() - started, 6), "library_loads": loads,
                          "candidate_reuses": sum(a["execution"]["candidate_set_reused"] for a in attempts)},
            "limitations": ["Ranked discovery hypotheses, not validated synthetic steps or condition recommendations.",
                            "Nonstereo hydrogen counts and omitted peripheral groups are unconstrained; specified stereo hydrogen counts are retained. Inspect source differences.",
                            "Unresolved atom correspondence cannot establish construction of the selected core.",
                            "Ranking uses at most ten hydrated hits per query, not every indexed observation.",
                            "Broad queries inspect a bounded product sample; their ranks and counts are incomplete.",
                            "Zero returned hits do not establish absence of routes or related literature."]}

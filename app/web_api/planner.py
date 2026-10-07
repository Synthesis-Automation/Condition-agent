"""Thin HTTP composition for deterministic interactive route planning."""

from __future__ import annotations

import json
from typing import Any, Literal

from fastapi import APIRouter, HTTPException, Request
from pydantic import Field

from core_retrosynthesis.interactive_planning import (
    PlanningChoice,
    candidate_for_choice,
    edit_session,
    find_node,
    normalize_settings,
    record_conditions,
    record_search,
    record_stock,
    restore_session,
    session_response,
    start_session,
)
from .contracts import (
    StrictRequest,
    RetrosynthesisRequest,
    RetrosynthesisConditionsRequest,
    envelope,
)
from .runtime import error_payload

router = APIRouter(prefix="/api/v1/retrosynthesis/planner")


class PlanningRequest(StrictRequest):
    """One explicit user action on a browser-owned planning session."""

    action: Literal[
        "start",
        "restore",
        "search",
        "select",
        "clear",
        "remove",
        "stop",
        "reopen",
        "undo",
        "redo",
        "conditions",
        "stock",
    ]
    session: dict[str, Any] | None = None
    target_smiles: str | None = Field(default=None, max_length=20_000)
    node_id: str = Field(default="root", min_length=1, max_length=200)
    settings: dict[str, Any] = Field(default_factory=dict)
    search_id: str | None = Field(default=None, max_length=200)
    strategy_index: int = Field(default=0, ge=0, strict=True)
    realization_index: int = Field(default=0, ge=0, strict=True)


@router.post("")
def planning_action(payload: PlanningRequest, request: Request) -> dict[str, Any]:
    """Compose session edits with the same single-step and condition services."""
    try:
        if payload.action == "start":
            session = start_session(payload.target_smiles or "")
        else:
            if payload.session is None:
                raise ValueError("A planning session is required")
            if len(json.dumps(payload.session)) > 20_000_000:
                raise ValueError("Planning session exceeds the 20 MB import limit")
            session = restore_session(payload.session)
            if payload.action == "search":
                node = find_node(session.root, payload.node_id)
                settings = normalize_settings(payload.settings)
                result = request.app.state.runtime.retrosynthesize(
                    RetrosynthesisRequest(target_smiles=node.smiles, **settings)
                )
                session = record_search(session, node.node_id, settings, result)
            elif payload.action == "stock":
                node = find_node(session.root, payload.node_id)
                evidence = request.app.state.runtime.planning_stock(node.smiles)
                session = record_stock(session, node.node_id, evidence)
            elif payload.action == "conditions":
                node = find_node(session.root, payload.node_id)
                if node.choice is None:
                    raise ValueError("Choose a step before looking up conditions")
                search, candidate, _ = candidate_for_choice(session, node.choice)
                evidence = request.app.state.runtime.retrosynthesis_conditions(
                    RetrosynthesisConditionsRequest(
                        reaction_smiles=candidate.get("condition_query_reaction_smiles")
                        or candidate["proposed_reaction_smiles"],
                        library_mode=search.settings["library_mode"],
                        top_k=3,
                        preferred_reaction_ids=candidate.get(
                            "condition_precedent_reaction_ids", []
                        ),
                        starting_materials=candidate["precursor_smiles"],
                        intended_product=node.smiles,
                        operator_hint=candidate.get("operator_id"),
                        use_forward_validation=search.settings[
                            "use_forward_validation"
                        ],
                        include_l0=search.settings["include_l0"],
                    )
                )
                session = record_conditions(session, node.node_id, evidence)
            elif payload.action != "restore":
                choice = None
                if payload.action == "select":
                    if not payload.search_id:
                        raise ValueError("Select a retained search result")
                    choice = PlanningChoice(
                        payload.search_id,
                        payload.strategy_index,
                        payload.realization_index,
                    )
                session = edit_session(session, payload.action, payload.node_id, choice)
        return envelope(session_response(session))
    except (ValueError, LookupError, TypeError) as exc:
        raise HTTPException(status_code=422, detail=error_payload(exc)) from exc
    except (FileNotFoundError, RuntimeError) as exc:
        raise HTTPException(status_code=503, detail=error_payload(exc)) from exc

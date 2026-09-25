"""Direct scientific operations over existing public domain implementations."""

from __future__ import annotations

from dataclasses import asdict
import gzip
import inspect
import json
from pathlib import Path
from typing import Any, Mapping

from .store import canonical_bytes, InvestigationStore


class ScientificOperations:
    """Lazy application adapters; no LLM, scientific rule, or network dependency."""

    NAMES = (
        "analyze_reaction", "analyze_molecule", "recommend_conditions",
        "get_precedents", "get_procedures", "resolve_recipe", "assess_recipe",
        "inspect_condition_precedents", "propose_condition_adaptation",
        "plan_routes", "revise_routes",
    )

    def __init__(self, store: InvestigationStore) -> None:
        self.store = store
        self._recommender: Any = None
        self._planner: Any = None
        self._route_runs: dict[bytes, Any] = {}

    def catalog(self) -> list[dict[str, str]]:
        """Describe the explicit operations available to a workspace client."""
        return [{"name": name, "signature": str(inspect.signature(getattr(self, name))),
                 "description": inspect.getdoc(getattr(self, name)) or ""}
                for name in self.NAMES]

    def invoke(self, operation: str, arguments: Mapping[str, Any]) -> Any:
        """Dispatch only named application operations, never caller-supplied imports."""
        if operation not in self.NAMES:
            raise ValueError(f"Unknown scientific operation: {operation}")
        return getattr(self, operation)(**arguments)

    def _path(self, name: str) -> Path:
        entry = self.store.manifest["baseline"]["artifacts"].get(name)
        if not entry or entry["status"] != "present":
            raise FileNotFoundError(f"Capability requires configured artifact: {name}")
        return Path(entry["path"])

    def _conditions(self) -> Any:
        if self._recommender is None:
            from condition_recommender import GenericConditionRecommender

            self._recommender = GenericConditionRecommender.from_path(
                self._path("condition_index"), shared_core_path=self._path("shared_core_index"),
            )
        return self._recommender

    def analyze_reaction(self, reaction_smiles: str) -> Any:
        """Observe and interpret structures, preserving ambiguity and supplied mapping evidence."""
        from reactive_taxonomy import featurize_reaction

        return featurize_reaction(reaction_smiles)

    def analyze_molecule(self, smiles: str) -> Any:
        """Return the canonical graph-derived target audit, including unassessed properties."""
        from reactive_taxonomy import audit_target

        return audit_target(smiles)

    def recommend_conditions(
        self, reaction_smiles: str, top_k: int = 5, search_scope: str = "automatic",
    ) -> Any:
        """Run canonical shared-core retrieval, compatibility, aggregation, and ranking."""
        if type(top_k) is not int or not 1 <= top_k <= 50:
            raise ValueError("top_k must be an integer between 1 and 50")
        return self._conditions().recommend(
            reaction_smiles, top_k=top_k, search_scope=search_scope,
        )

    def get_precedents(
        self, reaction_ids: list[str], offset: int = 0, limit: int = 20,
    ) -> dict[str, Any]:
        """Retrieve complete indexed records by ID; these are reduced records, not raw source files."""
        if not isinstance(reaction_ids, list) or not all(isinstance(item, str) for item in reaction_ids):
            raise ValueError("reaction_ids must be a list of strings")
        if len(reaction_ids) > 100 or type(limit) is not int or not 1 <= limit <= 100:
            raise ValueError("At most 100 reaction IDs and limit 1..100 are supported")
        if type(offset) is not int or offset < 0:
            raise ValueError("offset must be a nonnegative integer")
        index = self._conditions().index
        positions = tuple(dict.fromkeys(
            position for identity in reaction_ids
            for position in index.reaction_ids.get(identity, ())
        ))
        return {
            "records": [asdict(row) for row in index.select(positions[offset:offset + limit])],
            "missing_reaction_ids": [identity for identity in reaction_ids if not index.reaction_ids.get(identity)],
            "total": len(positions), "offset": offset,
            "next_offset": offset + limit if offset + limit < len(positions) else None,
            "record_scope": "indexed_fields_only",
            "index_scope": index.precedent_scope.value,
        }

    def get_procedures(self, reaction_ids: list[str]) -> dict[str, Any]:
        """Read all matching observed procedure records; missing text remains missing."""
        if not isinstance(reaction_ids, list) or not all(isinstance(item, str) for item in reaction_ids):
            raise ValueError("reaction_ids must be a list of strings")
        if len(reaction_ids) > 100:
            raise ValueError("At most 100 reaction IDs are supported")
        selected = set(reaction_ids)
        path = self._path("procedure_catalog")
        records = []
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                value = json.loads(line)
                if value.get("reaction_id") in selected:
                    records.append(value)
        found = {record["reaction_id"] for record in records}
        return {"records": records, "missing_reaction_ids": sorted(selected - found),
                "source_path": str(path), "origin": "source_report"}

    def resolve_recipe(
        self, components: list[dict[str, Any]], temperature_c: float | None = None,
        time_h: float | None = None, concentration_m: float | None = None,
        atmosphere: str | None = None,
    ) -> Any:
        """Normalize typed identifiers and contextual roles through condition_registry."""
        from condition_registry import ConditionComponentInput, build_resolved_recipe_from_inputs

        return build_resolved_recipe_from_inputs(
            [ConditionComponentInput(**item) for item in components],
            temperature_c=temperature_c, time_h=time_h, concentration_m=concentration_m,
            atmosphere=atmosphere,
        )

    def assess_recipe(self, reaction_smiles: str, recipe: dict[str, Any]) -> Any:
        """Assess compatibility; unknown/invalid_input are not chemical conflicts or success."""
        from condition_recommender import assess_reaction_recipe

        return assess_reaction_recipe(reaction_smiles, recipe)

    def inspect_condition_precedents(
        self, reaction_smiles: str, reaction_ids: list[str], offset: int = 0, limit: int = 20,
    ) -> dict[str, Any]:
        """Compare selected indexed observations and link source procedures without merging them."""
        from condition_recommender import compare_condition_evidence
        from condition_recommender.generic_indexing import GenericIndexedReaction

        page = self.get_precedents(reaction_ids, offset=offset, limit=limit)
        result = asdict(compare_condition_evidence(
            reaction_smiles, [GenericIndexedReaction(**row) for row in page["records"]],
        ))
        try:
            procedures = self.get_procedures([row["reaction_id"] for row in page["records"]])
        except FileNotFoundError:
            procedures = {"records": [], "availability": "catalog_unavailable"}
        else:
            procedures["availability"] = "catalog_available"
        for item in result["precedents"]:
            observation = item["observation"]
            candidates = [record for record in procedures["records"]
                          if record["reaction_id"] == observation["reaction_id"]]
            item["procedure_observations"] = [record for record in candidates
                                               if record.get("observation_id") == observation["observation_id"]]
            item["reaction_level_procedures"] = [record for record in candidates if not record.get("observation_id")]
            item["procedure_link_scope"] = "exact_observation_id_or_explicitly_unassigned_reaction_record"
        result["procedure_catalog"] = procedures
        result["page"] = {key: value for key, value in page.items() if key != "records"}
        return result

    def propose_condition_adaptation(
        self, source_ref: str, observation_id: str, components: list[dict[str, Any]],
        operating_conditions: dict[str, Any], change_reasons: dict[str, str],
        evidence_refs: list[str], assumptions: list[str], risks: list[str],
    ) -> dict[str, Any]:
        """Record a proposed complete recipe, attributed changes and canonical assessment.

        Use a completed inspect_condition_precedents call. Supply the entire new
        component list and operating values; omitted operating values remain unknown.
        Reasons must cover every changed bucket/operating field. This records an
        agent hypothesis, not a recommendation admitted to a dataset or proof of transfer.
        """
        from .adaptation import record_adaptation

        return record_adaptation(self, source_ref, observation_id, components,
                                 operating_conditions, change_reasons, evidence_refs, assumptions, risks)

    def _route_planner(self, include_conditions: bool) -> Any:
        from dataclasses import replace
        from chem_coworker.retrosynthesis import RetrosynthesisCoworker
        from chem_coworker.multistep import MultistepRetrosynthesisCoworker

        if self._planner is None:
            retro = RetrosynthesisCoworker.from_path(self._path("retro_library"))
            self._planner = MultistepRetrosynthesisCoworker.from_retrosynthesis_coworker(
                retro, stock_path=self._path("stock_index"),
            )
        return replace(self._planner, condition_recommender=self._conditions() if include_conditions else None)

    def _plan(self, settings: Mapping[str, Any], exclusions: tuple[Any, ...] = ()) -> Any:
        from chem_coworker.contracts import MultistepRetrosynthesisRequest

        if "review" in settings:
            raise ValueError("External workspace operations do not invoke an internal LLM review")
        request = MultistepRetrosynthesisRequest(**settings)
        return self._route_planner(request.include_conditions).plan(request, candidate_exclusions=exclusions)

    def _route_payload(self, response: Any) -> dict[str, Any]:
        from core_retrosynthesis import collect_route_refinement_issues, enumerate_route_repair_proposals, verify_planned_route

        routes = (*response.result.routes, *response.result.partial_routes) if response.result else ()
        assessments = []
        for route in routes:
            issues = collect_route_refinement_issues(route)
            assessments.append({
                "route_id": route.route_id, "verification": verify_planned_route(route).to_dict(),
                "issues": [issue.to_dict() for issue in issues],
                "repair_proposals": [proposal.to_dict() for issue in issues
                                     for proposal in enumerate_route_repair_proposals(route, issue)],
            })
        return {"response": response.to_dict(), "assessments": assessments}

    def plan_routes(self, settings: dict[str, Any]) -> dict[str, Any]:
        """Search bounded routes with the existing planner, then expose issues and repair choices."""
        response = self._plan(settings)
        payload = self._route_payload(response)
        self._route_runs[canonical_bytes(payload)] = (response, settings, ())
        return payload

    def revise_routes(self, source_ref: str, intent: dict[str, Any]) -> dict[str, Any]:
        """Revise an evidenced route choice; retain the original and recheck the new search.

        A fresh session replays the saved source call to recover typed domain objects,
        and requires exact result parity before applying the requested revision.
        """
        from core_retrosynthesis import (
            RouteRefinementIntent, build_route_refinement_plan,
            collect_route_refinement_issues, summarize_route_refinement,
        )

        source = self.store.read_artifact(source_ref)
        if source.get("operation") not in {"plan_routes", "revise_routes"} or source.get("execution_status") != "completed":
            raise ValueError("source_ref must identify a completed route call")
        key = canonical_bytes(source["result"])
        if key not in self._route_runs:
            replay = self.invoke(source["operation"], source["arguments"])
            if canonical_bytes(replay) != key:
                raise ValueError("Source route replay differs from saved evidence")
        response, settings, exclusions = self._route_runs[key]
        if response.result is None:
            raise ValueError("Source search has no route result")
        value = dict(intent)
        value["issue_ids"] = tuple(value.get("issue_ids", ()))
        request = RouteRefinementIntent(**value)
        route = next((item for item in (*response.result.routes, *response.result.partial_routes)
                      if item.route_id == request.source_route_id), None)
        if route is None:
            raise ValueError("Unknown source route ID")
        plan = build_route_refinement_plan(route, request, collect_route_refinement_issues(route))
        updated_exclusions = (*exclusions, plan.exclusion)
        revised = self._plan(settings, updated_exclusions)
        payload = self._route_payload(revised)
        payload["source_ref"] = source_ref
        payload["intent"] = request.to_dict()
        payload["refinement"] = (
            summarize_route_refinement(route, request, revised.result).to_dict()
            if revised.result is not None else None
        )
        self._route_runs[canonical_bytes(payload)] = (revised, settings, updated_exclusions)
        return payload

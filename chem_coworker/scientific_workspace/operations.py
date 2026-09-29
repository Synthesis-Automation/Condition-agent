"""Direct scientific operations over existing public domain implementations."""

from __future__ import annotations

from dataclasses import asdict
import gzip
import inspect
import json
from pathlib import Path
from typing import Any, Mapping

from .store import InvestigationStore


class ScientificOperations:
    """Lazy application adapters; no LLM, scientific rule, or network dependency."""

    NAMES = (
        "analyze_reaction", "analyze_molecule", "recommend_conditions",
        "compare_molecules", "inspect_reactive_sites",
        "get_precedents", "get_procedures", "resolve_recipe", "assess_recipe",
        "search_fragment_precedents", "suggest_search_fragments",
        "inspect_condition_precedents", "propose_condition_adaptation",
        "disconnect_target",
        "assess_route_step", "assess_route_proposal",
        "assess_route_step_forward",
        "inspect_route_step", "inspect_step_precedents", "revise_route_branch", "compare_route_proposals",
    )

    def __init__(self, store: InvestigationStore) -> None:
        self.store = store
        self._recommender: Any = None
        self._proposal_library: Any = None

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

    def compare_molecules(
        self, left_smiles: str, right_smiles: str, core_smiles: str | None = None,
        timeout_seconds: int = 2,
    ) -> Any:
        """Compare target/precedent cores, substituents and stereo; alignments are not reaction maps.

        Optional core_smiles anchors the comparison; otherwise MCS has a bounded
        search. Inspect ambiguity, partial coverage and stereo warnings before
        transferring a precedent. No index, forward model or network is needed.
        """
        from reactive_taxonomy import compare_molecules

        return compare_molecules(left_smiles, right_smiles, core_smiles, timeout_seconds)

    def inspect_reactive_sites(
        self, smiles: str, selected_atom_ids: list[int] | None = None, radius: int = 2,
    ) -> Any:
        """Inspect local motifs, site descriptors and other possible sites; no selectivity prediction.

        IDs use returned canonical SMILES, also used by suggest_search_fragments.
        Omit selection to inspect all sites and get atom IDs; choose a region for
        focused follow-up. Stereo gaps and descriptor provenance remain explicit.
        """
        from reactive_taxonomy import inspect_reactive_sites

        return inspect_reactive_sites(smiles, selected_atom_ids, radius)

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

    def suggest_search_fragments(
        self, target_smiles: str, limit: int = 5, selected_atom_ids: list[int] | None = None,
    ) -> dict[str, Any]:
        """Optionally suggest up to five target-derived search regions, without searching.

        Candidates overlap and are not precursors or synthesis recommendations.
        Inspect reasons, boundaries and cautions, then choose your own query for
        search_fragment_precedents. Explicit selections use the returned canonical
        target's zero-based atom IDs and are never silently expanded.
        """
        from reactive_taxonomy.search_fragments import suggest_search_fragments

        return suggest_search_fragments(target_smiles, limit, selected_atom_ids).to_dict()

    def search_fragment_precedents(
        self, query: str, query_format: str = "smiles", topology: str = "preserve_rings",
        limit: int = 5, timeout_seconds: int = 10, target_smiles: str | None = None,
    ) -> dict[str, Any]:
        """Find product cores and local construction evidence in a prebuilt fragment_index.

        Supply one connected core. SMILES permits peripheral substitution while
        preserving rings; explicit SMARTS/subgraph permits deliberate broadening.
        Supply target_smiles when the query represents a core of that target;
        mismatches are recorded as errors before scanning, not zero-hit results.
        Inspect saved hits for source records and procedure chunks. No automatic
        mapping, forward check, route expansion, or index rebuild is performed.
        """
        from .fragment_search import run_fragment_search

        return run_fragment_search(self, {"query": query, "query_format": query_format,
                                         "topology": topology, "limit": limit,
                                         "timeout_seconds": timeout_seconds,
                                         "target_smiles": target_smiles})

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

    def disconnect_target(
        self, target_smiles: str, top_k: int = 3, max_realizations_per_strategy: int = 3,
        max_templates_to_apply: int = 40, max_candidates_to_validate: int = 10,
        use_context: bool = True, include_l0: bool = True, include_conditions: bool = False,
        condition_top_k: int = 3, condition_minimum_pool_size: int | None = None,
        unrestricted_condition_fallback: bool = False,
    ) -> Any:
        """Generate single-step strategies for one target; never expand precursor branches.

        The agent chooses realizations, subsequent targets and stopping decisions.
        Requires retro_library only unless conditions are requested. Preserves the
        existing engine's evidence and warnings; internal LLM review is disabled.
        """
        from chem_coworker.contracts import RetrosynthesisRequest
        from chem_coworker.retrosynthesis import RetrosynthesisCoworker

        request = RetrosynthesisRequest(
            target_smiles=target_smiles, top_k=top_k,
            max_realizations_per_strategy=max_realizations_per_strategy,
            max_templates_to_apply=max_templates_to_apply,
            max_candidates_to_validate=max_candidates_to_validate,
            use_context=use_context, include_l0=include_l0,
            include_conditions=include_conditions, condition_top_k=condition_top_k,
            condition_minimum_pool_size=condition_minimum_pool_size,
            unrestricted_condition_fallback=unrestricted_condition_fallback,
        )
        coworker = RetrosynthesisCoworker(
            library=self._external_route_library(), library_path=self._path("retro_library"),
            condition_recommender=self._conditions() if include_conditions else None,
        )
        return coworker.disconnect(request)

    def _external_route_library(self) -> Any:
        """Load the recorded operator library without requiring a stock index."""
        from core_retrosynthesis import load_generic_library

        if self._proposal_library is None:
            self._proposal_library = load_generic_library(self._path("retro_library"))
        return self._proposal_library

    def assess_route_step(
        self, proposal: dict[str, Any], include_conditions: bool = False,
        include_forward: bool = False, evidence_refs: list[str] | None = None,
    ) -> dict[str, Any]:
        """Assess target_smiles/precursor_smiles through canonical external-proposal gates.

        Optional mapped_reaction_smiles, proposed_conditions (resolved recipe),
        and sources remain untrusted input. Unknown is not impossible. Conditions
        retrieval is opt-in; supplied recipes are assessed separately. Forward
        prediction is a separate bounded assess_route_step_forward call after
        a saved route assessment. include_forward must remain False.
        """
        from .route_investigation import assess_step

        return assess_step(self, proposal, include_conditions, include_forward, evidence_refs)

    def assess_route_proposal(
        self, proposal: dict[str, Any], unavailable_starting_materials: list[str] | None = None,
        include_conditions: bool = False, include_forward: bool = False,
        evidence_refs: list[str] | None = None,
    ) -> dict[str, Any]:
        """Assess a target_smiles and steps graph; each step needs external_step_id and step proposal fields.

        Preserve unsupported steps as hypotheses. Declared unavailable starting
        materials are graph-matched against route leaves; this is not a stock lookup.
        The returned artifact is the source_ref for inspection, revision and comparison.
        include_forward must remain False; use assess_route_step_forward only
        for a consequential uncertainty in one eligible saved step.
        """
        from .route_investigation import assess_route

        return assess_route(self, proposal, unavailable_starting_materials, include_conditions, include_forward, evidence_refs)

    def assess_route_step_forward(
        self, source_ref: str, step_id: str, question: str, timeout_seconds: int = 30,
    ) -> dict[str, Any]:
        """Optionally challenge one eligible step of a saved route within 1..30 seconds.

        State a question about competing products that could change the route
        decision. Requires a baseline-pinned prebuilt forward_library. Loading,
        source compatibility checks and prediction share the deadline; timed-out
        or failed attempts retain diagnostics without changing the source route.
        This is separate from the mandatory structural checks in retrosynthesis.
        """
        from .forward_check import assess_step_forward

        return assess_step_forward(self, source_ref, step_id, question, timeout_seconds)

    def inspect_route_step(self, source_ref: str, step_id: str) -> dict[str, Any]:
        """Inspect a recorded proposal step's gates, neighboring steps, molecule audits and recipe assessment."""
        from .route_investigation import inspect_step

        return inspect_step(self, source_ref, step_id)

    def inspect_step_precedents(
        self, source_ref: str, step_id: str | None = None, realization_id: str | None = None,
        offset: int = 0, limit: int = 3,
    ) -> dict[str, Any]:
        """Inspect actual supporting reactions for one saved realization or assessed step.

        Use realization_id for disconnect_target, step_id for a route assessment
        or revision, and neither for assess_route_step. Follow page.next_offset.
        Source reactions, product comparison, scoped counts, and available source
        conditions remain evidence, not proof that the proposed step will work.
        """
        from .step_precedents import inspect_step_precedents

        return inspect_step_precedents(self, source_ref, step_id, realization_id, offset, limit)

    def revise_route_branch(
        self, source_ref: str, remove_step_ids: list[str], replacement_steps: list[dict[str, Any]],
        reason: str, evidence_refs: list[str] | None = None, assumptions: list[str] | None = None,
        risks: list[str] | None = None,
    ) -> dict[str, Any]:
        """Explicitly replace/add steps, preserve the source, and reassess every step and route topology.

        Use a completed proposal assessment/revision source_ref. Reusing an ID
        requires explicitly removing it; remove_step_ids=[] extends a leaf branch.
        Supply reason and nonempty risks. Condition settings and declared material
        constraints are inherited. Optional forward checks remain separate.
        Invalid revisions remain inspectable, not accepted.
        """
        from .route_investigation import revise_branch

        return revise_branch(self, source_ref, remove_step_ids, replacement_steps, reason, evidence_refs, assumptions, risks)

    def compare_route_proposals(self, source_refs: list[str]) -> dict[str, Any]:
        """Compare 2–5 recorded alternatives for the same target, checks and constraints; no automatic winner."""
        from .route_investigation import compare_routes

        return compare_routes(self, source_refs)

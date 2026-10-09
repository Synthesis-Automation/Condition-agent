"""Direct scientific operations over existing public domain implementations."""

from __future__ import annotations

import gzip
import inspect
import json
from dataclasses import asdict
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, Mapping

from ..core.operation_contracts import OperationDefinition
from ..core.store import InvestigationStore

if TYPE_CHECKING:
    from .fragment_search import FragmentWorker


def _fragment_replay_result(result: Any) -> Any:
    """Retain fragment chemistry while excluding per-attempt execution telemetry."""
    if isinstance(result, dict):
        return {key: value for key, value in result.items() if key != "execution"}
    return result


def _discovery_replay_result(result: Any) -> Any:
    """Exclude attempt timings while preserving decisions, scope and chemistry."""
    if isinstance(result, dict):
        return {key: _discovery_replay_result(value) for key, value in result.items() if key != "execution"}
    if isinstance(result, list):
        return [_discovery_replay_result(value) for value in result]
    return result


def _forward_replay_result(result: Any) -> Any:
    """Compare forward stage outcomes without attempt-specific logs and timings."""
    if not isinstance(result, dict):
        return result
    execution = dict(result.get("execution", {}))
    execution.pop("diagnostics", None)
    execution["timings"] = {
        key: value for key, value in execution.get("timings", {}).items()
        if key != "elapsed_seconds"
    }
    execution["stages"] = [
        {key: value for key, value in stage.items()
         if key not in {"elapsed_seconds", "duration_seconds"}}
        for stage in execution.get("stages", [])
    ]
    return {**result, "execution": execution}


class ScientificOperations:
    """Lazy application adapters; no LLM, scientific rule, or network dependency."""

    DEFINITIONS = (
        OperationDefinition("analyze_reaction"),
        OperationDefinition("analyze_molecule"),
        OperationDefinition("prepare_literature_reaction", contract_version="3", evidence_arguments=("source_ref", "indexed_ref"),
                            result_evidence_field="evidence_refs"),
        OperationDefinition("recommend_conditions", required_artifacts=("condition_index", "shared_core_index")),
        OperationDefinition("generate_weak_label_screening_array",
                            required_artifacts=("weak_label_records", "weak_label_recipe_catalog")),
        OperationDefinition("compare_molecules"),
        OperationDefinition("inspect_reactive_sites"),
        OperationDefinition("assess_starting_material"),
        OperationDefinition("get_precedents", contract_version="2", required_artifacts=("condition_index", "shared_core_index")),
        OperationDefinition("get_observations", required_artifacts=("evidence_catalog",)),
        OperationDefinition("get_routes", required_artifacts=("route_catalog",)),
        OperationDefinition("get_procedures", required_artifacts=("procedure_catalog",)),
        OperationDefinition("resolve_recipe"),
        OperationDefinition("assess_recipe", contract_version="2"),
        OperationDefinition("assess_proposed_recipe", contract_version="2", evidence_arguments=("evidence_refs",)),
        OperationDefinition("search_captured_sources", evidence_arguments=("source_refs",),
                            result_evidence_field="evidence_refs"),
        OperationDefinition("inspect_route_inputs", contract_version="2", evidence_arguments=("source_ref", "source_refs"),
                            result_evidence_field="evidence_refs"),
        OperationDefinition("search_fragment_precedents", required_artifacts=("fragment_index",),
                            evidence_arguments=("query_variant_ref",),
                            execution_status_field="execution_status",
                            replay_comparison="scientific_result_excluding_fragment_execution_telemetry",
                            replay_projection=_fragment_replay_result),
        OperationDefinition("find_synthesis_precedents", required_artifacts=("fragment_index",),
                            execution_status_field="execution_status",
                            replay_comparison="scientific_result_excluding_discovery_execution_telemetry",
                            replay_projection=_discovery_replay_result),
        OperationDefinition("suggest_search_fragments"),
        OperationDefinition("propose_fragment_queries"),
        OperationDefinition("investigate_fragment_precedent", evidence_arguments=("source_ref",)),
        OperationDefinition("inspect_condition_precedents", contract_version="2",
                            required_artifacts=("condition_index", "shared_core_index")),
        OperationDefinition("propose_condition_adaptation", evidence_arguments=("source_ref", "evidence_refs")),
        OperationDefinition(
            "disconnect_target", contract_version="2", required_artifacts=("retro_library",),
            usage_policy=(
                "Retrosynthesis in this workspace uses single-step calls; the agent owns "
                "multi-step planning. Do not invoke the built-in multistep planner, "
                "including through custom Python scripts."
            ),
        ),
        OperationDefinition("assess_route_step", contract_version="2", required_artifacts=("retro_library",),
                            evidence_arguments=("evidence_refs",), result_evidence_field="evidence_refs"),
        OperationDefinition("assess_retro_validity", contract_version="3", required_artifacts=("retro_library",),
                            evidence_arguments=("source_ref", "forward_ref", "evidence_refs"),
                            usage_policy="Assess concrete realizations; ordinal evidence is not success probability. "
                            "Use a saved bounded forward check only for consequential uncertainties."),
        OperationDefinition(
            "disconnect_composite", required_artifacts=("composite_library", "composite_catalog"),
            usage_policy=(
                "Return bounded composite actions, each retaining two physical reactions. "
                "The agent owns route expansion and stopping decisions. Count two physical "
                "steps; inspect both reactions and conditions before selecting a route."
            ),
        ),
        OperationDefinition("assess_route_proposal", contract_version="2", required_artifacts=("retro_library",),
                            evidence_arguments=("evidence_refs",), result_evidence_field="evidence_refs"),
        OperationDefinition("assess_route_step_forward", required_artifacts=("forward_library",),
                            evidence_arguments=("source_ref",), execution_status_field="execution_status",
                            replay_comparison="scientific_result_and_stage_outcomes_excluding_forward_execution_telemetry",
                            replay_projection=_forward_replay_result),
        OperationDefinition("inspect_route_step", evidence_arguments=("source_ref",)),
        OperationDefinition("inspect_step_precedents", contract_version="2", evidence_arguments=("source_ref",)),
        OperationDefinition("revise_route_branch", contract_version="2", required_artifacts=("retro_library",),
                            evidence_arguments=("source_ref", "evidence_refs"), result_evidence_field="evidence_refs"),
        OperationDefinition("compare_route_proposals", evidence_arguments=("source_refs",)),
    )

    def __init__(self, store: InvestigationStore, *, fragment_worker: FragmentWorker | None = None) -> None:
        names = tuple(item.name for item in self.DEFINITIONS)
        if len(set(names)) != len(names):
            raise ValueError("Scientific operation names must be unique")
        if any(not callable(getattr(self, name, None)) for name in names):
            raise ValueError("Every scientific operation requires a registered implementation")
        self.store = store
        self.fragment_worker = fragment_worker
        self._recommender: Any = None
        self._proposal_library: Any = None

    def definition(self, operation: str) -> OperationDefinition:
        """Resolve a reviewed capability declaration without importing caller code."""
        for definition in self.DEFINITIONS:
            if definition.name == operation:
                return definition
        raise ValueError(f"Unknown scientific operation: {operation}")

    def catalog(self) -> list[dict[str, Any]]:
        """Describe the explicit operations available to a workspace client."""
        return [{"name": item.name, "contract_version": item.contract_version,
                 "signature": str(inspect.signature(getattr(self, item.name))),
                 "description": inspect.getdoc(getattr(self, item.name)) or "",
                 "required_artifacts": list(item.required_artifacts),
                 "evidence_arguments": list(item.evidence_arguments),
                 "usage_policy": item.usage_policy,
                 "execution_status_field": item.execution_status_field,
                 "replay_comparison": item.replay_comparison}
                for item in self.DEFINITIONS]

    def invoke(self, operation: str, arguments: Mapping[str, Any]) -> Any:
        """Dispatch only named application operations, never caller-supplied imports."""
        definition = self.definition(operation)
        return getattr(self, definition.name)(**arguments)

    def prepare_literature_reaction(
        self, source_ref: str, source_id: str, title: str, locator: str, structure_evidence: str,
        reactants: list[dict[str, Any]], products: list[dict[str, Any]],
        structure_origin: Literal["source_explicit", "reconstructed_from_description"] = "reconstructed_from_description",
        indexed_ref: str | None = None,
        reaction_id: str | None = None, conditions: list[dict[str, Any]] | None = None,
        yield_info: dict[str, Any] | None = None, source_conflicts: list[dict[str, Any]] | None = None,
        limitations: list[str] | None = None,
        scheme_refs: list[str] | None = None,
    ) -> dict[str, Any]:
        """Prepare source drawings with exact participant passages and optional indexed graph reuse.

        Each side has 1..20 participants with name, compound_id, evidence_ref and
        material_form. Supply smiles for reconstructions/source-explicit graphs;
        omit smiles with indexed_ref/reaction_id. reported_formula is optional.
        Source discrepancies use description and at least two exact excerpt_refs.
        compound_id permits a literal name up to 300 characters. scheme_refs may
        cite up to five w.capture_source_image records from the same source.
        Read w.prepared_literature_reaction(ref) into the answer without dumping it.
        """
        from .literature_reactions import prepare_literature_reaction

        return prepare_literature_reaction(
            self.store, source_ref=source_ref, source_id=source_id, title=title, locator=locator,
            structure_evidence=structure_evidence, reactants=reactants, products=products,
            structure_origin=structure_origin, indexed_ref=indexed_ref, reaction_id=reaction_id,
            conditions=conditions, yield_info=yield_info, source_conflicts=source_conflicts, limitations=limitations,
            scheme_refs=scheme_refs,
        )

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

    def assess_starting_material(
        self, smiles: str, mw_threshold: float | None = None,
        allow_registry_stop: bool = True, allow_literature_stop: bool = True,
        allow_mw_stop: bool = True, unavailable_starting_materials: list[str] | None = None,
    ) -> Any:
        """Assess a route leaf through registry, exact product lookup, then MW.

        Default MW cutoff is 200 g/mol (strictly below), from the versioned policy.
        fragment_index is optional; unavailable/error checks remain explicit.
        Stops are planning assumptions, never verified commercial availability.
        Pass explicit unavailable materials; preserve leaf assumptions in answers.
        """
        from condition_recommender.starting_materials import assess_starting_material

        entry = self.store.manifest["baseline"]["artifacts"].get("fragment_index")
        return assess_starting_material(
            smiles, fragment_index=entry["path"] if entry else None,
            mw_threshold=mw_threshold, allow_registry_stop=allow_registry_stop,
            allow_literature_stop=allow_literature_stop, allow_mw_stop=allow_mw_stop,
            unavailable_starting_materials=(
                [] if unavailable_starting_materials is None else unavailable_starting_materials
            ),
        )

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

    def generate_weak_label_screening_array(
        self, reaction_smiles: str, array_size: int = 24,
        source_reaction_type_hint: str | None = None,
    ) -> Any:
        """Select diverse intact weak-label recipes for screening, with unverified-source warnings.

        Requires baseline-pinned weak_label_records and weak_label_recipe_catalog,
        but no structural condition index. The query must have verified graph edits.
        Returns up to array_size recipes, not predicted yields or invented mixtures.
        An optional source type may narrow retrieval only when consistent with the graph.
        """
        from condition_recommender import (
            generate_weak_label_screening_array,
            load_weak_label_retrieval_rules,
            weak_label_recipe_catalog_path,
        )

        limit = int(load_weak_label_retrieval_rules()["screening_candidate_limit"])
        if type(array_size) is not int or not 1 <= array_size <= limit:
            raise ValueError(f"array_size must be an integer between 1 and {limit}")
        records = self._path("weak_label_records")
        catalog = self._path("weak_label_recipe_catalog")
        if catalog != weak_label_recipe_catalog_path(records).resolve():
            raise ValueError("weak_label_recipe_catalog must be the catalog beside weak_label_records")
        return generate_weak_label_screening_array(
            reaction_smiles, records_path=records, array_size=array_size,
            source_reaction_type_hint=source_reaction_type_hint,
        )

    def get_precedents(
        self, reaction_ids: list[str], offset: int = 0, limit: int = 20,
        view: Literal["summary", "full"] = "summary",
    ) -> dict[str, Any]:
        """Page precedent summaries; request view='full' for complete indexed chemistry."""
        if view not in {"summary", "full"}:
            raise ValueError("view must be summary or full")
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
        records = [asdict(row) for row in index.select(positions[offset:offset + limit])]
        if view == "summary":
            keys = ("reaction_id", "observation_id", "canonical_reaction_id", "reaction_smiles",
                    "yield_pct", "source_dataset", "reference_id", "publication_year", "recipe_id",
                    "resolved_recipe", "condition_uncertain", "chemistry_status", "condition_status",
                    "condition_stage_status", "outcome_status", "precedent_tier", "reaction_label")
            records = [{**{key: row[key] for key in keys if key in row},
                        "chemistry_warnings": (row.get("reaction_core") or {}).get("warnings", []),
                        "details": {"operation": "get_observations", "observation_ids": [row["observation_id"]]},
                        "omitted_sections": ["signature", "reaction_core", "molecular_features", "fallback_descriptor"]}
                       for row in records]
        return {
            "records": records,
            "missing_reaction_ids": [identity for identity in reaction_ids if not index.reaction_ids.get(identity)],
            "total": len(positions), "offset": offset,
            "next_offset": offset + limit if offset + limit < len(positions) else None,
            "record_scope": "indexed_summary" if view == "summary" else "indexed_fields_only",
            "index_scope": index.precedent_scope.value,
        }

    def get_observations(
        self, observation_ids: list[str], fields: list[str] | None = None,
        offset: int = 0, limit: int = 10, max_bytes: int = 524288,
    ) -> dict[str, Any]:
        """Fetch canonical evidence by exact observation ID, with field selection and a byte budget.

        Default fields expose structures, conditions, outcome and admission status.
        Request source, signature or reaction_core explicitly for detailed evidence.
        Oversized records retain identity and list omitted fields for a narrower query.
        """
        from condition_recommender.processed_catalog import ProcessedCatalog
        if not isinstance(observation_ids, list) or len(observation_ids) > 100 or not all(isinstance(i, str) for i in observation_ids):
            raise ValueError("Supply at most 100 observation IDs")
        if offset < 0 or not 1 <= limit <= 20 or not 1024 <= max_bytes <= 4 * 1024 * 1024:
            raise ValueError("Invalid observation page or byte budget")
        if fields is not None and (len(fields) > 100 or not all(isinstance(i, str) for i in fields)):
            raise ValueError("fields must contain at most 100 field names")
        selected = fields if fields is not None else ["reaction_smiles", "source_dataset", "reference_id",
                    "resolved_recipe", "yield_pct", "chemistry_status", "condition_status", "admission_tier",
                    "admission_reasons", "warnings"]
        ids = list(dict.fromkeys(observation_ids))
        catalog = ProcessedCatalog(self._path("evidence_catalog"))
        records, missing, consumed, size = [], [], 0, 0
        for identity in ids[offset:offset + limit]:
            try:
                row = catalog.observation(identity, selected)
            except KeyError:
                missing.append(identity)
                consumed += 1
                continue
            encoded_size = len(json.dumps(row, ensure_ascii=False).encode("utf-8"))
            if encoded_size > max_bytes:
                row = {"observation_id": identity, "reaction_id": row.get("reaction_id"),
                       "omitted_fields": sorted(set(row) - {"observation_id", "reaction_id"}),
                       "status": "record_exceeds_budget_select_fewer_fields", "requested_bytes": encoded_size}
                encoded_size = len(json.dumps(row).encode())
            if records and size + encoded_size > max_bytes:
                break
            records.append(row)
            consumed += 1
            size += encoded_size
        next_offset = offset + consumed
        return {"records": records, "missing_observation_ids": missing, "total": len(ids),
                "offset": offset, "next_offset": next_offset if next_offset < len(ids) else None,
                "returned_bytes": size, "selected_fields": selected, "record_scope": "canonical_evidence"}

    def get_routes(self, route_ids: list[str], offset: int = 0, limit: int = 10,
                   include_tree: bool = False, include_source: bool = False) -> dict[str, Any]:
        """Page route summaries and step-observation joins; fetch source or typed tree explicitly."""
        from core_retrosynthesis.processed_routes import ProcessedRouteCatalog
        if not isinstance(route_ids, list) or len(route_ids) > 100 or not all(isinstance(i, str) for i in route_ids):
            raise ValueError("Supply at most 100 route IDs")
        if offset < 0 or not 1 <= limit <= 20:
            raise ValueError("Invalid route page")
        catalog = ProcessedRouteCatalog(self._path("route_catalog"))
        ids = list(dict.fromkeys(route_ids))
        records, missing = [], []
        for identity in ids[offset:offset + limit]:
            try:
                records.append(catalog.route(identity, include_tree=include_tree, include_source=include_source))
            except KeyError:
                missing.append(identity)
        return {"records": records, "missing_route_ids": missing, "total": len(ids), "offset": offset,
                "next_offset": offset + limit if offset + limit < len(ids) else None}

    def find_synthesis_precedents(self, target_smiles: str, limit: int = 10,
                                 timeout_seconds: int = 90) -> dict[str, Any]:
        """Automatically search cores and add context; inspect source hits before proposing a route."""
        from .fragment_search import run_fragment_search

        return run_fragment_search(self, {"target_smiles": target_smiles, "limit": limit,
                                         "timeout_seconds": timeout_seconds}, automatic=True)

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

    def propose_fragment_queries(
        self, target_smiles: str, query: str, query_format: str = "smiles",
        topology: str = "preserve_rings", aromatic_atom_ids: list[int] | None = None,
    ) -> dict[str, Any]:
        """Preview target-validated peripheral, ring-boundary and selected C/N edits.

        C/N positions use query atom IDs from this operation. Choose a returned
        variant for search_fragment_precedents; no search is executed here.
        """
        from reactive_taxonomy.fragment_broadening import propose_fragment_queries

        return propose_fragment_queries(target_smiles, query, query_format, topology, aromatic_atom_ids)

    def investigate_fragment_precedent(self, source_ref: str, observation_id: str) -> dict[str, Any]:
        """Inspect and test a selected observation from a saved target-derived search.

        Source conditions, procedures and comparisons survive compilation failure.
        Only admitted source operators are tested, in memory; no production library
        or global search is required. Partial searches remain inspectable.
        """
        from core_retrosynthesis.fragment_investigation import investigate_fragment_precedent
        from .step_selection import _call

        payload = _call(self.store, source_ref, {"search_fragment_precedents"})
        search = payload["result"]
        target = (search.get("target_validation") or {}).get("target_smiles")
        if not target:
            raise ValueError("Choose a saved fragment search with target_smiles")
        return {**investigate_fragment_precedent(search, observation_id, target), "source_ref": source_ref}

    def search_fragment_precedents(
        self, query: str, query_format: str = "smiles", topology: str = "preserve_rings",
        limit: int = 5, timeout_seconds: int = 30, target_smiles: str | None = None,
        search_side: str = "product",
        query_variant_ref: str | None = None, query_variant_id: str | None = None,
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

        if query_variant_ref is not None or query_variant_id is not None:
            from reactive_taxonomy.fragment_search import indexed_product
            from .step_selection import _call

            if not query_variant_ref or not query_variant_id or not target_smiles:
                raise ValueError("Supply query_variant_ref, query_variant_id and target_smiles together")
            preview = _call(self.store, query_variant_ref, {"propose_fragment_queries"})["result"]
            selected = [item for item in preview["variants"] if item["variant_id"] == query_variant_id]
            if (len(selected) != 1 or preview["target_smiles"] != indexed_product(target_smiles)[0]
                    or any(selected[0][key] != value for key, value in (
                        ("query", query), ("query_format", query_format), ("topology", topology)))):
                raise ValueError("Query request differs from the selected saved variant")
        return run_fragment_search(self, {"query": query, "query_format": query_format,
                                         "topology": topology, "limit": limit,
                                         "timeout_seconds": timeout_seconds,
                                         "target_smiles": target_smiles, "search_side": search_side})

    def get_procedures(self, reaction_ids: list[str], offset: int = 0, limit: int = 100) -> dict[str, Any]:
        """Read all matching observed procedure records; missing text remains missing."""
        if not isinstance(reaction_ids, list) or not all(isinstance(item, str) for item in reaction_ids):
            raise ValueError("reaction_ids must be a list of strings")
        if len(reaction_ids) > 100:
            raise ValueError("At most 100 reaction IDs are supported")
        selected = set(reaction_ids)
        path = self._path("procedure_catalog")
        if path.suffix == ".sqlite":
            from condition_recommender.processed_catalog import ProcessedCatalog
            result = ProcessedCatalog(path).procedures(reaction_ids, offset=offset, limit=limit)
            return {**result, "source_path": str(path), "origin": "source_report"}
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
        from condition_registry import (
            ConditionComponentInput,
            build_resolved_recipe_from_inputs,
        )

        return build_resolved_recipe_from_inputs(
            [ConditionComponentInput(**item) for item in components],
            temperature_c=temperature_c, time_h=time_h, concentration_m=concentration_m,
            atmosphere=atmosphere,
        )

    def assess_recipe(self, reaction_smiles: str, recipe: dict[str, Any]) -> Any:
        """Assess compatibility; unknown/invalid_input are not chemical conflicts or success."""
        from condition_recommender import assess_reaction_recipe

        return assess_reaction_recipe(reaction_smiles, recipe)

    def assess_proposed_recipe(
        self, reaction_smiles: str, components: list[dict[str, Any]],
        operating_conditions: dict[str, Any], evidence_refs: list[str] | None = None,
    ) -> dict[str, Any]:
        """Normalize and check the actual proposed recipe, with explicit process stages.

        Components use raw_identifier, source_field, identifier_type='name', optional
        source_role_hint, amount, amount_unit and provenance. Operating fields are
        temperature_c, time_h, concentration_m, atmosphere, stages and declared_absences.
        Stages use stage_index, temperature_c, time_h, atmosphere and provenance.
        Cite evidence_refs; an unknown check is not evidence of incompatibility.
        """
        from .investigation_checks import assess_proposed_recipe

        return assess_proposed_recipe(self, reaction_smiles, components, operating_conditions, evidence_refs)

    def search_captured_sources(
        self, queries: list[str], source_refs: list[str] | None = None, offset: int = 0, limit: int = 5,
    ) -> dict[str, Any]:
        """Search 1..10 literal names/labels in captured roots; page matches using next_offset.

        By default search all recorded source text, skipping failed fetches. Matches
        retain exact offsets for record_source_excerpt, not verified chemical claims.
        """
        from .investigation_checks import search_captured_sources

        return search_captured_sources(self, queries, source_refs, offset, limit)

    def inspect_route_inputs(
        self, source_ref: str, leaf_queries: list[dict[str, Any]] | None = None, source_refs: list[str] | None = None,
    ) -> dict[str, Any]:
        """Inspect all route leaves and search saved sources before deciding to stop.

        Use source_ref from assess_route_proposal/revise_route_branch. leaf_queries
        uses {smiles: actual_leaf, terms: [literal_source_name_or_label]}. Omit
        leaf_queries to inspect all leaves without text search. Missing or empty
        terms remain unsearched; stock assumptions are retained, not promoted.
        Follow source-search next_offset using search_captured_sources for details.
        """
        from .investigation_checks import inspect_route_inputs

        return inspect_route_inputs(self, source_ref, leaf_queries, source_refs)

    def inspect_condition_precedents(
        self, reaction_smiles: str, reaction_ids: list[str], offset: int = 0, limit: int = 20,
    ) -> dict[str, Any]:
        """Compare selected indexed observations and link source procedures without merging them."""
        from condition_recommender import compare_condition_evidence
        from condition_recommender.generic_indexing import GenericIndexedReaction

        from .source_catalogs import reference_records

        page = self.get_precedents(reaction_ids, offset=offset, limit=limit, view="full")
        result = asdict(compare_condition_evidence(
            reaction_smiles, [GenericIndexedReaction(**row) for row in page["records"]],
        ))
        references, reference_status = reference_records(
            self, {row["reference_id"] for row in page["records"] if row.get("reference_id")},
        )
        try:
            procedures = self.get_procedures([row["reaction_id"] for row in page["records"]])
        except FileNotFoundError:
            procedures = {"records": [], "availability": "catalog_unavailable"}
        else:
            procedures["availability"] = "catalog_available"
        for item in result["precedents"]:
            observation = item["observation"]
            item["reference_record"] = references.get(observation.get("reference_id"))
            candidates = [record for record in procedures["records"]
                          if record["reaction_id"] == observation["reaction_id"]]
            item["procedure_observations"] = [record for record in candidates
                                               if record.get("observation_id") == observation["observation_id"]]
            item["reaction_level_procedures"] = [record for record in candidates if not record.get("observation_id")]
            item["procedure_link_scope"] = "exact_observation_id_or_explicitly_unassigned_reaction_record"
        result["procedure_catalog"] = procedures
        result["reference_catalog_status"] = reference_status
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
        required_disconnection_bond: list[int] | None = None,
        focus_target_smiles: str | None = None,
    ) -> Any:
        """Generate single-step strategies for one target; never expand precursor branches.

        The agent chooses realizations, subsequent targets and stopping decisions.
        Requires retro_library only unless conditions are requested. Preserves the
        existing engine's evidence and warnings; internal LLM review is disabled.
        Optional required_disconnection_bond selects two zero-based canonical
        atom IDs. Supply focus_target_smiles from molecular inspection to bind
        the selection to that target. Only verified forward formation of that
        bond qualifies; bond-order changes alone do not. Search remains bounded.
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
            required_disconnection_bond=required_disconnection_bond,
            focus_target_smiles=focus_target_smiles,
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

    def disconnect_composite(
        self, target_smiles: str, top_k: int = 5, max_templates_to_apply: int = 50,
        max_candidates_to_validate: int = 12, use_context: bool = True,
        include_l0: bool = True, include_one_step_fallbacks: bool = True,
    ) -> dict[str, Any]:
        """Propose common coupled transformations as one action with two physical steps.

        Uses separately pinned composite_library and composite_catalog artifacts.
        Target-specific site coupling and both graph transformations are checked.
        Conditions, one-pot execution and material supply remain unassessed; these
        are predicted proposals with source precedents, not observed target routes.
        """
        from core_retrosynthesis import (
            load_composite_strategy_catalog, load_generic_library, search_composite_actions,
        )

        catalog = load_composite_strategy_catalog(self._path("composite_catalog"))
        result = search_composite_actions(
            target_smiles, load_generic_library(self._path("composite_library")), catalog.strategies,
            top_k=top_k, max_templates_to_apply=max_templates_to_apply,
            max_candidates_to_validate=max_candidates_to_validate, use_context=use_context,
            include_l0=include_l0, include_one_step_fallbacks=include_one_step_fallbacks,
        )
        return {**result.to_dict(), "catalog_id": catalog.catalog_id}

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
        Instead of copying a mapping, saved_candidate={source_ref, realization_id,
        strategy_id} reuses that exact disconnection's reconstruction and revalidates
        it. Keep target_smiles/precursor_smiles; they must match the saved candidate.
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
        Each step may include saved_candidate={source_ref, realization_id, strategy_id}
        to retain an exact disconnection's mapping and provenance during reassessment.
        """
        from .route_investigation import assess_route

        return assess_route(self, proposal, unavailable_starting_materials, include_conditions, include_forward, evidence_refs)

    def assess_retro_validity(
        self, proposal: dict[str, Any] | None = None, source_ref: str | None = None,
        step_id: str | None = None, realization_id: str | None = None,
        forward_ref: str | None = None, candidate_limit: int = 128, match_limit: int = 20,
        evidence_refs: list[str] | None = None,
        strategy_id: str | None = None, precursor_smiles: str | None = None,
    ) -> dict[str, Any]:
        """Grade concrete precursor-to-target evidence and suggest an advisory next action.

        Supply proposal or a saved source_ref: use realization_id for a
        disconnect_target result, step_id for a route, neither for assess_route_step.
        Optional pinned condition/shared-core indexes add whole-reaction and L0/L1/L2
        support. Optional forward_ref joins an existing bounded forward check; no
        prediction is run here. Inspect validity.status, cautions and unresolved_checks.
        Ranks 4..0 mean evidence strength, never experimental success probability.
        Template realization IDs can repeat. Supply strategy_id and precursor_smiles
        from the selected saved candidate; ambiguous IDs are rejected.
        """
        from .retro_validity import assess_validity

        return assess_validity(
            self, proposal, source_ref, step_id, realization_id, forward_ref,
            candidate_limit, match_limit, evidence_refs,
            strategy_id, precursor_smiles,
        )

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
        strategy_id: str | None = None, precursor_smiles: str | None = None,
    ) -> dict[str, Any]:
        """Inspect actual supporting reactions for one saved realization or assessed step.

        Use realization_id for disconnect_target, step_id for a route assessment
        or revision, and neither for assess_route_step. Follow page.next_offset.
        offset must be nonnegative; limit is 1..5 (default 3).
        Template realization IDs can repeat. Supply strategy_id and precursor_smiles
        from the selected saved candidate; ambiguous IDs are rejected.
        Source reactions, product comparison, scoped counts, and available source
        conditions remain evidence, not proof that the proposed step will work.
        """
        from .step_precedents import inspect_step_precedents

        return inspect_step_precedents(self, source_ref, step_id, realization_id, offset, limit,
                                       strategy_id, precursor_smiles)

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

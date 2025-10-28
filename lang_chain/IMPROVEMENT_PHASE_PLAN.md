# ChemTools AI Assistant Enhancement Roadmap

## Phase 1 – Structured Tool Interfaces

- Replace JSON-string returns with `StructuredTool` + Pydantic response models for key tools (`normalize_smiles_tool`, `recommend_conditions_tool`, `find_reagent_tool`, `add_reagent_tool`, etc.).
- Surface standardized error payloads so the agent can branch reliably on failures.
- Update `ChemToolsAgent` to lean on tool schemas (function-calling metadata, LLM provider configs) and adjust CLI/tool listings to display new signatures.
- Extend `tests/test_llm.py` (or add a LangChain-specific pytest module) to cover schema validation and backward compatibility with direct `.invoke()` calls.

**Exit criteria:** LLM providers can parse tool outputs without JSON decoding; unit tests assert schema compliance; CLI renders updated tool metadata.

**Status:** ✅ Completed — structured tool schemas, standardized responses, and CLI/test updates are in place.

## Phase 2 – Recommendation Result Reuse

- Introduce an internal cache/dependency injection in `chemtools_wrapper.py` so `recommend_conditions_tool` and `list_supported_cores_tool` share a normalized reaction fingerprint + constraint spec key.
- Return the raw precedent pack alongside simplified data (behind a flag) to avoid redundant DRFP retrieval across follow-up agent turns.
- Add cache invalidation hooks for constraint changes and explicit CLI command to clear session cache.
- Benchmark representative reactions to confirm latency reduction and ensure memory does not balloon in conversational use.

**Exit criteria:** Sequential tool calls reuse cached packs; measured latency drop (>30% on repeat queries); CLI exposes cache management command.

**Status:** 🟢 Cache sharing shipped (`chemtools_wrapper` now reuses DRFP results) and CLI exposes cache stats/clear commands; next capture latency benchmarks to quantify gains.

## Phase 3 – Rich Precedent Rationales

- Wrap `chemtools.explain.for_pack` as a dedicated LangChain tool and/or integrate its output into `recommend_conditions_tool` summaries.
- Provide structured fields (`core_summary`, `base_summary`, `temp_range`) so the agent can cite precedent rationale without extra prompting.
- Update system prompt and examples to encourage rationale-aware responses; expand direct-tool tests to verify narrative content.

**Exit criteria:** Recommendations include chemistry-native rationales sourced from deterministic explainers; new tests validate presence and formatting.

## Phase 4 – Tool Surface Expansion

- Audit ChemTools modules (e.g., `condition_core.py`, `dataset_analytics.py`, taxonomy utilities) and prioritize deterministic helpers for exposure as tools.
- Design typed wrappers with concise docstrings aligned to agent goals (substrate fingerprint breakdowns, ligand coverage stats, success-rate analytics).
- Refresh README/QUICKSTART with new capabilities and provide example prompts demonstrating advanced analytical workflows.

**Exit criteria:** At least three new high-value tools are available with documentation updates and usage demos.

## Phase 5 – Provider & Testing Flexibility

- Abstract `get_llm_client` to support provider plugins (Azure OpenAI, Anthropic, local mocks) and enable dependency injection for offline tests.
- Replace ad-hoc `lang_chain/test_tools.py` script with pytest suites leveraging fixtures/mocked LLM responses; integrate into CI.
- Document configuration matrix (env vars, .env templates) and add smoke-test targets for no-network environments.

**Exit criteria:** Tests run without live API keys; new providers can be registered via configuration only; CI verifies LangChain layer deterministically.

---

**Next Steps:** Validate scope with stakeholders, align resourcing per phase, and open tracking issues/epics before implementation.

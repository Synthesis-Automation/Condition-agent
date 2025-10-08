# Condition MCP — Implementation Plan

_Thin MCP façade over your Condition-Agent stack (tools, resources, prompts), with tests, CI, and deploy._

## 1) Goals & Scope

**Goal:** Make your reaction-condition recommender callable, testable, and UI‑agnostic by exposing it through a small **Model Context Protocol (MCP)** server (“**Condition MCP**”).  
**Scope (deliverables):**

- **Tools** (callable functions): normalize → detect_family → featurize → recommend_core → rank_tails (bases/ligands/solvents) → assemble_conditions → apply_constraints → explain_choice.
- **Resources** (read-only): `Condition_Rule_Library` (Buchwald C–N, Suzuki, Ullmann), small benchmark datasets, I/O schemas.
- **Prompts** (workflow recipes): `playbook.buchwald_cn`, `playbook.suzuki`, `playbook.ullmann`.
- **Client integration**: agent/CLI/VSCode can auto‑discover & call.
- **Quality guardrails**: golden tests, determinism tests, schema contract tests.
- **Packaging & deploy**: containerized server; single endpoint callable by any MCP‑aware client.

---

## 2) Architecture Overview

```
[Agent / UI (Chat, CLI, VSCode)]
        │
        │  (MCP transport: stdio or socket)
        ▼
  ┌─────────────────────────────┐
  │        Condition MCP        │
  │  - tools (function calls)   │
  │  - resources (read-only)    │
  │  - prompts (templates)      │
  └─────────────────────────────┘
        │               ▲
        │ calls         │ reads
        ▼               │
  ┌───────────────────────────────┐
  │  Your Code & Data (“Core”)    │
  │  chemtools, Condition_Rule_…  │
  │  KNN/DRFP models, taxonomies  │
  └───────────────────────────────┘
```

- The MCP server is a **thin adapter**. It doesn’t change your internals; it just defines **stable JSON contracts**.

### Auto‑discovery (how clients find things)

- On connect, clients enumerate what the server offers by sending **list messages** for _tools_, _resources_, and _prompts_.
- The server returns **manifests** (names, parameter schemas, descriptions, optional metadata). Clients render these as callable “capabilities” (like slash‑commands or tool cards).
- When users or agents choose one, the client sends a **call** for a tool, a **read** for a resource, or a **get/instantiate** for a prompt, all with JSON arguments validated by the server.
- This means **no hard‑coded glue** in the client; new tools/resources/prompts become instantly visible after you ship them.

---

## 3) Tooling Contracts (first‑class APIs)

Each tool must have **clear, versioned JSON schemas** for inputs/outputs.

| Tool                   | Purpose                                              | Minimal Input                                             | Output                                                  |
| ---------------------- | ---------------------------------------------------- | --------------------------------------------------------- | ------------------------------------------------------- |
| `normalize_reaction`   | Split/canonicalize reaction SMILES, optional mapping | `smiles_rxn: str`                                         | `{reactants: str[], products: str[], mapped_rxn?: str}` |
| `detect_family`        | Route to Buchwald‑C–N, Suzuki, Ullmann, etc.         | `reactants: str[]`                                        | `{family: str, confidence: number}`                     |
| `featurize_substrates` | Steric/electronic/flags                              | `reactants: str[]`                                        | `{descriptors: object}`                                 |
| `recommend_core`       | Pick “ConditionCore”                                 | `family: str, descriptors: object, k?: int`               | `{core_candidates: [{core_id, score, rationale?}]}`     |
| `rank_bases`           | Tail ranking (bases)                                 | `core_id: str, descriptors: object, constraints?: object` | `[{base_id, score, why?}]`                              |
| `rank_ligands`         | Tail ranking (ligands)                               | same                                                      | `[{ligand_id, score, why?}]`                            |
| `rank_solvents`        | Tail ranking (solvents)                              | same                                                      | `[{solvent_id, score, why?}]`                           |
| `assemble_conditions`  | Compose core + tails                                 | `{core_id, base_id?, ligand_id?, solvent_id?}`            | `{ConditionSet}`                                        |
| `apply_constraints`    | Filter by flags/rules                                | `{ConditionSet, flags}`                                   | `{ConditionSet, violations: []}`                        |
| `explain_choice`       | Human‑readable rationale + rule cites                | `{ConditionSet, inputs}`                                  | `{narrative, citations: []}`                            |

**Schema versioning:** Include `schema_version` in every output; maintain a changelog.

---

## 4) Resources (read-only registries)

- `cond://rules/buchwald_cn.json` — Condition_Rule_Library entry (if‑then heuristics, precedence, known pitfalls).
- `cond://rules/suzuki.json`, `cond://rules/ullmann.json` — ditto.
- `dataset://benchmarks/buchwald_cn.jsonl` — 50 curated reactions with ground‑truth cores/tails for smoke tests.
- `doc://schemas/condition_set.json` — machine‑readable schema for ConditionSet.
- `tax://ligands.json`, `tax://bases.json`, `tax://solvents.json` — your taxonomies (with CAS, synonyms, embeddings optional).

**Metadata to add:** `lastModified`, `source`, `license`, `checksum`, `tags` (e.g., `role=base`).

---

## 5) Prompts (workflow templates)

- `playbook.buchwald_cn(aryl_halide, amine, constraints?)`  
  Guides an agent to call tools in sequence: detect → featurize → recommend_core → rank tails → assemble → apply_constraints → explain.
- `playbook.suzuki(aryl_halide, boronate, constraints?)`
- `playbook.ullmann(aryl_halide, amine/phenol, constraints?)`

Keep them **short, deterministic**, and reference the **tool names** rather than verbose instructions.

---

## 6) Repository Layout

```
condition-mcp/
├─ server/                  # MCP façade
│  ├─ main.py               # entrypoint
│  ├─ tools/                # thin wrappers calling your core
│  │  ├─ normalize.py
│  │  ├─ detect_family.py
│  │  └─ ...
│  ├─ resources/
│  │  ├─ rules/...
│  │  ├─ datasets/benchmarks/...
│  │  └─ schemas/condition_set.json
│  ├─ prompts/
│  │  ├─ playbook_buchwald_cn.txt
│  │  └─ ...
│  └─ mcp_config/           # server manifest/config if needed
├─ core/                    # your existing code (submodule or path)
├─ tests/
│  ├─ golden/benchmarks/... # truth sets for top‑k recovery
│  ├─ contract/             # schema validation tests
│  └─ smoke/                # end‑to‑end prompt-to-tools
├─ docker/
│  └─ Dockerfile
├─ pyproject.toml / uv.lock
└─ README.md
```

---

## 7) Phased Work Plan (WBS)

### Phase A — Foundation

- [ ] Pin environment (Python, RDKit, chemtools); write `uv.lock`/`requirements.txt`.
- [x] Define **ConditionSet** JSON schema (IDs, loadings, units, rationale, citations). ✅ See `condition_mcp/resources/schemas/condition_set.json`

### Phase B — Tools (thin wrappers)

- [x] Implement `normalize_reaction`, `detect_family`, `featurize_substrates` wrappers. ✅ Delivered via `condition_mcp/tools/` and exposed in the Gradio MCP tab
- [ ] Implement `recommend_core`, `rank_bases`, `rank_ligands`, `rank_solvents` (call your current models/rules).
- [ ] Implement `assemble_conditions`, `apply_constraints`, `explain_choice`.
- [x] Add runtime **input validation** + output `schema_version`. ✅ Pydantic-backed request/response models wrap each tool

### Phase C — Resources

- [ ] Import `Condition_Rule_Library` JSON (Buchwald, Suzuki, Ullmann).
- [ ] Create **benchmark mini‑packs** (≈50/ family) with expected core/tails.
- [x] Publish `doc://schemas/condition_set.json` and taxonomy snapshots. ✅ Schema published; taxonomy snapshots still pending

### Phase D — Prompts

- [ ] Author the three playbooks; keep them short and tool‑oriented.
- [ ] Add parameter docs + usage examples.

### Phase E — Tests & CI

- [x] **Contract tests**: JSON schema validation for every tool output. ✅ `tests/condition_mcp/test_tools_basic.py`
- [ ] **Determinism**: same input → identical ranked list (within tie‑break rules).
- [ ] **Golden**: top‑k recovery on benchmarks; enforce minimal thresholds.
- [ ] **Safety**: forbidden pairs (e.g., solvent/base constraints) always filtered.
- [ ] CI pipeline: lint, type‑check, test, container build, artifact publish.

### Phase F — Packaging & Deployment

- [ ] Build Docker image (non‑root user; read‑only resources; env vars for model paths).
- [ ] Expose MCP over stdio (for local dev) and TCP (for remote)—secure as needed.
- [ ] Publish versioned container tags (e.g., `v0.1.0`).

### Phase G — Client Integration

- [ ] Test with an **MCP Inspector/CLI** to verify auto‑discovery (list tools/resources/prompts).
- [ ] Wire to your chat UI / VS Code extension that speaks MCP.
- [ ] Provide quickstart snippets for calling tools (copy‑paste for teammates).

### Phase H — Observability & Ops

- [ ] Structured logging (tool name, latency, input hash, output hash, status).
- [ ] Stats: success rate, p95 latency, rule‑coverage per family.
- [ ] Error taxonomy + actionable messages for clients.

---

## 8) Data & Rule Strategy

- Start with a **minimal** Condition_Rule_Library; include rule IDs and citations so `explain_choice` can reference them.
- Maintain **taxonomies** (bases, ligands, solvents) with CAS, synonyms, role flags, and optional **embeddings** for future vector search.
- Keep **benchmarks small but curated**; expand after the first pass.

---

## 9) Security & Compliance

- Expose only necessary tools/resources.
- Validate inputs strictly (types, length, allowed ranges).
- Sanitize file/resource paths; resources are **read‑only**.
- Consider rate limiting / concurrency caps; protect paid model keys via server‑side env vars.

---

## 10) Versioning & Compatibility

- Semantic versioning for the server (`MAJOR.MINOR.PATCH`).
- Each tool output includes `schema_version`.
- Deprecate by supporting N‑1 schema versions with clear warnings.

---

## 11) Acceptance Criteria (Definition of Done)

- Clients can **auto‑discover** tools/resources/prompts without hard‑coded names.
- All tools pass **contract tests** and **golden tests** at configured thresholds.
- A user can run `playbook.buchwald_cn` end‑to‑end and receive a valid `ConditionSet` with rationale.
- Container image published and runnable with a one‑liner.
- README includes: quickstart, examples, schemas, and troubleshooting.

---

## 12) Nice‑to‑Have (Post‑v0.1)

- Add **`simulate_htE_grid`** tool to generate DoE grids.
- Add **`rank_additives`** tool using your additive/acid taxonomy.
- Add **confidence calibration** (e.g., isotonic or conformal) for scores.
- Emit **ConditionID** stable identifiers for tracking.
- Optional **embedding search** resource for reagent similarity lookup.

---

## 13) Risk Register & Mitigations

- **Non‑determinism** in ranking → set seeds, stable tie‑breakers, snapshot model weights.
- **Schema drift** → strict contract tests, changelog discipline.
- **Performance regressions** → p95 latency budget + CI perf gate.
- **Data licensing** → store only metadata or licensed subsets in resources.

---

## 14) Quick Usage Examples (client view)

- **List capabilities**: client asks for the server’s tools/resources/prompts; a UI shows them as callable actions.
- **Call a tool**: invoke `recommend_core` with `{family:"Buchwald-CN", descriptors:{...}, k:5}`; receive ranked `core_candidates`.
- **Read a resource**: fetch `cond://rules/buchwald_cn.json` to display rules alongside the recommendation.
- **Run a playbook**: call `playbook.buchwald_cn` with `{aryl_halide:"Cl‑Ar...", amine:"H2N‑Ar...", constraints:{no_chlorinated:true}}` to orchestrate the full chain.

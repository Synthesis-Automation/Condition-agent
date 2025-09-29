# MCP Paper Server Integration Plan (for Conditions-Recommendation Codebase)

**Owner:** Condition-Agent team  
**Target systems:** ChemTools / Condition-Agent orchestrator, precedent retrieval, HTE loop  
**Goal:** Integrate the MCP rules server (Buchwald C–N today; Ullmann C–N next) as a first-class tool the orchestrator can call to **apply**, **merge**, and **audit** rules; leave hooks for a future `extract_from_paper` parser.

---

## 1) Scope & Success Criteria

### In-scope
- Vendor the MCP rules server stub (`server.py`) into the repo and make it runnable from our monorepo tooling.
- Provide a **Python Client Adapter** that the orchestrator can call with a simple function:
  - `mcp_rules.apply(reaction_smiles: str, features: dict) -> Suggestions`
  - `mcp_rules.merge(candidates: dict, strategy="append_as_candidate") -> MergeResult`
  - `mcp_rules.audit(rule_id: str) -> RuleBody`
- Add **config** to select the ruleset JSON (Buchwald vs Ullmann).
- Add **feature-mapping** utilities that convert our substrate analyzer output → the server’s `features` schema.
- Add **tests** (unit + integration) and **CI** checks.
- Provide a **graceful fallback** when the MCP server is not present.

### Definition of Done
- CI green on Linux/macOS/Windows.
- A demo notebook/CLI shows: (1) apply on a few reactions, (2) merge a candidate rule, (3) audit an id.
- Configurable ruleset path via env/CLI/JSON config.
- Telemetry/logs for each MCP call (latency, errors, payload size).

---

## 2) Repository Layout Changes

Add a new top-level package and a client adapter:

```
/agents/
/chemtools/
/condition_agent/
/data/rules/
  ├─ buchwald_cn.json           # existing
  ├─ ullmann_cn.json            # new (from us; candidate status entries)
/mcp_rules_server/
  ├─ server.py                  # provided stub
  ├─ requirements.txt
  ├─ README.md
  └─ .env.example
/clients/
  └─ mcp_rules_client.py        # new: subprocess-based adapter
/tests/
  ├─ test_mcp_rules_client.py   # unit tests with a mocked server
  └─ test_rules_integration.py  # starts a real server on stdio and exercises flows
/examples/
  └─ demo_rules_cli.py          # convenience CLI for developers
```

Notes:
- We keep the MCP server **vendored** to avoid cross-repo drift. Can replace with a submodule later.
- `data/rules/ullmann_cn.json` is optional now; included for symmetry.

---

## 3) Dependencies

Server (vendored):  
```
modelcontextprotocol>=0.1.0
anyio>=3.7
jsonpatch>=1.33
pydantic>=2.5
```

Client/Orchestrator (new): no extra deps required; we’ll use `subprocess` + stdio JSON-RPC to talk MCP. Optionally add a lightweight client helper later.

---

## 4) Configuration

Priority order: **CLI flag > ENV > default**.
- `RULES_PATH` — absolute path to rules JSON, defaults to `data/rules/buchwald_cn.json`.
- `MCP_SERVER_BIN` — command to run the server, default `python mcp_rules_server/server.py`.
- `MCP_STARTUP_TIMEOUT_S` — default 10s.
- `MCP_LOG_LEVEL` — info|debug.

Example `.env`:
```env
RULES_PATH=./data/rules/buchwald_cn.json
MCP_SERVER_BIN=python ./mcp_rules_server/server.py
MCP_STARTUP_TIMEOUT_S=10
MCP_LOG_LEVEL=info
```

---

## 5) Feature Mapping (Condition-Agent → rules.apply)

We standardize the features dict the server expects (minimal matcher, nested keys allowed).

**Expected keys (today):**
```json
{
  "electrophile": { "class": "aryl chloride|aryl bromide|aryl iodide|heteroaryl bromide", "subclass": "optional e.g., sterically hindered pyrazine" },
  "nucleophile": { "class": "primary aniline|secondary aliphatic amine|amide|sulfonamide|ammonia" },
  "mode": "batch|flow|hto",
  "temperature_C":  optional number,
  "solvent":        optional string,
  "base":           optional string
}
```

**Implement mappers:**
- `map_electrophile(reactants) -> {"class": "...", "subclass": "..."}` using existing substrate typing (RDKit + our heuristics).
- `map_nucleophile(reactants) -> {"class": "..."}` using amine classifier (primary/secondary; aniline vs aliphatic; amide/sulfonamide).
- `infer_context(options) -> {"mode": "...", ...}` from user constraints/HTE settings.
- Compose into a single `features` dict for `rules.apply`.

Add robust unit tests for these mappers with hand-picked edge cases (ortho-blocked anilines, electron-poor aryls, heteroaryls).

---

## 6) Client Adapter (Subprocess stdio)

Create `clients/mcp_rules_client.py`:

```python
import os, json, subprocess, threading, queue, time, sys, shlex, atexit

class McpRulesClient:
    def __init__(self, server_cmd: str, rules_path: str, startup_timeout_s: float = 10.0):
        self.server_cmd = server_cmd
        self.rules_path = rules_path
        self.startup_timeout_s = startup_timeout_s
        self.p = None
        self._rq = queue.Queue()
        self._wq = queue.Queue()
        self._id = 0

    def start(self):
        env = os.environ.copy()
        cmd = f"{self.server_cmd} --rules {shlex.quote(self.rules_path)}"
        self.p = subprocess.Popen(shlex.split(cmd), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1)
        atexit.register(self.stop)

    def stop(self):
        if self.p and self.p.poll() is None:
            try:
                self.p.terminate()
                self.p.wait(timeout=3)
            except Exception:
                self.p.kill()
        self.p = None

    def _rpc(self, method: str, params: dict):
        self._id += 1
        req = {"jsonrpc":"2.0","id":self._id,"method":method,"params":params}
        self.p.stdin.write(json.dumps(req) + "\n")
        self.p.stdin.flush()
        line = self.p.stdout.readline()
        if not line:
            raise RuntimeError("MCP server returned no data.")
        resp = json.loads(line)
        if "error" in resp:
            raise RuntimeError(resp["error"])
        return resp.get("result")

    def apply(self, reaction_smiles: str, features: dict, max_suggestions: int = 5) -> dict:
        return self._rpc("rules.apply", {"reaction_smiles": reaction_smiles, "features": features, "max_suggestions": max_suggestions})

    def merge(self, candidates: dict, strategy: str = "append_as_candidate") -> dict:
        return self._rpc("rules.merge", {"candidates": candidates, "strategy": strategy})

    def audit(self, rule_id: str) -> dict:
        return self._rpc("rules.audit", {"rule_id": rule_id})
```

**Graceful fallback:** if the server fails to start, raise a typed error and allow the orchestrator to switch to a pure-local rule engine or default heuristics.

---

## 7) Orchestrator Wiring

Create a thin service inside `condition_agent/`:

```python
# condition_agent/services/rules_service.py
from clients.mcp_rules_client import McpRulesClient
from condition_agent.features.mapping import build_features  # you implement

class RulesService:
    def __init__(self, cfg):
        self.client = McpRulesClient(cfg.mcp_server_bin, cfg.rules_path, cfg.mcp_startup_timeout_s)
        self.client.start()

    def suggest_conditions(self, reaction_smiles, substrate_objects, context_opts):
        features = build_features(substrate_objects, context_opts)
        return self.client.apply(reaction_smiles, features)

    def merge_candidates(self, candidates, strategy="append_as_candidate"):
        return self.client.merge(candidates, strategy=strategy)

    def audit_rule(self, rule_id: str):
        return self.client.audit(rule_id)
```

Update the agent’s tool registry to expose a `rules.suggest` tool that calls `RulesService.suggest_conditions(...)`.

---

## 8) Example CLI

Add `/examples/demo_rules_cli.py` for developers:

```python
import json, os, sys
from condition_agent.config import load_config
from condition_agent.services.rules_service import RulesService

if __name__ == "__main__":
    cfg = load_config()
    svc = RulesService(cfg)
    rxn = "FC(F)(F)c1ccc(Cl)cc1.Cc1ccc(N)cc1>>Cc1ccc(Nc2ccc(C(F)(F)F)cc2)cc1"
    features = {
        "electrophile": {"class":"aryl chloride"},
        "nucleophile": {"class":"primary aniline"},
        "mode": "batch"
    }
    result = svc.client.apply(rxn, features)
    print(json.dumps(result, indent=2))
```

---

## 9) Tests

### Unit
- `test_mcp_rules_client.py`:
  - mocks server process, verifies JSON-RPC formatting and error handling.
  - golden tests for request/response encoding.

- Feature mapping tests:
  - given SMILES, ensure correct `electrophile/nucleophile` class labels.

### Integration
- `test_rules_integration.py`:
  - spin up real server with `buchwald_cn.json` (fixture).
  - call `apply` on 3–5 canonical cases (ArCl + aniline, ArBr + aniline, heteroaryl + amide).
  - `merge` a tiny candidate, then `audit` it; verify it appears in file & in-memory.

Add a small test dataset under `/tests/data/` with 3 inputs and expected keys (not exact scores).

---

## 10) CI & Dev Scripts

- Update CI to run server integration tests (Linux + Windows).
- Add `make` or `tox` tasks:
```
make rules-test          # runs only rules-related tests
make rules-demo          # starts server & runs demo CLI
```

- Pre-commit hook to ensure `data/rules/*.json` are valid JSON (and optionally sorted keys).

---

## 11) Observability & Logging

- Log each MCP call with: tool name, duration, payload bytes, top rule ids in response.
- Add a sampling logger to dump full responses when `MCP_LOG_LEVEL=debug`.
- Emit counters: `mcp_apply_success_total`, `mcp_apply_failure_total`, `mcp_latency_ms`.

---

## 12) Multi-Ruleset Support (Buchwald + Ullmann)

- Add `--rules` flag/environment for path; ship both JSON files:
  - `data/rules/buchwald_cn.json` (existing)
  - `data/rules/ullmann_cn.json` (candidate; you can promote later)
- Optionally run **two** server instances if you want isolation, or **one** server with a `rules.switch` tool (future). For now: one server per ruleset path is simplest.

---

## 13) Security & Safety

- Server runs **locally** via stdio; no network listener by default.
- Only read/write the configured rules file; create timestamped backups before merge.
- Limit payload size to 2 MB; reject larger `merge` candidates.
- Sanitize any downstream strings that could be used in shell contexts (already safe with stdio, but be defensive).

---

## 14) Rollout Plan

1. **Dev**: integrate client adapter + simple demo; run locally with Buchwald rules.  
2. **QA**: run integration tests, compare agent outputs before/after (no regression on existing workflows).  
3. **Pilot**: enable for a subset of agent runs (e.g., behind feature flag `rules_via_mcp=true`).  
4. **Full**: make MCP the default rules backend; keep fallback for 1–2 releases.  
5. **Next**: hook up Ullmann rules and start merging BO-derived candidates with `rules.merge` + manual review.

---

## 15) Work Items (for Codex)

- [ ] Add folder structure & vendor `/mcp_rules_server`.
- [ ] Implement `clients/mcp_rules_client.py` (subprocess JSON-RPC).
- [ ] Implement `condition_agent/services/rules_service.py`.
- [ ] Implement feature mappers in `condition_agent/features/mapping.py`:
  - [ ] `map_electrophile(...)`
  - [ ] `map_nucleophile(...)`
  - [ ] `build_features(...)`
- [ ] Add config plumbing (`load_config`, env overrides).
- [ ] Write unit tests (client + mapping) and integration tests (server round-trip).
- [ ] Add example CLI in `/examples/demo_rules_cli.py`.
- [ ] Wire into orchestrator tool registry (expose `rules.suggest`).
- [ ] CI scripts + pre-commit.
- [ ] Documentation updates in `docs/rules_mcp.md` (copy sections from this plan).

---

## 16) Acceptance Demo Script

```bash
# 1) Start demo (uses default rules: buchwald)
python ./examples/demo_rules_cli.py

# 2) Merge a candidate rule
python - <<'PY'
from clients.mcp_rules_client import McpRulesClient
cli = McpRulesClient("python ./mcp_rules_server/server.py", "./data/rules/buchwald_cn.json")
cli.start()
cand = {"playbooks":[{"id":"BH-CAND-DEMO-1","status":"candidate","when":{"electrophile":{"class":"aryl chloride"},"nucleophile":{"class":"aniline"}},"recipe":{"ligands":["XPhos"],"base":["K3PO4"],"solvent":["1,4-dioxane"],"temperature_C":[100]}}]}
print(cli.merge(cand))
print(cli.audit("BH-CAND-DEMO-1"))
PY

# 3) Apply on a sample reaction to see suggestions
python - <<'PY'
from clients.mcp_rules_client import McpRulesClient
cli = McpRulesClient("python ./mcp_rules_server/server.py", "./data/rules/buchwald_cn.json")
cli.start()
rxn = "FC(F)(F)c1ccc(Cl)cc1.Cc1ccc(N)cc1>>Cc1ccc(Nc2ccc(C(F)(F)F)cc2)cc1"
features = {"electrophile":{"class":"aryl chloride"},"nucleophile":{"class":"primary aniline"},"mode":"batch"}
print(cli.apply(rxn, features))
PY
```

---

## 17) Future Hooks

- Replace placeholder `rules.extract_from_paper` with a PDF→DSL extractor (Grobid/ScienceParse + prompt templates + section heuristics).
- Add `rules.evaluate` + `rules.promote` gates tied to HTE/precedent back-tests.
- Introduce a simple schema validator for rules JSON (pydantic models).

---

## Appendix: Data Contracts

**Apply → Suggestions (shape)**:
```json
{
  "input": {...},
  "matched_playbooks": ["BH-..."],
  "suggestions": [
    { "conditions": {"from_playbook":"BH-...","ligands":"XPhos","base":"K3PO4","solvent":"1,4-dioxane","temperature_C":100}, "score": 0.7, "notes": ["Guard ..."] }
  ]
}
```

**Merge → Result**:
```json
{ "message":"merge_ok", "rules_path":".../buchwald_cn.json", "version":"2025-09-29T12:34:56Z" }
```

**Audit → Rule body**:
```json
{ "section": "playbooks", "rule": { "...": "..." } }
```

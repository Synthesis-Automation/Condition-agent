# Iterative fragment-search development POC

This harness lets an agent prepare a graph query, search its saved reference,
inspect a returned record, and choose the next query from the evidence. It uses
the existing scientific workspace and canonical fragment matcher. No automatic
query policy, chemistry definitions, admission, ranking, or production APIs change.

Run one request at a time from the repository root:

```powershell
python -m examples.ai_native.iterative_search_poc results/ai_native/my_search request.json --output response.json
```

Initialize a new workspace with an existing index and its evidence catalog:

```json
{
  "action": "init",
  "objective": "Investigate preparation of the selected target core",
  "target_smiles": "TARGET_SMILES",
  "budget_seconds": 90,
  "artifacts": {
    "fragment_index": "PATH_TO_FRAGMENT_INDEX.sqlite",
    "evidence_catalog": "PATH_TO_LINKED_CATALOG.sqlite"
  }
}
```

Prepare a query using `action: prepare`, `query`, `query_format` (smiles or smarts),
`topology` (preserve_rings or subgraph), and `reason`. Every query is checked against
the target before searching. Revisions include `parent_search_ref` and describe
what changed. A timed-out parent remains usable as provenance, never as evidence
of absence. The returned `query_ref` identifies the immutable query preview.

Search with `action: search`, `query_ref`, `reason`, and `timeout_seconds` (1–30).
An optional `variant_id` selects an exact alternative from the preview. The caller
does not recopy its expression. A compact card preserves reaction SMILES, edit
witnesses, unresolved evidence, warnings, and source IDs. Witness omissions are
explicit. Full source data remain in the workspace artifact.

Inspect with `action: inspect`, `search_ref`, `observation_id`, and `reason`.
This opens the saved record, witnesses, and linked procedures; it does not run
transfer, establish usefulness, or authenticate the original publication.

`action: baseline` runs the existing automatic discovery, with the same cumulative
search-operation allowance as the iterative arm. `action: summary` reads the
recorded runs. Preparation, baseline verification, reasoning, and inspection are
outside that search allowance; workspace durations are separately recorded.
The pilot is a single-writer local tool, not a production session service.

## Performance control

`iterative_search_replay.py` is a worker for `ScientificWorkspace.run_python`.
Copy it into the investigation, then supply `index_path` and `search_refs`, with
the same references in `evidence_refs`. It loads one request-owned
`FragmentSearchSession` and replays already selected queries. It compares complete
scientific payloads, excluding execution telemetry, against the original calls.
This is a fixed-query timing control, not another adaptive search trial. The
private session hook is experimental and does not establish a supported public API.

Ordinary single-query workspace calls start a worker and reload the library on
each call. The persistent console below now supports live reuse. Compact cards
reduce agent-facing output but do not defer hydration of
full records in the existing matcher. Query editing still requires SMILES/SMARTS
or existing preview variants; a general natural-language query editor is not built.

## Initial development run

Local artifacts: `results/ai_native/iterative_search_poc_20261008/`.
The two targets came from the already inspected workflow examples. The root agent
selected each next query after reading the previous result; this was not an
independent blinded model trial. No web search or prior local pilot findings were
used. The precursor direction was informed by the known workflow and recorded as
such. Automatic discovery ran first for the aza target and second for the amino
target. Source-linked structures can support a useful lead without establishing
an executable preparation or the complete target synthesis. Counts of graph
construction hits are not counts of successful routes or independent publications.

The local HTML/JSON report records query traces, timed-out calls, source
inspections, unchanged-result replay checks, and limitations. These development
results do not satisfy independent chemistry review or untouched evaluation gates.

## Live console and agent pilot

The follow-up integration is
`python -u -m chem_coworker.scientific_workspace.console WORKSPACE`.
Use one persistent terminal and send one JSON line per request, reading its
response before choosing the next. This is the canonical workspace adapter;
the earlier request-file harness above is retained for its recorded experiment.

```json
{"operation":"search_fragment_precedents","arguments":{"query":"AGENT_CHOSEN_CORE","target_smiles":"TARGET","timeout_seconds":30},"reason":"Chemical question being tested"}
{"action":"inspect","source_ref":"SEARCH_REF","path":["result","hits",0,"record"],"reason":"Check actual source reaction"}
{"operation":"search_fragment_precedents","arguments":{"query":"REVISED_CORE","target_smiles":"TARGET","timeout_seconds":30},"parent_ref":"SEARCH_REF","reason":"Revision justified by previous evidence"}
{"action":"quit"}
```

`event.artifact_ref` identifies a saved call. Optional `propose_fragment_queries`
followed by `action: search_prepared` searches a saved parent or selected
`variant_id`. See the [agent guide](../../chem_coworker/scientific_workspace/task_playbooks/retrosynthesis_fragments.md)
for shapes. Inspection records mean opened fields, not adjudicated chemistry.
The console owns one child process, validates each request through existing
operations, and clears candidate subsets between queries. Quit, EOF and deadlines
release the child. A deadline retains partial/unknown evidence; a later query
reloads the library. Automatic discovery and ordinary workspace calls are unchanged.

Run the target-only live development pilot with the configured agent runtime:

```powershell
python -m examples.ai_native.fragment_loop_pilot --run-live --output results/ai_native/my_loop_pilot --fragment-index PATH_TO_FRAGMENT_INDEX.sqlite --timeout 600
```

The two new cases supply no fragment or route. The agent chooses queries,
inspects results and decides when to stop; automatic fragment selection and web
search are excluded for this integration test. Saved prompts, runtime logs,
answers, call evidence and `pilot.json` expose failures and timings. This is an
unpaired development pilot, so it cannot establish a quality improvement over
automatic search. Concurrent machine workloads can affect its timings.

Workspace API usage:

Use a saved Python script with the recorded workspace:
```python
from chem_coworker.scientific_workspace import ScientificWorkspace
w = ScientificWorkspace(@INVESTIGATION_ROOT@)
print(w.describe())  # Selected inputs, baseline identity and recent event references.
print([entry for entry in w.operations.catalog() if entry['name'] in selected_names])
print(w.run_summary(operation_name, parameters))
```

Discover exact arguments for unfamiliar operations from the catalogue, choosing
relevant names yourself. Batch independent signature lookups; do not print the
whole catalogue or read implementation files solely to discover public arguments.
Reuse saved calls for unchanged inputs. Put executable script bodies inside
`if __name__ == '__main__':` so importing helpers does not repeat scientific calls.
Prefer saved Python files for structured inputs over nested shell quoting.
Use w.help(['capture_source_file', 'attach_literature_reaction']) for public helper
signatures and callable examples. Helpers are not in w.operations.catalog().
Do not guess old module paths or read implementation files to learn these APIs.
For several results, collect event.artifact_ref and print w.batch_summary(refs)
once. Its 16 KiB budget applies to the entire batch, with every omitted ref retained.
`w.run(...)` returns an InvestigationEvent, not a JSON object. Save only
`event.artifact_ref` in JSON; `w.call_summary(event)` reads the recorded result.
For repeated fragment queries, reuse one persistent terminal process:
`python -u -m chem_coworker.scientific_workspace.console WORKSPACE`, using the
recorded interpreter and investigation directory. Send one JSON line at a time:
`{"operation":"search_fragment_precedents","arguments":{"target_smiles":"TARGET","query":"CORE","timeout_seconds":30},"reason":"Construction question"}`.
Read each response before choosing the next query; allow time for the first index
load. Use the catalogue for query arguments and
`w.task_guide('retrosynthesis_fragments')` for additional console examples.
If a persistent terminal is unavailable, batch independent calls in one saved
Python process and reuse recorded results.

Manifest investigation.json identifies selected datasets, versions and limitations.
Events and artifacts retain complete calls and notes. Start with w.run_summary(...)
or w.call_summary(ref) for a saved result. No custom summary script is needed.
Console summaries are bounded to 16 KiB; any omitted decision detail is explicit.
Fragment summaries retain publication metadata and source reaction structures.
Inspect only a decision-critical field when needed:
```python
print(w.inspect_artifact(event.artifact_ref, path=('result', 'warnings'), limit=5))
```
Paths are literal JSON keys/list indices, not expressions; next_offset continues a
selected page. Read warnings and non-passing gates before relying on results.
Counts refer to saved results; displayed lists may contain only a preview. Do not
dump whole analyses, manifests or route trees to stdout or rerun a tool for detail.
Use w.store.read_artifact(ref) for complete data in custom analysis. The store
verifies and reconstructs linked compact sections; raw JSON may have a storage
envelope rather than the expanded result. Cite the parent scientific-call ref.

Recipe assessments report a versioned coverage object separately from admission
and score. Inspect coverage.capability_status and evaluated rule IDs: an empty
checked_requirements list does not mean conflict rules were skipped. A score of
one means no encoded penalty, not experimental confidence. No applicable capability
requirement, unresolved identities, and an unassessed reaction remain explicit.
Attribute each recipe check to its saved recipe; changing that recipe requires a
new assessment. Structural precedent support does not establish condition support.
Before recommending concrete operating conditions, use assess_proposed_recipe
with the actual complete condition components and operating_conditions, citing the
source evidence. It resolves identities through condition_registry and checks the
actual recipe. Preserve temperature/time stages and atmosphere when reported;
do not replace an addition/activation/workup sequence with an invented single stage.
If the signature or recipe is unresolved, retain that status; never invent mapping
or use an unrelated analogue check as the target recipe's validation.

Saved disconnect_target realization_id identifies a template realization and can
repeat across different concrete precursor sets. Select with strategy_id AND
precursor_smiles from the chosen saved candidate for inspect_step_precedents or
assess_retro_validity. Ambiguous IDs are rejected, not silently selected.

For route-finding tasks, apply these decision checks without needing to load a guide:

- Identify the shared synthesis bottleneck and a credible starting-material entry
  before expanding finishing variants. Spend the next search on the unresolved
  question most likely to change the route choice. Distinguish evidence for the
  scaffold, substitution position and partner compatibility from evidence for
  their combination in the proposed step.
- Before stopping at advanced route inputs, use inspect_route_inputs on the saved
  route assessment, supplying leaf_queries=[{'smiles': leaf, 'terms': [source_name]}].
  For unresolved preparation or supply, search exact structures/names in captured
  sources and accessible primary literature. Inspect the passages behind a lead;
  a name match, molecular-weight stop or vendor listing does not establish access.
  If a plausible upstream preparation or alternative entry addresses the gap,
  assign explicit structures, assess it and revise the branch. A limitation note
  alone should not replace an available decision-changing check.
- If discovery is partial or transfer reports incomplete_search, keep the source
  as a lead. When transfer could change the route choice, run an explicit
  search_fragment_precedents query using the promising hit's discovery.query,
  query_format and topology, with the target. Narrow a still-broad query to one
  relevant attached region while preserving the core and substitution position.
  Investigate the new search result; inspecting an old artifact does not complete
  a search. Never clear a partial flag or infer absence from a bounded empty result.
- Compare precursor complexity and preparation burden, including protection and
  deprotection costs, alongside structural support. Reassess a proposed change;
  more references or structural admission do not establish experimental feasibility.
- Keep follow-up work bounded to a question that could change the decision. Stop
  when evidence resolves it or distinct reasonable searches leave an explicit gap;
  record the attempted check, remaining uncertainty and reason for stopping. If
  upstream access remains unresolved, label the proposal as a route from advanced
  inputs and compare alternatives on that same starting scope. State the preferred
  route and its shared bottleneck in the concise answer summary.

Notes use w.store.note(kind, text, evidence_refs=(ref,)); kinds are 'hypothesis',
'decision', 'question', 'limitation' and 'review'. Record branch choices as 'decision'.
This prompt, catalogue and optional guides are the starting reference. For unresolved
usage questions, read only the relevant section of @README@. Do not read the README
at startup.

Custom analysis is a normal capability. Read repository source, compose public APIs,
and save scripts and outputs inside this investigation. Attach useful files with
w.store.attach_file and label derived analysis. For recorded calculations:
```python
event = w.run_python('analysis.py', parameters,
                     evidence_refs=(source_ref,), timeout_seconds=60)
```
The script receives input.json and output.json paths as sys.argv[1:3]. Read parameters
and evidence from the input and write JSON output. Code, inputs, exit status and
output are recorded. Execution records a calculation, not its correctness.

Literature access is an application capability, independent of deterministic chemistry:
- w.fetch_source(url, title=title) saves a raw HTML/text/PDF snapshot and extracted
  text, including recorded retrieval failures.
- w.inspect_source(source.artifact_ref, query='Example', limit=4000) reads bounded
  passages with character/line/page locations. Follow next_offset for details.
- w.record_source_excerpt(source.artifact_ref, start=START, end=END,
  locator='Example 1') saves an exact passage. Alternatively supply excerpt=VERBATIM_TEXT.
- w.capture_source(text, url=url, title=title, locator=locator) preserves text already
  inspected with browser/search tools. It is agent_supplied_excerpt, not an independent
  fetch of the URL. Keep raw external sources separate from local-corpus evidence.
- w.capture_source_file('browser.txt', url=url, locator=locator, reference_id=REF1)
  imports an exact saved UTF-8 browser export without retyping it or printing it.
  JSON strings are supported; structured JSON requires text_path=[literal keys].
  reference_id is optional agent-attributed publication identity, not verification.
  reference_id must be REF1:<64 lowercase hex> from an indexed record, not a DOI;
  put a DOI in url/title or omit reference_id for a new external publication.
- w.search_captured_sources is not a helper: use w.run('search_captured_sources',
  {'queries': [literal_name], 'limit': 3}). Read matches, record exact offsets, and
  follow next_offset. It searches captured roots, including earlier upstream text.
- w.capture_source_image('scheme.png', source_ref=source_ref, locator='Scheme 1, p. 2')
  preserves a local PNG/JPEG inside the investigation. Pass its artifact_ref as
  scheme_refs when preparing a literature drawing. A stored screenshot supports
  visual comparison; it does not certify graph assignment or experimental stereo.
- If direct access reports network_permission_denied, use an available browser and
  capture passages for subsequent URLs while the same denial applies. Preserve the
  failure/acquisition scope. Do not repeatedly retry or relax the sandbox.
  Subsequent direct fetches are skipped after a recorded permission denial.
  Use retry_network=True only after a concrete network-permission change.
- w.capabilities() checks local imports and data presence. Requested web/model settings
  do not establish working tools; disclose unavailable search, full text or PDF parsing.

Search snippets are leads, not inspected experiments. Source snapshots can omit
drawings, tables or scans; inspect the original when these matter. An unspecified
stereoisomer does not establish exact stereo. Literature identity needs comparison
with the actual graph, chemical form, protecting groups and specified stereochemistry;
computed CIP labels do not establish experimental enantiopurity.
After exact/close searches, state the missing decision-relevant evidence. Repeat a
search only with a distinct structure, source or question that could change the answer.
Usually try one alternative access path after a failed publisher/SI request, then use
accessible primary evidence or disclose the gap. This is advice, not a hard quota.
Once the answer is supported or its limits are clear, finalize rather than keep collecting.

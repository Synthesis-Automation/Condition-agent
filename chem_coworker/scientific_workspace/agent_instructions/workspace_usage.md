Workspace API usage:

Use a saved Python script with the recorded workspace:
```python
from chem_coworker.scientific_workspace import ScientificWorkspace
w = ScientificWorkspace(@INVESTIGATION_ROOT@)
print(w.store.summary())  # Prior calls and notes for continuation.
print([entry for entry in w.operations.catalog() if entry['name'] in selected_names])
event = w.run(operation_name, parameters)
print(w.call_summary(event))
```

Discover exact arguments for unfamiliar operations from the catalogue, choosing
relevant names yourself. Batch independent signature lookups; do not print the
whole catalogue or read implementation files solely to discover public arguments.
Reuse saved calls for unchanged inputs. Put executable script bodies inside
`if __name__ == '__main__':` so importing helpers does not repeat scientific calls.
Prefer saved Python files for structured inputs over nested shell quoting.

Manifest investigation.json identifies selected datasets, versions and limitations.
Events and artifacts retain complete calls and notes. Start with w.call_summary(event).
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
- If direct access reports network_permission_denied, use an available browser and
  capture passages for subsequent URLs while the same denial applies. Preserve the
  failure/acquisition scope. Do not repeatedly retry or relax the sandbox.
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

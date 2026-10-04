# Literature preparation and bounded agent output

Status: development integration, 2026-10-04. Independent chemistry review and
untouched-evaluation release gates remain unchanged.

## Problem and resulting behavior

The October 4 agent investigation found a paper through local scaffold matches,
then authored source structures and several JSON-summary scripts manually. The
saved command output exceeded 250,000 characters. SMILES validity was checked,
but the source's compound assignments were not independently verified.

New investigations can compose exact source excerpts and saved corpus structures
through `prepare_literature_reaction`. It produces `literature_reaction.v2`, with
an immutable preparation reference, evidence for each participant, explicit
material-form notes, graph checks and separate source-acquisition provenance.
The answer validator rejects changes to the prepared block. Existing v1 drawings
and their historical self-review identities remain readable; they do not acquire
invented verification. The v2 preparation workflow is the authoring path for new
drawings. Saved evidence is not rewritten or admitted into an index.

Graph auditing belongs to `reactive_taxonomy.structure_audit`; the workspace
adapter records and composes results. The audit preserves isotope identity,
specified stereo, charges and all disconnected salt components. It reports
invalid graphs, formula disagreements and unspecified stereo. It does not
translate names, recognize scheme images, verify compound numbering or tautomer
assignments, infer atom correspondence, or establish feasibility.

## Efficient inspection

```python
from chem_coworker.scientific_workspace import ScientificWorkspace

w = ScientificWorkspace(investigation_directory)
print(w.describe())
print(w.run_summary("analyze_molecule", {"smiles": target_smiles}))
print(w.call_summary(saved_call_ref))
print(w.inspect_artifact(saved_call_ref, path=("result", "hits"), offset=3, limit=3))
```

Complete scientific inputs and results remain saved. Console call summaries are
bounded to 16 KiB and explicitly identify omitted detail. Fragment-search briefs
include source reaction structures, IDs, source yields and available publication
metadata, including nested source bibliography. Read the indicated saved fields
before relying on truncated structures, omitted warnings or non-passing gates.
Saved-result inspection does not rerun chemistry or fabricate new call events.

CLI equivalents:

```powershell
python -m chem_coworker.scientific_workspace describe WORKSPACE
python -m chem_coworker.scientific_workspace catalog WORKSPACE --operation prepare_literature_reaction
python -m chem_coworker.scientific_workspace call-summary WORKSPACE sha256:REF
python -m chem_coworker.scientific_workspace show WORKSPACE sha256:REF --path '["result","warnings"]'
```

`show` now defaults to bounded, paged inspection. Use `show ... --full` explicitly
for a complete expanded JSON export. This is a CLI output-contract change;
scientific records, chemistry results and web API routes are unchanged.

## Capture source text without transcription

Save browser-obtained text in a UTF-8 file inside the investigation, then import it:

```python
source = w.capture_source_file(
    "browser.txt", url=paper_url, locator="Experimental section",
    reference_id=corpus_reference_id,
)
print(w.call_summary(source))
```

A JSON string export is automatically unwrapped. Structured JSON requires a
literal `text_path`, for example `["document", "text"]`. Source bytes, selection
path and text are hashed; the captured source is still agent-supplied and not
HTTP-verified. `reference_id` is an explicit agent-attributed bibliography link.
It does not certify that the exported text came from that publication.

After a recorded network-permission denial, later direct fetches are skipped and
linked to that failure. Use the available browser/web tool instead. If the
transport permissions actually change, `fetch_source(..., retry_network=True)`
or CLI `fetch-source ... --retry-network` makes an explicit new attempt.

## Prepare source reactions

First save exact source passages with `record_source_excerpt`, using character
offsets or a unique literal passage. Every reactant and product needs an excerpt
from the same captured root, a source compound label or name, and material-form
notes. For reconstructed structures, supply explicit SMILES:

```python
event = w.run("prepare_literature_reaction", {
    "source_ref": source.artifact_ref,
    "source_id": "lit1", "title": "Source Example 2", "locator": "Example 2",
    "structure_evidence": "Ethanol was oxidized to acetaldehyde.",
    "reactants": [{
        "name": "Ethanol", "compound_id": "Ethanol", "smiles": "CCO",
        "evidence_ref": passage_ref, "material_form": "Neutral parent; purity unavailable",
    }],
    "products": [{
        "name": "Acetaldehyde", "compound_id": "acetaldehyde", "smiles": "CC=O",
        "evidence_ref": passage_ref, "material_form": "Source isolation form unverified",
    }],
})
print(w.call_summary(event))
step["literature_reactions"] = [w.prepared_literature_reaction(event.artifact_ref)]
```

Prefer recorded source graphs when available: supply `indexed_ref` from a saved
`search_fragment_precedents` or `get_precedents` call and `reaction_id`, and omit
all participant SMILES. Participants must match that record's side/component
counts, including salt components. The capture's `reported_reference_id` must
match the indexed publication identity. Ambiguous duplicate records are rejected.
The UI labels this as indexed source structures, separately from reconstruction
and source-explicit SMILES. It does not promote the bibliography join to verified.

Optional `reported_formula` must occur literally in the participant's excerpt;
disagreements remain visible rather than repairing the structure. Conditions and
yields must be reported claims quoting captured text and citing `source_id`.
Proposed or normalized target recipes stay on the proposed step.

For conflicting source details, supply `source_conflicts` entries with
`description` and at least two exact `excerpt_refs`. Both passages are linked and
the unresolved discrepancy remains visible. For example, preserve discussion
text naming DCM and experimental text naming DCE rather than silently merging
the solvents. The tool retains identified conflicts; it does not automatically
detect or resolve all textual chemistry disagreements.

Answer-finalization failures now enter the activity/debug log even if an outer
script catches the exception and exits with code zero. Failed drafts are not
published, and operational failure events are not citable scientific evidence.

## Deployment and validation scope

Restart the scientific-chat server and start a new investigation after updating
the scientific adapter, audit or answer validators. No corpus rebuild is required.
Old investigations remain readable. Prompt improvements cover both normal and
tools-only modes. This reduces avoidable output and script overhead; a faster
end-to-end agent run remains to be measured independently. No search quota,
chemistry filter, condition-ranking rule or release gate is relaxed.

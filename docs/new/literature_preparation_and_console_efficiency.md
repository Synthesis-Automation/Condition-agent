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

## Follow-up improvements from the two latest investigations

The October 4 evening investigations used the prepared-drawing and brief-output
paths. Their remaining gaps motivated the following additions. These do not
change chemistry definitions, rank external literature as verified structures,
or establish experimental feasibility.

### Select a concrete candidate

`realization_id` identifies a template realization; it can repeat for different
precursor sets. Inspection and retro-validity calls now reject ambiguous IDs.
Use the selected saved strategy and complete precursor SMILES:

```python
candidate = saved_strategy['representative']
inspection = w.run('inspect_step_precedents', {
    'source_ref': disconnection_ref,
    'realization_id': candidate['realization_id'],
    'strategy_id': saved_strategy['strategy_id'],
    'precursor_smiles': candidate['precursor_smiles'],
})
```

The same selectors are accepted by `assess_retro_validity`. Reactant ordering and
atom-map labels do not change graph identity; specified stereo and salt components
remain significant. Shared IDs for identical concrete candidates are supported.
No template or corpus IDs are rewritten. These two operation contracts are now v2;
old evidence remains readable, while replay requires its recorded implementation.

### Check unresolved route inputs against captured sources

```python
inputs = w.run('inspect_route_inputs', {
    'source_ref': route_assessment_ref,
    'leaf_queries': [{'smiles': actual_leaf_smiles, 'terms': [literal_source_name]}],
})
search = w.run('search_captured_sources', {
    'queries': [literal_source_name], 'limit': 3,
})
```

Searches cover recorded roots, including upstream sections omitted from a cited
excerpt. Matches retain exact source offsets for `record_source_excerpt`; follow
`next_offset`. Failed fetches are skipped by default. Searches are literal leads,
not name-to-structure assignments. The input inspection reports every actual leaf,
reuses saved starting-material assessments when available, and retains missing
terms as unsearched. Molecular-weight stops remain assumptions about planning,
not evidence of supply. Inspect and assess source-supported upstream reactions
before extending a route with `revise_route_branch`; otherwise disclose the route
as beginning with advanced inputs. This tool does not automatically expand routes.

### Assess and attach the actual proposed recipe

```python
check = w.run('assess_proposed_recipe', {
    'reaction_smiles': proposed_reaction_smiles,
    'components': [{'raw_identifier': 'ethanol', 'source_field': 'proposal',
                    'identifier_type': 'name', 'source_role_hint': 'solvent'}],
    'operating_conditions': {
        'stages': [{'stage_index': 0, 'temperature_c': 0, 'time_h': 0.25},
                   {'stage_index': 1, 'temperature_c': 25, 'time_h': 2}],
    },
    'evidence_refs': [source_excerpt_ref],
})
draft = w.attach_recipe_check(draft, 's1', check.artifact_ref)
```

Supply the complete actual proposed condition components, not only the solvent
in this API example. Registry identity resolution and canonical compatibility
rules remain authoritative. Unknown identities/signatures, missing capability
coverage and hard conflicts stay explicit. Quantities and stage provenance are
supported; numerical conditions are not inferred from prose. Preserve reported
activation/addition/workup order in the attributed procedure. Stage rules do not
constitute simulation of that procedure. The current canonical rules assess the
recipe's top-level operating fields; stage-specific mixtures and temperatures
remain explicitly `stages_recorded_not_evaluated`. This gap appears in preflight
and the browser. `recipe_assessment_refs` is an optional
answer-step field, with exact structure/stereo binding and browser coverage views.
It refers to the saved recipe; altered prose is not automatically reconciled.
Historical answers/review hashes without this field remain supported.

### Preserve source images separately from reconstructed graphs

```python
image = w.capture_source_image('scheme.png', source_ref=source_ref,
                               locator='Scheme 1, page 2')
# Include scheme_refs=[image.artifact_ref] in prepare_literature_reaction.
```

PNG/JPEG files must be inside the investigation, at most 2 MiB and 20 million
pixels. Bytes and source lineage are preserved. The browser shows captured images
inside source-evidence details, with a shared six-image/4 MiB source-byte display
budget and explicit omissions. Image capture does not verify the page attribution
or chemical assignment. The preparation remains assignment-unverified. Literal
compound names now permit 300 characters; ambiguous name suffixes are unnecessary.
Preparation's operation contract is v2; the saved literature drawing remains v2.

### Discover helpers and avoid assembly retries

```python
print(w.help(['capture_source_file', 'attach_literature_reaction', 'finalize_answer']))
print(w.batch_summary(saved_event_refs))  # One shared 16 KiB budget, 1..30 refs.
draft = w.answer_template(answer_markdown)
# Populate scientific content, explicit molecule/step IDs and attribution.
draft = w.attach_literature_reaction(draft, 's1', preparation_ref)
print(w.answer_preflight(draft))
receipt = w.finalize_answer(draft_path, draft, findings=findings)
```

`help` describes named public helpers as well as registered scientific operations.
`capture_source_file` documents the distinction between a corpus REF1 identity and
a DOI; put a DOI in URL/title or omit the optional corpus join. Drawing attachment
adds its exact required source reference and immutable block together. Conflicting
source IDs are rejected. Preflight validates the draft and reports missing recipe
or route-input checks without rerunning scientific operations; follow `next_offset`
to inspect additional warnings. Those warnings are
advisory; unresolved scientific work must remain disclosed. Historical answers
are not invalidated by new authoring advice. Use literal strings/f-strings for
scientific prose containing `%`, rather than percent interpolation.

## Deployment and validation scope

Restart the scientific-chat server and start a new investigation after updating
the scientific adapter, audit or answer validators. No corpus rebuild is required.
Old investigations remain readable. Prompt improvements cover both normal and
tools-only modes. This reduces avoidable output and script overhead; a faster
end-to-end agent run remains to be measured independently. No search quota,
chemistry filter, condition-ranking rule or release gate is relaxed.

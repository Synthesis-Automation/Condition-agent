# Fragment precedent adviser v4

Optional supporting advice for `retrosynthesis`. Use it when constructing an unfamiliar
core is a decision-changing uncertainty. A supported route or clear disconnection may
already answer the question.

## Agent-chosen fragment loop

Propose a connected key fragment for a concrete question: formation of a ring or
junction, preparation of an intermediate, or installation of a functional group.
The fragment is a search hypothesis, not necessarily a precursor. Explain which
heteroatom positions, ring boundaries, valence and stereo constraints matter.
Search once, inspect actual source reactions, then decide whether to narrow,
broaden, change the selected region, follow a precursor, or stop. Do not execute
a predetermined ladder or require a minimum number of queries when evidence is
already sufficient. `suggest_search_fragments` is optional assistance.

For multiple searches, use one persistent terminal session:

```text
python -u -m chem_coworker.scientific_workspace.console WORKSPACE
```

Start with a persistent terminal (`tty=true` if the execution tool requires it),
then send one JSON line through that terminal's input tool and read the response
before choosing another query. Reuse this process; restarting it reloads the
index. The first search includes library loading; allow 30 seconds. Later calls
reuse the same library but search the full eligible index. Example request shapes:

```json
{"operation":"search_fragment_precedents","arguments":{"target_smiles":"TARGET","query":"AGENT_CHOSEN_CORE","timeout_seconds":30},"reason":"Find construction of the selected ring junction"}
{"operation":"propose_fragment_queries","arguments":{"target_smiles":"TARGET","query":"REVISED_CORE","query_format":"smarts","topology":"subgraph"},"parent_ref":"SEARCH_REF","reason":"Previous hits retained the ring; select the region around its formation instead"}
{"action":"search_prepared","query_ref":"PREVIEW_REF","reason":"Search the validated revised core"}
{"action":"inspect","source_ref":"SEARCH_REF","path":["result","hits",0,"matches"],"reason":"Check which bonds were actually formed"}
{"action":"quit"}
```

Response `event.artifact_ref` is the saved call ID. Inspect `record`, `procedures`
and other hit paths similarly; bounded inspection supports `offset` and `limit`.
`search_prepared` can select a saved `variant_id` without recopying SMARTS.
Errors, deadlines and empty results remain distinct. A killed worker restarts on
the next request. Quit when finished. Ordinary `w.run` calls remain available
when a persistent terminal is unavailable, with a fresh worker per search.

SMILES substructure queries do not automatically constrain carbon hydrogen count
or prohibit extra carbonyl substitution. If returned lactams differ from a target
amine, inspect that difference and add an explicit valence constraint where useful.
Omitting stereo permits discovery without establishing stereochemical transfer.
Changing a ring size or heteroatom position creates an analogue hypothesis: do not
claim it is a matching target fragment or silently replace the target.

## Optional automatic target-first entry

Use `find_synthesis_precedents(target_smiles=...)` as an optional aid when the user wants
relevant construction precedents without choosing query syntax. It selects
informative ring subsets, samples broad matches, adds context, and retains earlier
construction leads. It loads the fragment index once and reuses only completely
enumerated candidate sets for nested refinements. The default budget is 90 seconds.

Inspect `hits[].discovery`, references, source records and procedures. Explicit
relaxations include peripheral substitution, hydrogen counts and additional ring
fusion; N-H analogues can therefore appear for N-acetyl targets. Broad samples
remain partial, counts retain their precision, and unresolved correspondence is
not construction evidence. Ranking is a heuristic over retrieved candidates,
not independent chemistry validation or proof of a feasible route.

Use `call_summary`/`run_summary` for a compact view. For a promising source, copy
its `discovery.query`, `query_format` and `topology` into a recorded
`search_fragment_precedents` call with the target, then investigate the returned
observation using the detailed workflow below. No production retro library is
needed for discovery. The agent may choose its own fragments and query edits
directly; running automatic discovery first is not required.

## Purpose and context

Name the construction question and choose a connected core that represents it.
`suggest_search_fragments` offers graph-validated regions; suggestions overlap and
are not precursors. Inspect boundaries, distinctive ring topology and heteroatom
positions before choosing a query for `search_fragment_precedents`.

## Explicit query and precedent investigation

Use `propose_fragment_queries(target_smiles, query, query_format, topology)` when
an explicit relaxation can resolve the construction question. It previews removal
of peripheral context and relaxation of ring boundaries without searching. For
C/N analogues, select eligible input-query atom IDs from `query_atoms` and call
again with `aromatic_atom_ids`. Pyrrole [nH], charge and isotope constraints are
not generalized. Custom SMARTS remains available for other hypotheses.

Choose one variant, inspect `reason`, `relaxations` and target alignments, then
call `search_fragment_precedents` with its exact query/format/topology and target.
Include `query_variant_ref` and `query_variant_id` to retain and validate the saved
preview link. A mismatched parent query is an error, not a reason to relax silently.
Try a small number of informative alternatives; do not exhaust every variant.

For an interesting returned observation, call
`investigate_fragment_precedent(source_ref=search_ref, observation_id=...)`.
It preserves the actual source record/procedures, compares its product with the
target, and runs a bounded source-only transfer using the existing strict compiler.
It needs no production retro library. Inspect `source`, `comparison` and `transfer`
in the saved result. Compilation rejection must not erase the source's chemical
value: use its reference/procedure to investigate applicability or missing evidence.
Incomplete searches remain inspectable but cannot seed transfer. A complete empty
bounded transfer is not chemical impossibility. The operation does not change
library admission, infer experimental feasibility or validate transferred conditions.

## Scientific questions

- Which core matters chemically, rather than merely looking complex or scoring highly?
- Would one or two informative queries resolve its construction question?
- If a whole-target search finds only retention, would removing peripheral substituents
  reveal useful construction evidence while preserving the distinctive core?
- Which observed operation can transfer, and which handles or stereochemistry differ?

## Evidence distinctions

Supply the target for target-derived queries, including broadened ones, so the tool
checks the actual query semantics. A target mismatch is an error, not zero hits.
Changing ring topology or aromaticity is not peripheral substitution; disclose an
explicit SMARTS/subgraph relaxation and inspect whether it still matches the target.
`too_broad`, `partial` and zero hits retain their bounded index meaning.

Inspect local bond-change witnesses, exact observation/reference IDs, procedures and
source scope. A product containing the core does not establish construction. Separate
construction, modification, boundary changes, retention and unresolved correspondence;
these may overlap. A construction witness may establish only part of the core.
Support claims with specific inspected records and fields, preserving missing data.
Explain what transfers to this target and what remains a proposed adaptation.

## Stopping criteria

Stop once evidence informs the route choice or the bounded search leaves a stated gap.
Refine or search another region only for a remaining question; reuse unchanged queries.
A hit is not an automatically validated route or a reason for forward prediction.

# Fragment precedent adviser v1

Optional supporting advice for `retrosynthesis`. Use it when constructing an unfamiliar
core is a decision-changing uncertainty. A supported route or clear disconnection may
already answer the question.

## Purpose and context

Name the construction question and choose a connected core that represents it.
`suggest_search_fragments` offers graph-validated regions; suggestions overlap and
are not precursors. Inspect boundaries, distinctive ring topology and heteroatom
positions before choosing a query for `search_fragment_precedents`.

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

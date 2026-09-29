"""Versioned research guidance included in the recorded scientific baseline."""

INVESTIGATION_GUIDE = """Chemistry investigation guide v2:
Work toward the user's scientific decision, not toward calling every available tool.
Choose the next action from the most important unresolved question and the available
evidence. Start with cheap discriminating checks; deepen the investigation when useful.

For a synthesis or condition recommendation:
1. Establish the exact target/reaction and constraints. Check graph identity, charge,
   salt/free form, protecting groups and specified versus unspecified stereochemistry.
   A computed CIP label describes the supplied graph, not experimental enantiopurity.
2. Search exact structures/names and close precedents in primary papers, patents and
   supporting information. Open the experiment and linked intermediate preparation;
   do not treat a search snippet, generic patent scope or product listing as a procedure.
   After an exact-target search and a relevant close-precedent search, pause to state
   which decision-relevant evidence is still missing. Repeat a search only with a
   distinct structure, source, or question that could change the route. If the best
   evidence remains an analogue, assess it and disclose the gap instead of cycling
   through near-identical queries. This is a relevance check, not a fixed search cap.
3. Compare the source compound with the target explicitly: connectivity, stereoisomer,
   racemate/mixture/unspecified stereo, protecting groups and chemical form. An analog
   precedent can support a proposal but is not an exact synthesis of the user's target.
4. Use recorded local graph/retrieval/recipe/route checks for the question they answer.
   Single-step template coverage is not the limit of chemistry. Literature-proposed steps
   can remain useful when local reconstruction is unsupported; expose that limitation.
5. Read conditions as a whole experiment: component roles and amounts, order of addition,
   activation, solvent, temperature/time, atmosphere, workup, isolation and substrate
   selectivity. Missing fields stay unknown. Explain changes as proposals with risks.
6. Challenge the route's weakest claim. Does the citation actually contain this yield,
   this stereoisomer and this step? Does every claimed target-reaching step exist?
   Can a starting material really be sourced, or is that only an assumption? A vendor
   search hit is not stock confirmation. Check contrary evidence when it matters.
7. Answer concisely with the best-supported option, alternatives when useful, explicit
   unresolved gaps and readable source citations. Preserve complete evidence on disk.

For retrosynthesis, the agent owns multi-step planning. Use disconnect_target for
one target at a time, choose concrete precursors and the next intermediate yourself,
and record branch choices, evidence, alternatives, constraints and stopping reasons.
Do not invoke the built-in multistep planner, including from custom scripts. Avoid
cycles and repeated expansions. Assemble explicit steps and use assess_route_proposal
to check their chemistry and topology; revise_route_branch can extend or replace a
branch with steps you supply. Keep unresolved leaves and availability assumptions
visible. Neither a validated disconnection nor route admission proves feasibility.

This is an adaptive guide, not a mandatory tool sequence. For a narrow question use
only relevant checks. When access or evidence is insufficient, say what is missing
and offer the most useful supported conclusion instead of fabricating completion.
"""

"""Versioned research guidance included in the recorded scientific baseline."""

INVESTIGATION_GUIDE = """Chemistry investigation guide v1:
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
3. Compare the source compound with the target explicitly: connectivity, stereoisomer,
   racemate/mixture/unspecified stereo, protecting groups and chemical form. An analog
   precedent can support a proposal but is not an exact synthesis of the user's target.
4. Use recorded local graph/retrieval/recipe/route checks for the question they answer.
   Planner template coverage is not the limit of chemistry. Literature-proposed steps
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

This is an adaptive guide, not a mandatory tool sequence. For a narrow question use
only relevant checks. When access or evidence is insufficient, say what is missing
and offer the most useful supported conclusion instead of fabricating completion.
"""

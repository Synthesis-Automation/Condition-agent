# Chemistry investigation guide

The investigator uses the existing Codex tool loop and recorded scientific workspace.
The executable prompt guidance lives in
[`investigation_guide.py`](../../chem_coworker/scientific_workspace/investigation_guide.py),
which participates in the code baseline. This guide describes its scientific intent.

1. **Identify the question and target.** Check connectivity, charge, chemical form,
   protecting groups and specified stereochemistry. Calculated CIP labels refer to
   the supplied graph; they do not establish experimental enantiopurity.
2. **Find the experiment.** Search primary papers, patents and supporting information.
   Read the relevant example and precursor preparations. Search snippets, broad
   patent claims and vendor listings are leads; they are not experimental procedures.
3. **Compare source and target.** State exact versus analog evidence. Compare
   connectivity, stereochemistry, salt/free form and protecting groups. A cis
   racemate precedent does not establish a single-enantiomer route. Unspecified
   stereochemistry does not establish that a sample was racemic.
4. **Use the appropriate local checks.** Record graph analysis, retrieval, recipe
   compatibility, route assessment and custom scripts. Unsupported reconstruction
   means the tool has not established the step, not that the chemistry is impossible.
5. **Read conditions in context.** Preserve component identities/roles, quantities,
   addition order, activation, solvent, temperature, duration, atmosphere, workup
   and isolation. Leave missing information unknown. Separate proposed adaptations
   from conditions actually reported for that substrate.
6. **Challenge the weakest claim.** Check yield attribution, stereoisomer identity,
   missing steps, starting-material availability and conflicting precedents. Do not
   promote a procurement assumption into confirmed stock or an analog into an exact
   experimental synthesis. Record an honest self-review of the final draft.
7. **Present the decision clearly.** Lead with the best-supported conclusion, show
   attributed reaction steps, link original sources, and retain material limitations
   in the main answer. Keep detailed calls and raw snapshots available for inspection.

This is an adaptive guide. A short identity question does not require a full route
search; a difficult synthesis may need repeated literature/tool cycles. The agent
chooses actions that resolve uncertainty and can stop with a useful partial result.

Evidence capture and deterministic attribution checks preserve traceability. They
do not replace chemical judgment, independent review or experimental validation.

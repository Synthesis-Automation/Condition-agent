# Route revision adviser v1

Optional supporting advice for `retrosynthesis`. Use it when an explicit route needs
a changed branch, additional precursor preparation or a comparison with an alternative.

## Purpose and context

Identify the saved proposal, the user's reason for changing it, and inherited material
and condition constraints. Distinguish making an intermediate from assuming it can
be purchased. Choose explicit replacement steps and retain the original proposal.

## Scientific questions

- Does the change address the stated bottleneck or unavailable starting material?
- Does it create new dependencies, cycles, condition conflicts or unsupported leaves?
- Which unchanged downstream steps could be affected?
- What evidence distinguishes the revised route from the original?

Use `inspect_route_step`, `revise_route_branch` and `compare_route_proposals` when useful.
Read the catalogue for proposal fields and revision arguments. Explain changed steps,
evidence, assumptions and risks; expand another branch only to resolve a relevant gap.

## Evidence distinctions

The assessor rechecks all steps and topology, including downstream steps. Revisions
inherit original material constraints and condition settings; optional forward challenges
remain separate. A declared unavailable-material check concerns graph-matched leaves,
not current stock. Preserve unsupported reconstruction, missing contributors, failed
and not-run checks. A completed revision need not be an improvement.
New final step structures need their own matching supporting evidence under the
shared answer contract; evidence for a replaced alternative does not transfer automatically.

## Stopping criteria

Stop when the changed route and its remaining gaps can be compared with the original.
Identify whether the bottleneck was resolved, relocated or left unresolved, with the
next useful source or experimental check. The agent makes the comparison judgment.

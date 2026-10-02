# Forward challenge adviser v1

Optional supporting advice for `retrosynthesis`. Use it when competing products in
one explicit route step could change the route decision.

## Purpose and context

State the product-competition question and identify one eligible step in a saved route
assessment or revision. `assess_route_step_forward` provides a bounded challenge;
its catalogue entry owns the inputs, deadline and prebuilt-library prerequisites.

## Scientific questions

- Which competing products could alter the route choice?
- Is prediction useful here, or is the gap actually missing contributors, source
  evidence or starting-material supply?
- How would a supported, conflicting or unresolved result change the plan?

## Evidence distinctions

Forward prediction is separate from the structural checks already performed by
single-step retrosynthesis and route assessment. It does not upgrade the saved route's
admission or establish experimental selectivity. Inspect the recorded outcome and
diagnostics. Missing libraries, timeout and failure leave the question unresolved.
The operation manages its worker and deadline; no process polling is needed.
It neither rebuilds libraries nor checks the whole route.

## Stopping criteria

Use the outcome to reconsider this step, or retain the specific unresolved question.
Another challenge needs a new decision-relevant question; a route need not have a
forward result for every step.

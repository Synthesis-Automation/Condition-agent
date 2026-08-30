# Transition-provider orchestration design and status

## Purpose

This layer lets an agent select among deterministic single-step transition
providers without allowing the agent to write or edit molecular structures.
It adapts the useful control pattern from multi-model retrosynthesis agents to
the chemistry-first contracts already present in this repository.

The first implementation is deliberately shadow-only: it compares provider
outputs but does not mutate a route tree or alter the public planner policy.

## Ownership and dependency direction

`core_retrosynthesis.transition_orchestration` owns the search-control contract.
Providers may call existing `core_retrosynthesis` generators. Candidate
observation and chemistry remain owned by `reactive_taxonomy`; condition
identity and recommendation remain owned by their existing standalone packages.
`chem_coworker` may later expose identifier-only actions over this contract, but
must not implement provider chemistry.

## Agent action

```python
ExpandLeafAction(
    leaf_id="...",
    provider_id="...",
    proposal_budget=3,
)
```

The action intentionally contains no SMILES, SMARTS, reaction edits, condition
strings, or provider-specific parameters. `ExpansionState` resolves the leaf ID
to an environment-owned canonical molecule. The provider ID must be registered,
and the requested budget must not exceed its declared capability.

## Deterministic boundary

Each provider returns `TransitionProviderBatch`. The orchestrator then checks:

1. the provider-declared forward-validation status;
2. target and precursor structural identity;
3. agreement between candidate fields and reaction notation;
4. hard precursor and reaction compatibility dispositions;
5. independent reproduction of a generic reaction signature from retained
   mapped validation evidence when available;
6. duplicate transition identity.

Mapped validation evidence is used only after its map-stripped precursor and
product structures exactly match the candidate's display structures. This
avoids false rejection when correspondence inference from an unmapped display
reaction is ambiguous, while rejecting contradictory mapped evidence.

Rejected candidates remain visible with typed reasons. Admitted candidates carry
the provider ID and provider-local rank. Scores from different providers are not
compared because they are not assumed to be calibrated.

## Provider adapters

The current implementation includes:

- `CallableTransitionProvider`, for adapting an existing deterministic callable;
- `OperatorLadderTransitionProvider`, for the canonical generic operator ladder.

Adding another provider requires a stable metadata record and an `expand`
implementation returning ordinary `GenericDisconnectionCandidate` objects. A
provider does not receive permission to bypass candidate admission.

## Shadow comparison

`TransitionProviderOrchestrator.compare_shadow` runs a bounded ordered set of
providers on the same leaf and reports:

- raw, admitted, rejected, and duplicate counts per provider;
- provider-local diagnostics;
- the union of deterministic transition identities;
- provider attribution and local rank for every transition; and
- the number of transitions shared by more than one provider.

The comparison is read-only. Its immediate purpose is to determine whether
candidate providers are genuinely complementary before a router is allowed to
affect search.

## Current status

Implemented:

- identifier-only bounded action;
- immutable state, provider, result, rejection, and shadow-report contracts;
- canonical generic-operator adapter;
- deterministic provider registry and budget validation;
- target/reaction/signature/compatibility admission checks;
- within-provider and cross-provider deduplication;
- public package exports and focused regression tests.

Not yet implemented:

- conversion from a live `ReactionRouteTree` into an agent observation;
- `chem_coworker` tool exposure;
- persistent provider telemetry and replay artifacts;
- a provider router with planner influence;
- external single-step model adapters; or
- reinforcement learning.

External LLM, web, literature, chemist, and third-party SSR proposals require a
separate evidence-reconstruction boundary before they can become provider
candidates. See
[`external_retrosynthesis_proposal_admission_design.md`](external_retrosynthesis_proposal_admission_design.md).

## Next release gate

Before adding an agent router, register at least two meaningfully different
existing providers and run them in shadow mode on the untouched route-action
panel. Record signature admission, compatibility rejection, unique transition
coverage, observed-action recovery, solved-route contribution, calls, and wall
time. Promote routing only if complementarity survives scaffold- and
patent-disjoint evaluation. Structural-core observations may stratify that
evaluation, but they do not establish disconnection validity.

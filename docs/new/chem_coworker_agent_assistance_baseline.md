# Chem Coworker Assistance Baseline

**Frozen:** 2026-08-20  
**Purpose:** migration and replay baseline; not an activation approval

Before the assistance implementation, the complete repository suite passed:

```text
1106 passed in 210.65s
```

The existing condition, one-step, and multistep reviewer tests captured these
compatibility properties:

- schema-constrained provider responses;
- one retry for incomplete output;
- token and provider-attempt accounting;
- fail-closed unknown evidence references;
- unchanged deterministic results and stored ranks;
- visible provider failures with secret-like values redacted; and
- solved-versus-partial route preservation.

The implementation retains those reviewer tests as golden transport fixtures.
The three compatibility reviewers now delegate provider requests and retry
behavior to `chem_coworker.assistance.transport`; their domain-specific result
adapters remain until the external activation and serialized-parity gates permit
removal.

Current assistance definitions at the time of this baseline:

```text
assistance_prompt.v1.txt sha256:66838461a47b0f2260e171bcd860f03d99c82b2f3517772db99d1992ca9fba3c
assistance_policy.v1.json sha256:149cb18987ea6d43743feeccab33ae869c2d30b62056d160fe8875ec085774f6
assistance state schema: chem_coworker_assistance.v1
condition constraint schema: 1.0
route search policy schema: 1.0
generic recommendation result schema: 3.5
```

Activation remains blocked on the primary roadmap's blind chemist review,
untouched evaluation, and full-corpus gates. These require external review data
and are intentionally not inferred from unit-test success.

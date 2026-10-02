# BINAP registry identity correction

Date: 2026-10-02

## Problem and evidence

The registry audit reported two `UNDECLARED_AMBIGUOUS_IDENTIFIER` issues for
`cas:76189-55-4` and `cas:98327-87-8`. Both public audit tests correctly detected
this inconsistent curation.

[TCI's R-BINAP entry](https://www.tcichemicals.com/US/en/p/B1406) identifies CAS
76189-55-4 as the R enantiomer;
[TCI's racemic entry](https://www.tcichemicals.com/JP/en/p/B2383) identifies CAS
98327-87-8 as the racemate. The former registry record had a corrupted canonical
name and incorrectly included `rac-BINAP` and a `RAC-` long-name alias.

The unqualified `BINAP` alias also collided with the normalized racemate canonical
name `BINAP(+/-)`. That ambiguity had not been explicitly declared.

## Correction and contract impact

- Correct the R record's canonical name to `(R)-BINAP`.
- Remove racemate aliases from that record. Retain `rac-BINAP` on the racemate
  and move the `RAC-` long-name alias there.
- Declare the unqualified `BINAP` alias shared on both records. A bare name
  continues to return ambiguous candidates; no enantiomer is guessed.

Only these two records' names/aliases change. All 27,434 substance records remain
accepted, with 55,605 identifiers and no identifier audit issues. No CAS numbers,
substance IDs, roles, molecular graphs, normalization rules or schemas change.
The stored graphs do not encode the axial distinction; structure-only resolution
therefore remains ambiguous. Regression tests cover explicit names/CAS, the
previously misassigned racemate long name, and unqualified name/structure ambiguity.
The CLI test uses `has_errors`, which accounts for both substance and identifier issues.
The audit tests are retained; validation is unchanged.

Definition schema remains `substances.v2.jsonl`. Its content identity changes:

| Snapshot | SHA-256 |
| --- | --- |
| Before | `038f2bc8a0744e93869efb9ecb56b6c4f7189f677e8b4528a0a8688bc8d588cd` |
| Corrected | `1ea051ed3082d176c9cc5b0617ef9733ab591b4286c60a7c75b0de97f4d4c94e` |

Name-only inputs using the two racemate aliases now resolve to the racemate;
CAS-based identity and the unqualified `BINAP` ambiguity retain their behavior.
The corrupted former canonical name is removed.
Existing converted corpora and saved investigations are not rewritten. A new
conversion consumes the corrected registry; scientific baseline hashes detect the
definition change. Start a new investigation when using the new scientific baseline.
Long-running applications must reload their cached registry, normally by restarting.

## Validation

- `python -m condition_registry.cli validate --format json`: 27,434 accepted
  substances, 55,605 accepted identifiers, zero errors.
- `pytest -q tests/condition_registry`: 77 passed.
- Complete `pytest -q`: 2,279 passed, 2 skipped, zero failures in 611.28 seconds.
- Ruff on changed tests and `git diff --check`: passed.
- Definition comparison against HEAD: exactly two records changed, limited to
  names/aliases; all other fields and records are identical.

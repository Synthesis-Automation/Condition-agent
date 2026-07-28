# Connectivity-First Reaction Grammar: Phase 4 Implementation

## Status

Implemented on 2026-07-28.

Phase 4 migrates parity-complete production reconstruction from
operator-specific Python branches to versioned, normalized connectivity
rewrites. Authority is granted per grammar, not per reaction name, operator
class, or rewrite file.

## 1. Expanded normalized interfaces

The existing immutable site interfaces now also adapt:

- activated `electrophilic_center` sites to a center-leaving reactive link;
- carbonyl centers to polarized double-bond capacity;
- alcohol sites to localized C-O single-bond capacity;
- `eliminable_pair` sites to a departing link, backbone bond capacity, and
  one hydrogen-loss endpoint; and
- `nucleophile_anion` sites to a connection endpoint requiring the observed
  charge transition from -1 to 0.

These remain views over curated detector output. No arbitrary molecule atom is
promoted to a reactive endpoint.

## 2. Canonical role bindings

Connectivity rewrite schema version `1.2` adds validated grammar role
bindings. A shared program may use canonical roles such as
`leaving_source`, `joining_partner`, and `site`, while each grammar binds those
roles to its existing assignment names.

Execution retains the original grammar role in compatibility `BondChange`
labels. Role binding therefore removes duplicated chemistry programs without
changing public edit provenance.

## 3. Migrated rewrite families

The versioned definitions now implement:

- release and connect for neutral X-H and aromatic C-H partners;
- release and connect for anionic partners with charge normalization;
- two-anchor transfer coupling;
- explicit A-B and implicit A-H addition to alkenes and alkynes;
- beta elimination; and
- carbonyl reduction, alcohol oxidation, alkene hydrogenation, partial alkyne
  hydrogenation, and complete alkyne hydrogenation.

The same bounded instruction vocabulary remains in force: localized bond
state changes, schema-level hydrogen changes, observed formal-charge changes,
endpoint permutations, product seeds, and explicitly authorized projection.

## 4. Grammar-scoped authority

Each compiled rewrite exposes `authoritative_grammar_ids`.
`enumerate_operator_outcomes` uses a connectivity rewrite only when the
specific grammar ID is listed. A registered shadow rewrite is not
automatically production-authoritative.

Twenty-eight grammars pass the current authority gate:

- C-N, C-O, and C-S substitution in sp2 and sp3 contexts, including neutral
  acyl-sulfur and anionic sulfur partners;
- amide, ester, thioester, sulfonamide, sulfonate, Friedel-Crafts acylation,
  Sonogashira, activated-carbon substitution, Suzuki, and other-metal transfer
  coupling;
- explicit and X-H addition to alkenes and alkynes;
- beta-halo elimination; and
- the five simple bond-order-change grammars.

The generic executor now also represents a chemistry construction failure as
an `OperatorOutcome` with a null product when the edit program was resolved
but sanitization or projection failed. If an authoritative rewrite produces
no outcome at all, production uses an explicit legacy fallback so malformed
or not-yet-representable execution does not silently erase a candidate.

## 5. Deliberate legacy boundary

The following remain legacy-backed:

- `sp2_c_aromatic_ch_substitution`: its rewrite is registered for shadow
  evaluation, but two same-component candidate projections in the regression
  corpus do not yet have exact failed-outcome parity;
- `chan_lam_heteroatom_coupling`: its compatibility broken-bond role and
  transfer-link projection need one unambiguous normalized contract;
- `carbonyl_amine_reductive_coupling`: carbonyl-oxygen extrusion and paired
  hydrogen deltas remain a bespoke projection;
- `terminal_alkene_heck_coupling`: terminal-endpoint selection and alkene
  product topology remain bespoke; and
- multi-event center replacements: simultaneous projection and atom-index
  remapping remain in `apply_operator_sequence`.

No legacy branch is removed while any authoritative grammar still uses the
construction-absence fallback or a non-migrated grammar shares that branch.

## 6. Parity gates

The migration corpus runs every authoritative candidate assignment, not only
the ultimately selected product match. Exact parity requires:

- outcome IDs and ordering;
- predicted canonical isomeric product;
- ordered compatibility bond changes;
- warnings;
- ambiguity count; and
- stereochemical reconstruction.

The corpus includes neutral and anionic substitution, activated acyl and
sulfonyl centers, transfer coupling, explicit and implicit addition, both
alkyne reduction depths, ambiguous beta elimination, spectators, and
multi-site candidate enumeration. Separate tests prove that successful
authoritative execution does not enter the legacy function, that an absent
rewrite outcome uses the declared fallback, and that a shadow-only grammar
does not enter the production rewrite.

Validation commands:

```powershell
pytest -q tests/reactive_taxonomy/test_reaction_site_interfaces.py
pytest -q tests/reactive_taxonomy/test_connectivity_rewrite_executor.py
pytest -q tests/reactive_taxonomy
pytest -q
```

At implementation handoff, the focused Phase 4 gates pass 83 tests and the
complete repository suite passes 576 tests.

## 7. Handoff to Phase 5

Phase 5 should derive edit shapes from the normalized connectivity edit graph
rather than treating a grammar-declared archetype as authority. Before
removing more legacy code, it should also:

1. make projection-failure reasons typed;
2. close direct aromatic C-H and Chan-Lam parity gaps;
3. represent Heck and reductive-amination projections declaratively;
4. add a composite rewrite contract for multi-event transformations; and
5. remove legacy branches only when no production or fallback path references
   them.

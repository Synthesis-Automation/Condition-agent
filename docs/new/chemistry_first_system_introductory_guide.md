# Chemistry-First Reaction-System Introductory Guide

**Status:** Introductory guide to the current implementation  
**Audience:** Chemists, cheminformatics developers, and new contributors  
**Reviewed against the repository:** 2026-08-25  
**Scope:** Reaction observation, condition recommendation, forward synthesis,
and retrosynthesis

This guide introduces the vocabulary, data flow, mathematical basis, and
practical limits of the repository's chemistry-first reaction system. It is a
conceptual entry point, not the primary architecture authority. Detailed design
and implementation status remain in:

- [`type_agnostic_reaction_recommendation_implementation.md`](type_agnostic_reaction_recommendation_implementation.md)
- [`type_agnostic_forward_synthesis_design_and_status.md`](type_agnostic_forward_synthesis_design_and_status.md)
- [`type_agnostic_retrosynthesis_design_and_status.md`](type_agnostic_retrosynthesis_design_and_status.md)

## 1. The short version

Reactions often share the same transformation core, while differing mainly in substituents, unaffected functional groups, steric environment, and electronic context.

The central mental model is:

```text
molecule = labeled graph
reaction = differences between molecular graphs
edits = individual graph differences
reaction core = observed region containing the differences
operator = reusable graph-rewrite identity derived from observations
template = one executable realization of an operator
Rewrite = applying that rule to a molecule
round trip = check that the executable rule reproduces known chemistry
recommendation = rank compatible, evidence-backed alternatives
```

The molecular graph is the source of truth. Reaction names such as
*Suzuki-Miyaura* are useful annotations when supported, but they do not define
graph identity or override contradictory structural evidence.

All three application paths use the same broad pattern:

```text
hard structural and chemistry gates
                ↓
generation or precedent retrieval
                ↓
validation against graph evidence
                ↓
interpretable evidence-based ranking
```

The primary paths are deterministic. They combine graph matching, set and
fingerprint similarity, empirical support counts, and weighted ranking. They
are not one end-to-end neural model, and their scores are not probabilities of
experimental success.

## 2. System boundaries

The standalone packages have distinct responsibilities:

| Package | Responsibility |
| --- | --- |
| `reactive_taxonomy` | Molecular parsing, atom correspondence, edits, reaction cores, signatures, operators, and structural chemistry |
| `condition_registry` | Canonical identities, aliases, roles, and recipes for condition substances |
| `condition_recommender` | Dataset admission, precedent retrieval, compatibility, recipe aggregation, ranking, and explanations |
| `core_retrosynthesis` | Generic operator compilation, reverse application, forward validation, disconnection ranking, and route planning |
| `forward_synthesis` | Independent forward admission, operator application, product generation, and product ranking |

Application code composes these packages but should not own core chemistry
rules.

## 3. Molecular graphs

A molecule is represented as a labeled graph:

\[
G=(V,E,\ell_V,\ell_E)
\]

where:

- \(V\) is the set of atoms;
- \(E\) is the set of bonds;
- \(\ell_V\) contains atom properties such as element, formal charge,
  aromaticity, isotope, hydrogens, hybridization, radical state, and stereo;
- \(\ell_E\) contains bond properties such as order, aromaticity, and stereo.

For biphenyl, the two phenyl-ring subgraphs are joined by one C-C edge:

```text
Ph-Ph
```

One SMILES representation is:

```text
c1ccccc1-c2ccccc2
```

SMILES is a serialization of the graph, not the graph identity itself.
Different valid SMILES can describe the same molecule, so the system
canonicalizes structures before identity comparisons.

## 4. Atom correspondence and mapping

Atom correspondence determines which precursor atom became which product atom.
For a simplified Suzuki coupling:

```text
[C:1]-Br + [C:2]-B(OH)2 → [C:1]-[C:2]
```

the map numbers `1` and `2` are correspondence labels. They are not chemical
properties. They allow the system to determine:

- which product carbon came from each precursor;
- which atoms are retained in the product;
- which fragments depart from the recorded main product;
- which bond or atom states changed.

A supplied map is accepted only when it passes validation. When mapping must be
materialized, the mapping provider may use learned machinery, but the resulting
observation, identity, graph checks, and primary ranking remain explicit and
deterministic.

## 5. Graph edits

An edit is one normalized difference between precursor and product graphs. The
main edit types are:

| Edit | Meaning |
| --- | --- |
| Formed bond | No corresponding bond before; a bond exists after |
| Broken bond | A corresponding bond exists before; none exists after |
| Order change | The corresponding bond changes order |
| Hydrogen change | A schema-level hydrogen is gained or lost |
| Atom-state change | Charge, radical, isotope, aromaticity, hybridization, or stereo changes |

With an atom map \(M\), the bond-edit intuition is:

\[
\Delta_{\mathrm{formed}}=E_P\setminus M(E_R)
\]

\[
\Delta_{\mathrm{broken}}=M(E_R)\setminus E_P
\]

and the complete edit set is approximately:

\[
\Delta=
\Delta_{\mathrm{formed}}\cup
\Delta_{\mathrm{broken}}\cup
\Delta_{\mathrm{order}}\cup
\Delta_{\mathrm{state}}.
\]

For the Suzuki example, the central product-retained edit is formation of a
single C(sp2)-C(sp2) bond. Precursor-only C-Br and C-B handles are important for
constructing an executable realization, even though Br and B are absent from
the recorded main product.

Edits describe what changed. They do not claim a mechanism, transition state,
rate, yield, or optimal conditions.

## 6. Reaction observation, core, and signature

These contracts describe related but different views of a reaction.

### 6.1 Reaction observation

The observation owns parsed components, atom correspondence, normalized edits,
topology, completeness, spectators, evidence, warnings, and uncertainty. It is
the structure-derived record of what the system can justify.

### 6.2 Reaction core

The `ReactionCoreProjection` is a minimized, scaffold-abstracted view of the
observed change region. It contains:

- participating atom states before and after;
- bond and atom-state changes;
- connected edit events and their relationships;
- attachment ports to omitted subgraphs;
- retained, departing, appearing, or unresolved remote subgraphs;
- quality status, confidence, evidence, and warnings.

The reaction core says:

> This is the part of the observed molecular graphs that changed.

It is template-free evidence. It is not itself an instruction for generating
new molecules.

### 6.3 Reaction signature

The generic `ReactionSignature` is a normalized, versioned identity and feature
representation derived from the observation. It supports indexing,
compatibility checks, comparison, and explanation. Its identity is based on
chemistry and definition versions, not display labels or source reaction names.

## 7. Operators, realizations, and templates

An operator is a generalized graph transformation. It answers:

> If a molecule contains the required structural pattern, which graph change
> can be executed in the forward or reverse direction?

The repository separates several identities:

| Namespace | Meaning |
| --- | --- |
| `OP1` | Handle-independent normalized graph transformation |
| `SITE1` | Target-specific disconnection site |
| `SYN1` | Normalized synthon concept |
| `STRAT1` | `OP1 + SITE1 + SYN1`, one target-specific strategy |
| `REAL2` | Concrete precursor-handle realization |
| `GRT3` | Stored executable template at an abstraction level |
| `COMP2` | Evidence group for related handle completions |

The relationship is:

```text
one OP1 graph transformation
  ├─ one or more target sites
  ├─ one or more synthon concepts
  └─ one or more REAL2 handle realizations
       └─ one or more GRT3 executable templates
            ├─ L2
            ├─ L1
            └─ L0
```

The same handle-independent C(sp2)-C(sp2) coupling operator can have several
realizations:

```text
Suzuki-like:  C-B  + C-halide
Stille-like:  C-Sn + C-halide
Negishi-like: C-Zn + C-halide
Kumada-like:  C-Mg + C-halide
```

Consequently, *Suzuki* is not necessarily a unique `OP1`. It is commonly a
particular C-B/C-halide realization of a more general C-C coupling operator.

## 8. SMARTS and graph rewrites

SMILES describes a concrete molecular graph. SMARTS describes a graph query
pattern. Reaction SMARTS connects precursor and product patterns:

```text
precursor SMARTS >> product SMARTS
```

The retrosynthetic direction is stored explicitly as:

```text
product SMARTS >> precursor SMARTS
```

Conceptually, applying a graph rewrite means:

1. find a graph embedding of the template pattern in the input molecule;
2. preserve the explicit atom correspondence;
3. remove, form, or change the designated graph elements;
4. construct the output-side components;
5. sanitize and canonicalize the generated molecular graphs.

If \(L\) is the matched template side and \(R\) is the generated side, the
rewrite can be pictured as:

\[
L\rightarrow R,
\qquad
\phi:L\hookrightarrow G,
\qquad
G\mapsto G'.
\]

The implementation uses RDKit and RDChiral reaction execution; it should not be
read as a full mechanistic or quantum-chemical simulation.

## 9. Round trips and query-time validation

### 9.1 Source round trip

During operator extraction, a template must reproduce its own source reaction:

```text
known source precursors
          ↓ observed forward reaction
known source product
          ↓ extracted reverse template
recovered precursors
          ↓ canonical comparison
must agree with the source precursors
```

The comparison follows the admitted exact or explicitly relaxed stereo policy.
A failed round trip rejects the template. Passing establishes fidelity to the
source example, not universal experimental applicability.

### 9.2 Retrosynthesis forward validation

After reverse application proposes precursors, the proposed forward reaction is
reanalyzed:

```text
proposed precursors → target
```

The candidate is retained only when its reconstructed normalized operator
signature agrees with the operator that generated it:

\[
\operatorname{signature}
(\text{proposed precursors}\rightarrow\text{target})
=
\operatorname{signature}(\text{operator}).
\]

Retained retrosynthesis candidates therefore have
`forward_validation_status="verified_signature"`.

### 9.3 Forward-synthesis reverse validation

Forward-generated products are independently checked in the other direction:
the reverse operator must recover the participating precursor set. This avoids
accepting a product merely because an unconstrained reaction rule emitted it.

## 10. Synthons, handles, sites, and strategies

A **synthon** is an idealized fragment concept produced by disconnecting a
target before choosing isolable reagents. For biphenyl:

```text
Ph-Ph → Ph• + •Ph
```

A **handle** makes a synthon concretely realizable, for example halide,
boronic acid, boronate, organotin, organozinc, Grignard, or an activated acyl
group.

A **disconnection site** identifies the particular target bond or atom state
selected for retrosynthesis. One operator can match several sites in a complex
target, so operator recall and correct-site recall are different measurements.

A **strategy** groups candidates with the same operator, target site, and
synthon concept while retaining alternate handle realizations.

## 11. L2, L1, and L0 abstraction

Operators are compiled with different amounts of structural context:

| Level | Interpretation |
| --- | --- |
| L2 | Context-rich local environment; most specific |
| L1 | Reaction core and directly relevant handle context |
| L0 | Most generalized fallback |

Retrosynthesis searches in the order:

```text
L2 → L1 → L0
```

Specific templates improve selectivity but transfer less broadly. General
templates improve coverage but can match more sites and generate more competing
candidates. The returned result records the level used so fallback is explicit.

## 12. What operator extraction covers

The GUI's full-scale operator workflow uses the current default configuration:

- all completed rows in the selected combined corpus are attempted;
- `admission_mode="data_driven"`;
- L0, L1, and L2 are requested;
- `max_rows_per_shard=None`;
- there is no minimum frequency threshold for retaining a successfully
  admitted operator.

For each source row, compilation requires:

1. an available source reaction;
2. valid or materialized atom correspondence;
3. a valid recomputed reaction analysis;
4. an admitted, passing recomputed reaction core;
5. verified product completeness;
6. exactly one product;
7. canonicalizable participants;
8. mapped active atoms and a generic operator signature;
9. compilable templates;
10. successful source-precursor round trip.

Every successfully admitted row contributes its templates. Templates are
deduplicated exactly, and every resulting distinct `OP1` is retained. Thus the
precise coverage statement is:

> A completed build is exhaustive over successfully admitted observations in
> the selected corpus, not over every reaction a chemist might regard as valid
> and not over every chemically possible transformation.

The current Compact artifact illustrates the distinction:

| Build quantity | Count |
| --- | ---: |
| Source rows | 117,232 |
| Accepted observations | 76,542 |
| Acceptance fraction | 65.29% |
| Unique `OP1` operators | 3,710 |
| Handle realizations | 9,463 |
| Executable templates | 31,977 |

Major recorded losses include 20,359 atom-mapping failures, 16,244 unverified
cores, 2,377 source-round-trip failures, and 726 multi-product reactions. See
the current
[`build_report.json`](../../results/operator_retrosynthesis_poc/full_scale_v3/compact/build_report.json).

The compiler creates one normalized operator for the complete admitted edit set
of a source reaction. It does not enumerate every counterfactual subset of a
multi-event reaction.

## 13. Does a similar target guarantee retrieval?

No. It improves the chance of retrieval but is not sufficient. Retrosynthesis
uses the following stages:

```text
necessary product-feature retrieval
          ↓
product SMARTS substructure match
          ↓
bounded reverse template application
          ↓
precursor generation
          ↓
bounded forward validation
          ↓
ranking and top-k selection
```

Global structural similarity is primarily ranking evidence after structural
applicability. A globally similar molecule can differ exactly at the reactive
center and fail SMARTS matching. A globally dissimilar molecule can still match
when it contains the required local reactive structure.

The interactive web path deliberately bounds work. For a request with
`top_k=10`, it currently applies at most 100 templates and forward-validates at
most 30 generated candidates. The bounds grow with `top_k` but are capped at
300 template applications and 100 validations. These limits keep interaction
responsive, but mean that existence in the library is not identical to recall
in one bounded request.

Several different claims must therefore be separated:

| Claim | What is required |
| --- | --- |
| Dataset operator exists | At least one supporting source reaction passed admission |
| Operator applies to target | A stored product SMARTS matches and executes |
| Correct site is found | The intended target site survives generation, validation, and ranking |
| Exact precursors are found | The required handle realization is present and generated |
| Candidate is displayed | It survives the template, validation, and `top_k` budgets |

For the exact source product of an admitted template, source round trip proves
that the template can regenerate its source precursors. A bounded interactive
search can still omit it if too many higher-ranked applicable templates consume
the application or validation budget.

## 14. Worked Suzuki example

Suppose the dataset contains many admitted reactions of the form:

```text
Ar-Cl + Ar'-B(OH)2 → Ar-Ar'
```

Compilation observes the C(sp2)-C(sp2) formation, constructs a shared C-C
coupling `OP1`, and retains C-B/C-Cl handle realizations and L2/L1/L0 executable
templates.

For the target biphenyl:

```text
c1ccccc1-c2ccccc2
```

the intended sequence is:

```text
match a product-side aryl-aryl bond
          ↓
apply a reverse C-C coupling template
          ↓
propose Clc1ccccc1.OB(O)c1ccccc1
          ↓
reanalyze the proposed forward reaction
          ↓
verify the C-C coupling operator signature
```

In a direct check against the current Compact library on 2026-08-25, using the
web-equivalent `top_k=10` budgets of 100 template applications and 30 forward
validations, the search returned:

```text
rank 1: Stille-like C-Sn/C-Br realization
rank 2: Suzuki-like C-B/C-Cl realization
```

The Suzuki-like candidate was:

```text
Clc1ccccc1.OB(O)c1ccccc1
```

and passed `verified_signature` validation. Both candidates shared the same
handle-independent C-C coupling `OP1`. This demonstrates three points:

1. dataset-derived Suzuki-like chemistry can transfer to a compatible target;
2. the result may be annotated generically as `c_c_coupling`, because a named
   family is not the routing identity;
3. finding the Suzuki realization does not guarantee ranking it first. A
   `top_k=1` query would have hidden it in this particular run.

More numerous and more diverse independent Suzuki precedents improve local
environment coverage, handle coverage, support, and completion priors. They do
not turn graph applicability or ranking into a guarantee.

## 15. Common mathematical basis

The primary mathematical tools are:

1. labeled-graph representation and substructure matching;
2. normalized graph-edit identities;
3. set and multiset similarities;
4. molecular fingerprints and Tanimoto similarity;
5. hard Boolean compatibility constraints;
6. logarithmically saturated evidence counts;
7. weighted linear scores with explicit missing-data handling;
8. smoothed empirical completion frequencies;
9. deterministic score bands, grouping, and tie-breaking.

In general, candidate selection has the form:

\[
\operatorname{rank}
\left(
\{c\mid H_1(c)=H_2(c)=\cdots=1\}
\right),
\]

where each \(H_i(c)\in\{0,1\}\) is a mandatory chemistry or evidence gate.
Similarity cannot compensate for a failed hard gate.

### 15.1 Jaccard similarity

For feature sets \(A\) and \(B\):

\[
J(A,B)=\frac{|A\cap B|}{|A\cup B|}.
\]

For multiplicity-preserving event collections:

\[
J_m(A,B)=
\frac{\sum_x\min(n_A(x),n_B(x))}
{\sum_x\max(n_A(x),n_B(x))}.
\]

### 15.2 Morgan-Tanimoto similarity

Forward and retrosynthetic molecular comparisons use 2,048-bit, radius-2
Morgan fingerprints. For bit sets \(A\) and \(B\):

\[
T(A,B)=
\frac{|A\cap B|}
{|A|+|B|-|A\cap B|}.
\]

The score used for a query is generally the maximum Tanimoto similarity to a
stored precedent for that template.

### 15.3 Saturated support

Evidence counts use a diminishing-return transformation:

\[
U(n;k)=
\min\left(
1,
\frac{\log(1+n)}{\log(1+k)}
\right),
\]

where \(n\) is support and \(k\) is the configured saturation count.

### 15.4 Missing evidence

When a ranking component is genuinely unavailable, available weights are
renormalized:

\[
S=\frac{\sum_{i\in A}w_i x_i}{\sum_{i\in A}w_i},
\]

where \(A\) is the set of available components. Missing evidence is therefore
not silently treated as a perfect match or as zero when the policy explicitly
permits renormalization.

## 16. Mathematics of condition recommendation

Condition recommendation starts from a known or proposed forward reaction. It
does not generate molecular products or precursors.

### 16.1 Retrieval ladder

Precedents are retrieved through increasingly general structural tiers:

```text
exact signature
→ handle signature
→ high-confidence named family, when available
→ transformation signature
→ environment neighbors
→ bond-edit signature
→ reaction-core exact/local/context tiers
→ edit-graph neighbors
```

Net-edit and reactive-center compatibility are enforced before similarity.
The default ladder seeks at least two independent evidence units before
accepting a tier, while retaining an explicitly labelled limited-support result
when appropriate.

### 16.2 Reaction-signature similarity

Condition precedent similarity is a weighted sum:

\[
S_{\mathrm{reaction}}=\sum_i w_i s_i.
\]

Current components and weights are:

| Component | Weight |
| --- | ---: |
| Edit topology | 0.20 |
| Reaction events | 0.12 |
| Reaction topology | 0.10 |
| Reactive handles | 0.16 |
| Attachment contexts | 0.12 |
| Local environment | 0.11 |
| Spectators | 0.06 |
| Transformation class | 0.07 |
| Named family | 0.06 |

Named-family identity therefore contributes only 6% and cannot override an
incompatible edit graph.

### 16.3 Recipe compatibility

A hard conflict gives:

\[
C_{\mathrm{recipe}}=0
\]

and excludes the precedent. Otherwise:

\[
C_{\mathrm{recipe}}=
\max\left(0,1-\sum_j p_j\right),
\]

where \(p_j\) are matched soft penalties. Related penalties can share a group,
in which case only the strongest matched penalty in that group is used.

### 16.4 Similarity-weighted yield shrinkage

For precedent \(i\), define:

\[
q_i=\max(S_i,0.05).
\]

The current recipe-level expected-yield evidence is:

\[
\hat y=
\frac{\sum_i q_i y_i+\lambda\mu_{\mathrm{pool}}}
{\sum_i q_i+\lambda},
\qquad \lambda=3,
\]

where \(y_i\) is reported yield and \(\mu_{\mathrm{pool}}\) is the mean yield
of independent yield-bearing precedents in the retrieved pool. This is an
empirical shrinkage estimate, not a universally calibrated yield prediction.

### 16.5 Default condition score

After precedents are grouped by canonical resolved recipe, the current balanced
score is:

\[
\begin{aligned}
S_{\mathrm{condition}}={}&
0.45S_{\mathrm{reaction}}
+0.15Y
+0.15U_{\mathrm{independent}}\\
&+0.08U_{\mathrm{reaction\ breadth}}
+0.05U_{\mathrm{dataset\ diversity}}\\
&+0.08C_{\mathrm{recipe}}
+0.04C_{\mathrm{certainty}}.
\end{aligned}
\]

Chemist-selectable profiles can instead emphasize reactant category,
functional-group tolerance, independent evidence, expected yield, or procedure
completeness. Hard chemistry and compatibility gates remain in force.

Definitions:

- [`generic_similarity.v1.json`](../../condition_recommender/definitions/generic_similarity.v1.json)
- [`generic_ranking.v1.json`](../../condition_recommender/definitions/generic_ranking.v1.json)
- [`chemist_ranking_profiles.v1.json`](../../condition_recommender/definitions/chemist_ranking_profiles.v1.json)
- [`compatibility.v1.json`](../../condition_recommender/definitions/compatibility.v1.json)

## 17. Mathematics of forward synthesis

Forward synthesis starts with precursor graphs and enumerates structurally
possible products from independently admitted bidirectional operators.

### 17.1 Retrieval and hard admission

The conservative index checks precursor component count and required atomic
numbers. Forward SMARTS application then assigns components to operator
reactant patterns. A candidate must pass:

1. product sanitization;
2. reverse recovery of participating precursors;
3. valid reaction analysis and materialized core/signature;
4. operator-edit agreement;
5. no hard supplied-recipe or coarse-profile conflict.

### 17.2 Forward ranking

For a valid product candidate:

\[
\begin{aligned}
S_{\mathrm{forward}}={}&
0.35T_{\mathrm{precursors}}
+0.20Q_{\mathrm{specificity}}
+0.20U_{\mathrm{support}}\\
&+0.15A_{\mathrm{edits}}
+0.10C_{\mathrm{recipe}}.
\end{aligned}
\]

Here:

\[
Q_{\mathrm{specificity}}=
\min\left(
1,
\frac{N_{\mathrm{template\ atoms}}}
{N_{\mathrm{starting-material\ heavy\ atoms}}}
\right),
\]

\(A_{\mathrm{edits}}=1\) for an admitted candidate, and
\(C_{\mathrm{recipe}}\) is optional. Missing optional evidence causes explicit
weight renormalization.

The final score can receive a 0.05 virtual-copy penalty and a bounded condition
profile adjustment:

\[
S'=\operatorname{clip}_{[0,1]}
(S_{\mathrm{forward}}-P_{\mathrm{virtual}}+A_{\mathrm{profile}}).
\]

See
[`forward_ranking.v3.json`](../../forward_synthesis/definitions/forward_ranking.v3.json).

## 18. Mathematics of retrosynthesis

Retrosynthesis starts with a target product graph and applies operators in
reverse.

### 18.1 Applicability

A product template is applicable only when its pattern embeds in the target:

\[
\phi:G_{\mathrm{template\ product}}
\hookrightarrow G_{\mathrm{target}}.
\]

Necessary product tokens prune the library, but product SMARTS matching and
graph application remain authoritative.

### 18.2 Preliminary score

For each generated precursor set:

\[
P=
0.50T_{\mathrm{product}}
+0.28T_{\mathrm{precursors}}
+0.12U_{\mathrm{support}}
+0.10Q_{\mathrm{specificity}}.
\]

When contextual evidence is enabled:

\[
S_{\mathrm{retro}}=0.85P+0.15C_{\mathrm{context}}.
\]

The context term is:

\[
C_{\mathrm{context}}=
0.30J_{\mathrm{spectators}}
+0.35J_{\mathrm{substituent\ features}}
+0.35S_{\mathrm{site\ profiles}},
\]

where the site-profile score combines categorical agreement with normalized
steric-accessibility and electronic-activation distances.

Only candidates that pass forward operator-signature validation reach the
final ranking stages.

### 18.3 Hierarchical site, synthon, and realization scores

Validated candidates are ranked as sites, then synthons, then concrete
realizations:

\[
S_{\mathrm{site}}=
0.65S_{\mathrm{structural}}
+0.10U_{\mathrm{support}}
+0.25S_{\mathrm{complexity\ reduction}},
\]

\[
S_{\mathrm{synthon}}=
0.80S_{\mathrm{structural}}
+0.20U_{\mathrm{support}},
\]

\[
S_{\mathrm{realization}}=
0.65S_{\mathrm{structural}}
+0.15U_{\mathrm{support}}
+0.20P_{\mathrm{completion}}.
\]

### 18.4 Completion priors

Handle-completion evidence uses Laplace smoothing:

\[
P(g\mid c)=
\frac{n_g+\alpha}{N+\alpha K},
\qquad \alpha=1,
\]

where \(n_g\) is independent support for completion group \(g\), \(N\) is
total support in the context, and \(K\) is the number of alternatives. Evidence
backs off through:

```text
operator + synthon → operator → global
```

This term helps order Suzuki-like, Stille-like, Negishi-like, and other
realizations of a shared disconnection. It is an empirical completion frequency,
not the probability that the laboratory reaction will succeed.

### 18.5 Score bands and diversity

Structural scores are partitioned into bands of width 0.05:

\[
b(c)=
\left\lfloor
\frac{1-S_{\mathrm{structural}}}{0.05}
\right\rfloor.
\]

Hierarchy, completion evidence, precursor realism, condition evidence, and
diversity are guarded so they do not freely promote a substantially worse
structural candidate. Diversity groups candidates by operator, target site, and
synthon identity instead of returning only minor handle variants.

Definitions:

- [`retrosynthesis_ranking.v3.json`](../../core_retrosynthesis/definitions/retrosynthesis_ranking.v3.json)
- [`hierarchical_ranking.v2.json`](../../core_retrosynthesis/definitions/hierarchical_ranking.v2.json)

## 19. How to interpret the outputs

A retained retrosynthetic or forward candidate establishes that:

- a stored structural operator was applicable;
- graph execution generated sanitized molecules;
- required round-trip and signature checks passed;
- the result has traceable operator, template, and precedent identities;
- its ranking components and cautions can be inspected.

A condition recommendation establishes that:

- compatible precedents were retrieved through an explicit tier;
- raw conditions were resolved into canonical recipes where possible;
- hard conflicts were excluded;
- independent support and similarity were aggregated transparently.

None of these outputs alone establishes:

- experimental success;
- major-product selectivity;
- reaction mechanism;
- calibrated yield;
- optimal catalyst, solvent, temperature, or time;
- precursor availability or stability;
- scale-up safety.

Forward and retrosynthesis scores are explicitly deterministic priorities, not
probabilities. Condition similarity uses chemistry-prior weights, while the
default recipe ranking has sample-level calibration but is not a universal
statistical model. Human review and experimental validation remain necessary.

## 20. Compact glossary

| Term | Short definition |
| --- | --- |
| Molecular graph | Atoms as labeled vertices and bonds as labeled edges |
| Atom mapping | Correspondence between precursor and product atoms |
| Edit | One normalized bond or atom-state difference |
| Reaction observation | Complete structure-derived evidence and uncertainty |
| Reaction core | Minimized observed change region and its attachments |
| Reaction signature | Versioned normalized identity/features for indexing and comparison |
| Operator | Reusable handle-independent graph transformation |
| Realization | Concrete precursor-handle implementation of an operator |
| Template | Executable forward/reverse SMARTS at one abstraction level |
| Graph rewrite | Pattern match followed by explicit graph changes |
| Round trip | Reapplication check against known source chemistry |
| Synthon | Idealized fragment concept before selecting real handles |
| Handle | Functional group that realizes a synthon experimentally |
| Disconnection site | Specific target bond or state selected for reversal |
| Strategy | Operator, site, and synthon identity grouped together |
| Precedent | Dataset reaction supporting an operator or recipe |
| Independent support | Support counted by distinct evidence units/references |
| Compatibility | Hard conflicts plus explicit soft penalties |
| Retrieval level | Specificity tier at which evidence was found |
| Ranking score | Deterministic ordering priority, not success probability |

## 21. Further reading

- [`reaction_featurization_workflow.md`](reaction_featurization_workflow.md)
- [`reaction_condition_recommendation_design_for_chemists.md`](reaction_condition_recommendation_design_for_chemists.md)
- [`route_core_projection_design_and_status.md`](route_core_projection_design_and_status.md)
- [`core_retrosynthesis/README.md`](../../core_retrosynthesis/README.md)
- [`forward_synthesis/README.md`](../../forward_synthesis/README.md)

# Synthesis-route search workflow: fused aminoketone case study

**Date:** 2026-10-08  
**Purpose:** Document the method behind the two proposed routes and SVG so that it can be compared with a chemistry agent using dedicated tools.  
**Target:** `CC(C)C[C@@H]1CN2[C@@H](c3c(C2)cccc3)CC1=O`  
**Outcome:** Two literature-informed proposals; no experimental synthesis of this exact stereochemical target was verified.

## 1. Scope and provenance

This is a stage-level account of the delivered answer, its source evidence, structural checks, and limitations—not a complete historical tool-call transcript. The final SVG and companion JSON identify which routes were actually delivered. They take precedence over any intermediate draft.

Three kinds of information are separated throughout:

- **Delivered result:** Route structures, references, and qualifications present in the original answer and files.
- **Documentation audit:** Source and structure checks performed while preparing this workflow document.
- **Replay or improvement:** A reproducible implementation or comparison requirement, not a claim that the original run executed that exact command.

The original exhaustive search-query list, original final drawing script, per-stage timings, token usage, and cost are not available in the delivered files. Query examples below are explicitly a replay recipe, not a reconstructed execution log. The chat displayed approximately 24 min 55 s for the earlier response; that is a user-interface elapsed time, not instrumented tool or compute time.

**Method classification:** This was already a tool-assisted workflow: LLM-guided chemical interpretation and literature adaptation, general web retrieval, and local cheminformatics. The useful comparison is therefore **general-purpose tools versus dedicated chemistry retrieval/planning/validation tools**, not “LLM without tools versus agent with tools.”

## 2. Workflow at a glance

```text
Input SMILES
    ↓
Minimal input repair → parse → lock target connectivity and stereochemistry
    ↓
Exact-target / identifier retrieval, where possible
    ↓  if an exact preparation is not verified
Scaffold and intermediate retrieval + related synthesis literature
    ↓
Read experimental sections → separate exact precursor evidence from analogies
    ↓
Construct explicit candidate intermediates and reaction steps
    ↓
Check connectivity, ring size, functional groups, and stereochemical completeness
    ↓
Rank proposals qualitatively; identify the decisive unverified transformation
    ↓
Validate molecular representations with RDKit
    ↓
Compose route SVG from validated structures and evidence-labelled arrows
    ↓
Cross-check JSON ↔ SVG ↔ written answer; disclose unresolved chemistry
```

The chemical hypotheses and the computational checks serve different purposes. A parser can confirm the intended molecule is represented; it cannot confirm that a proposed reaction will produce it.

## 3. Division of work between model and tools

| Task | Method behind the delivered result | What it does not establish |
|---|---|---|
| Interpret the target | Chemical interpretation, supported by RDKit graph and stereochemistry checks | A published identity or a known preparation |
| Find precedents | General web/patent retrieval; retained primary references R1–R3 | Exhaustive coverage of structure-indexed literature |
| Propose routes | LLM-guided adaptation of documented chemistry | A retrosynthesis-engine prediction or an experimentally validated route |
| Represent intermediates | Explicit SMILES followed by RDKit parsing and descriptor checks | Product selectivity, reaction yield, or intermediate stability |
| Rank routes | Qualitative comparison of operations, evidence, and unresolved steps | A calibrated success probability or optimized route cost |
| Present routes | Vector molecular drawings with SVG layout, annotations, and embedded JSON | Independent chemical validation from the appearance of the scheme |
| Audit this documentation | Read the final files; reopen primary sources; inspect patent-page images; rerun structural checks | A new experimental validation of either route |

No results from a proprietary structure-search database, dedicated retrosynthesis planner, forward-reaction predictor, reaction atom-mapper, or supplier inventory API are established in the retained deliverables. Such capabilities must not be credited to this baseline without an actual result record.

## 4. Stage 1 — Normalize and lock the target

### Input correction

The user supplied:

```text
CC(C)C[C@@H]1CN2[C@@H]\(c3c(C2)cccc3)CC1=O
```

The backslash immediately before the opening parenthesis was treated as a formatting escape. Removing that single character gives the target stated above. Both `[C@@H]` annotations were preserved.

This is a **case-specific repair**, not a rule to delete backslashes generally: SMILES can use directional bond notation to encode alkene stereochemistry. Store the raw input, the accepted normalized form, and an edit description separately.

### Recorded and independently rechecked identity

| Field | Result |
|---|---|
| Canonical isomeric SMILES | `CC(C)C[C@@H]1CN2Cc3ccccc3[C@H]2CC1=O` |
| Molecular formula | C₁₆H₂₁NO |
| Calculated average molecular weight | 243.35 g/mol |
| InChIKey | `DSHSQHJCJXXGIL-UKRRQHHQSA-N` |
| Specified carbon stereocenters | R,R |
| Ring sizes returned for this graph | 5, 6, 6 |
| Chemical interpretation used | Isoindoline fused to a six-membered aminoketone ring; isobutyl substituent |

These are calculations or interpretations of the supplied structure, not evidence of experimental synthesis. The atom indices in the JSON are software indices, not IUPAC locants. A change from `@@` to `@` in an equivalent canonical SMILES is not, by itself, a stereochemical inversion: atom traversal changes must be considered.

**Critical target check:** Do not silently replace this 6–5–6 framework with the related 6–6–6 framework used in the main analogue references. The ring-size difference remains a transferability gap.

## 5. Stage 2 — Search iteratively, but preserve search provenance

The final answer is based on a documented precursor and related ring-construction chemistry. An exact-target experimental route was not verified. That statement does **not** establish that the molecule is novel or absent from the literature.

For reproduction, use the following query ladder. These are representative queries, not certified historical queries:

| Search level | Example query | Purpose |
|---|---|---|
| Exact identifier | `"DSHSQHJCJXXGIL-UKRRQHHQSA-N"` | Find an indexed exact structure record |
| Connectivity-level identifier | `"DSHSQHJCJXXGIL"` | Look beyond a particular stereochemical identifier; verify any hit structurally |
| Exact structure string | `"CC(C)C[C@@H]1CN2Cc3ccccc3[C@H]2CC1=O"` | Search literal indexed SMILES; a negative result is weak evidence |
| Scaffold and motif | `isoindoline fused piperidinone isobutyl synthesis` | Retrieve a relevant chemical family without claiming an exact name |
| Advanced precursor | `"isoindol-1-yl" "acetic acid methyl ester"` | Locate preparations of the isoindoline-acetate building block |
| Reaction-family analogue | `tetrabenazine malonate Mannich Dieckmann synthesis` | Find precedent for Route A's construction strategy |
| Alternative ring closure | `tetrabenazine aldehyde alkenyl iodide amino cyclization` | Find precedent for Route B's assembly strategy |

Once a useful document is found, follow its preparation numbers, compound references, and experimental cross-references. Do not replace a source read with the search-result snippet.

For a benchmark, log every actual query, timestamp, tool, returned document identifier, retrieved range, and access failure. Separate **no relevant hit**, **not indexed**, **inaccessible source**, and **search not performed**.

## 6. Stage 3 — Build an evidence register

### Evidence categories

Use two independent fields: the evidence supporting a precedent and the status of the proposed target step.

| Precedent evidence | Meaning |
|---|---|
| Exact experimental | The same relevant structures and transformation are reported; stereochemistry and isolation status still require explicit checking |
| Experimental analogue | An experimentally reported transformation differs in substrate, scaffold, substitution, or stereochemistry |
| Reaction-class precedent | Support exists for a general transformation, but the substrate match is weak |
| Unvalidated proposal | The proposed step lacks an adequately matched experimental example |

A target step can have `precedent = experimental analogue` and simultaneously `target_step_status = proposed`. A reference attached to an arrow is not proof of an exact-substrate experiment.

### Primary sources retained in the answer

| ID | Source and locator | What it supports | Boundary |
|---|---|---|---|
| R1 | US20040082590A1, Preparation 1C and Preparation 3C, Step A | Preparation of the isoindoline-acetate building-block family | Not synthesis of the final target or assignment of its requested absolute configuration |
| R2 | US2830993A, Example 11, referring to operations in Example 1 | Isobutylmalonate-based Mannich/Dieckmann/decarboxylation strategy | Experimental scaffold is the 6–6–6 analogue, not this 6–5–6 target |
| R3 | US8008500B2, Examples 3, 4, 7, and 10 | Aldehyde preparation, alkenyl coupling, oxidation, and amino cyclization | Related tetrabenazine-series substrates, not the proposed isoindoline intermediates |

### A concrete extraction failure found during this audit

The Google Patents HTML for R2 associates an **n-amyl** passage with “Example 11.” The scanned original shows **Example 11 is the isobutyl case**, while the n-amyl case is **Example 12**. The relevant location is PDF page 4, printed column 7. The page image therefore confirms the isobutyl citation used in the delivered answer. [R2]

This is a useful benchmark test: preserve page/column provenance and compare the original page when an extracted example number or substituent is inconsistent. Text extraction is not the chemical ground truth.

### Avoid treating a telescoped operation as multiple isolated products

R1 Preparation 1C, Step C combines deprotection/work-up and re-protection before isolation/resolution. The SVG subdivides that operation to show free amine A1. Preparation 3C, Step A separately reports access to the free amine by deprotecting the protected building block. Thus, the scheme is a structural route representation, not proof that every drawn node was independently isolated in the same experimental passage. [R1]

Record `isolated_intermediate`, `telescoped`, and `stereochemistry_reported` explicitly rather than relying on a single “reported” label.

## 7. Stage 4 — Construct the two final candidate routes

The node identifiers below match the original JSON. They are proposed target intermediates, not the source patents' compound numbers.

### Common entry

```text
P0  N-Boc-2-bromobenzylamine
    → P1  ortho-aminomethyl cinnamate derivative
    → A1  methyl isoindolin-1-ylacetate
    → C0  its N-Boc derivative
```

R1 supplies the experimental building-block basis. The graphic's subdivisions and unspecified stereochemistry require the qualifications in Section 6.

### Route A — Mannich/Dieckmann design

```text
A1 + M (dimethyl isobutylmalonate) + formaldehyde source
    → A2  amino triester
    → A3  cyclic keto diester
    → Tmix  target connectivity, stereochemistry unspecified
    → T  requested R,R target, conditional on successful separation/assignment
```

The target-specific proposal adapts the R2 strategy. The key unverified step is **A2 → A3**. The original ranking preferred this as the shorter proposed sequence from the free-amine entry, not as a proven higher-yielding route.

The practical gates proposed in the answer were: obtain the tethered precursor, establish the intended ring closure, then test survival through hydrolysis/decarboxylation. Resolution is not worth developing before the required connectivity has been identified.

### Route B — Aldehyde/enone design

```text
C0
    → B1  aldehyde
    + V (2-iodo-4-methylpent-1-ene)
    → B2  allylic alcohol
    → B3  enone
    → Tmix  target connectivity, stereochemistry unspecified
    → T  requested R,R target, conditional on successful separation/assignment
```

This adapts R3 with a different ring-closing design. The proposed decision points include stopping reduction at the aldehyde and establishing the final intramolecular aza-Michael closure. The transient alcohol stereocenter does not solve the final stereochemical requirement.

### Stereochemistry and ranking limits

Neither candidate is an established asymmetric synthesis. `Tmix` is a representation with unspecified stereochemistry, **not evidence that all stereoisomers form**, that R,R is present in useful proportion, or that chromatography can resolve it.

The final separation arrow is a development task, not a validated reaction. The requested stereoisomer must be assigned independently; elution order is not an absolute-configuration assignment.

The ranking was qualitative. There was no calibrated route score, measured yield prediction, verified purchasing plan, or end-to-end cost calculation. Both routes share the isoindoline precursor family and unresolved stereochemical delivery; they are chemically different proposals, but not independent solutions to every risk.

## 8. Stage 5 — Validate representations without overstating validation

### Checks completed in the documentation audit

The audit used Python **3.13.5** and RDKit **2025.09.4**. These are the audit environment versions; the original versions were not recorded in the delivered artifacts.

| Check | Result |
|---|---|
| Parse all intermediate/reactant SMILES | 13 of 13 passed |
| Recompute stored canonical isomeric SMILES | All matched |
| Recompute stored molecular formulas | All matched |
| Compare final target to normalized input | Matched including represented stereochemistry |
| Check target CIP assignments | R,R |
| Check route node references | Every referenced node exists |
| Compare JSON to SVG-embedded metadata | Identical |
| Parse the SVG as XML | Passed |
| Inspect the existing PNG preview | Inspected for agreement with the route labels and structures |
| Check for embedded raster image elements in SVG | None; molecular drawings are vector paths |

### Checks not established by those results

- No atom-mapped, fully balanced reaction set was produced. Separate valid reactant/product SMILES do not prove a valid transformation.
- No forward-reaction model verified the products, and no experimental data validate the target-specific steps.
- No target yields, diastereomer ratios, enantiomeric purities, separation recoveries, or supplier availability were established.

For a more capable comparison workflow, add bond-change and atom-provenance checks, explicitly including reagent-derived atoms and departing groups. Do not require identical molecular formulas across a reaction; test a chemically appropriate balance. Keep experimental feasibility separate from computational consistency.

## 9. Stage 6 — Generate and check the SVG

A reproducible implementation is:

1. Load one structured route object containing node SMILES, edge definitions, conditions, evidence classes, and source locators.
2. Parse each molecule and generate 2D coordinates using RDKit; use `rdMolDraw2D.MolDraw2DSVG` for molecular fragments.
3. Place the fragments on a parent SVG canvas, then add arrows, conditions, node labels, evidence legends, and reference links.
4. Embed the same route data in SVG metadata; export the companion JSON from that object.
5. Render a PNG preview with an SVG renderer, inspect it, and check that diagram annotations agree with the structured data.

These are replay implementation details, not an assertion that an archived original script is available. The retained SVG does establish vector molecular drawings, text annotations, source links, and embedded route metadata. A particular original rendering library or layout algorithm should not be inferred without its script or log.

The original legend distinguishes reported building-block chemistry, proposed target adaptations, and unresolved stereoisomer separation. Maintain those distinctions in exported diagrams. A polished drawing must not visually promote an analogue-supported step to an experimentally verified one.

## 10. Comparison framework for your chemistry agent

### Hold the task constant

Give both workflows the same raw target, required stereochemistry, route count, starting-material policy, output schema, and time budget. Record differences in access to subscription literature, structure indexing, source PDFs, and inventory data. Otherwise, retrieval coverage and planning quality become confounded.

Score routes after normalizing structures and grouping genuinely equivalent strategies. Changes to solvent or protecting-group notation should not automatically count as a different synthetic strategy.

| Dimension | Baseline demonstrated here | What a dedicated-tool workflow should demonstrate |
|---|---|---|
| Target fidelity | Input repair plus RDKit structure/stereo checks | Same identity; explicit and reviewable normalization |
| Exact-route retrieval | No exact-target synthesis verified | Structure-linked experimental hit, or honest search/access limitations |
| Precedent retrieval | Three retained primary sources | Relevant examples with validated substrate/product structures and locators |
| Evidence quality | Precursor evidence separated from target adaptations | Step-level exact/analogue/proposed status without citation inflation |
| Reaction consistency | Valid individual molecular representations | Atom/bond-change checks and documented handling of reagents/byproducts |
| Stereochemical delivery | Target identity checked; delivery unresolved | Supported stereochemical sequence or an explicit unresolved branch |
| Operational completeness | Grouped steps and proposed conditions | Isolation/telescoping, work-up, purification, missing operations, and compatibility checks |
| Starting-material access | No live supplier/inventory validation | Exact purchasable identity, stereoisomer, source, and availability timestamp |
| Route diversity | Two different construction strategies sharing an entry family | Distinct disconnections, with correlated risks identified |
| Traceability and efficiency | Final artifacts retained; incomplete execution telemetry | Versioned tool calls, retrieved evidence, failed calls, elapsed time, and cost |

### Hard failures versus remaining uncertainty

Treat wrong connectivity, wrong ring size, loss of required stereochemistry in the claimed final product, and an unrelated cited example as correctness failures. Treat a clearly labelled unsupported cyclization or unresolved resolution as a **route-development gap** rather than hiding it inside a passing parser result.

Useful metrics include source-locator accuracy, exact-step coverage, analogue-step coverage, unresolved critical steps, stereochemical completeness, independently isolated operations, and validated leaf-material coverage. Define the operation-count denominator before comparing routes: one SVG arrow can contain several operations.

Do not combine these metrics into a success probability unless the scoring model has been calibrated against appropriate experimental outcomes. A longer route with stronger evidence may be more actionable than a shorter route with an unsupported central transformation.

## 11. Suggested machine-readable step record

The following is an **improved schema example**, not a record emitted during the original answer. It deliberately leaves unsupported measurements and unextracted source structures null.

```json
{
  "route_id": "A",
  "step_id": "A2_to_A3",
  "input_nodes": ["A2"],
  "output_nodes": ["A3"],
  "transformation": "intramolecular Dieckmann condensation",
  "target_step_status": "proposed",
  "precedent": {
    "evidence_class": "experimental_analogue",
    "document_id": "US2830993A",
    "locator": "Example 11; operations cross-referenced to Example 1",
    "url": "https://patents.google.com/patent/US2830993A/en",
    "source_reaction_smiles": null,
    "source_atom_mapping": null,
    "example_identity_checked_against_pdf": true,
    "verification_phase": "documentation_audit"
  },
  "applicability_gaps": [
    "6-6-6 source framework versus 6-5-6 target framework",
    "different aromatic substitution",
    "target-specific cyclization and stereochemical outcome unverified"
  ],
  "target_conditions_status": "analogue-derived proposal, not validated",
  "target_isolated_yield_percent": null,
  "target_dr": null,
  "target_ee_percent": null,
  "experimental_validation": false
}
```

For a complete implementation, also store condition fields separately from their provenance, source passages, structure-extraction confidence, individual reagent roles, work-up, isolation status, tool versions, and external atom sources. An atom-mapping algorithm's output requires validation; it is not experimental evidence.

## 12. Minimal reproducible structure audit

Save this code as `audit_route_structures.py` next to the original JSON, or pass the JSON path as an argument. It checks the delivered representations; it does not predict reaction feasibility. RDKit must be installed in the selected Python environment.

```python
from __future__ import annotations

import argparse
import json
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, rdMolDescriptors


def parse_smiles(smiles: str, label: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES for {label}: {smiles}")
    return mol


def canonical(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def audit(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    # Only the documented formatting escape is repaired for this case.
    normalized = data["target_input"].replace("\\(", "(")
    target = parse_smiles(normalized, "input target")
    target_key = canonical(target)
    if target_key != data["target_canonical_smiles"]:
        raise ValueError("Stored target differs from the normalized input")

    for node_id, node in data["structures"].items():
        mol = parse_smiles(node["smiles"], node_id)
        if canonical(mol) != node["canonical_smiles"]:
            raise ValueError(f"Canonical structure mismatch: {node_id}")
        if rdMolDescriptors.CalcMolFormula(mol) != node["formula"]:
            raise ValueError(f"Formula mismatch: {node_id}")

    final_target = parse_smiles(data["structures"]["T"]["smiles"], "T")
    if canonical(final_target) != target_key:
        raise ValueError("Final target identity/stereochemistry mismatch")
    for route_id, nodes in data["routes"].items():
        missing = set(nodes) - set(data["structures"])
        if missing:
            raise ValueError(f"Unknown nodes in route {route_id}: {missing}")

    return {
        "rdkit_version": rdBase.rdkitVersion,
        "canonical_target": target_key,
        "formula": rdMolDescriptors.CalcMolFormula(target),
        "molecular_weight": round(Descriptors.MolWt(target), 2),
        "chiral_centers": Chem.FindMolChiralCenters(
            target, includeUnassigned=True
        ),
        "checked_structures": len(data["structures"]),
        "experimental_feasibility_checked": False,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "json_path", nargs="?", type=Path,
        default=Path("fused_aminoketone_route_structures.json")
    )
    args = parser.parse_args()
    print(json.dumps(audit(args.json_path), ensure_ascii=False, indent=2))
```

The code was executed against the delivered JSON during preparation of this document. It returned 13 checked structures, C₁₆H₂₁NO, 243.35 g/mol, and R,R target assignments.

## 13. Intermediate SMILES for direct comparison

These are copied from the delivered JSON. Except for T, potential stereocenters are intentionally not assigned. This notation must not be interpreted as a measured stereoisomer distribution.

| Node | Role | SMILES |
|---|---|---|
| P0 | Boc-protected bromobenzylamine | `CC(C)(C)OC(=O)NCc1ccccc1Br` |
| P1 | Cinnamate intermediate | `COC(=O)/C=C/c1ccccc1CNC(=O)OC(C)(C)C` |
| A1 | Free isoindoline-acetate | `COC(=O)CC1NCc2ccccc21` |
| C0 | N-Boc isoindoline-acetate | `COC(=O)CC1N(C(=O)OC(C)(C)C)Cc2ccccc21` |
| M | Dimethyl isobutylmalonate | `COC(=O)C(CC(C)C)C(=O)OC` |
| A2 | Amino triester | `COC(=O)CC1N(CC(CC(C)C)(C(=O)OC)C(=O)OC)Cc2ccccc21` |
| A3 | Cyclic keto diester | `CC(C)CC1(C(=O)OC)CN2Cc3ccccc3C2C(C(=O)OC)C1=O` |
| B1 | Aldehyde | `O=CCC1N(C(=O)OC(C)(C)C)Cc2ccccc21` |
| V | Alkenyl iodide | `C=C(I)CC(C)C` |
| B2 | Allylic alcohol | `C=C(CC(C)C)C(O)CC1N(C(=O)OC(C)(C)C)Cc2ccccc21` |
| B3 | Enone | `C=C(CC(C)C)C(=O)CC1N(C(=O)OC(C)(C)C)Cc2ccccc21` |
| Tmix | Target connectivity, stereo unspecified | `CC(C)CC1CN2Cc3ccccc3C2CC1=O` |
| T | Requested stereochemical target | `CC(C)C[C@@H]1CN2[C@@H](c3c(C2)cccc3)CC1=O` |

## 14. Source and artifact references

### Primary literature/patents

**R1.** *Piperazine- and piperidine-derivatives as melanocortin receptor agonists*, US20040082590A1. Relevant locators: Preparation 1C, Steps A–C; Preparation 3C, Step A; paragraphs [0837]–[0844].  
Source: `https://patents.google.com/patent/US20040082590A1/en`

**R2.** *Quinolizine derivatives*, US2830993A. Relevant locators: Example 11 and its cross-reference to Example 1. For the example-number/substituent check, inspect PDF page 4, printed column 7.  
Source: `https://patents.google.com/patent/US2830993A/en`  
Original PDF: `https://patentimages.storage.googleapis.com/69/13/22/a2b4a05e618bd3/US2830993.pdf`

**R3.** *Intermediates useful for making tetrabenazine compounds*, US8008500B2. Relevant locators: Examples 3, 4, 7, and 10.  
Source: `https://patents.google.com/patent/US8008500B2/en`

These patents are chemical evidence sources here. No conclusions about present legal status, freedom to operate, or commercial rights are made.

### Delivered artifacts audited

- `fused_aminoketone_synthetic_routes.svg` — route drawings, evidence legend, reference links, and embedded route data.
- `fused_aminoketone_synthetic_routes.png` — preview of the SVG.
- `fused_aminoketone_route_structures.json` — 13 molecular records and the two route-node sequences.

Artifact checksums identify the exact files audited; they do not authenticate the chemical claims:

```text
SVG SHA-256:
766f1ea0eaa729c2dd0ac2b529a85ce1a24702f869d4737521ca0e1b0e1f2bd4

JSON SHA-256:
1421b2721436d30bd867aa637f3a28a8baf3d986ffd42e4b7ffb2ebd39f967a0
```

## Bottom line for comparison

This baseline produced chemically explicit, literature-linked route proposals and structurally consistent graphics. It did **not** establish a working synthesis of the specified enantiomer. The most meaningful improvement from dedicated tools would be stronger substrate-specific evidence, explicit reaction and stereochemical validation, reliable starting-material access, and reproducible provenance—not simply more routes or more tool calls.

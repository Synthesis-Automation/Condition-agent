# Synthesis route-search workflow: worked example and agent comparison

Date: 2026-10-07  
Purpose: document the workflow actually used in this conversation, so it can be compared with a chemistry agent equipped with dedicated tools.

This is a record of the methods, evidence, outputs, and limitations. It is not an experimentally validated synthesis or a full transcript of every search.

## 1. Target and result

**User-supplied SMILES**

```text
CC1(CC(C)(C)C(c(c1n2C(C)=O)c3c2ccnc3)=O)C
```

**RDKit canonical SMILES**

```text
CC(=O)n1c2c(c3cnccc31)C(=O)C(C)(C)CC2(C)C
```

| Property | Calculated result |
|---|---|
| Molecular formula | C17H20N2O2 |
| Neutral monoisotopic mass | 284.15247788 Da |
| InChIKey | ZUGALAZVIKDDDM-UHFFFAOYSA-N |
| Structural features used in the search | N-acetyl azaindole core; fused cyclohexanone; two gem-dimethyl groups |

**Outcome:** no verified preparation of the exact target was located in the web searches performed. Two literature-inspired routes were proposed. The preferred route uses an aminopyridine–diketone condensation and Pd-mediated ring closure. The alternative uses C-arylation of a diketone followed by nitro reduction and cyclization.

“Not located” describes the search outcome; it does not establish that the compound or its synthesis is absent from the literature.

## 2. Tools actually used

| Tool or capability | Actual use | What it did not establish |
|---|---|---|
| Python + RDKit | Parse and canonicalize SMILES; calculate formula, exact mass, and InChIKey; inspect connectivity; construct intermediate structures; check final graph identity | Reaction feasibility, selectivity, yield, or commercial availability |
| Web search: `system2_search_query` | Search scaffold terms, substituent patterns, building-block names, identifiers, and reaction-method terms | Exhaustive chemical structure or reaction-database retrieval |
| Web retrieval: `open` and `find` | Read indexed patent examples, publisher pages, and accessible article text; locate relevant procedures | Universal access to full papers, supporting information, or structure images |
| Assistant chemical assessment | Propose disconnections, assess analogue relevance, rank candidate routes, and identify missing evidence | An independent predictive score or experimental confirmation |
| RDKit `MolDraw2DSVG` | Draw molecular structures in the route scheme | Validation of the transformations represented by arrows |
| Python SVG composition + CairoSVG | Add arrows, conditions, evidence notes, and source links; render a PNG preview | Automated chemistry checking |
| Image inspection | Check the rendered scheme for readable labels and visible structural/layout problems | Proof of synthetic correctness |
| File-saving connector | Save the finished SVG and this Markdown document | Chemical validation |

No ASKCOS, AiZynthFinder, Syntheseus, commercial retrosynthesis service, reaction-template engine, or trained single-step retrosynthesis model was used. No SciFinder/Reaxys search, dedicated chemical substructure search, atom-mapping service, reaction-feasibility model, or live supplier-inventory API was used.

The workflow was **web literature retrieval + chemical assessment + deterministic structure handling**.

## 3. Workflow actually followed

### Stage A — Normalize and inspect the structure

1. Parse the supplied SMILES with `Chem.MolFromSmiles`.
2. Generate canonical SMILES with `Chem.MolToSmiles`.
3. Calculate formula and exact mass with `rdMolDescriptors.CalcMolFormula` and `Descriptors.ExactMolWt`.
4. Generate the InChIKey with `Chem.MolToInchiKey`.
5. Inspect atom connectivity and a 2D structure drawing.
6. Identify the features that a route must preserve: pyridine nitrogen position, ring fusion, ketone position, both gem-dimethyl groups, and acetylation of the pyrrolic nitrogen.

These operations checked the input structure and provided search terms. They did not generate synthesis routes automatically.

### Stage B — Search for the target and nearby chemistry

Searches were iterative, combining scaffold names, substitution patterns, identifiers, and transformation terms. Representative queries actually used include:

```text
"ZUGALAZVIKDDDM"
"tetramethyl" "azacarbazol" synthesis
"tetramethyl" "pyrrolo[3,2-c]"
"tetramethyl" "carbazol-4-one"
"tetramethyl" "cyclohexane-1,3-dione"
"carbazolones" "palladium" synthesis diones
"3-bromo-4-aminopyridine" "cyclohex"
"3-iodo-4-nitropyridine" synthesis
"60681-10-9" synthesis
```

Several scaffold-name queries returned irrelevant structures or other heterocyclic isomers. Search terms were adjusted using the target connectivity. A matching name fragment or molecular formula was not treated as an exact structure match.

The useful search directions were:

- Halogenated aminopyridine + cyclic 1,3-diketone → enaminone → fused azaindole.
- Ortho-nitroaryl-substituted 1,3-diketone → reductive cyclization.
- Preparation of the correctly substituted tetramethyl diketone.

### Stage C — Read and classify the precedents

For promising results, relevant example numbers and experimental text were opened and searched. The main sources were:

| Source | Evidence recovered | Transfer limitation |
|---|---|---|
| JP2010500961A, Example 37(1)–(2) | A brominated aminopyridine condenses with cyclohexane-1,3-dione, followed by Pd-catalyzed cyclization | Different aminopyridine substitution and an unmethylated diketone |
| Janreddy et al., *Eur. J. Org. Chem.* 2011, 2360–2365 | C-arylated cyclic diketones undergo Fe/AcOH-mediated reductive cyclization to carbazolones | Does not verify the proposed nitropyridine/tetramethyl substrate combination |
| US10160695B2, Figure 1B | Halonitropyridines can undergo NO2 displacement by fluoride | Flags a competing pathway; does not prove the same outcome with a diketone enolate |
| ACS chapter on triketone HPPD herbicides | Describes preparation of the tetramethyl diketone by dianion methylation of a trimethyl precursor | Read through accessible reproduced chapter text; a complete experimental procedure was not verified |

Evidence came from patent HTML, publisher information, and accessible reproduced article/chapter text. Supporting information and all graphical substrate scopes were not exhaustively inspected. A direct Python attempt to retrieve patent HTML returned HTTP 503; the needed example text was available through the web retrieval tool.

### Stage D — Construct and rank route proposals

The following considerations were used qualitatively, without numerical scores:

- Correct molecular connectivity and substituent placement.
- Relevance of actual experimental examples.
- Number of new ring-forming steps.
- Regioselectivity and chemoselectivity risks.
- Availability or preparability of required fragments.
- Clearly identified gaps between the precedent and the proposed substrate.

**Route 1 was ranked first** because it starts with the complete carbocyclic skeleton, uses a halogen to define ring closure, and has an experimental aminopyridine/enaminone/Pd-cyclization precedent. Its hindered condensation and exact cyclization remain unverified.

**Route 2 was ranked lower** because the initial C-arylation lacks an established procedure for this exact combination and must distinguish productive C–C coupling from O-arylation or competing substitution. The downstream reductive cyclization has reaction-class precedent.

These rankings are chemical judgments based on retrieved evidence, not measured success probabilities.

### Stage E — Validate structures and prepare the scheme

1. Construct SMILES for the starting fragments and proposed intermediates.
2. Parse them with RDKit and calculate their formulas.
3. Programmatically add an acetyl group to the proposed N–H tricycle.
4. Confirm that its canonical molecular graph matches both the normalized target and the user's original SMILES.
5. Draw the molecules using RDKit and assemble the route SVG with arrows and annotations.
6. Render the SVG to PNG with CairoSVG and visually inspect it.
7. Save the SVG with explicit labels distinguishing proposed transformations from literature precedents.

The final graph-identity check verifies the destination of the proposed N-acetylation. It does not constitute atom-mapped validation of every reaction or prove that either route works.

## 4. Concrete route outputs

### Route 1 — Preferred proposal

```text
3-Bromopyridin-4-amine + tetramethyl diketone D
    → pyridyl enaminone E
    → N–H fused tricycle B
    → N-acetyl target
```

| Step | Candidate conditions | Evidence status |
|---|---|---|
| Condensation | Catalytic p-TsOH, toluene, reflux, water removal | Adapted from an experimental analogue |
| Pd ring closure | PdCl2(PPh3)2, Cs2CO3, toluene, reflux | Adapted from an experimental analogue |
| N-acetylation | Ac2O/DMAP/Et3N | Proposed; exact substrate procedure not located |

### Route 2 — Higher-risk proposal

```text
3-Iodo-4-nitropyridine + tetramethyl diketone D
    → C-arylated diketone C
    → N–H fused tricycle B by nitro reduction/cyclization
    → N-acetyl target
```

| Step | Candidate operation | Evidence status |
|---|---|---|
| C-arylation | Develop a coupling at the pyridine C–I bond | No validated catalyst/condition set supplied for this combination |
| Reductive cyclization | Fe/AcOH, reflux | Reaction-class precedent; substrate transfer unverified |
| N-acetylation | Same proposed final step as Route 1 | Unverified for this substrate |

### Structures supplied for independent checking

| Identifier | Description | SMILES |
|---|---|---|
| D | 4,4,6,6-Tetramethylcyclohexane-1,3-dione | `CC1(C)CC(C)(C)C(=O)CC1=O` |
| A | 3-Bromopyridin-4-amine | `Nc1ccncc1Br` |
| E | Proposed pyridyl enaminone | `CC1(C)CC(C)(C)C(=O)C=C1Nc1ccncc1Br` |
| B | Proposed N–H tricycle | `CC1(C)CC(C)(C)C(=O)c2c1[nH]c1ccncc21` |
| I | 3-Iodo-4-nitropyridine | `O=[N+]([O-])c1ccncc1I` |
| C | Proposed C-arylated diketone, drawn in keto form | `CC1(C)CC(C)(C)C(=O)C(c2cnccc2[N+](=O)[O-])C1=O` |

D is not dimedone. Using dimedone would place the methyl groups incorrectly. No live stock check was performed for D or the other building blocks.

## 5. Main gaps and potential failure points

1. **Exact-target retrieval:** text searches can miss compounds represented only as structure images or database records.
2. **Analogue transfer:** a different pyridine nitrogen position and increased steric hindrance can substantially change reactivity.
3. **Starting-material readiness:** the route depends on obtaining D; its full preparative sequence and current availability remain unresolved.
4. **Route 2 completeness:** its first step is a development proposal, not an executable literature procedure.
5. **Validation depth:** valid SMILES, plausible atom connectivity, and a clean diagram do not validate reaction chemistry.
6. **Performance measurement:** elapsed time, search cost, and model/tool token usage were not systematically benchmarked.
7. **Search efficiency:** broad scaffold-name searches produced substantial noise. Dedicated chemical structure/reaction retrieval could make this stage more targeted.

## 6. How to compare with your tool-equipped agent

For an independent comparison, give your agent only the original target SMILES and the same task. Keep the proposed routes hidden until it has produced its own result. Record its tool access and resource budget separately.

| Evaluation dimension | This workflow's baseline | What to ask your agent to demonstrate |
|---|---|---|
| Exact structure recognition | RDKit normalization and connectivity inspection | Preserve ring fusion, nitrogen positions, and all methyl groups |
| Exact precedent retrieval | Exact synthesis not located | Provide a matching structure, source, and example identifier if found |
| Reaction evidence | Related experimental precedents | Show the actual precedent substrate and explain the structural differences |
| Route completeness | Route 1 proposed; Route 2 has an unresolved first step | Specify every transformation and identify unresolved steps |
| Building-block access | Literature lead for D; no stock check | Verified purchase options or a sourced preparative sequence |
| Chemical validation | Structure/formula checks and final graph identity | Add atom mapping, balance checks, selectivity assessment, and model outputs where useful |
| Condition provenance | Analogue-derived and proposed conditions labeled | Distinguish reported conditions from predicted or manually suggested ones |
| Confidence calibration | No target yields predicted | Avoid transferring analogue yields to untested substrates |
| Citation auditability | Patent example numbers and paper DOI | Make every claimed experimental step traceable |
| Efficiency | Not formally measured | Record elapsed time, tool calls, retrieval volume, and cost |

A stronger answer need not contain more routes. Finding a closer precedent, resolving the preparation of D, or eliminating an implausible step would be more meaningful improvements.

Suggested evidence labels for both systems:

- **Exact experimental:** the same substrate/product transformation is reported.
- **Close experimental analogue:** a similar transformation is reported, with differences explicitly recorded.
- **Reaction-class precedent:** the chemistry is established, but the substrate match is weak.
- **Unvalidated proposal:** a chemically motivated suggestion without adequate matching experimental evidence.

## 7. Sources and artifact

1. **JP2010500961A**, Example 37(1)–(2): aminopyridine condensation and Pd cyclization.  
   https://patents.google.com/patent/JP2010500961A/en
2. **Janreddy et al.**, “An Easy Access to Carbazolones and 2,3-Disubstituted Indoles,” *European Journal of Organic Chemistry* 2011, 2360–2365.  
   https://doi.org/10.1002/ejoc.201001357
3. **US10160695B2**, Figure 1B and related experimental discussion: halonitropyridine substitution behavior.  
   https://patents.google.com/patent/US10160695B2/en
4. **“The Synthesis and Structure–Activity Relationships of the Triketone HPPD Herbicides,”** ACS Symposium Series chapter: upstream diketone preparation lead.  
   https://doi.org/10.1021/bk-2001-0774.ch002  
   Accessible reproduced text used during the search: https://dokumen.pub/agrochemical-discovery-insect-weed-and-fungal-control-9780841237247-9780841218338-0-8412-3724-7.html
5. Route diagram previously generated in this conversation: **azaindole_ketone_routes.svg**.

No experimental yields for the exact target or its proposed intermediates were established in this work.

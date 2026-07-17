# Reactive-Handle Taxonomy and Reaction-Labeling System

## Implementation Design Specification

**Status:** Draft v1  
**Target implementation:** Python 3.12 + RDKit  
**Primary consumer:** Codex or another coding agent  
**Encoding:** UTF-8

---

## 1. Purpose

This document specifies a compact, extensible system that converts:

- compound SMILES; and
- reaction SMILES

into chemist-readable reactive-site labels such as:

- `Ar–Br`
- `Ar–OMs`
- `Ar–B(OH)2`
- `Ar–NH2`
- `R1R2–NH`
- `Ar–OH`
- `R–SH`
- `R–C(O)OH`

The system should capture the **minimum structural requirements of a synthetic transformation** while treating the rest of each molecule as context.

For example, a Buchwald–Hartwig C–N coupling is represented by the minimal grammar:

```text
Ar–X + N–H → Ar–N
```

The identities of the amine substituents, aryl substituents, steric environment, heteroatoms, and other functional groups are retained as context features. They are not required to define the core reaction family.

This design deliberately does **not** attempt to classify every structure or every reaction in organic chemistry.

---

## 2. Core Design Principle

### 2.1 Minimal reactive requirement

Each reaction family is defined by the smallest set of reactive-site requirements needed to explain the key bond changes.

Examples:

```text
C–N coupling:
sp2-C–leaving-group + N–H → sp2-C–N

C–O coupling:
sp2-C–leaving-group + O–H → sp2-C–O

C–S coupling:
sp2-C–leaving-group + S–H → sp2-C–S

Suzuki coupling:
sp2-C–leaving-group + C–B → C–C

Amide formation:
acyl electrophilic center + N–H → acyl-C–N
```

Everything outside these minimal sites is represented separately as context.

### 2.2 Reactive-site taxonomy, not whole-molecule taxonomy

The system classifies **reactive sites**, not entire compounds.

A multifunctional molecule may contain multiple sites:

```text
4-bromoacetophenone
├── leaving-group site: Ar–Br
└── retained motif/context: Ar–C(O)R
```

A compound may also contain several candidate sites of the same family. Each site must have a unique atom-index-based identifier.

### 2.3 Structured representation first, label second

The canonical machine representation is a structured object derived from the molecular graph.

Chemist-facing labels such as `Ar–Br` and `R1R2–NH` are generated from that object.

Do not use display labels as the only machine identity.

---

## 3. Version 1 Scope

Version 1 supports six site families.

```text
reactive_site
├── leaving_group
├── pronucleophile_XH
├── transfer_group
├── electrophilic_center
├── aromatic_CH
└── unsaturated_bond
```

### 3.1 Leaving-group sites

General forms:

```text
Context–X
```

Examples:

```text
Ar–Br
HeteroAr–Cl
Alkenyl–OTf
Alkyl–OMs
```

### 3.2 Pronucleophilic X–H sites

General form:

```text
reactive atom + retained contexts + one or more H
```

Examples:

```text
Ar–NH2
R1R2–NH
Ar–OH
R–SH
RSO2–NHR
RC(O)–NHR
R–NH–NH2
```

The initial supported center atoms are:

```text
N
O
S
C(sp) for terminal alkynes
```

### 3.3 Transfer-group sites

General form:

```text
Context–M
```

`M` is a transferable organometallic or organometalloid handle.

Examples:

```text
Ar–B(OH)2
HeteroAr–B(OR)2
Alkenyl–BF3K
Ar–ZnX
Ar–MgX
Ar–SnR3
Ar–SiR3
```

### 3.4 Electrophilic functional centers

General form:

```text
retained context–electrophilic center–leaving/activatable group
retained context–polarized multiple bond
```

Initial focus:

```text
acyl centers
sulfonyl centers
aldehyde and ketone carbonyl centers
```

Examples:

```text
R–C(O)OH
R–C(O)Cl
R–C(O)OR*
R–S(O)2Cl
R–CH=O
R2C=O
```

This family supports both substitution at acyl/sulfonyl centers and addition to aldehyde/ketone carbonyls. The `reaction_mode` field distinguishes those mechanisms.

### 3.5 Aromatic C–H sites

An aromatic C–H site is localized on an aromatic carbon bearing an implicit or explicit hydrogen. The detector emits `ArH` for carbocyclic aromatic systems and `HetArH` for heteroaromatic systems. Regioselectivity belongs to reaction-family descriptors rather than the generic SMARTS definition.

### 3.6 Unsaturated carbon–carbon bonds

Alkenes and alkynes are bond-localized sites with two retained endpoints. Their `bond` topology is distinct from an `edge`, where an anchor is connected to a removable or transferred handle. A terminal alkyne may validly emit both an `unsaturated_bond` site and a `pronucleophile_XH` site because the two records represent different reaction modes.

---

## 4. Vocabulary

The vocabulary should be small, stable, and versioned.

## 4.1 Context tokens

Initial context tokens:

```text
Ar
HeteroAr
Alkenyl
Alkynyl
Alkyl
C(O)R
C(O)OR
C(O)NR
SO2R
P(O)R2
N
O
S
H
Other
```

### Definitions

| Token | Meaning |
|---|---|
| `Ar` | Aromatic carbocyclic context |
| `HeteroAr` | Aromatic heterocyclic context |
| `Alkenyl` | Non-aromatic sp2 carbon context |
| `Alkynyl` | sp carbon context |
| `Alkyl` | Non-aromatic saturated carbon context |
| `C(O)R` | Acyl/carbonyl-carbon context |
| `C(O)OR` | Carbamate/ester-like carbonyl context |
| `C(O)NR` | Urea/amidic carbonyl context |
| `SO2R` | Sulfonyl context |
| `P(O)R2` | Phosphoryl context |
| `N`, `O`, `S` | Direct heteroatom context |
| `H` | Retained hydrogen context when needed |
| `Other` | Unsupported or unresolved context |

Context tokens describe the local environment near the reactive site. Exact molecular fragments remain available through the original molecule and atom indices.

## 4.2 Leaving-group tokens

Initial leaving groups:

```text
F
Cl
Br
I
OTf
OTs
OMs
```

Optional later additions:

```text
ONf
OPO(OR)2
N2+
SR
SeR
iodonium
```

## 4.3 Transfer-group tokens

Initial transfer handles:

```text
B(OH)2
B(OR)2
BF3K
BPin
BMIDA
ZnX
MgX
SnR3
SiR3
```

The canonical token should identify the chemically relevant transfer family. Exact substituents can be recorded separately.

## 4.4 Pronucleophile center tokens

Initial centers:

```text
N
O
S
Csp
```

Do not include arbitrary C–H sites in v1. Carbon pronucleophiles should be added by explicit supported classes.

## 4.5 Electrophilic-center families

Initial center families:

```text
Acyl
Sulfonyl
```

Optional later additions:

```text
Phosphoryl
Carbonyl
Imidoyl
ActivatedAlkene
```

---

## 5. Canonical Site Data Model

All detected sites use a common envelope.

```json
{
  "schema_version": "1.0",
  "site_id": "mol0:atom5:site0",
  "site_type": "pronucleophile_XH",
  "topology": "atom",
  "component_index": 0,
  "atom_indices": [5],
  "bond_indices": [],
  "canonical_signature": "XH|N|H2|Ar",
  "chemist_label": "Ar–NH2",
  "context_features": {},
  "confidence": 1.0,
  "warnings": []
}
```

Required fields:

| Field | Description |
|---|---|
| `schema_version` | Site schema version |
| `site_id` | Stable identifier using the reactive center/anchor atom or bond locus |
| `site_type` | One of the supported site families |
| `topology` | `edge`, `atom`, `center`, or `bond` |
| `component_index` | Reactant/product component index |
| `atom_indices` | Atoms defining the site |
| `bond_indices` | Bonds defining the site |
| `canonical_signature` | Stable machine-readable signature |
| `chemist_label` | Generated display label |
| `context_features` | Optional overlays |
| `confidence` | Detection confidence from 0 to 1 |
| `warnings` | Ambiguities or unsupported features |

---

## 6. Site-Type Schemas

## 6.1 Leaving-group site

```json
{
  "site_type": "leaving_group",
  "topology": "edge",
  "anchor_atom_index": 3,
  "handle_atom_indices": [2],
  "anchor_context": "Ar",
  "handle_token": "Br",
  "connector_bond_index": 2,
  "connector_bond_type": "SINGLE",
  "canonical_signature": "LG|Ar|Br",
  "chemist_label": "Ar–Br"
}
```

Canonical format:

```text
LG|<anchor_context>|<handle_token>
```

Examples:

```text
LG|Ar|Br
LG|HeteroAr|Cl
LG|Alkenyl|OTf
LG|Alkyl|OMs
```

## 6.2 Pronucleophile X–H site

```json
{
  "site_type": "pronucleophile_XH",
  "topology": "atom",
  "center_atom_index": 8,
  "center_element": "N",
  "formal_charge": 0,
  "aromatic": false,
  "h_count": 2,
  "contexts": ["Ar"],
  "canonical_signature": "XH|N|H2|Ar",
  "chemist_label": "Ar–NH2",
  "derived_family": "primary_arylamine"
}
```

Canonical format:

```text
XH|<center>|H<h_count>|<sorted_context_tokens>
```

Examples:

```text
XH|N|H2|Alkyl
XH|N|H2|Ar
XH|N|H1|Alkyl,Alkyl
XH|N|H1|Ar,Alkyl
XH|N|H1|SO2R,Alkyl
XH|O|H1|Ar
XH|S|H1|Alkyl
XH|Csp|H1|Alkynyl
```

## 6.3 Transfer-group site

```json
{
  "site_type": "transfer_group",
  "topology": "edge",
  "anchor_atom_index": 4,
  "handle_atom_indices": [7, 8, 9],
  "anchor_context": "Ar",
  "handle_token": "B(OH)2",
  "connector_bond_index": 4,
  "canonical_signature": "TM|Ar|B(OH)2",
  "chemist_label": "Ar–B(OH)2"
}
```

Canonical format:

```text
TM|<anchor_context>|<transfer_handle>
```

## 6.4 Electrophilic-center site

```json
{
  "site_type": "electrophilic_center",
  "topology": "center",
  "center_family": "Acyl",
  "center_atom_index": 4,
  "multiple_bond_atom_indices": [5],
  "retained_context": "Ar",
  "leaving_or_activatable_group": "OH",
  "activation_state": "latent",
  "canonical_signature": "EC|Acyl|Ar|OH|latent",
  "chemist_label": "Ar–C(O)OH"
}
```

Canonical format:

```text
EC|<center_family>|<retained_context>|<leaving_group>|<activation_state>
```

Examples:

```text
EC|Acyl|Alkyl|OH|latent
EC|Acyl|Ar|Cl|activated
EC|Acyl|Ar|OR|ester
EC|Sulfonyl|Ar|Cl|activated
EC|Carbonyl|aldehyde|Alkyl|addition
EC|Carbonyl|ketone|Alkyl,Alkyl|addition
```

Acyl and sulfonyl centers use `reaction_mode: substitution`. Aldehyde and ketone carbonyl centers use `reaction_mode: addition` and expose `center` and `heteroatom` atom roles.

## 6.5 Aromatic C–H site

```json
{
  "site_type": "aromatic_CH",
  "topology": "atom",
  "atom_roles": {"center": [3]},
  "handle_token": "HetArH",
  "ring_context": "HeteroAr",
  "h_count": 1,
  "canonical_signature": "CH|HetArH",
  "chemist_label": "HeteroAr–H"
}
```

## 6.6 Unsaturated-bond site

```json
{
  "site_type": "unsaturated_bond",
  "topology": "bond",
  "atom_indices": [2, 3],
  "bond_indices": [2],
  "atom_roles": {"endpoint_a": [2], "endpoint_b": [3]},
  "handle_token": "Alkene",
  "bond_order": 2,
  "canonical_signature": "PI|Alkene",
  "chemist_label": "C=C"
}
```

---

## 7. Context Classification

Context classification is a core service used by all site detectors.

## 7.1 Input

```text
molecule
anchor atom
neighbor atom or retained fragment
excluded site atoms
```

## 7.2 Output

```json
{
  "token": "Ar",
  "subtype": "phenyl",
  "atom_indices": [3, 4, 5, 6, 7, 8],
  "features": {
    "heteroatom_count": 0,
    "substituted": true
  }
}
```

## 7.3 Precedence

Use the following precedence when a neighboring atom could match multiple broad categories:

```text
SO2R
> P(O)R2
> C(O)NR
> C(O)OR
> C(O)R
> HeteroAr
> Ar
> Alkenyl
> Alkynyl
> Alkyl
> N
> O
> S
> Other
```

The most reactivity-defining local context must win.

Examples:

```text
RC(O)–NHR
```

The N-attached carbonyl carbon must be `C(O)R`, not `Alkyl`.

```text
RSO2–NHR
```

The N-attached sulfur must be `SO2R`, not `S`.

## 7.4 Aromatic context

Classify an attached aromatic carbon as:

- `Ar` if its aromatic ring contains only carbon;
- `HeteroAr` if the relevant aromatic ring contains at least one heteroatom.

For fused or ambiguous systems, use the ring system containing the attachment atom.

## 7.5 Context depth

Default context depth:

```text
first shell:
directly bonded retained neighbors

second shell:
only reactivity-dominant patterns needed to identify a composite token
```

Examples requiring limited second-shell inspection:

- carbonyl carbon → inspect O/N/O substitution;
- sulfonyl sulfur → inspect two oxo atoms;
- neighboring N in hydrazines → inspect its substitution and H count.

Do not recursively encode the full molecule into the core site signature.

---

## 8. Pronucleophile X–H Detection

## 8.1 Candidate atoms

Version 1 candidates:

```text
neutral N with total H count >= 1
neutral O with total H count >= 1
neutral S with total H count >= 1
terminal alkyne carbon with total H count = 1
```

Formal charge and valence restrictions should be configurable.

## 8.2 Hydrogen count

Use total hydrogen count after RDKit sanitization:

```python
atom.GetTotalNumHs(includeNeighbors=True)
```

Do not require explicit `[H]` atoms in SMILES.

## 8.3 Nitrogen precedence

Classify N–H sites using precedence:

```text
aromatic [nH]
> sulfonyl-attached N
> phosphoryl-attached N
> urea-like N
> carbamate-like N
> amide-like N
> hydrazine/hydrazide N
> free amine
```

The site signature remains compositional. The derived family is an additional field.

Examples:

```text
N{Ar; H2}
→ XH|N|H2|Ar
→ Ar–NH2
→ primary_arylamine

N{Alkyl,Alkyl; H1}
→ XH|N|H1|Alkyl,Alkyl
→ R1R2–NH
→ secondary_dialkylamine

N{SO2R,Alkyl; H1}
→ XH|N|H1|SO2R,Alkyl
→ RSO2–NHR
→ secondary_sulfonamide
```

## 8.4 Hydrazine support

For:

```text
R–NH–NH2
```

detect two separate N–H sites.

Substituted nitrogen:

```text
XH|N|H1|Alkyl,N
```

Terminal nitrogen:

```text
XH|N|H2|N
```

Add second-shell context to distinguish:

```text
neighbor_N_h_count
neighbor_N_contexts
neighbor_N_attached_to_acyl
neighbor_N_attached_to_sulfonyl
```

Example structured output:

```json
{
  "canonical_signature": "XH|N|H2|N",
  "chemist_label": "R–NH–NH2",
  "derived_family": "terminal_hydrazine_NH2",
  "context_features": {
    "neighbor_N_h_count": 1,
    "neighbor_N_contexts": ["Alkyl"]
  }
}
```

---

## 9. Label Rendering

Label rendering must be deterministic and separate from site detection.

## 9.1 General rules

- Use Unicode en dash `–` for chemist-facing labels.
- Use ASCII-safe canonical signatures internally.
- Use `R1`, `R2`, `Ar1`, and `Ar2` where distinguishable substituent placeholders are useful.
- Use compact conventional labels rather than exposing internal graph syntax.

## 9.2 Leaving-group labels

```text
LG|Ar|Br                 → Ar–Br
LG|HeteroAr|Cl           → HeteroAr–Cl
LG|Alkenyl|OTf           → Alkenyl–OTf
LG|Alkyl|OMs             → R–OMs
```

## 9.3 Transfer-group labels

```text
TM|Ar|B(OH)2             → Ar–B(OH)2
TM|HeteroAr|BPin         → HeteroAr–BPin
TM|Alkenyl|ZnX           → Alkenyl–ZnX
```

## 9.4 N–H labels

```text
XH|N|H2|Alkyl            → R–NH2
XH|N|H2|Ar               → Ar–NH2
XH|N|H1|Alkyl,Alkyl      → R1R2–NH
XH|N|H1|Ar,Alkyl         → Ar–NH–R
XH|N|H1|Ar,Ar            → Ar1Ar2–NH
XH|N|H1|C(O)R,Alkyl      → RC(O)–NHR
XH|N|H1|SO2R,Alkyl       → RSO2–NHR
XH|N|H2|N                → N–NH2
```

## 9.5 O–H and S–H labels

```text
XH|O|H1|Alkyl            → R–OH
XH|O|H1|Ar               → Ar–OH
XH|S|H1|Alkyl            → R–SH
XH|S|H1|Ar               → Ar–SH
```

## 9.6 Electrophilic-center labels

```text
EC|Acyl|Alkyl|OH|latent  → R–C(O)OH
EC|Acyl|Ar|Cl|activated  → Ar–C(O)Cl
EC|Sulfonyl|Ar|Cl|activated
                           → Ar–S(O)2Cl
```

---

## 10. Compound Analysis Pipeline

```text
SMILES
↓
RDKit parse and sanitize
↓
standardize molecule
↓
split components if required
↓
detect candidate reactive sites
↓
classify contexts
↓
generate canonical signatures
↓
generate chemist labels
↓
calculate optional context overlays
↓
return structured result
```

### 10.1 Standardization

Recommended initial operations:

- parse SMILES;
- sanitize;
- normalize aromaticity;
- preserve atom maps when present;
- optionally normalize common salts and counterions;
- do not neutralize structures by default if charge is chemically relevant;
- retain original SMILES and canonical SMILES.

### 10.2 Output

```json
{
  "input_smiles": "Nc1ccccc1",
  "canonical_smiles": "Nc1ccccc1",
  "sites": [
    {
      "site_type": "pronucleophile_XH",
      "canonical_signature": "XH|N|H2|Ar",
      "chemist_label": "Ar–NH2",
      "derived_family": "primary_arylamine"
    }
  ]
}
```

---

## 11. Reaction Analysis Pipeline

```text
reaction SMILES
↓
parse reactants, agents, products
↓
standardize each component
↓
retain or generate atom mapping
↓
detect all candidate sites on both sides
↓
calculate mapped bond changes
↓
associate changed atoms/bonds with detected sites
↓
apply reaction-family grammars
↓
assign reactant roles
↓
generate concise reaction label
↓
return contexts and confidence
```

## 11.1 Atom mapping

Atom mapping is strongly recommended.

If the input reaction is not mapped:

1. call a configured atom mapper;
2. store the mapping method and version;
3. calculate bond changes from mapped atoms;
4. lower confidence if mapping is ambiguous.

The site detector itself should not depend on mapping. Mapping is needed for reaction-event assignment.

## 11.2 Bond changes

Represent graph differences explicitly:

```json
{
  "broken_bonds": [
    {"a1": 1, "a2": 2, "bond_type": "SINGLE"}
  ],
  "formed_bonds": [
    {"a1": 1, "a2": 8, "bond_type": "SINGLE"}
  ],
  "changed_bonds": []
}
```

Hydrogen changes may be inferred when hydrogens are implicit.

For Buchwald coupling, the required heavy-atom evidence is:

```text
C–X bond absent in product
C–N bond newly present in product
```

Loss of N–H is inferred from the change in N substitution/H count.

---

## 12. Reaction Grammar Registry

Reaction families are defined by minimal site requirements plus expected bond changes.

## 12.1 Common grammar schema

Component organization is declared independently of the chemical roles:

```json
{
  "role_relationships": [
    {
      "roles": ["electrophile", "nucleophile"],
      "component_relation": "same_or_different"
    }
  ]
}
```

`same` restricts a grammar to intramolecular assignments, `different` requires
separate reactant components, and `same_or_different` lets one grammar support
both. Operators must merge removals safely when roles share a component and
must validate the resulting graph normally. `distinct_components` is retained
only as a legacy loader contract and is not used by current definitions.

```json
{
  "id": "buchwald_hartwig_cn",
  "display_name": "Buchwald–Hartwig C–N coupling",
  "required_roles": {
    "electrophile": {
      "site_type": "leaving_group",
      "anchor_context_any": ["Ar", "HeteroAr", "Alkenyl"]
    },
    "nucleophile": {
      "site_type": "pronucleophile_XH",
      "center_any": ["N"],
      "h_count_min": 1
    }
  },
  "required_bond_changes": {
    "break": ["anchor–leaving_group"],
    "form": ["anchor–N"]
  },
  "product_expectation": {
    "n_h_count_delta": -1
  }
}
```

## 12.2 C–N coupling

Minimal grammar:

```text
Ar/HeteroAr/Alkenyl–X + N–H → C–N
```

Applicable reaction families may include:

```text
Buchwald–Hartwig
Ullmann/Goldberg
SNAr
other catalytic or non-catalytic C–N substitutions
```

The minimal site grammar identifies the transformation class. Catalyst/reagent evidence may distinguish named mechanisms or reaction families.

## 12.3 C–O coupling

```text
Ar/HeteroAr/Alkenyl–X + O–H → C–O
```

## 12.4 C–S coupling

```text
Ar/HeteroAr/Alkenyl–X + S–H → C–S
```

## 12.5 Suzuki-type C–C coupling

```text
Ar/HeteroAr/Alkenyl–X + C–B → C–C
```

Required site types:

```text
leaving_group + transfer_group(B)
```

## 12.6 Other transfer couplings

```text
leaving_group + transfer_group(Zn/Mg/Sn/Si) → C–C
```

Named family classification may use the transfer handle:

```text
B → Suzuki
Zn → Negishi
Mg → Kumada
Sn → Stille
Si → Hiyama
```

## 12.7 Chan–Lam-type coupling

```text
C–B + N–H/O–H/S–H → C–N/O/S
```

Required site types:

```text
transfer_group(B) + pronucleophile_XH
```

## 12.8 Sonogashira-type coupling

```text
Ar/HeteroAr/Alkenyl–X + C(sp)–H → C–C
```

Required site types:

```text
leaving_group + pronucleophile_XH(center=Csp)
```

## 12.9 Amide formation

Minimal grammar:

```text
acyl electrophilic center + N–H → acyl-C–N
```

Examples:

```text
R–C(O)OH + R1R2NH → R–C(O)NR1R2
R–C(O)Cl + RNH2 → R–C(O)NHR
```

Required site types:

```text
electrophilic_center(center_family=Acyl)
+
pronucleophile_XH(center=N)
```

## 12.10 Ester formation

```text
acyl electrophilic center + O–H → acyl-C–O
```

## 12.11 Thioester formation

```text
acyl electrophilic center + S–H → acyl-C–S
```

## 12.12 Sulfonamide formation

```text
sulfonyl electrophilic center + N–H → S–N
```

## 12.13 Alkyl C–N/O/S substitution

```text
Alkyl–X + N–H/O–H/S–H → Alkyl–N/O/S
```

Required site types are `leaving_group(anchor_context=Alkyl)` and
`pronucleophile_XH`. This grammar identifies a structural substitution class;
it must not select SN1, SN2, Mitsunobu, or a protection family without separate
condition or mechanistic evidence.

## 12.14 Terminal-alkene Heck coupling

```text
Ar/HeteroAr–X + H2C=CHR → Ar/HeteroAr–CH=CHR
```

Required site types are `leaving_group` and
`unsaturated_bond(handle_token=Alkene)`. At least one alkene endpoint must have
two hydrogens in the initial conservative grammar. The operator attaches the
aryl anchor at that endpoint and records loss of one endpoint hydrogen. Exact
stereochemical verification remains separate from constitutional matching; the
operator must not invent E/Z geometry.

## 12.15 Friedel–Crafts acylation

```text
Ar/HeteroAr–H + R–C(O)X → Ar/HeteroAr–C(O)–R
```

Required site types are
`electrophilic_center(center_family=Acyl, activation_state=activated)` and
`aromatic_CH`. The operator removes the acyl leaving fragment, joins the acyl
center to the aromatic carbon, and records loss of aromatic hydrogen. Latent
carboxylic acids are excluded. For substituted rings, exact observed-product
reconstruction must resolve the participating aromatic site before the named
family is assigned.

---

## 13. Example: Buchwald C–N Reaction

Input:

```text
Brc1ccccc1C(=O)C.Nc1ccccc1>>CC(=O)c1ccccc1Nc1ccccc1
```

## 13.1 Reactant site detection

Electrophile:

```json
{
  "site_type": "leaving_group",
  "anchor_context": "Ar",
  "handle_token": "Br",
  "canonical_signature": "LG|Ar|Br",
  "chemist_label": "Ar–Br"
}
```

Nucleophile:

```json
{
  "site_type": "pronucleophile_XH",
  "center_element": "N",
  "h_count": 2,
  "contexts": ["Ar"],
  "canonical_signature": "XH|N|H2|Ar",
  "chemist_label": "Ar–NH2",
  "derived_family": "primary_arylamine"
}
```

Retained electrophile context:

```json
{
  "motif": "Ar–C(O)R",
  "relationship_to_reactive_site": "ortho",
  "context_feature": "ortho_acyl"
}
```

## 13.2 Product site

```json
{
  "site_type": "pronucleophile_XH",
  "center_element": "N",
  "h_count": 1,
  "contexts": ["Ar", "Ar"],
  "canonical_signature": "XH|N|H1|Ar,Ar",
  "chemist_label": "Ar1Ar2–NH",
  "derived_family": "secondary_diarylamine"
}
```

## 13.3 Reaction output

```json
{
  "reaction_class": "C_N_coupling",
  "reaction_family": "buchwald_hartwig_cn",
  "reaction_label": "Ar–Br + Ar–NH2 → Ar1Ar2–NH",
  "roles": {
    "electrophile": "LG|Ar|Br",
    "nucleophile": "XH|N|H2|Ar"
  },
  "bond_changes": {
    "broken": ["Ar–Br"],
    "formed": ["Ar–N"]
  },
  "contexts": {
    "electrophile": ["ortho_acyl"],
    "nucleophile": ["primary_arylamine"]
  },
  "confidence": 0.98
}
```

---

## 14. Context Overlays for Condition Recommendation

Core site signatures should remain small. Additional condition-relevant features are overlays.

## 14.1 Electrophile overlays

Examples:

```text
leaving_group_identity
aryl_vs_heteroaryl
heteroaryl_ring_type
ortho_substitution_count
electron_poor
electron_rich
additional_halide
additional_coordination_site
sterically_hindered
```

## 14.2 Nucleophile overlays

Examples:

```text
primary_vs_secondary
alkyl_vs_aryl
cyclic
ring_size
alpha_branching_count
aryl_ortho_substitution_count
additional_basic_N
multiple_NH_sites
hydrazine_site_type
amide_or_sulfonamide
sterically_hindered
```

## 14.3 Acyl-center overlays

Examples:

```text
carboxylic_acid
acyl_halide
activated_ester
anhydride
ester
sterically_hindered
alpha_stereocenter
additional_acid
additional_amine
```

## 14.4 Separation rule

Do not create a new stable site signature for every overlay combination.

Preferred:

```json
{
  "canonical_signature": "XH|N|H2|Ar",
  "context_features": {
    "aryl_ortho_substitution_count": 2,
    "electron_poor": true
  }
}
```

Avoid:

```text
XH|N|H2|ortho_disubstituted_electron_poor_Ar
```

unless later validation proves that the distinction must be part of the stable vocabulary.

---

## 15. Condition-Recommendation Feature Object

```json
{
  "reaction_family": "buchwald_hartwig_cn",
  "electrophile_signature": "LG|Ar|Br",
  "partner_signature": "XH|N|H2|Ar",
  "formed_bond": "C–N",
  "electrophile_context": {
    "ortho_acyl": true,
    "heteroaryl": false,
    "steric_class": "moderate"
  },
  "partner_context": {
    "derived_family": "primary_arylamine",
    "aryl_ortho_substitution_count": 0,
    "steric_class": "low"
  }
}
```

Recommended hierarchical fallback:

```text
exact site signatures + overlays
↓
exact site signatures
↓
site-family pair
↓
reaction-family prior
```

Example:

```text
LG|Ar|Br + XH|N|H2|Ar + ortho_acyl
↓
LG|Ar|Br + XH|N|H2|Ar
↓
aryl leaving group + primary arylamine
↓
C–N coupling global prior
```

---

## 16. Proposed Python Package Structure

```text
reactive_taxonomy/
├── __init__.py
├── models.py
├── vocabulary/
│   ├── context_tokens.json
│   ├── handle_tokens.json
│   ├── reaction_grammars.json
│   └── rendering_rules.json
├── standardize.py
├── context.py
├── sites/
│   ├── __init__.py
│   ├── leaving_groups.py
│   ├── pronucleophiles.py
│   ├── transfer_groups.py
│   └── electrophilic_centers.py
├── labels.py
├── mapping.py
├── bond_changes.py
├── reaction_classifier.py
├── features.py
├── api.py
└── tests/
    ├── test_context.py
    ├── test_leaving_groups.py
    ├── test_pronucleophiles.py
    ├── test_transfer_groups.py
    ├── test_electrophilic_centers.py
    ├── test_labels.py
    └── test_reactions.py
```

---

## 17. Suggested Python Models

Use dataclasses or Pydantic models.

```python
from typing import Literal
from pydantic import BaseModel, Field


SiteType = Literal[
    "leaving_group",
    "pronucleophile_XH",
    "transfer_group",
    "electrophilic_center",
    "aromatic_CH",
    "unsaturated_bond",
]


class ReactiveSite(BaseModel):
    schema_version: str = "1.0"
    site_id: str
    site_type: SiteType
    topology: Literal["edge", "atom", "center", "bond"]
    component_index: int
    atom_indices: list[int]
    bond_indices: list[int]
    canonical_signature: str
    chemist_label: str
    context_features: dict = Field(default_factory=dict)
    confidence: float = 1.0
    warnings: list[str] = Field(default_factory=list)


class ReactionClassification(BaseModel):
    schema_version: str = "1.0"
    reaction_class: str
    reaction_family: str | None = None
    reaction_label: str
    reactant_sites: list[ReactiveSite]
    product_sites: list[ReactiveSite]
    roles: dict[str, str]
    bond_changes: dict
    context_features: dict = Field(default_factory=dict)
    confidence: float
    warnings: list[str] = Field(default_factory=list)
```

---

## 18. Public API

## 18.1 Analyze a compound

```python
analyze_compound(
    smiles: str,
    *,
    include_context_features: bool = True,
) -> CompoundAnalysis
```

Example:

```python
result = analyze_compound("Nc1ccccc1")
```

Expected label:

```text
Ar–NH2
```

## 18.2 Analyze a reaction

```python
analyze_reaction(
    reaction_smiles: str,
    *,
    atom_map_if_missing: bool = True,
    include_context_features: bool = True,
) -> ReactionClassification
```

## 18.3 Detect only specific site families

```python
detect_sites(
    mol,
    site_types={
        "leaving_group",
        "pronucleophile_XH",
    },
) -> list[ReactiveSite]
```

## 18.4 Render a 2D motif

Optional:

```python
render_site_svg(
    site: ReactiveSite,
    *,
    width: int = 300,
    height: int = 180,
) -> str
```

Use RDKit dummy atoms and custom display labels for `R1`, `R2`, `Ar`, etc.

---

## 19. Detection Strategy

Use modular detectors rather than a single large SMARTS catalog.

Each detector may combine:

- small SMARTS patterns;
- atom properties;
- local graph traversal;
- context classification;
- explicit precedence rules.

Avoid a single giant SMARTS expression for free amines or other broad classes.

Recommended pattern:

```python
def detect_pronucleophile_sites(mol):
    sites = []

    for atom in mol.GetAtoms():
        if not is_supported_protic_center(atom):
            continue

        h_count = get_total_h_count(atom)
        if h_count < 1:
            continue

        contexts = classify_retained_contexts(mol, atom)
        derived_family = derive_protic_family(atom, contexts)
        signature = build_xh_signature(atom, h_count, contexts)
        label = render_xh_label(atom, h_count, contexts)

        sites.append(...)

    return sites
```

---

## 20. Confidence and Ambiguity

Potential ambiguity cases:

```text
multiple equivalent leaving groups
multiple N–H sites
hydrazines
tautomeric heterocycles
unmapped reaction SMILES
ambiguous atom mapping
salts and protonation states
unsupported transfer groups
```

Return warnings rather than silently forcing a result.

Example:

```json
{
  "confidence": 0.72,
  "warnings": [
    "Two candidate N–H sites are present.",
    "Reaction atom mapping was generated automatically.",
    "Site assignment depends on mapper output."
  ]
}
```

---

## 21. Non-Goals for Version 1

Do not attempt to fully support:

- every C–H activation;
- pericyclic reactions;
- radical rearrangements;
- photochemical mechanisms;
- all protecting-group transformations;
- all oxidation and reduction reactions;
- all carbonyl additions;
- complete mechanistic classification;
- complete IUPAC functional-group taxonomy.

The system should be extensible, but v1 should prioritize common synthetic handles.

---

## 22. Initial Acceptance Tests

## 22.1 Leaving groups

| SMILES | Expected signature | Expected label |
|---|---|---|
| `Brc1ccccc1` | `LG|Ar|Br` | `Ar–Br` |
| `Clc1ncccc1` | `LG|HeteroAr|Cl` | `HeteroAr–Cl` |
| `CCOS(=O)(=O)C` | `LG|Alkyl|OMs` when detected as mesylate | `R–OMs` |

## 22.2 Pronucleophiles

| SMILES | Expected signature | Expected label |
|---|---|---|
| `CN` | `XH|N|H2|Alkyl` | `R–NH2` |
| `Nc1ccccc1` | `XH|N|H2|Ar` | `Ar–NH2` |
| `CCNCC` | `XH|N|H1|Alkyl,Alkyl` | `R1R2–NH` |
| `CNc1ccccc1` | `XH|N|H1|Ar,Alkyl` | `Ar–NH–R` |
| `CC(=O)NC` | `XH|N|H1|C(O)R,Alkyl` | `RC(O)–NHR` |
| `CS(=O)(=O)NC` | `XH|N|H1|SO2R,Alkyl` | `RSO2–NHR` |
| `Oc1ccccc1` | `XH|O|H1|Ar` | `Ar–OH` |
| `CCS` | `XH|S|H1|Alkyl` | `R–SH` |
| `CNN` | two N–H sites | hydrazine-site labels |

## 22.3 Transfer groups

| SMILES | Expected signature | Expected label |
|---|---|---|
| `OB(O)c1ccccc1` | `TM|Ar|B(OH)2` | `Ar–B(OH)2` |
| aryl-BPin test | `TM|Ar|BPin` | `Ar–BPin` |

## 22.4 Electrophilic centers

| SMILES | Expected signature | Expected label |
|---|---|---|
| `CC(=O)O` | `EC|Acyl|Alkyl|OH|latent` | `R–C(O)OH` |
| `CC(=O)Cl` | `EC|Acyl|Alkyl|Cl|activated` | `R–C(O)Cl` |
| `O=S(=O)(Cl)c1ccccc1` | `EC|Sulfonyl|Ar|Cl|activated` | `Ar–S(O)2Cl` |

## 22.5 Reaction tests

### Buchwald C–N

Input:

```text
Brc1ccccc1C(=O)C.Nc1ccccc1>>CC(=O)c1ccccc1Nc1ccccc1
```

Expected:

```text
Ar–Br + Ar–NH2 → Ar1Ar2–NH
```

### Suzuki

Expected minimal pattern:

```text
Ar–Br + Ar–B(OH)2 → Ar–Ar
```

### Amide formation

Expected minimal pattern:

```text
R–C(O)OH + R–NH2 → R–C(O)NHR
```

---

## 23. Development Milestones

### Milestone 1: Compound-site detection

Implement:

- molecule standardization;
- context classifier;
- leaving-group detection;
- N/O/S–H detection;
- organoboron detection;
- acyl-center detection;
- canonical signatures;
- chemist labels;
- unit tests.

### Milestone 2: Reaction mapping and bond changes

Implement:

- reaction parsing;
- mapped-reaction support;
- optional mapper adapter;
- broken/formed bond extraction;
- site-to-event association.

### Milestone 3: Reaction grammar classification

Implement:

- C–N;
- C–O;
- C–S;
- Suzuki;
- other transfer couplings;
- Chan–Lam;
- Sonogashira;
- amide formation;
- ester formation;
- sulfonamide formation.

### Milestone 4: Context overlays

Implement high-value steric and electronic descriptors.

### Milestone 5: Condition-recommendation integration

Produce hierarchical reaction feature objects suitable for rule-based retrieval, nearest-neighbor retrieval, or machine-learning models.

---

## 24. Final Architectural Rule

The core system should answer three questions:

```text
1. What reactive handles are present?
2. Which handles participated in the reaction?
3. What local contexts modify their reactivity?
```

The minimal reaction grammar defines the transformation:

```text
Ar–X + N–H/O–H/S–H → Ar–N/O/S
```

The site signatures define the substrates:

```text
LG|Ar|Br
XH|N|H2|Ar
```

The context overlays support condition selection:

```text
ortho_acyl
primary_arylamine
sterically_hindered
heteroaryl_coordination_site
```

This separation should remain the central design constraint:

> **Core reactive requirements define the reaction.  
> Local context explains substrate behavior.  
> Exact molecular structure remains available but is not encoded into the core taxonomy.**

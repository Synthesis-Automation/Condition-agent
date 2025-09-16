# Plan: Automatic Role Identification → Type‑Specific Enrichment → Taxonomy Persistence

**Goal:** Given a reagent with a known CAS number, (1) assign its role/type (ligand, metal precursor, base, solvent, additive, etc.), (2) fetch/compute type‑specific facts for the relevant taxonomy, and (3) persist the result with provenance to the correct taxonomy files and the knowledge graph (KG).

---

## Architecture (pipeline per reagent)

**ingest (CAS)** → **resolve identity** → **role inference** → **type‑specific enrichment** → **validate** → **persist to taxonomies** → **emit KG edges** → **(optional) human review**

---

## 1) Identify the reagent role/type

### 1.1 Identity resolution (deterministic first)

- **Normalize CAS**: regex check, strip spaces, zero‑pad segments.
- **Lookups (ordered)**: internal registry → cached resolver DB → (when online) PubChem/CompTox/Cactus → vendor dumps.
- **Persist**: `cas`, `inchikey`, `smiles`, `name`, `aliases[]`, `sources[]`, `resolution_confidence`.

### 1.2 Signals for role inference (priority order)

**A. Structure/name rules (hard)**

- d‑block metal + anions (`Cl`, `Br`, `OAc`, `acac`, `BF4`, `PF6`, `OTf`) → `metal_precursor`.
- P(III) trivalent phosphorus or known phosphine name (`XPhos`, `SPhos`, `tBuBrettPhos`, `PPh3`, …) → `ligand:phosphine`.
- `bpy|phen|terpy|dtbbpy` in name or N–N chelate SMARTS → `ligand:diimine`.
- `IMes|IPr|SIPr|SIMes|PEPPSI` or imidazolium salt motif → `ligand:NHC (or precursor)`.
- Base dictionaries: alkoxides, amides, carbonates, phosphates, hydrides → `base`.
- Solvent whitelist (DN/AN/ET(30) table) **and** amount in mL/µL is largest → `solvent`.
- Coupling reagents dictionary (HATU, DCC, EDC, T3P, BTFFH…) → `coupling_reagent`.
- Supported catalysts (`Pd/C`, `Pd(OH)₂/C`, `SiliaCat…`) → `supported_precatalyst`.

**B. Context signals (from reaction row, if available)**

- `mol% ≤ 20` & co‑used with a metal salt → **ligand**.
- `eq ≥ 0.5` & not volumetric → **base**.
- Appears in “solvent” column or the **largest volume** → **solvent**.
- Co‑occurs strongly with _known cores/families_ (KG lift/support) → boost that role.

**C. ML disambiguation (multi‑label)**

- Features: name text, Morgan FP, descriptors (donors, charge, metal flags), context (eq, mol%, family), KG stats.
- Model: one‑vs‑rest gradient boosting + **conformal prediction** to abstain OOD.

**Output example**

```json
{
  "primary_role": "ligand",
  "subrole": "phosphine_monodentate",
  "confidence": 0.89,
  "evidence": {
    "rules": ["R_P_MONO"],
    "context": ["mol%≤10 with Pd"],
    "stats": { "lift_BH_N_1.9": true }
  }
}
```

---

## 2) Type‑specific enrichment (taxonomy fields)

### 2.1 Ligands

- **Structural**: donor type (P/N/C‑carbene), denticity, chelate ring size, **bite angle** (if bidentate), **cone angle** (Tolman; table lookup), **%V_bur** (Sterimol/%V_bur pipeline), heteroatoms, MW.
- **Families & aliases**: Buchwald class (XPhos/SPhos/etc.), biaryl vs chelating diphosphine vs NHC/diimine; common abbreviations.
- **Electronic**: approximate TEP (if available), HOMO/LUMO (if you keep a computed table).
- **Handling**: air/moisture sensitivity, storage notes.
- **Usage priors**: reaction families where prevalent; **requires metal** (always true).
- **Cross‑links**: known pre‑ligated complexes (e.g., `PdCl2(dppf)`).

### 2.2 Metal precursors

- **Metal/oxidation state**, counter‑anions, **pre‑ligated?**, **requires_added_ligand?**, activation protocol hint (e.g., Pd(II)‑phosphine needs base; `Pd₂(dba)₃` → add phosphine).
- **Solubility class** (qualitative), air/moisture sensitivity, support (if any).
- **Preferred families** and typical ligand partners (from KG co‑occurrence).

### 2.3 Bases

- **Category**: alkoxide / amide / carbonate / phosphate / fluoride / hydride / organic superbase.
- **Basicity**: conjugate‑acid **pKₐ (DMSO)** and water; Brønsted vs Lewis; **nucleophilicity class**.
- **Coordination tendency** (binds metals low/moderate/high), **hygroscopicity**, **bp/mp**, **solubility hints**.
- **Compatibility flags**: “avoid with acid chlorides in DMSO”, “forms hemiaminal with aldehydes”, etc.

### 2.4 Solvents

- **Physical**: DN/AN, **E_T(30)/E_T(N)**, dielectric, Hansen (δD/δP/δH), viscosity, bp, mp, water miscibility.
- **Protic/aprotic**, coordinating/non‑coordinating, peroxide risk, **ESG/safety** class.

### coupling reagents

- e.g. DCC for amine formation

### 2.5 Additives // salts

- **Function** (halide scavenger, fluoride source, phase‑transfer, oxidant/reductant).
- **Strength/scale**: typical eq range; common pairings.

**Populating the fields**

- Compute from `SMILES` (RDKit descriptors, donor counts, metal presence).
- Join to **reference tables** (cone angles, %V_bur, pKₐ, DN/AN, E_T(30)).
- Attach **provenance** for each field: `computed|table|literature` with `source_id`.

---

## 3) Persist to the corresponding compound taxonomy

### 3.1 Data contracts (pydantic)

```python
class BaseEntry(BaseModel):
    cas: str; inchikey: str|None; smiles: str|None; name: str; aliases: list[str]=[]
    base_class: str                    # "alkoxide","amide","carbonate",...
    pka_dmso: float|None = None
    pka_water: float|None = None
    nucleophilicity: str|None = None
    coordination_tendency: str|None = None  # "low|moderate|high"
    hygroscopic: bool|None = None
    compat_flags: list[str] = []
    sources: list[dict] = []
    status: str = "provisional"; confidence: float = 0.0

class LigandEntry(BaseModel):
    cas: str; inchikey: str|None; smiles: str|None; name: str; aliases: list[str]=[]
    ligand_class: str                  # "phosphine_mono","diphosphine","NHC","diimine"
    denticity: int
    donor_type: str                    # "P","N","Ccarbene"
    bite_angle: float|None = None
    cone_angle: float|None = None
    vbur_percent: float|None = None
    families: list[str] = []           # "Buchwald","Josiphos","Xantphos"
    handling: list[str] = []
    sources: list[dict] = []
    status: str = "provisional"; confidence: float = 0.0

class MetalPrecursorEntry(BaseModel):
    cas: str; inchikey: str|None; smiles: str|None; name: str; aliases: list[str]=[]
    metal: str; ox_state: int|None = None
    anions: list[str] = []
    preligated: bool = False
    requires_added_ligand: bool|None = None
    support: str|None = None           # "C","silica",...
    handling: list[str] = []
    partner_ligands: list[str] = []    # common ligand names
    sources: list[dict] = []
    status: str = "provisional"; confidence: float = 0.0

class RegistryEntry(BaseModel):
    cas: str; inchikey: str|None; smiles: str|None; name: str; aliases: list[str]=[]
    primary_role: str                  # "ligand","metal_precursor","base","solvent","additive"
    subrole: str|None = None
    sub_taxonomy_ref: str|None = None  # e.g., "ligands:CAS"
    evidence: dict = {}
    confidence: float = 0.0
    sources: list[dict] = []
    status: str = "provisional"
```

### 3.2 Stores & file layout

- `registry/registry.jsonl` (all compounds, 1 row per CAS)
- `taxonomies/ligands.jsonl`
- `taxonomies/bases.jsonl`
- `taxonomies/metal_precursors.jsonl`
- `taxonomies/solvents.jsonl`
- `taxonomies/additives.jsonl`

All as **JSONL or Parquet**, versioned with CI schema checks.

### 3.3 Merge policy

- If CAS exists → **update** fields only if:
  - the new value has higher confidence, or
  - the field was previously missing.
- Keep **provenance & diffs** (append new `sources[]`, stamp `data_version`).
- Conflict (e.g., two distinct InChIKeys for same CAS) → set `status="needs_review"` and open a review ticket; do **not** overwrite.

### 3.4 KG emission

On write, create/update nodes and relations:

- `(:Compound {cas})-[:HAS_ROLE {role, subrole, confidence}]->(:ReactionFamily)` (aggregated later)
- If ligand: `(:Ligand {cas})-[:IS]->(:Compound {cas})`
- If metal precursor: `(:MetalPrecursor {cas})-[:IS]->(:Compound {cas})`
- Link **co‑occurrence** edges to cores/tails from recent datasets with weights.

---

## ETL pipeline (modules & flow)

1. `resolve.py`

   - `normalize_cas`, `resolve_identity(cas) → Identity` (with caching)

2. `role_rules.yml` + `role_infer.py`

   - Stage‑1 rules (SMARTS + regex + lists)
   - Stage‑2 heuristic scorer
   - Stage‑3 ML (optional) + conformal wrapper

3. `enrich_<role>.py`

   - `enrich_ligand(identity) → LigandEntry`
   - `enrich_base(identity) → BaseEntry`
   - `enrich_metal_precursor(identity) → MetalPrecursorEntry`
   - Uses RDKit + your reference tables (cone angle, %V_bur, pKₐ, solvent properties)

4. `validate.py`

   - Pydantic/SHACL checks; mandatory fields per role

5. `persist.py`

   - Upserts into `registry.jsonl` and the specific `taxonomies/*.jsonl`
   - Writes provenance; bumps dataset/model versions

6. `kg_writer.py`

   - Upserts nodes/edges in Neo4j (or your graph store)

7. `qa_eval.py`
   - Batch metrics (precision@1, abstention rate), sample reports for human QA

**Batch driver sketch**

```python
for cas in new_cas_list:
    ident = resolve_identity(cas)
    role = infer_role(ident, context=row_context)   # returns primary_role, subrole, confidence, evidence
    if role.primary_role == "ligand":
        entry = enrich_ligand(ident); entry.confidence = role.confidence
        persist_registry_and_taxonomy(ident, entry, role)
    elif role.primary_role == "base":
        entry = enrich_base(ident); ...
    elif role.primary_role == "metal_precursor":
        entry = enrich_metal_precursor(ident); ...
    else:
        entry = enrich_other(ident); ...
    write_kg_edges(ident, role, cooccurrence_from_dataset)
```

---

## Guardrails & promotion

- **Abstain** if confidence < 0.5 → stash in `staging/` for review.
- **Promote** `provisional → curated` when: consistent role in ≥N precedents, no conflicts for a period, passes spot QA.
- **Auto‑learn**: retrain role classifier monthly using curated entries; refresh co‑occurrence stats.

---

## Deliverables checklist

- [ ] `role_rules.yml` (≈60 starter rules: Pd/Ni/Cu + phosphines/NHC/bpy + common bases/solvents/additives)
- [ ] `role_infer.py` (rules → heuristic → ML + conformal)
- [ ] `schemas.py` (pydantic models above)
- [ ] `resolve.py` (CAS→Identity with cache)
- [ ] `enrich_ligand.py`, `enrich_base.py`, `enrich_metal_precursor.py`, `enrich_solvent.py`
- [ ] `persist.py` (JSONL/Parquet upserts, provenance merge)
- [ ] `kg_writer.py` (Cypher upserts for nodes/edges)
- [ ] `validate.py` + CI tests (pytest + SHACL shapes)
- [ ] `qa_eval.py` (metrics + sampled confusion/abstention)
- [ ] Reference tables: cone angles, %V_bur, pKₐ(DMSO/water), solvent DN/AN/ET(30), Hansen

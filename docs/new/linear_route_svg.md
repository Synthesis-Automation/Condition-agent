# Three-column linear reaction schemes

Implemented 2026-10-03. This is presentation of supplied structures, not a
reaction-feasibility assessment or a chemistry-release gate.

The reusable `visualization.render_route_scheme_svg` returns standalone SVG bytes.
It lays out the sequence `starting material, arrow, intermediate, arrow, ...`
in three columns, reading left to right and then down. Intermediates appear once;
the final row may be incomplete. Compound names and numbers are omitted.

## Usage

```python
from pathlib import Path
from visualization import render_route_scheme_svg

svg = render_route_scheme_svg([
    "CCO>>CC=O",
    "CC=O.N>>CCN",
    "CCN.CC(=O)Cl>>CCNC(C)=O",
])
Path("route.svg").write_bytes(svg)
```

These are illustrative structures, not an experimentally validated synthesis.

```powershell
python -m visualization.route_cli examples/linear_route_scheme.json results/route_scheme/linear_route.svg
```

CLI input is either a JSON list of reaction SMILES or an object with `steps`,
optional `title`, optional zero-based `main_reactant_index`, and optional boolean
`reagents_only` (default false for explicitly annotated CLI schemes). Each step is a
string or an object:

```json
{
  "reaction_smiles": "CC=O.N>>CCN",
  "basis": "proposed",
  "reagents": [{"text": "Supplied reagent name", "basis": "proposed"}],
  "above": [{"text": "Supplied reagent text", "basis": "proposed"}],
  "below": [{"text": "Supplied solvent, temperature and time", "basis": "proposed"}],
  "yield_info": null,
  "product_index": null
}
```

Annotations may also be plain strings, with basis `supplied`. The Python
equivalent uses immutable `RouteSchemeStep` and `SchemeAnnotation` objects.
Nothing infers conditions, yield, compound names, or evidence status. Original
reaction SMILES and annotation attribution are retained in SVG JSON metadata.

## Geometry and component policy

- Versioned geometry: `visualization/definitions/route_scheme.v1.json`, schema 1.2.
- Arrows occupy 80% of their block width, centered under the annotations. This
  shortens their visible span without changing molecular scale. Arrow blocks start
  at 160 px; adjacent blocks have a 24 px gap. Long reagent names wrap at 20 characters.
- Three content-sized slots per row. Each slot fits its own block, so a large
  molecule does not widen an adjacent arrow. Row heights adapt
  to their contents. Structures reuse the existing `web_consistent` bond scale
  and colors; large structures increase canvas dimensions rather than shrink.
- The first block contains all first-step reactants with plus signs. Supply
  `main_reactant_index` to move the remaining reactants above the first arrow.
  Subsequent carried intermediates are determined by exact connectivity.
- Additional reactants and supplied middle-section agents are drawn above the
  arrow. Text above/below and yield are optional. No annotation is truncated.
  With `reagents_only=True`, only explicit `reagents` annotations are displayed;
  above/below procedure text and yield remain in SVG metadata. No reagent names
  are extracted or guessed from prose.
- Additional products are retained below their arrow under “Other products.”
  A unique next-step match selects the route product. Ambiguous selections,
  including a multiproduct last step, require `product_index`.
- Dot-separated fragments are treated as component occurrences. Salt grouping
  is not inferred. Specify product selection when needed; no counterion is
  silently removed. Stereo, protonation, isotope and tautomer differences are
  not automatically reconciled.

`reactive_taxonomy.reaction_sequence.connect_reaction_sequence` owns structural
matching. Canonical isomeric identities ignore map labels and SMILES traversal
order, but preserve charge, isotope and stereochemistry. A selected product must
match exactly one reactant of the next step. This does not validate atom mapping,
reaction edits, or chemical feasibility. It never forces a disconnected link.

## Workspace integration

Structured scientific answers use this renderer for declared linear routes.
The complete scheme precedes per-step evidence, rationale, limitations, and
experimental details. A visible caption explains the row reading order. Only
explicit reagent names and supplied partner structures appear on route arrows.
Full conditions and yields remain in experimental details; reagent attribution
is shown beside the step explanation and citations in the details. Structured
answer steps add an attributed `reagents` list (default empty for saved v2 answers).
Supporting source reaction drawings are expandable while observed conditions,
source identity and cautions remain visible.
Individual step drawings remain available in the response. Agent presentation
instructions request complete route objects rather than repeated target schemes
or prose-only SMILES sequences.
Branched dependencies, invalid structures or ambiguous connections retain the
dependency overview and individual reaction drawings, with a display warning.
The browser preserves the scheme scale using horizontal scrolling when needed.

Restart the workspace server and reload saved conversations to use the updated
presentation. Old procedure paragraphs no longer appear on route arrows; reagent
names must be supplied explicitly in new or revised answers. Existing saved answers
remain readable without inferred reagent labels. API routes are unchanged; route
presentation objects add `drawing_status`, `scheme_width` when drawn, and
`drawing_warning` on fallback. No corpus conversion or index rebuild is required.

This initial layout does not align common scaffolds across steps, abbreviate
substructures, combine operations under one arrow, or draw convergent routes as
a continuous molecular scheme. These can be added without replacing the public
route renderer.

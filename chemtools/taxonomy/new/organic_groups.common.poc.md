# Common Organic Chemistry “Groups” for Attachment-Point SMARTS Templating (POC)

This document defines **chemist-style group labels** like **Ar-**, **Vinyl-**, **-CHO**, **-CN**, **-B(OH)2**, etc., to support
**attachment-point SMARTS templating**.

## Design idea

We separate “groups” into:

- **Context groups** (left side): e.g., **Ar-**, **R-**, **Vinyl-**  
  These represent the *attachment environment* (what you are attaching *from*).
- **Core groups** (right side): e.g., **-CHO**, **-CN**, **-NO2**, **-Bpin**  
  These represent the *functional group core* (what you are attaching *to*).
- **Tail groups**: e.g., **-OTf**, **-OTs**, **-OMs**  
  These are used via an **O-linker** template (Ar–O–SO2R).

Each group has:

- `id` (stable, code-friendly)
- `name` (chemist-facing label, e.g. “Ar-” or “-CHO”)
- `kind`: `context` / `core` / `tail`
- `smarts`: the SMARTS fragment (often includes an attachment label like `:1` or `:2`)

## Templates

- `single_bond`: `{A}{B}`
  - Use for **Ar–X**, **Ar–CHO**, **Ar–CN**, **Ar–B(OH)2**, etc.
- `via_oxygen`: `{A}O{B}`
  - Use for **Ar–OTf**, **Ar–OTs**, **Ar–OMs**, etc. (pseudohalides)

## Examples (generated SMARTS)

- **Aromatic aldehyde (Ar–CHO)**  
  `Ar-` + `-CHO` via `single_bond` → `[c:1][CX3H1:2](=O)`

- **Aryl bromide (Ar–Br)**  
  `Ar-` + `-Br` via `single_bond` → `[c:1][Br]`

- **Aryl triflate (Ar–OTf)**  
  `Ar-` + `-OTf` via `via_oxygen` → `[c:1]OS(=O)(=O)C(F)(F)F`

- **Aryl boronic acid (Ar–B(OH)2)**  
  `Ar-` + `-B(OH)2` via `single_bond` → `[c:1][B:2](O)O`

## Notes / limitations

- The chemist label **Ar-** is implemented as `[c:1]` (aryl carbon attachment).
  If you need “any aromatic atom attachment”, use **Arom-** (`[a:1]`).
- Some functional classes (e.g., ketones, epoxides, rearrangement motifs) are better handled as **atomic SMARTS features**
  rather than generated from two-piece templates.
- Tail groups like **-OTf** are not meant to be used alone; they are attached through oxygen with `via_oxygen`.

## Machine-readable file

- `organic_groups.common.poc.json`

# Reaction Drawing POC

A local, browser-only proof of concept for drawing a chemical reaction with
Ketcher and exporting the drawing as reaction SMILES.

## Run locally

```powershell
npm install
npm run dev
```

Open `http://127.0.0.1:5173/`.

## Current scope

- Draw or edit a reaction in Ketcher.
- Load a small example reaction.
- Export `reactants>>products` reaction SMILES.
- Copy the generated SMILES.
- Perform a lightweight format check in the browser.

The POC does not yet call RDKit or the repository's recommendation packages.
Those integrations belong in a thin local Python API in a later phase.

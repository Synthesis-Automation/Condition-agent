# HTE Recommender API (FastAPI)

This service exposes a JSON endpoint that returns the same output format as the HTE GUI export.

## Start the service (PowerShell)

```powershell
# From repo root
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt

# Run the API
uvicorn app.main:app --reload --port 8000
```

Open docs at `http://127.0.0.1:8000/docs`.

## Test the endpoint (PowerShell)

```powershell
$payload = @{
  input = @{
    reaction_smiles = "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
    reaction_type_filter = "C_N_Coupling"
    catalyst_filter = "Pd"
  }
} | ConvertTo-Json -Depth 6


$payload = @{
  input = @{
    reaction_smiles = "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
  }
} | ConvertTo-Json -Depth 6


$response = Invoke-RestMethod `
  -Method Post `
  -Uri "http://127.0.0.1:8000/hte/recommendations" `
  -ContentType "application/json" `
  -Body $payload

# Print the first recommendation steps (adds + conditions)
$first = $response.recommended_conditions | Select-Object -First 1
if ($null -eq $first) {
  "No recommendations returned."
} else {
  $first.steps |
    ForEach-Object {
      if ($_.action -eq "add") {
        $chem = $_.chemical
        "{0}: {1} (role={2}, eq={3})" -f $_.action, $chem.name, $chem.role, $chem.equivalents
      } else {
        "{0}: {1} at {2} C for {3} h" -f $_.action, $_.name, $_.temperature_C, $_.time_h
      }
    }
}

# Or print the full JSON (pretty)
# $response | ConvertTo-Json -Depth 10
```

## Test with HTML (browser)

Open `docs/HTE_RECOMMENDER_TESTER.html` in a browser and click **Send Request**.

Note: If the API runs on a different host/port, update the endpoint field in the page.

Expected: a JSON payload with `meta`, `input`, `detection`, and `recommended_conditions`, where
`recommended_conditions[0].steps` interleaves `add` actions and `condition` actions.

## Default equivalents and solvent concentration

When the recommender output does not include `equivalents`, use these defaults:

- First reactant: `1.0`
- Other reactants: `1.5`
- Catalyst or ligand: `0.05`
- Base or acid: `2.0`
- Other reagents: `1.5`

When the output does not include solvent concentration, use:

- `concentration_M`: `0.2`

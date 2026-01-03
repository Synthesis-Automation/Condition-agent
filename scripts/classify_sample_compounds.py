import csv
import sys
from pathlib import Path

# Add the current directory to sys.path to import chemtools
sys.path.append(str(Path.cwd()))

try:
    from chemtools.featurizers.molecule import featurize_molecule
except ImportError:
    print("Error: Could not import featurize_molecule from chemtools.featurizers.molecule")
    sys.exit(1)

input_file = "examples/sample_compounds.csv"
output_file = "examples/sample_compounds_classified.csv"

if not Path(input_file).exists():
    print(f"Error: Input file {input_file} not found.")
    sys.exit(1)

# Options for featurization
options = {
    "discovery_mode": True,
    "include_ar_h": True
}

results = []

print(f"Classifying compounds from {input_file}...")

with open(input_file, mode="r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames + ["detected_motifs", "undocumented_motifs"]
    
    rows = list(reader)
    for i, row in enumerate(rows):
        if i % 20 == 0:
            print(f"Processing {i}/{len(rows)}...")
        smiles = row["smiles"]
        try:
            payload = featurize_molecule(smiles, options=options)
            motifs = payload.get("motifs", [])
            
            # Extract names and flag undocumented
            documented = []
            undocumented = []
            
            for m in motifs:
                cid = m.get("compound_id", "unknown")
                if m.get("undocumented", False):
                    undocumented.append(cid)
                else:
                    documented.append(cid)
            
            # Unique the lists for the summary CSV to reduce noise
            row["detected_motifs"] = "; ".join(sorted(list(set(documented))))
            row["undocumented_motifs"] = "; ".join(sorted(list(set(undocumented))))
            
        except Exception as e:
            row["detected_motifs"] = f"ERROR: {str(e)}"
            row["undocumented_motifs"] = ""
            
        results.append(row)

with open(output_file, mode="w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)

print(f"Successfully classified {len(results)} compounds and saved to {output_file}")

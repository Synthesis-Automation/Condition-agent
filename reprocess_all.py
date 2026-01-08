import os
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Import the processing function
from scripts.A_convert_to_hte_format import process_reaction_dataset

def main():
    dataset_dir = PROJECT_ROOT / "data" / "reaction_dataset"
    output_dir = PROJECT_ROOT / "data" / "HTE_db"
    
    # We will reprocess the datasets that already have canonical CSVs
    # Plus any other key ones we want to ensure are updated
    target_stems = [
        "C_N_Coupling",
        "Sulfonamide_formation",
        "SNAr-CN"
    ]
    
    for stem in target_stems:
        input_path = dataset_dir / f"{stem}.jsonl"
        output_path = output_dir / f"{stem}_canonical.csv"
        
        if input_path.exists():
            print(f"\n>>> Reprocessing {stem}...")
            try:
                process_reaction_dataset(
                    str(input_path), 
                    str(output_path), 
                    drop_no_catalyst=False  # Keep all for now to see transformation variety
                )
            except Exception as e:
                print(f"Error processing {stem}: {e}")
        else:
            print(f"Skipping {stem}, input file not found.")

if __name__ == "__main__":
    main()

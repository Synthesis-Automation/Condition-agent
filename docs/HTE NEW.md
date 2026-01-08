I have completed the requested improvements to the motif detection system and built the recommendation infrastructure for C-S couplings.

Key Enhancements
Resolved Methyl Specificity: Added explicit methyl-based motifs (e.g., CH3-NH2, CH3-Cl, CH3-OH, CH3-CN) to the organic_compounds.v1.3.json registry. This ensures that molecules like methylamine are detected with high specificity rather than being lumped into generic Alkyl categories.
Structural Complexity Prioritization: Implemented an automated tie-breaking system in motif_registry.py. It calculates a "complexity score" for every SMARTS pattern based on atom weights and bond density. Specific motifs like CH3-NH2 now naturally outrank generic ones like Alkyl-NH2 because they represent a more precise chemical environment.
C-S Coupling Dataset Integration:
Converted the 11,420-reaction C_S_Coupling.jsonl into a standardized C_S_Coupling_canonical.csv.
Updated A_convert_to_hte_format.py to handle system paths and large-scale data processing correctly.
Unified Recommendation Index:
Enhanced index_builder.py to extract high-fidelity tags and SMILES from HTE-style CSVs.
Successfully built a unified index of 418,855 records, providing structural and tag-based recommendations across 9 major reaction families.
Verification Results
I tested the system with a methylamine C-N coupling and a C-S coupling query:

Methylamine Reaction: Correctly detected CH3-NH2 as a specific motif and prioritized it over generic Any-NH2.
C-S Coupling: Successfully retrieved high-yield conditions (e.g., NiBr2/Zn/K2CO3 and Cu2O/Py) from the newly indexed dataset.
Discovery Mode: Confirmed that the "Discovery Mode" (enabled by default in the CLI) correctly identifies undocumented scaffolds and substituents without crashing.
The system is now fully equipped to provide expert-level condition recommendations for Carbon-Sulfur cross

-couplings.

Gemini 3 Flash (Preview) • 0.33x
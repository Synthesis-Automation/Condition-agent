pls follow the example of C_N coupling rule json, make json rule for ring closing metathesis



it could inlcude pd, cu and other metal catalyzed reactions. part of data could be use table 5 in the attached pdf file. 






for Ullman C-N bond formation reaction (copper catalyzed), please list several most common conditions and when they should be used (e.g. chiral acids, reaction of steric substrate)
Reductive amination — Convert carbonyls + amines → secondary/tertiary amines in one pot.

for SnAr reaction, please find several most common conditions and when they should be used (e.g. steric substrate). and generate a rule json file using the attached C-N formation rule as a template. 


Instruction (single-shot, file + JSON)



Please work on "SnAr" reaction

You will generate one reaction rule file in my schema and save it as a UTF-8 JSON file for download.
Reply format:

First line: a single markdown link to the saved file, exactly [Download JSON](sandbox:/mnt/data/<<SAFE_FILENAME>>.json)

Blank line

The raw JSON object (no code fences, no extra prose).

Schema to replicate exactly (same as my rule files):
Top-level keys: name, reaction_type, version (YYYY-MM-DD), evaluation (multi-line string), applies_if, default_rule, base_rules, modifiers.
applies_if uses all / any arrays of boolean feature tokens.
Each item in base_rules has name, id, description, optional reactant_features (all / any / not), and conditions (domain keys).
modifiers is an array of { id, when: [tokens], suggest }. Tokens may include route:<method> and symptom:<observation>.
Key naming: snake_case; include units in key names where helpful (temperature_C, pressure_bar, solvent_volume_mL_per_mmol, catalyst_loading_molpct).
Numbers are numeric; ranges may be strings like "40–60". IDs must be unique (BR_… for base rules, MOD_… for modifiers).

Target reaction (fill these):

<REACTION_NAME> = <<PUT NAME HERE>>

<REACTION_DESCRIPTION> = <<1–2 line overview>>

Token vocabulary (fill):
Provide a comma-separated list of snake_case boolean tokens relevant to this class:
<<TOKEN_VOCAB>>
(Examples to adapt: aldehyde_present, ketone_present, aryl_chloride_present, heteroaryl_electrophile_present, steric_hindrance_high, electron_poor_ring_present, hydrogenation_sensitive_present, polarity_high, base_sensitive_present, scale_up, route:H2, route:photoredox, symptom:slow reaction.)

applies_if (fill):

all:  <<APPLIES_IF_ALL>>   # tokens that must all be true
any:  <<APPLIES_IF_ANY>>   # at least one must be true


default_rule (requirements):
A conservative, broadly applicable fallback with clean conditions and concise notes on setup/order-of-addition and common pitfalls.

base_rules (6–10 buckets, easy→hard) (fill):
Provide diverse rules covering substrate classes, difficulty, and typical lab choices. Include realistic conditions (equiv/loadings/temps/time/solvent/additives).
Planned buckets:
<<BASE_RULE_BUCKETS>>

conditions keys to choose from (adapt):
<<CONDITION_KEYS>>
(e.g., catalyst, ligand, metal_precatalyst, reducing_agent, oxidant, base, acid_or_buffer, additives, lewis_acid, hydrogen_source, pressure_bar, light_source_nm, solvent, co_solvent, water_equiv, temperature_C, time_h, solvent_volume_mL_per_mmol, atmosphere, water_control, notes)

modifiers (include 6–12):
Cross-cutting tips triggered by tokens/symptoms (e.g., water removal, pH windows, protecting-group cautions, scale-up swaps, deactivation fixes).

evaluation:
Short, actionable playbook: when each bucket applies; key tradeoffs (sterics/electronics/polarity); FG cautions; optimization levers. Include safety notes where relevant.

Output & file-save contract (must do):

Generate the JSON object in memory.

Save it as UTF-8 to /mnt/data/<<SAFE_FILENAME>>.json.

Reply with exactly:

Line 1: [Download JSON](sandbox:/mnt/data/<<SAFE_FILENAME>>.json)

Blank line

The raw JSON object (no code fences, no extra text).
If a code-execution/file tool is unavailable, skip saving and just output the raw JSON object.
Use today’s date for version.

Fill for this run:

<<SAFE_FILENAME>> = <<your_file_basename>>

<<TOKEN_VOCAB>> = <<comma-separated tokens>>

<<APPLIES_IF_ALL>> = [<<tokens>>]

<<APPLIES_IF_ANY>> = [<<tokens>>]

<<BASE_RULE_BUCKETS>> = <<numbered list with one-line “when to use”>>

<<CONDITION_KEYS>> = <<keys relevant to this reaction>>

Begin.
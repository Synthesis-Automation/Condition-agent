1. **Family Information**:

   - **family**: 'Mitsunobu'
   - **metadata_id**: auto
   - **metadata_name**: auto
   - **metadata_version**: auto

2. **Creation Date**:

   - **created_date**: The creation date in YYYY-MM-DD format (optional, defaults to today).

3. **Status and Tags**:

   - **status**: auto.
   - **tags**: auto.

4. **Reference Reactions**:

   - **reference_reactions**: ['CCO.OC1=CC=CC=C1>>CCOC2=CC=CC=C2'].

5. **Protocol Text**:
   - **protocol_text**: A solution of benzyl alcohol (0.200 g, 1.85 mmol), 4-hydroxybenzaldehyde (0.226 g, 1.85 mmol), and PPh3(0.582 g, 2.22 mmol) was stirred in dry THF (20 mL) at 0 °C under a nitrogen atmosphere. To this mixture was added dropwise DIAD (0.44 mL, 2.22 mmol) over a period of 5 min, and the reaction was monitored by TLC. After complete disappearance of starting material (1 h), the solvent was evaporated under reduced pressure and the resulting oil purified by flash column chromatography (hexane/AcOEt, 8/2). Phenyl ether (0.297 g, 76%) was finally obtained as a white powder after precipitation from CH2Cl2/petroleum ether.

## Desktop Rule-Builder UI

- Launch the Qt interface via `python -m chem_assistant.gui.app`.
- The window mirrors Copilot-style experiences:
  - Left panel: LangGraph chat history driven by `ChemToolsAgent`.
  - Right panel: constraint summary, cache stats, and a system log.
  - Bottom controls: send messages, clear chat, manage constraints/cache, list tools, open the rule editor, or run the autofill workflow.
- **Rule Builder Editor** opens a form-based interface for metadata, references, scope, applies-if, default rule, base rules, and modifiers. You can load/save any `data/rule_db_v2/*.json` file, validate before saving, and stay inside the GUI.
- **Autofill Draft** wraps `rule_builder_autofill_tool`: provide metadata + protocol text, review the generated JSON/validation issues, then push the draft into the editor for polishing.
- Constraint, cache, and tool-list functionality match the CLI commands, so you can manage the entire chem_assistant workflow without the terminal.

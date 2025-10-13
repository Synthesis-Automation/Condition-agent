# CLI Visual Demo: Reaction Type Routing

## 🎯 Quick Start

```bash
# Test mode - validate input without API call
python app/cli_recommend.py --test
```

## 📋 Example Flow

### Step 1: Launch CLI
```
══════════════════════════════════════════════════════════════════════
  🧪 CHEMISTRY CONDITION RECOMMENDATION CLI
  Powered by LLM-assisted natural language parsing
  ⚠️  TEST MODE: Will stop before actual API submission
══════════════════════════════════════════════════════════════════════
```

### Step 2: Input Reaction
```
Please provide your reaction details:

Reaction SMILES (reactants>>products format):
Example: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1
> Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1
```

### Step 3: Describe Requirements
```
📝 Describe your requirements in natural language (optional):
Examples:
  • no strong base, room temperature
  • avoid expensive catalysts, prefer DMF solvent
  • mild conditions, no air-sensitive reagents
  • (press Enter to skip if no specific requirements)

> use copper catalyst, temperature below 120C
```

### Step 4: LLM Parsing
```
🤖 Parsing input with LLM...
✅ LLM parsing complete
```

### Step 5: Review Parsed Request
```
──────────────────────────────────────────────────────────────────────
📋 PARSED REQUEST:

✅ Reaction: Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1
   Type: C_N_Coupling

📌 Constraints:
  • metal_preference: Cu
  • temperature:
      max: 120

✅ Request is VALID and ready to submit!
```

### Step 6: Reaction Type Determination ⭐ NEW
```
──────────────────────────────────────────────────────────────────────

🔬 REACTION TYPE DETERMINATION

Initial/Provided Type: C_N_Coupling
Detected Family: Ullmann_CN (confidence: 0.90)
Detected Metal Catalyst: Cu

✓ Final Reaction Type: C_N_Coupling_Cu
Routing: metal-specific: Ullmann (rule-based)
🎯 Recommendation System: RULE-BASED (deterministic constraints)

ℹ️  C-N Coupling Routing:
   → Ullmann reaction (copper-catalyzed)
   → Uses rule-based constraint matching
```

### Step 7: Final Confirmation
```
──────────────────────────────────────────────────────────────────────
📋 SUBMISSION SUMMARY

Reaction: Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1
Type: C_N_Coupling_Cu
Constraints: 2 specified

🔍 API Request Preview:
{
  "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
  "reaction_type": "C_N_Coupling_Cu",
  "k": 10,
  "limit": 5,
  "constraints": {
    "metal_preference": "Cu",
    "temperature": {"max": 120}
  },
  "relax": {}
}
```

### Step 8: Test Mode Exit
```
══════════════════════════════════════════════════════════════════════
⚠️  TEST MODE: Stopping before actual submission
✅ Input process validation complete!
📋 The request is valid and ready to be submitted in production mode.
══════════════════════════════════════════════════════════════════════
```

## 🔄 Comparison: Copper vs Palladium

### Ullmann (Copper) - Rule-Based
```
Requirements: "use copper catalyst"

🔬 REACTION TYPE DETERMINATION
Detected Metal Catalyst: Cu
✓ Final Reaction Type: C_N_Coupling_Cu
🎯 Recommendation System: RULE-BASED (deterministic constraints)

ℹ️  C-N Coupling Routing:
   → Ullmann reaction (copper-catalyzed)
   → Uses rule-based constraint matching
```

### Buchwald (Palladium) - Machine Learning
```
Requirements: "prefer palladium catalyst"

🔬 REACTION TYPE DETERMINATION
Detected Metal Catalyst: Pd
✓ Final Reaction Type: C_N_Coupling_Pd
🎯 Recommendation System: MACHINE LEARNING (similarity-based)

ℹ️  C-N Coupling Routing:
   → Buchwald-Hartwig reaction (palladium-catalyzed)
   → Uses ML-based similarity search
```

## 🎨 Color Legend

When running with colors enabled (default):

- 🟢 **Green**: Success states, valid inputs, confirmation prompts
- 🔴 **Red**: Errors, validation issues, blocking problems
- 🟡 **Yellow**: Warnings, examples, test mode, routing info
- 🔵 **Cyan**: Headers, informational text, prompts
- 🟣 **Magenta**: Main titles, section headers

## 🧪 Test Cases

### Test Case 1: Ullmann Reaction
```bash
Reaction: Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1
Requirements: use copper catalyst, no strong base

Expected:
- Detected Metal: Cu
- Final Type: C_N_Coupling_Cu
- System: RULE-BASED
```

### Test Case 2: Buchwald Reaction
```bash
Reaction: Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1
Requirements: prefer palladium catalyst, mild conditions

Expected:
- Detected Metal: Pd
- Final Type: C_N_Coupling_Pd
- System: MACHINE LEARNING
```

### Test Case 3: Suzuki Coupling
```bash
Reaction: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1
Requirements: (leave empty)

Expected:
- Detected Family: Suzuki_CC
- Final Type: Suzuki
- System: AUTO-DETECT
```

### Test Case 4: Auto-Detection
```bash
Reaction: Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1
Requirements: (leave empty)

Expected:
- Detected Family: Ullmann_CN
- Final Type: C_N_Coupling_Cu (default for C-N)
- System: RULE-BASED
```

## 📊 Routing Decision Tree

```
Input Reaction SMILES + Constraints
        ↓
Detect Reaction Family
        ↓
Is C-N Coupling?
    Yes → Check Metal Catalyst
        ├─ Cu? → C_N_Coupling_Cu (Ullmann, Rule-Based)
        ├─ Pd? → C_N_Coupling_Pd (Buchwald, ML)
        ├─ Ni? → C_N_Coupling_Ni (Nickel, ML)
        └─ None → Default to detected family
    No → Use detected family (Suzuki, Sonogashira, etc.)
        ↓
Submit to Appropriate Recommendation Engine
```

## 🚀 Quick Commands

```bash
# Test Ullmann routing
python app/cli_recommend.py --test
# Enter: Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1
# Enter: use copper catalyst

# Test Buchwald routing
python app/cli_recommend.py --test
# Enter: Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1
# Enter: prefer palladium catalyst

# Test without colors (for logs)
python app/cli_recommend.py --test --no-color

# Production mode (actual API call)
python app/cli_recommend.py
```

## 💡 Key Features

✅ **Intelligent Routing**: Automatically routes to rule-based or ML system
✅ **Metal Detection**: Extracts catalyst from natural language
✅ **Visual Feedback**: Color-coded output for easy reading
✅ **Safe Testing**: Test mode validates without API calls
✅ **Transparent**: Shows exactly which system will be used
✅ **User-Friendly**: Natural language input, structured output

## 📚 Related Docs

- **CLI_REACTION_TYPE_ROUTING.md**: Comprehensive technical documentation
- **CLI_ENHANCEMENT_SUMMARY.md**: Summary of all changes
- **CLI_QUICKSTART.md**: General CLI usage guide

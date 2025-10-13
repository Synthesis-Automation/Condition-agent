# CLI Quick Reference Card

## 🚀 Launch Commands

```bash
# Test mode (no API call)
python app/cli_recommend.py --test

# Production mode
python app/cli_recommend.py

# Without colors
python app/cli_recommend.py --test --no-color

# With debug logging
python app/cli_recommend.py --test --debug
```

## 🎯 Reaction Type Routing

| User Says | Detected Metal | Final Type | System |
|-----------|----------------|------------|--------|
| "use copper catalyst" | Cu | C_N_Coupling_Cu | RULE-BASED |
| "prefer palladium" | Pd | C_N_Coupling_Pd | ML-BASED |
| "with CuI" | Cu | C_N_Coupling_Cu | RULE-BASED |
| "nickel catalyst" | Ni | C_N_Coupling_Ni | ML-BASED |
| (no metal specified) | None | Auto-detected | AUTO |

## 🎨 Color Codes

- 🟢 **Green**: Valid ✅ Success ✓
- 🔴 **Red**: Error ❌ Invalid ✗
- 🟡 **Yellow**: Warning ⚠️ Info ℹ️
- 🔵 **Cyan**: Headers, Prompts
- 🟣 **Magenta**: Main Titles

## 📋 Input Examples

### Ullmann (Copper)
```
Reaction: Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1
Requirements: use copper catalyst, no strong base
→ Routes to: C_N_Coupling_Cu (rule-based)
```

### Buchwald (Palladium)
```
Reaction: Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1
Requirements: prefer palladium catalyst
→ Routes to: C_N_Coupling_Pd (ML-based)
```

### Suzuki
```
Reaction: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1
Requirements: (empty)
→ Routes to: Suzuki (auto-detected)
```

## 🔍 What Gets Detected

### From Natural Language
- Temperature limits → `temperature: {max: X, min: Y}`
- Base strength → `base_strength: "weak" | "moderate" | "strong"`
- Metal preference → `metal_preference: "Cu" | "Pd" | "Ni"`
- Required reagents → `required_reagents: [...]`
- Excluded reagents → `exclude_reagents: [...]`
- Solvent preference → `solvent_preference: "polar_aprotic" | ...`

### From Reaction Structure
- Aryl halide + amine → C-N coupling
- Aryl halide + boronic acid → Suzuki
- Aryl halide + terminal alkyne → Sonogashira
- Carboxylic acid + amine → Amide coupling

## ⚡ Quick Test

```bash
python test_cli_routing.py
```

Shows routing for 5 test cases:
1. ✓ C-N + Cu → Ullmann (rule)
2. ✓ C-N + Pd → Buchwald (ML)
3. ✓ C-N auto → Detected
4. ✓ Suzuki → Auto
5. ✓ User type + Cu → Routed

## 🎯 Decision Flow

```
Input → LLM Parse → Validate → Detect Type → Route → Submit
         ↓           ↓          ↓             ↓
      Constraints  Refine   Metal Check   System Choice
```

## 📊 Systems

| System | Reaction Types | Method |
|--------|----------------|--------|
| RULE-BASED | C_N_Coupling_Cu (Ullmann) | Constraint matching |
| ML-BASED | C_N_Coupling_Pd (Buchwald), Suzuki, etc. | Similarity search |

## 🔧 Flags

| Flag | Description |
|------|-------------|
| `--test` | Stop before API call |
| `--no-color` | Disable colors |
| `--debug` | Show debug logs |
| `--provider` | LLM provider (openai/aliyun) |
| `--model` | LLM model name |
| `--api-key` | API key override |

## 📝 Files

- `app/cli_recommend.py` - Main CLI
- `test_cli_routing.py` - Test script
- `CLI_REACTION_TYPE_ROUTING.md` - Full docs
- `CLI_ENHANCEMENT_SUMMARY.md` - Summary
- `CLI_VISUAL_DEMO.md` - Visual guide

## 💡 Tips

- Use `--test` when validating inputs
- Metal preference overrides detection
- Empty requirements are OK
- Colors improve readability
- Test mode shows full API request

## 🐛 Troubleshooting

| Issue | Solution |
|-------|----------|
| No colors | Check terminal supports ANSI |
| Wrong routing | Check metal preference |
| LLM error | Verify API key set |
| Import error | Run from project root |

## 📞 Help

```bash
python app/cli_recommend.py --help
```

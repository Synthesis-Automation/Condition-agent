# Quick Start Guide: New Rule Databases

## Available Databases

Three new rule databases are now available in `data/rule_db/`:

1. **C_O_coupling_db.json** - C–O Coupling (Aryl–O bond formation)
2. **RCM_db.json** - Ring-Closing Metathesis  
3. **sonogashira_db.json** - Sonogashira Coupling

## Quick Test

Run the test script to verify everything works:

```bash
python scripts/test_new_rules.py
```

Expected output:
```
C-O Coupling         ✅ PASSED
RCM                  ✅ PASSED
Sonogashira          ✅ PASSED
```

## Example Usage

### C-O Coupling Example
```python
from chemtools.rule.engine import RuleEngine

engine = RuleEngine.from_file("data/rule_db/C_O_coupling_db.json")
reaction = "Oc1ccccc1.Brc1ccccc1>>c1ccc(Oc2ccccc2)cc1"
rec = engine.recommend(reaction)
print(rec.format_summary())
```

### RCM Example
```python
engine = RuleEngine.from_file("data/rule_db/RCM_db.json")
reaction = "C=CCCC=C>>C1=CCCC1"
rec = engine.recommend(reaction)
print(rec.format_summary())
```

### Sonogashira Example
```python
engine = RuleEngine.from_file("data/rule_db/sonogashira_db.json")
reaction = "Brc1ccccc1.C#C>>C#Cc1ccccc1"
rec = engine.recommend(reaction)
print(rec.format_summary())
```

## Database Contents

### C-O Coupling (5 base rules, 7 modifiers)
- Phenols / ArOH (Pd-based)
- Aliphatic alcohols (Cu-based)
- Hindered/tertiary alcohols
- Aryl chlorides / heteroaryl halides
- Pd precatalyst-free route

**Top condition families:**
- RockPhos–Pd G3 / Cs2CO3 / toluene
- BippyPhos systems
- CuI / diamine / KOtBu sets

### RCM (6 base rules, 10 modifiers)
- Small/medium rings (5-7 member)
- Medium rings (8-12 member)
- Macrocycles (≥13 member)
- Hindered targets (tri/tetrasubstituted alkenes)
- Amine/sulfur-containing substrates
- Stereoselective variants

**Catalysts:**
- Hoveyda–Grubbs II (HG-II) - general
- Nitro-Grela / Grela - fast-initiating
- Various Grubbs generations

### Sonogashira (5 base rules, 11 modifiers)
- Aryl iodides/bromides (standard)
- Aryl chlorides (hard electrophiles)
- Vinyl halides/triflates
- Heteroaryl halides
- Copper-free variants

**Key systems:**
- PdCl2(PPh3)2 / CuI / Et3N (classic)
- XPhos–Pd G3 (for Ar-Cl)
- Copper-free tracks (to avoid Glaser)

## CLI Usage

```bash
# C-O Coupling
python -m chemtools.rule.cli "Oc1ccccc1.Brc1ccccc1>>c1ccc(Oc2ccccc2)cc1" \
    --database C_O_coupling_db

# RCM
python -m chemtools.rule.cli "C=CCCC=C>>C1=CCCC1" \
    --database RCM_db

# Sonogashira
python -m chemtools.rule.cli "Brc1ccccc1.C#C>>C#Cc1ccccc1" \
    --database sonogashira_db
```

## Integration with ChemTools

These databases are automatically integrated:

```python
import chemtools as chem

# List all available databases
databases = chem.rules.list_databases()
print(databases)  # Includes C_O_coupling_db.json, RCM_db.json, sonogashira_db.json

# Load and use
db = chem.rules.load_database("C_O_coupling_db")
result = chem.rules.match(db, "Oc1ccccc1.Brc1ccccc1>>...")
```

## Troubleshooting

If you encounter issues:

1. **Validate structure**: `python scripts/validate_new_rules.py`
2. **Check feature detection**: Enable debug logging
3. **Verify reaction type**: Use detection API

```python
from chemtools.detection import detect_reaction

result = detect_reaction("Brc1ccccc1.C#C>>C#Cc1ccccc1")
print(result['family'])  # Should show 'sonogashira'
```

## Notes

- All databases use HTE-informed condition families
- Modifiers apply automatically based on detected features
- Each database includes detailed evaluation notes
- Compatible with existing chemtools workflow

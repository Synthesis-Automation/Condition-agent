# Protocol Database Tools

## CLI Commands

### Build/Rebuild Index

```bash
# Build index with DRFP fingerprints (stored in separate NPZ file for efficiency)
python -m chemtools.protocol.cli build --force

# Build with custom directory
python -m chemtools.protocol.cli build --force --protocol-dir "c:\Git-softwares\Condition-agent\data\protocol_db"

# View index statistics
python -m chemtools.protocol.cli stats
```

**Note**: DRFP fingerprints are now stored in a separate `.protocol_drfp.npz` file (compressed binary format) instead of being embedded in the JSON index. This reduces file size by ~98% and improves loading speed by ~10x.

### Recommend Protocols

```bash
# Find matching protocols for a reaction
python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1"

# Specify number of results
python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --k 5

# Filter by reaction family
python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --family Suzuki

# Filter by tags
python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --tags "coupling,Pd"

# Save results to JSON
python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --output results.json --pretty

# Disable SMARTS filtering (use DRFP similarity only)
python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --no-smarts-filter
```

### Validate Protocols

```bash
# Validate all protocols (check SMARTS matching)
python -m chemtools.protocol.validate_protocols

# Validate specific file
python -m chemtools.protocol.validate_protocols --file "Suzuki_Protocol.json"

# Show only errors
python -m chemtools.protocol.validate_protocols --errors-only

# Export validation report
python -m chemtools.protocol.validate_protocols --output validation_report.json

# Verbose output with details
python -m chemtools.protocol.validate_protocols --verbose

# Fail with exit code 1 if any protocols are invalid (useful for CI/CD)
python -m chemtools.protocol.validate_protocols --fail-on-error
```

## Common SMARTS Issues and Best Practices

### ❌ Common Problems

#### 1. RDKit H-count Specification Errors

**Problem**: Patterns like `[CH]`, `[CH2]`, `[CH3]`, `[NH]`, `!H0` cause RDKit errors:

```
Pre-condition Violation: getNumImplicitHs() called without preceding call to calcImplicitValence()
```

**Solution**:

```smarts
# BAD - causes RDKit error
Br[CH].[c,C,n,o,s]B>>[c,C,n,o,s][CH]

# GOOD - simplified
BrC.[c,C,n,o,s]B>>[c,C,n,o,s]C

# OR - explicit H count if critical
Br[C;H1].[c,C,n,o,s]B>>[c,C,n,o,s][C;H1]
```

#### 2. Invalid Chemical Abbreviations

**Problem**: `B2pin2`, `Opin`, `Bpin` are not valid SMARTS (they're common chemistry abbreviations)

**Solution**:

```smarts
# BAD - not valid SMARTS
IC.B2pin2>>CB(Opin)

# GOOD - actual SMARTS structure
IC.BB>>CB
# OR simplified
IC.B>>CB
```

#### 3. Overly Restrictive Patterns

**Problem**: Patterns like `[!#1]C#C[H]` require specific terminal alkyne with explicit hydrogen

**Solution**:

```smarts
# BAD - too restrictive
[c,C,n,o,s]I.[!#1]C#C[H]>>[c,C,n,o,s]C#C[!#1]

# GOOD - flexible
[c,C,n,o,s]I.C#C>>[c,C,n,o,s]C#C
```

#### 4. Explicit Hydrogen Specifications

**Problem**: Using `[H]` or `/C([H])=C([H])/` can cause parsing issues

**Solution**:

```smarts
# BAD - explicit hydrogens
[c,C,n,o,s]/C([H])=C([H])/[H].O=C([c,C,n,o,s])F>>...

# GOOD - simplified without H
C=C.O=C(F)>>C(C)C(=O)
```

#### 5. Negated H-count Queries

**Problem**: `!H0`, `!H1` trigger implicit valence calculation errors

**Solution**:

```smarts
# BAD - negated H count
[C;!a;!H0]CI.[C;!a;!H0]C(Cl)=O>>...

# GOOD - simplified
CI.C(Cl)=O>>CC(C)=O
```

### ✅ Best Practices

#### DO:

- ✅ Use simple, general patterns: `C`, `N`, `O`, `[c,C,n,o,s]`
- ✅ Test SMARTS in RDKit before adding to protocol
- ✅ Use validation tool before committing new protocols:
  ```bash
  python -m chemtools.protocol.validate_protocols --file "YourProtocol.json" --verbose
  ```
- ✅ Keep patterns flexible to match similar reactions
- ✅ Use `[C;H1]`, `[C;H2]` format if H-count is critical
- ✅ Replace `[a]` with `[c,C,n,o,s]` for aromatic/aliphatic matching

#### DON'T:

- ❌ Use `[CH]`, `[CH2]`, `[CH3]` - causes RDKit errors
- ❌ Use `!H0`, `!H1` - triggers implicit valence issues
- ❌ Use `[H]` explicitly unless absolutely necessary
- ❌ Use chemical abbreviations like `B2pin2`, `Opin`, `Bpin` in SMARTS
- ❌ Add trailing commas in JSON (not valid in strict JSON)
- ❌ Make patterns overly specific with stereochemistry unless required

### Example: Good Protocol SMARTS

```json
{
  "reaction": {
    "reaction_smiles": "CCBr.c1ccccc1B(O)O>>CCc1ccccc1",
    "reaction_SMARTS": [
      "[c,C,n,o,s]Br.OB(O)[c,C,n,o,s]>>[c,C,n,o,s][c,C,n,o,s]"
    ],
    "family": "Suzuki_Miyaura",
    "tags": "Suzuki; coupling; Pd"
  }
}
```

### Validation Workflow

```bash
# 1. Create or modify protocol JSON file

# 2. Validate the protocol
python -m chemtools.protocol.validate_protocols --file "YourProtocol.json" --verbose

# 3. If valid, rebuild index
python -m chemtools.protocol.cli build --force

# 4. Verify index statistics
python -m chemtools.protocol.cli stats
```

python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --k 3

### Testing SMARTS Patterns

You can test SMARTS patterns in Python before adding them to protocols:

```python
from rdkit import Chem
from rdkit.Chem import AllChem

# Test reaction SMARTS
rxn_smiles = "CCBr.c1ccccc1B(O)O>>CCc1ccccc1"
pattern = "[c,C,n,o,s]Br.OB(O)[c,C,n,o,s]>>[c,C,n,o,s][c,C,n,o,s]"

# Parse reaction
rxn = AllChem.ReactionFromSmarts(rxn_smiles, useSmiles=True)
rxn_pattern = AllChem.ReactionFromSmarts(pattern)

# Check if pattern is valid
if rxn_pattern:
    print("✅ Pattern is valid SMARTS")
else:
    print("❌ Pattern failed to parse")
```

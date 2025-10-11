# Protocol Recommendation System - Quick Start

## ✨ What's New

The **Protocol Recommendation Module** helps you find the most suitable chemical reaction protocols from a curated database using DRFP (Differential Reaction Fingerprint) similarity.

## 🚀 Quick Start for Windows Users

### Option 1: PowerShell Script (Recommended)

```powershell
# Interactive mode
.\test_protocol.ps1

# Quick query
.\test_protocol.ps1 "BrC1CCCCC1.c1ccccc1B(O)O>>c1ccccc1C1CCCCC1" -k 3
```

### Option 2: Batch Script

```cmd
test_protocol.bat "BrC1CCCCC1.c1ccccc1B(O)O>>c1ccccc1C1CCCCC1" -k 3
```

## 📚 Documentation

- **User Guide**: `docs/PROTOCOL_CLI_GUIDE.md` - Complete command reference and examples
- **Quick Start**: `docs/PROTOCOL_QUICKSTART.md` - 2-minute tutorial
- **Technical Docs**: `docs/PROTOCOL_MODULE.md` - Full API documentation

## 🧪 Example Queries

### Suzuki Coupling
```powershell
.\test_protocol.ps1 "BrC1CCCCC1.c1ccccc1B(O)O>>c1ccccc1C1CCCCC1"
```

### Borylation
```powershell
.\test_protocol.ps1 "CCCCCI.B(B1OC(C)(C)C(C)(C)O1)B2OC(C)(C)C(C)(C)O2>>CCCCB1OC(C)(C)C(C)(C)O1"
```

### With Filters
```powershell
.\test_protocol.ps1 "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --tags suzuki,palladium -k 5
```

## 📊 What You Get

For each matching protocol:
- **Similarity Score**: DRFP-based similarity (0-100%)
- **Protocol Details**: Title, journal, DOI, URL
- **Reaction Information**: SMILES, family, tags
- **Extracted Conditions**: Catalyst, ligand, base, solvent, temperature, time, atmosphere

## 🔧 Building the Index

Before first use, build the protocol index:

```powershell
python -m chemtools.protocol.cli build
```

The index is automatically updated when protocol files change.

## 💡 Interactive Mode

For testing multiple reactions:

```powershell
.\test_protocol.ps1
```

Then use commands:
- Type reaction SMILES to search
- `set k 5` - Show top 5 results
- `set family <name>` - Filter by reaction family
- `set tags tag1,tag2` - Filter by tags
- `settings` - View current settings
- `clear` - Clear all filters
- `help` - Show help
- `quit` - Exit

## 🎯 Current Database

- **16 protocols** from Organic Syntheses
- **16 reaction families** (Suzuki, Buchwald-Hartwig, borylation, etc.)
- **71 unique tags** for filtering

## 🛠️ For Developers

See `docs/PROTOCOL_MODULE.md` for:
- API reference
- Python usage examples
- Index management
- Integration with other modules

---

**Happy experimenting!** 🧪

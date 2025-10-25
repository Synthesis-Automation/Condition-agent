# Quick Start: LLM Workflow

## 🚀 What Changed

**Old workflow**: Messy output with 7+ JSON keys, unclear final entry
**New workflow**: Clean 3-key output, pure LLM, auto-detect role

## 📦 Installation Complete

✅ New prompt templates added to `llmtools/prompts.py`
✅ New classifier module created: `llmtools/reagent_classifier.py`
✅ New workflow function: `generate_taxonomy_entry_llm()` in `reagent_taxonomy_qt.py`
✅ Test script: `data-processor/test_llm_workflow.py`

## 🧪 Test Now

```bash
# 1. Set API key
export ALIYUN_API_KEY=your-key     # DeepSeek (recommended)
# OR
export OPENAI_API_KEY=your-key     # OpenAI

# 2. Run test
cd data-processor
python test_llm_workflow.py
```

**Expected output**:
```
🧪 Testing LLM Workflow Implementation

============================================================
Test: Triethylamine (CAS 121-44-8)
============================================================
✓ Using DeepSeek V3 (Aliyun)
✓ Registry directory: .../data/reagent_db

------------------------------------------------------------
Running LLM workflow...
------------------------------------------------------------

Status: ready_to_save

============================================================
Workflow Results
============================================================

✓ Step 1: Identity Resolved
  Name: Triethylamine
  CAS: 121-44-8
  SMILES: CCN(CC)CC
  Formula: C6H15N

✓ Step 2: Role Classification
  Role: base (confidence: 95%)
  Reasoning: "Tertiary amine with basicity"
  Model: deepseek-v3 (250 tokens, 1200ms)

✓ Step 3: Field Assignment
  Family: tertiary_amines_aliphatic
  Fields: {
      "basicity": "strong",
      "nucleophilicity": "moderate",
      "sterics": "unhindered"
  }
  Confidence: 92%
  Model: deepseek-v3 (450 tokens, 1500ms)

✓ Step 4: Verification
  Approved: ✓ YES
  Issues: 0
  Model: deepseek-v3 (350 tokens, 1100ms)

============================================================
✅ All tests passed!
============================================================
```

## 📖 Use in Your Code

```python
from pathlib import Path
from llmtools.clients import LLMClient
from reagent_taxonomy_qt import generate_taxonomy_entry_llm

# Initialize LLM client
client = LLMClient(provider="aliyun", model="deepseek-v3")
registry_dir = Path("data/reagent_db")

# Generate entry
result = generate_taxonomy_entry_llm(
    cas="121-44-8",                    # Triethylamine
    registry_dir=registry_dir,
    llm_client=client,
)

# Check status
if result["status"] == "ready_to_save":
    entry = result["entry"]
    print(f"✅ Ready to save: {entry['name']}")
    
elif result["status"] == "needs_review":
    print(f"⚠️ Needs review: {result['message']}")
    
else:  # "error"
    print(f"❌ Error: {result['error']}")
```

## 📊 Output Format

**Clean 3-key structure**:

```json
{
  "status": "ready_to_save",        // or "needs_review" or "error"
  
  "workflow": {
    "step1_identity": {...},        // CAS resolution
    "step2_role": {...},            // LLM role classification
    "step3_fields": {...},          // LLM field assignment
    "step4_verification": {...}     // LLM quality check
  },
  
  "entry": {                        // Final entry to save
    "name": "Triethylamine",
    "cas": "121-44-8",
    "molecular_formula": "C6H15N",
    "smiles": "CCN(CC)CC",
    "roles": {
      "base": {
        "families": ["tertiary_amines_aliphatic"],
        "basicity": "strong",
        "nucleophilicity": "moderate",
        "sterics": "unhindered"
      }
    },
    "abbreviations": ["TEA"]
  }
}
```

## 🎯 Next Steps

**Option 1: Test First** (recommended)
```bash
python data-processor/test_llm_workflow.py
```

**Option 2: Integrate with UI**
- Ask me to add UI mode toggle
- I'll connect the new workflow to the GUI

**Option 3: Request Features**
- Batch processing (multiple CAS)
- Confidence thresholds
- Multi-role support
- Auto-save when ready

## 📚 Documentation

**Quick Guides**:
- `LLM_WORKFLOW_SUMMARY.md` - Executive summary (this is most important!)
- `LLM_WORKFLOW_IMPLEMENTATION.md` - Technical deep dive
- `LLM_QUICKSTART.md` - Previous quick start (old workflow)

**Test Scripts**:
- `test_llm_workflow.py` - Test new pure LLM workflow
- `test_markdown_fix.py` - Test DeepSeek markdown parsing

## 🔑 Key Benefits

✅ **Auto-detect role** - No manual selection needed
✅ **Clean output** - 3 keys instead of 7+
✅ **Built-in QA** - Verification step catches errors
✅ **Transparent** - Full workflow trace with confidence scores
✅ **Pure LLM** - No mixed deterministic/LLM confusion

## ⚠️ Important Notes

1. **API Keys Required**: Set `ALIYUN_API_KEY` or `OPENAI_API_KEY`
2. **Registry Must Exist**: `data/reagent_db/` directory with schema
3. **Old Workflow Unchanged**: `generate_taxonomy_entry()` still works
4. **DeepSeek Compatible**: Handles markdown code fences automatically

## 🐛 Troubleshooting

**"LLM support not available"**:
```bash
pip install openai  # or your LLM provider SDK
```

**"Registry directory not found"**:
```bash
# Make sure you're in the right directory
cd Condition-agent
ls data/reagent_db  # Should exist
```

**"No data from PubChem"**:
- Check internet connection
- Verify CAS number is valid
- Try with `name_override` parameter

**DeepSeek returns markdown fences**:
- ✅ Already handled! The `_strip_markdown_fences()` function removes them

## 📞 Need Help?

Just ask! I can:
- Add UI integration
- Create more test cases
- Add batch processing
- Explain any part in detail
- Fix any issues you encounter

# Optional LLM Utilities

`llmtools` contains provider clients and optional review helpers. It does not
analyze reactions, select reaction families, retrieve precedents, score
recipes, or validate chemistry.

The condition coworker is deterministic by default. If natural-language
narration is added later, it must consume an existing typed recommendation
result and may not alter recipes, evidence, cautions, or confidence.

```python
from llmtools import LLMClient

client = LLMClient(provider="openai", model="gpt-5.6-terra")
```

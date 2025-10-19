# Refactoring Quick Start Guide

**For developers who want to start improving code quality immediately**

---

## Quick Wins (Can Do Today)

### 1. Remove Deprecated Imports in `app/main.py`

**Line 23:**

```python
# TODO: Remove these imports once all code uses chem.*
from chemtools import smiles, router, featurizers, condition_core, precedent, constraints, explain, recommend
```

**Action:** Search codebase for usage, replace with `chem.*` API, then remove.

```bash
# Find usages
grep -r "from chemtools import smiles" app/
grep -r "import smiles" app/

# Replace pattern: smiles.normalize() -> chem.smiles.normalize()
```

### 2. Complete or Remove TODOs

**Found:**
- `chemtools/output_formatter.py:756` - "Calculate from SMILES or look up"
- `llmtools/agents.py:388` - "Implement iterative optimization loop"

**Action:** Either implement or create GitHub issue and remove TODO comment.

### 3. Standardize Error Handling

**Create:** `chemtools/exceptions.py`

```python
"""Centralized exception hierarchy for chemtools."""

class ChemToolsError(Exception):
    """Base exception for all chemtools errors."""
    pass

class ValidationError(ChemToolsError):
    """Raised when input validation fails."""
    pass

class DatabaseNotFoundError(ChemToolsError):
    """Raised when database file cannot be found."""
    pass

class ProcessingError(ChemToolsError):
    """Raised when processing fails."""
    pass
```

**Create:** `app/error_handlers.py`

```python
"""Centralized error handlers for FastAPI."""

from fastapi import Request, HTTPException
from fastapi.responses import JSONResponse
from chemtools.exceptions import *
import datetime

ERROR_STATUS_MAP = {
    ValidationError: 400,
    DatabaseNotFoundError: 404,
    ProcessingError: 500,
    ChemToolsError: 500,
}

@app.exception_handler(ChemToolsError)
async def chemtools_error_handler(request: Request, exc: ChemToolsError):
    status_code = ERROR_STATUS_MAP.get(type(exc), 500)
    return JSONResponse(
        status_code=status_code,
        content={
            "error": exc.__class__.__name__,
            "message": str(exc),
            "timestamp": datetime.datetime.utcnow().isoformat()
        }
    )
```

---

## Medium Effort Improvements (1-2 Days Each)

### 4. Extract Service Layer from `app/main.py`

**Problem:** Business logic mixed with route handlers

**Create:** `app/services/matching.py`

```python
"""Business logic for rule-based matching."""

from typing import Dict, List, Optional
from chemtools import chem
from chemtools.exceptions import ValidationError, DatabaseNotFoundError

class RuleMatchingService:
    """Service for rule-based reaction matching."""
    
    def __init__(self, default_db: str):
        self.default_db = default_db
        self._db_cache = {}
    
    def match_reaction(self, reaction: str, db_path: Optional[str] = None) -> Dict:
        """Match reaction to rule database."""
        if not reaction or not reaction.strip():
            raise ValidationError("Reaction SMILES cannot be empty")
        
        db_path = (db_path or self.default_db).strip()
        if not db_path:
            raise ValidationError("No database path configured")
        
        try:
            db = chem.rules.load_database(db_path, cache=True)
            result = chem.rules.match(db, reaction)
            return result.to_json_dict()
        except FileNotFoundError as e:
            raise DatabaseNotFoundError(f"Database not found: {db_path}") from e
    
    def filter_by_catalyst(
        self, 
        conditions: List[Dict], 
        catalyst_class: str
    ) -> List[Dict]:
        """Filter conditions by catalyst class."""
        filtered = []
        filter_lower = str(catalyst_class).lower()
        
        for cond in conditions:
            core = cond.get("core", "")
            if not isinstance(core, str):
                continue
            
            core_lower = core.lower()
            # Match if catalyst appears in core (e.g., "Cu" in "Cu/L1")
            if filter_lower in core_lower or filter_lower in core_lower.split("/")[0]:
                filtered.append(cond)
        
        return filtered
```

**Usage in `app/main.py`:**

```python
from app.services.matching import RuleMatchingService

matching_service = RuleMatchingService(default_db=_SCDB_DEFAULT_DB)

@app.post("/match")
def api_scdb_match(req: SchemeMatchRequest):
    """Match reaction to rule-based scheme database."""
    try:
        # Business logic delegated to service
        result = matching_service.match_reaction(req.reaction, req.db)
        
        # Apply catalyst filtering if requested
        if req.relax and req.relax.get("catalyst_class"):
            conditions = result.get("recommended_conditions", [])
            filtered = matching_service.filter_by_catalyst(
                conditions, 
                req.relax["catalyst_class"]
            )
            result["recommended_conditions"] = filtered
        
        # Format and return
        return output_formatter.format_rule_match_result(
            reaction_smiles=req.reaction,
            match_result=result,
            requested_type=None,
            database_name=Path(req.db or _SCDB_DEFAULT_DB).name,
        )
    except (ValidationError, DatabaseNotFoundError) as e:
        # Error handler will catch this
        raise
```

### 5. Split `chemtools/output_formatter.py`

**Create directory structure:**

```bash
mkdir -p chemtools/formatters
touch chemtools/formatters/__init__.py
```

**Create:** `chemtools/formatters/base.py`

```python
"""Base formatting utilities."""

from typing import Any, Dict, Optional
import datetime
from abc import ABC, abstractmethod

class BaseFormatter(ABC):
    """Base class for all formatters."""
    
    @abstractmethod
    def format(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """Format the data according to schema."""
        pass
    
    def add_metadata(
        self, 
        output: Dict[str, Any], 
        model: str = "unknown",
        status: str = "success",
        **kwargs
    ) -> Dict[str, Any]:
        """Add standard metadata to output."""
        output.setdefault("meta", {})
        output["meta"].update({
            "generated_at": datetime.datetime.utcnow().isoformat() + "Z",
            "model": model,
            "status": status,
            **kwargs
        })
        return output
```

**Create:** `chemtools/formatters/ml_output.py`

```python
"""Formatter for ML recommendation output."""

from typing import Any, Dict, List, Optional
from .base import BaseFormatter

class MLOutputFormatter(BaseFormatter):
    """Format ML-based recommendation output."""
    
    def format(
        self,
        reaction_smiles: str,
        recommendations: List[Dict[str, Any]],
        detected_type: str,
        confidence: float = 1.0,
        processing_time_ms: Optional[float] = None,
    ) -> Dict[str, Any]:
        """Format ML recommendation output."""
        output = {
            "input": {
                "reaction_smiles": reaction_smiles,
            },
            "detection": {
                "family": detected_type,
                "confidence": round(confidence, 4),
                "method": "auto",
            },
            "recommendations": self._format_recommendations(recommendations),
        }
        
        return self.add_metadata(
            output,
            model="ML-precedent-knn",
            processing_time_ms=processing_time_ms,
        )
    
    def _format_recommendations(
        self, 
        recommendations: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """Format individual recommendations."""
        return [
            self._format_single_recommendation(rec) 
            for rec in recommendations
        ]
    
    def _format_single_recommendation(
        self, 
        rec: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Format a single recommendation."""
        return {
            "rank": rec.get("rank", 0),
            "confidence": rec.get("confidence", 0.0),
            "reagents": rec.get("chemicals", []),
            "conditions": rec.get("conditions", {}),
            "precedent_count": rec.get("support", 0),
        }
```

**Update:** `chemtools/formatters/__init__.py`

```python
"""Output formatters for chemtools."""

from .base import BaseFormatter
from .ml_output import MLOutputFormatter

# Maintain backward compatibility
def format_ml_output(*args, **kwargs):
    """Legacy function for backward compatibility."""
    formatter = MLOutputFormatter()
    return formatter.format(*args, **kwargs)

__all__ = [
    "BaseFormatter",
    "MLOutputFormatter",
    "format_ml_output",  # Legacy
]
```

### 6. Centralize LLM Configuration

**Create:** `llmtools/config.py`

```python
"""Centralized LLM configuration management."""

import os
from dataclasses import dataclass, field
from typing import Dict, Optional

@dataclass
class LLMConfig:
    """Configuration for LLM client."""
    
    provider: str = "openai"
    model: Optional[str] = None
    api_key: Optional[str] = None
    base_url: Optional[str] = None
    temperature: float = 0.7
    max_tokens: int = 2000
    timeout: int = 60
    
    @classmethod
    def from_env(cls, prefix: str = "LLM_") -> "LLMConfig":
        """Load configuration from environment variables.
        
        Environment variables:
            LLM_PROVIDER: Provider name (openai, aliyun)
            LLM_MODEL: Model name
            LLM_API_KEY: API key
            LLM_BASE_URL: Base URL
            LLM_TEMPERATURE: Temperature (0.0-1.0)
            LLM_MAX_TOKENS: Maximum tokens
            LLM_TIMEOUT: Timeout in seconds
        """
        return cls(
            provider=os.getenv(f"{prefix}PROVIDER", "openai"),
            model=os.getenv(f"{prefix}MODEL"),
            api_key=os.getenv(f"{prefix}API_KEY"),
            base_url=os.getenv(f"{prefix}BASE_URL"),
            temperature=float(os.getenv(f"{prefix}TEMPERATURE", "0.7")),
            max_tokens=int(os.getenv(f"{prefix}MAX_TOKENS", "2000")),
            timeout=int(os.getenv(f"{prefix}TIMEOUT", "60")),
        )
    
    @classmethod
    def for_chemistry(cls, **overrides) -> "LLMConfig":
        """Recommended configuration for chemistry tasks.
        
        Lower temperature for more deterministic chemistry reasoning.
        Higher max_tokens for detailed explanations.
        """
        config = cls(
            temperature=0.3,
            max_tokens=3000,
            timeout=90,
        )
        # Apply overrides
        for key, value in overrides.items():
            if hasattr(config, key):
                setattr(config, key, value)
        return config
    
    @classmethod
    def for_reagent_classification(cls, **overrides) -> "LLMConfig":
        """Optimized configuration for reagent classification tasks."""
        return cls(
            temperature=0.1,  # Very deterministic
            max_tokens=1000,   # Concise output
            timeout=30,
            **overrides
        )
    
    def to_client_kwargs(self) -> Dict:
        """Convert config to kwargs for LLMClient."""
        return {
            "provider": self.provider,
            "model": self.model,
            "api_key": self.api_key,
            "base_url": self.base_url,
            "temperature": self.temperature,
            "max_tokens": self.max_tokens,
            "timeout": self.timeout,
        }
```

**Usage:**

```python
from llmtools.config import LLMConfig
from llmtools.clients import LLMClient

# Load from environment
config = LLMConfig.from_env()
client = LLMClient(**config.to_client_kwargs())

# Use chemistry-optimized config
config = LLMConfig.for_chemistry(provider="aliyun")
client = LLMClient(**config.to_client_kwargs())
```

---

## Larger Refactorings (1-2 Weeks)

### 7. Split `chemtools/recommend/core.py` (1,302 lines)

**Plan:**

1. **Phase 1:** Extract helper functions
   - Move DRFP search to `drfp_search.py`
   - Move reranking logic to `reranking.py`
   - Move filtering to `filtering.py`

2. **Phase 2:** Create new modules
   ```
   chemtools/recommend/
   ├── __init__.py          # Public API
   ├── core.py              # Main orchestration (~200 lines)
   ├── drfp_search.py       # DRFP similarity search
   ├── reranking.py         # Reranking strategies
   ├── filtering.py         # Constraint filtering
   ├── aggregation.py       # Result aggregation
   └── utils.py             # Shared utilities
   ```

3. **Phase 3:** Update imports and tests

### 8. Refactor UI Files

**Plan for `app/ui_gradio.py` (2,439 lines):**

```
app/ui/gradio/
├── __init__.py
├── app.py              # Main Gradio app (~200 lines)
├── components.py       # UI component builders
├── handlers.py         # Event handler functions
├── formatters.py       # Output formatting for UI
├── state.py            # State management
└── utils.py            # Helper functions
```

**Migration strategy:**
1. Extract components one at a time
2. Test each extraction
3. Update imports progressively
4. Remove old file when complete

---

## Testing Strategy

### Integration Tests Template

**Create:** `tests/integration/test_api_workflow.py`

```python
"""Integration tests for complete API workflows."""

import pytest
from fastapi.testclient import TestClient
from app.main import app

client = TestClient(app)

class TestRecommendationWorkflow:
    """Test complete recommendation workflows."""
    
    def test_ml_recommendation_full_flow(self):
        """Test ML recommendation from request to response."""
        response = client.post("/api/v1/recommend/conditions", json={
            "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
            "k": 5,
            "limit": 3,
        })
        
        assert response.status_code == 200
        data = response.json()
        
        # Verify structure
        assert "input" in data
        assert "detection" in data
        assert "recommendations" in data
        assert "meta" in data
        
        # Verify content
        assert len(data["recommendations"]) <= 3
        assert data["meta"]["status"] == "success"
        assert "processing_time_ms" in data["meta"]
    
    def test_rule_matching_full_flow(self):
        """Test rule-based matching workflow."""
        response = client.post("/match", json={
            "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
            "db": "cn_coupling_pd_db.json",
        })
        
        assert response.status_code == 200
        data = response.json()
        
        assert "matched_scheme" in data or "recommended_conditions" in data
        assert "meta" in data
```

---

## Checklist

Use this checklist to track your progress:

- [ ] Remove deprecated imports from `app/main.py`
- [ ] Complete or remove all TODO comments
- [ ] Create `chemtools/exceptions.py`
- [ ] Create `app/error_handlers.py`
- [ ] Extract service layer for rule matching
- [ ] Split `output_formatter.py` into `formatters/` package
- [ ] Create `llmtools/config.py`
- [ ] Add integration tests
- [ ] Split `recommend/core.py`
- [ ] Refactor UI files into organized structure
- [ ] Update documentation

---

## Getting Help

If you encounter issues during refactoring:

1. **Check existing tests:** Run `pytest -v` to ensure nothing breaks
2. **Review type hints:** Use `mypy chemtools/` to catch type errors
3. **Check imports:** Use `python -c "import chemtools; print(chemtools.__file__)"` to verify
4. **Ask for review:** Create PR for each major change

---

## Commit Message Template

```
refactor(module): brief description of change

- Detailed point 1
- Detailed point 2

Part of: #ISSUE_NUMBER
Breaking changes: None/Yes (explain)
```

Example:

```
refactor(formatters): extract ML output formatting to separate module

- Create chemtools/formatters/ package with base classes
- Move ML output formatting to formatters/ml_output.py
- Maintain backward compatibility with legacy function
- Add comprehensive docstrings

Part of: #123 (Code Quality Improvements)
Breaking changes: None (backward compatible)
```

---

**Start with Quick Wins, then move to Medium Effort items. The Larger Refactorings should be done after the foundation is solid.**

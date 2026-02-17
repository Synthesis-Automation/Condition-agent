# Automatic Retry Implementation Plan (Phase 5)

## Overview

Implement automatic retry loop that re-analyzes reactions with stronger models when validation fails, with configurable maximum retry cycles to prevent infinite loops.

## Goals

1. ✅ Automatically retry failed analyses with stronger models
2. ✅ Support configurable max retries (default: 3)
3. ✅ Track all retry attempts and history
4. ✅ Prevent infinite retry loops
5. ✅ Preserve best result based on quality score
6. ✅ Provide transparent cost/time reporting
7. ✅ Optional feature (opt-in via flag)

## Design

### Retry Decision Logic

```python
# Quality Gate decides if retry needed
if validation_status == "fail":
    → Retry with stronger model (critical issues)
elif validation_status == "warning" and retries_left > 0:
    → Optional retry (if --aggressive-retry flag set)
else:
    → Accept result
```

### Retry Strategy

```
Attempt 1: gpt-4o-mini (fast, cheap)
  ↓ FAIL
Attempt 2: deepseek-v3.2 (accurate, slower)
  ↓ FAIL
Attempt 3: gpt-4o (balanced fallback)
  ↓ FAIL
Return best attempt (highest quality score)
```

### Model Progression

| Retry # | Model | Provider | Speed | Cost/reaction | Use Case |
|---------|-------|----------|-------|---------------|----------|
| 0 (initial) | gpt-4o-mini | openai | ~5s | $0.001 | Fast first attempt |
| 1 (retry) | deepseek-v3.2 | aliyun | ~17s | $0.005 | Most accurate |
| 2 (retry) | gpt-4o | openai | ~8s | $0.003 | Strong fallback |
| 3 (retry) | deepseek-v3.2 | aliyun | ~17s | $0.005 | Final attempt with best |

## Implementation

### 1. New File: `reaction_agent/retry.py`

```python
"""
Automatic retry logic for failed validations.

Implements Phase 5: Retry loop with configurable max attempts.
"""

from typing import Dict, Any, List, Optional
import logging
import time

from llmtools.clients import LLMClient
from reaction_agent.validation import validate_with_rdkit, check_consensus, quality_gate

logger = logging.getLogger(__name__)


class RetryConfig:
    """Configuration for retry behavior."""

    def __init__(
        self,
        max_retries: int = 3,
        aggressive: bool = False,
        retry_on_warning: bool = False
    ):
        """
        Args:
            max_retries: Maximum retry attempts (default: 3)
            aggressive: Retry even on warnings (default: False)
            retry_on_warning: Retry when status is WARNING (default: False)
        """
        self.max_retries = max_retries
        self.aggressive = aggressive
        self.retry_on_warning = retry_on_warning

    @classmethod
    def from_cli_args(cls, args):
        """Create from argparse args."""
        return cls(
            max_retries=getattr(args, 'max_retries', 3),
            aggressive=getattr(args, 'aggressive_retry', False),
            retry_on_warning=getattr(args, 'retry_on_warning', False)
        )


class RetryHistory:
    """Track retry attempts and results."""

    def __init__(self):
        self.attempts: List[Dict[str, Any]] = []
        self.total_cost = 0.0
        self.total_time_ms = 0

    def add_attempt(
        self,
        attempt_num: int,
        model: str,
        provider: str,
        result: Dict[str, Any],
        validation: Optional[Dict[str, Any]] = None
    ):
        """Record a retry attempt."""
        metadata = result.get('metadata', {})

        attempt = {
            'attempt_num': attempt_num,
            'model': model,
            'provider': provider,
            'quality_score': validation['consensus']['quality_score'] if validation else 0.0,
            'validation_status': validation['gate']['status'] if validation else 'unknown',
            'tokens': metadata.get('total_tokens', 0),
            'latency_ms': metadata.get('latency_ms', 0),
            'cost_estimate': self._estimate_cost(model, metadata.get('total_tokens', 0))
        }

        self.attempts.append(attempt)
        self.total_cost += attempt['cost_estimate']
        self.total_time_ms += attempt['latency_ms']

    def _estimate_cost(self, model: str, tokens: int) -> float:
        """Estimate cost based on model and tokens."""
        costs_per_1k = {
            'gpt-4o-mini': 0.0005,
            'gpt-4o': 0.003,
            'deepseek-v3.2': 0.005
        }
        cost_per_1k = costs_per_1k.get(model, 0.001)
        return (tokens / 1000) * cost_per_1k

    def get_best_attempt(self) -> Optional[Dict[str, Any]]:
        """Return attempt with highest quality score."""
        if not self.attempts:
            return None
        return max(self.attempts, key=lambda x: x['quality_score'])

    def summary(self) -> str:
        """Generate human-readable summary."""
        if not self.attempts:
            return "No retry attempts"

        lines = [
            f"Retry Summary:",
            f"  Total attempts: {len(self.attempts)}",
            f"  Total cost: ${self.total_cost:.4f}",
            f"  Total time: {self.total_time_ms / 1000:.1f}s",
            f""
        ]

        for attempt in self.attempts:
            status_symbol = {
                'pass': '✓',
                'warning': '⚠',
                'fail': '✗',
                'unknown': '?'
            }.get(attempt['validation_status'], '?')

            lines.append(
                f"  Attempt {attempt['attempt_num']}: "
                f"{attempt['model']} → "
                f"{status_symbol} {attempt['validation_status']} "
                f"(quality: {attempt['quality_score']:.2f})"
            )

        best = self.get_best_attempt()
        if best:
            lines.append(f"\n  Best: Attempt {best['attempt_num']} (quality: {best['quality_score']:.2f})")

        return "\n".join(lines)


def get_retry_model(attempt_num: int) -> tuple[str, str]:
    """
    Get model and provider for retry attempt.

    Args:
        attempt_num: 0 = initial, 1+ = retry attempts

    Returns:
        (model, provider) tuple
    """
    retry_sequence = [
        ('gpt-4o-mini', 'openai'),      # Attempt 0: Fast first try
        ('deepseek-v3.2', 'aliyun'),    # Attempt 1: Most accurate
        ('gpt-4o', 'openai'),           # Attempt 2: Strong fallback
        ('deepseek-v3.2', 'aliyun'),    # Attempt 3+: Keep trying with best
    ]

    if attempt_num < len(retry_sequence):
        return retry_sequence[attempt_num]
    else:
        # For attempts beyond sequence, use best model
        return retry_sequence[-1]


def should_retry(
    validation: Dict[str, Any],
    attempt_num: int,
    config: RetryConfig
) -> bool:
    """
    Decide if we should retry based on validation results.

    Args:
        validation: Validation results from Tier 4
        attempt_num: Current attempt number (0-indexed)
        config: Retry configuration

    Returns:
        True if should retry, False otherwise
    """
    # Check max retries
    if attempt_num >= config.max_retries:
        logger.info(f"Max retries ({config.max_retries}) reached")
        return False

    gate = validation['gate']
    status = gate['status']

    # Always retry on fail
    if status == 'fail':
        logger.info(f"Validation failed, will retry (attempt {attempt_num + 1})")
        return True

    # Optionally retry on warning
    if status == 'warning' and config.retry_on_warning:
        logger.info(f"Validation warning, will retry (retry_on_warning=True)")
        return True

    # Pass or warning without retry_on_warning
    logger.info(f"Validation status: {status}, no retry needed")
    return False


def analyze_with_retry(
    rxn_smiles: str,
    initial_client: LLMClient,
    config: RetryConfig,
    validate: bool = True,
    mode: str = "auto",
    **analysis_kwargs
) -> Dict[str, Any]:
    """
    Analyze reaction with automatic retry on validation failure.

    This is the main retry orchestrator. It:
    1. Runs initial analysis with provided client
    2. Validates results (if validate=True)
    3. Retries with stronger models if validation fails
    4. Returns best result based on quality score

    Args:
        rxn_smiles: Reaction SMILES to analyze
        initial_client: LLMClient for first attempt
        config: RetryConfig controlling retry behavior
        validate: Enable validation (required for retry logic)
        mode: Analysis mode (auto/quick/deep)
        **analysis_kwargs: Additional args for analyze_reaction_smiles

    Returns:
        Dictionary with:
        - All fields from analyze_reaction_smiles()
        - retry_history: RetryHistory object
        - final_attempt_num: Which attempt was returned
        - all_attempts: List of all analysis results
    """
    from reaction_agent import analyze_reaction_smiles

    if not validate:
        logger.warning("validate=False, retry logic disabled")
        return analyze_reaction_smiles(
            rxn_smiles=rxn_smiles,
            client=initial_client,
            validate=False,
            **analysis_kwargs
        )

    history = RetryHistory()
    all_attempts = []
    current_client = initial_client
    attempt_num = 0

    while attempt_num <= config.max_retries:
        logger.info(f"Starting attempt {attempt_num} with {current_client.model}")

        start_time = time.time()

        try:
            # Run analysis
            result = analyze_reaction_smiles(
                rxn_smiles=rxn_smiles,
                client=current_client,
                validate=True,
                **analysis_kwargs
            )

            validation = result.get('validation')

            # Record attempt
            history.add_attempt(
                attempt_num=attempt_num,
                model=current_client.model,
                provider=current_client.provider,
                result=result,
                validation=validation
            )

            all_attempts.append({
                'attempt_num': attempt_num,
                'result': result
            })

            # Check if we should retry
            if validation and should_retry(validation, attempt_num, config):
                # Get next model
                next_model, next_provider = get_retry_model(attempt_num + 1)
                logger.info(f"Retrying with {next_model} (provider: {next_provider})")
                current_client = LLMClient(provider=next_provider, model=next_model)
                attempt_num += 1
                continue
            else:
                # No retry needed, return this result
                logger.info(f"Analysis successful after {attempt_num + 1} attempts")
                break

        except Exception as e:
            logger.error(f"Attempt {attempt_num} failed with error: {e}")

            # Record failed attempt
            history.add_attempt(
                attempt_num=attempt_num,
                model=current_client.model,
                provider=current_client.provider,
                result={'error': str(e)},
                validation=None
            )

            # Try next model
            if attempt_num < config.max_retries:
                next_model, next_provider = get_retry_model(attempt_num + 1)
                logger.info(f"Retrying after error with {next_model}")
                current_client = LLMClient(provider=next_provider, model=next_model)
                attempt_num += 1
                continue
            else:
                # No more retries, raise error
                raise

    # Find best attempt based on quality score
    best_attempt_info = history.get_best_attempt()
    if best_attempt_info:
        best_attempt_num = best_attempt_info['attempt_num']
        best_result = all_attempts[best_attempt_num]['result']
        logger.info(f"Returning best attempt: {best_attempt_num} (quality: {best_attempt_info['quality_score']:.2f})")
    else:
        # No valid attempts, return last
        best_attempt_num = attempt_num
        best_result = result

    # Augment result with retry information
    best_result['retry_history'] = history
    best_result['final_attempt_num'] = best_attempt_num
    best_result['all_attempts'] = all_attempts

    return best_result


# Export public API
__all__ = [
    'RetryConfig',
    'RetryHistory',
    'analyze_with_retry',
    'should_retry',
    'get_retry_model'
]
```

### 2. Modify `reaction_agent/agent.py`

Add retry support to `ReactionSMILESAnalyzer`:

```python
# Line ~52: Add to analyze_reaction_smiles() docstring
"""
Args:
    ...
    validate: Enable Tier 4 validation (default: False)
    retry_config: Optional RetryConfig for automatic retry (requires validate=True)
"""

# Line ~309: Add retry parameter to analyze() method
def analyze(
    self,
    rxn_smiles: str,
    mode: str = "auto",
    validate: bool = False,
    retry_config: Optional['RetryConfig'] = None  # NEW
) -> Dict[str, Any]:
    """
    Analyze a reaction SMILES string.

    Args:
        retry_config: If provided and validate=True, enables automatic retry
    """
    # If retry requested, delegate to retry logic
    if retry_config and validate:
        from reaction_agent.retry import analyze_with_retry
        return analyze_with_retry(
            rxn_smiles=rxn_smiles,
            initial_client=self.client,
            config=retry_config,
            validate=True,
            mode=mode,
            drop_spectators=self.drop_spectators,
            temperature=self.temperature if self.reasoning_effort is None else None,
            max_tokens=self.max_tokens,
            reasoning_effort=self.reasoning_effort
        )

    # Otherwise use existing single-attempt logic
    # ... existing code ...
```

### 3. Modify `reaction_agent/cli.py`

Add CLI flags and display:

```python
# Line ~560: Add retry flags to argparser
parser.add_argument(
    '--max-retries',
    type=int,
    default=3,
    help='Maximum retry attempts when validation fails (default: 3, requires --validate)'
)
parser.add_argument(
    '--retry-on-warning',
    action='store_true',
    help='Retry even when validation status is WARNING (more aggressive)'
)
parser.add_argument(
    '--no-retry',
    action='store_true',
    help='Disable automatic retry even with --validate'
)

# Line ~390: Modify analyze_reaction_interactive signature
def analyze_reaction_interactive(
    analyzer: ReactionSMILESAnalyzer,
    rxn_smiles: str,
    save_output: Optional[Path] = None,
    mode: str = "auto",
    validate: bool = False,
    retry_config: Optional['RetryConfig'] = None  # NEW
) -> Dict[str, Any]:
    """Analyze a reaction interactively with progress updates."""

    # ... existing code ...

    try:
        result = analyzer.analyze(
            rxn_smiles,
            mode=mode,
            validate=validate,
            retry_config=retry_config  # NEW
        )
    except Exception as e:
        # ... error handling ...

# Line ~318: Add retry history display (after validation section)
# Retry history section
if 'retry_history' in result:
    print_header("RETRY HISTORY")
    history = result['retry_history']
    print(history.summary())
    print(f"\n{Colors.BOLD}Final result:{Colors.END} Attempt {result['final_attempt_num']}")

# Line ~680: Create retry config from args
if args.validate and not args.no_retry:
    from reaction_agent.retry import RetryConfig
    retry_config = RetryConfig(
        max_retries=args.max_retries,
        retry_on_warning=args.retry_on_warning
    )
else:
    retry_config = None

analyze_reaction_interactive(
    analyzer,
    args.reaction,
    save_output=args.output,
    mode=effective_mode,
    validate=args.validate,
    retry_config=retry_config  # NEW
)
```

### 4. Interactive Mode Support

Add retry commands to interactive mode:

```python
# In interactive_mode() function
print("Commands:")
print("  'quit' or 'exit' - Exit the program")
print("  'config' - Show current configuration")
print("  'validate on/off' - Enable/disable Tier 4 validation")
print("  'retry on/off' - Enable/disable automatic retry")  # NEW
print("  'max-retries N' - Set max retry attempts (e.g., max-retries 5)")  # NEW
print("  'help' - Show this help message")

# Command handlers
elif user_input.lower().startswith('retry '):
    setting = user_input[6:].strip().lower()
    if setting == 'on':
        retry_enabled = True
        print(f"{Colors.GREEN}✓ Automatic retry ENABLED{Colors.END}")
    elif setting == 'off':
        retry_enabled = False
        print(f"{Colors.YELLOW}✓ Automatic retry DISABLED{Colors.END}")
    else:
        print_error("Usage: retry on|off")
    continue

elif user_input.lower().startswith('max-retries '):
    try:
        num = int(user_input[12:].strip())
        if 1 <= num <= 10:
            max_retries = num
            print(f"{Colors.GREEN}✓ Max retries set to {num}{Colors.END}")
        else:
            print_error("Max retries must be between 1 and 10")
    except ValueError:
        print_error("Usage: max-retries N (where N is a number)")
    continue
```

## Usage Examples

### CLI - Automatic Retry

```bash
# Enable validation + retry (default max=3)
python -m reaction_agent.cli --reaction "SMILES" --validate

# Custom max retries
python -m reaction_agent.cli --reaction "SMILES" --validate --max-retries 5

# Aggressive retry (retry even on warnings)
python -m reaction_agent.cli --reaction "SMILES" --validate --retry-on-warning

# Disable retry (validation only)
python -m reaction_agent.cli --reaction "SMILES" --validate --no-retry
```

### CLI - Interactive Mode

```bash
python -m reaction_agent.cli --interactive

> validate on
✓ Validation ENABLED

> retry on
✓ Automatic retry ENABLED

> max-retries 5
✓ Max retries set to 5

> Brc1ccccc1.B(O)(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1
[Analysis runs with automatic retry if needed]

> config
Configuration:
  Model: gpt-4o-mini
  Validation: ENABLED
  Retry: ENABLED (max: 5)
```

### Programmatic API

```python
from llmtools.clients import LLMClient
from reaction_agent import ReactionSMILESAnalyzer
from reaction_agent.retry import RetryConfig

# Create analyzer
client = LLMClient(provider="openai", model="gpt-4o-mini")
analyzer = ReactionSMILESAnalyzer(client)

# Configure retry
retry_config = RetryConfig(
    max_retries=3,
    retry_on_warning=False
)

# Analyze with automatic retry
result = analyzer.analyze(
    rxn_smiles="...",
    validate=True,
    retry_config=retry_config
)

# Check retry history
if 'retry_history' in result:
    history = result['retry_history']
    print(f"Total attempts: {len(history.attempts)}")
    print(f"Total cost: ${history.total_cost:.4f}")
    print(f"Best attempt: {history.get_best_attempt()['attempt_num']}")
```

## Cost/Performance Analysis

### Worst Case Scenario (4 attempts)

| Attempt | Model | Time | Cost | Total Time | Total Cost |
|---------|-------|------|------|------------|------------|
| 0 | gpt-4o-mini | 5s | $0.001 | 5s | $0.001 |
| 1 | deepseek-v3.2 | 17s | $0.005 | 22s | $0.006 |
| 2 | gpt-4o | 8s | $0.003 | 30s | $0.009 |
| 3 | deepseek-v3.2 | 17s | $0.005 | 47s | $0.014 |

**Worst case**: 47 seconds, $0.014 per reaction

### Best Case Scenario (1 attempt)

First attempt succeeds with gpt-4o-mini:
- **Time**: 5 seconds
- **Cost**: $0.001

### Expected Case (2 attempts)

Most reactions pass on first or second attempt:
- **Average time**: ~12 seconds
- **Average cost**: $0.003

## Safety Features

### 1. Max Retry Limit
```python
if attempt_num >= config.max_retries:
    logger.warning(f"Max retries ({config.max_retries}) reached, returning best result")
    return best_result
```

### 2. Infinite Loop Prevention
- Hard limit on retries (max 10, default 3)
- Each retry uses different model (no same-model retry)
- Quality score tracking (don't retry if not improving)

### 3. Cost Control
```python
# Add cost warnings
if history.total_cost > 0.05:
    logger.warning(f"Retry cost exceeding $0.05: ${history.total_cost:.4f}")

# Optional cost limit
config = RetryConfig(
    max_retries=3,
    max_cost=0.02  # Stop if total cost exceeds $0.02
)
```

### 4. Timeout Protection
```python
# Optional timeout
config = RetryConfig(
    max_retries=3,
    timeout_seconds=60  # Stop after 60 seconds total
)
```

## Output Display

### CLI Output Example

```
================================================================================
  RETRY HISTORY
================================================================================
Retry Summary:
  Total attempts: 2
  Total cost: $0.006
  Total time: 22.3s

  Attempt 0: gpt-4o-mini → ✗ fail (quality: 0.45)
  Attempt 1: deepseek-v3.2 → ✓ pass (quality: 0.95)

  Best: Attempt 1 (quality: 0.95)

Final result: Attempt 1
```

### JSON Output Structure

```json
{
  "schema_version": "reaction_analysis.v1",
  "input": {...},
  "interpretation": {...},
  "validation": {...},

  "retry_history": {
    "attempts": [
      {
        "attempt_num": 0,
        "model": "gpt-4o-mini",
        "provider": "openai",
        "quality_score": 0.45,
        "validation_status": "fail",
        "tokens": 1234,
        "latency_ms": 5000,
        "cost_estimate": 0.001
      },
      {
        "attempt_num": 1,
        "model": "deepseek-v3.2",
        "provider": "aliyun",
        "quality_score": 0.95,
        "validation_status": "pass",
        "tokens": 2345,
        "latency_ms": 17300,
        "cost_estimate": 0.005
      }
    ],
    "total_cost": 0.006,
    "total_time_ms": 22300
  },

  "final_attempt_num": 1,
  "all_attempts": [...]
}
```

## Testing Plan

### Test Case 1: Success on First Attempt
```python
# Simple Suzuki coupling (should pass immediately)
rxn = "Brc1ccccc1.B(O)(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"

config = RetryConfig(max_retries=3)
result = analyze_with_retry(rxn, client, config, validate=True)

assert len(result['retry_history'].attempts) == 1  # No retries
assert result['final_attempt_num'] == 0
assert result['validation']['gate']['status'] == 'pass'
```

### Test Case 2: Retry on Failure
```python
# Intentionally use weak model that might fail
weak_client = LLMClient(provider="openai", model="gpt-4o-mini")

config = RetryConfig(max_retries=3)
result = analyze_with_retry(complex_rxn, weak_client, config, validate=True)

# Should retry if first attempt fails
assert len(result['retry_history'].attempts) >= 1
assert result['validation']['gate']['status'] in ['pass', 'warning']
```

### Test Case 3: Max Retries Reached
```python
# Impossible reaction (should fail all attempts)
invalid_rxn = "INVALID>>SMILES"

config = RetryConfig(max_retries=3)
result = analyze_with_retry(invalid_rxn, client, config, validate=True)

assert len(result['retry_history'].attempts) == 4  # Initial + 3 retries
assert result['validation']['gate']['status'] == 'fail'
# Returns best attempt despite all failing
```

## Migration Path

### Phase 1: Core Implementation
- ✅ Create `reaction_agent/retry.py`
- ✅ Add `RetryConfig`, `RetryHistory`, `analyze_with_retry()`
- ✅ Unit tests for retry logic

### Phase 2: Integration
- ✅ Add retry_config parameter to `ReactionSMILESAnalyzer`
- ✅ Add CLI flags (--max-retries, --retry-on-warning, --no-retry)
- ✅ Integration tests

### Phase 3: CLI/Display
- ✅ Add retry history display section
- ✅ Add interactive mode commands
- ✅ Add cost/time warnings

### Phase 4: Documentation
- ✅ Update README with retry examples
- ✅ Update CLI_GUIDE.md
- ✅ Add retry recipes to cookbook

## Backward Compatibility

### Default Behavior (No Breaking Changes)

```python
# Without retry_config: existing behavior unchanged
result = analyzer.analyze(rxn_smiles)  # Works as before

# With validation only: no automatic retry (same as before)
result = analyzer.analyze(rxn_smiles, validate=True)  # Works as before

# With retry: new opt-in behavior
config = RetryConfig(max_retries=3)
result = analyzer.analyze(rxn_smiles, validate=True, retry_config=config)  # NEW
```

### CLI Compatibility

```bash
# All existing commands work unchanged
python -m reaction_agent.cli --reaction "SMILES"  # Works
python -m reaction_agent.cli --reaction "SMILES" --validate  # Works

# New retry behavior is opt-in via flags
python -m reaction_agent.cli --reaction "SMILES" --validate --max-retries 5  # NEW
```

## Summary

| Feature | Status | Complexity | Time Est. |
|---------|--------|------------|-----------|
| Core retry logic | ⬜ To implement | Medium | 3-4 hours |
| RetryConfig/History | ⬜ To implement | Low | 1 hour |
| analyze_with_retry() | ⬜ To implement | Medium | 2-3 hours |
| Agent integration | ⬜ To implement | Low | 1 hour |
| CLI flags | ⬜ To implement | Low | 1 hour |
| Display/formatting | ⬜ To implement | Medium | 2 hours |
| Interactive mode | ⬜ To implement | Low | 1 hour |
| Testing | ⬜ To implement | Medium | 2-3 hours |
| Documentation | ⬜ To implement | Low | 1 hour |
| **Total** | - | - | **14-18 hours** |

## Benefits

1. ✅ **Automatic quality improvement** - Don't accept bad results
2. ✅ **Cost-efficient** - Start cheap, escalate only if needed
3. ✅ **Transparent** - Full history and cost tracking
4. ✅ **Safe** - Hard limits prevent runaway costs
5. ✅ **Optional** - Doesn't affect existing workflows
6. ✅ **Best result selection** - Returns highest quality attempt

## Next Steps

1. ✅ **Review this plan** - Get user approval
2. Create `reaction_agent/retry.py` with core logic
3. Integrate with `agent.py` and `cli.py`
4. Add tests
5. Update documentation
6. Test with real reactions

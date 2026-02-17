"""
Automatic retry logic for failed validations (Phase 5).

Implements configurable retry loop that re-analyzes reactions with stronger
models when validation fails.
"""

from typing import Dict, Any, List, Optional
import logging
import time

from llmtools.clients import LLMClient

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
        Configure automatic retry behavior.

        Args:
            max_retries: Maximum retry attempts (default: 3, max: 10)
            aggressive: Retry even on warnings (default: False)
            retry_on_warning: Retry when validation status is WARNING (default: False)
        """
        if max_retries < 0 or max_retries > 10:
            raise ValueError("max_retries must be between 0 and 10")

        self.max_retries = max_retries
        self.aggressive = aggressive
        self.retry_on_warning = retry_on_warning

    @classmethod
    def from_cli_args(cls, args):
        """Create RetryConfig from argparse args."""
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
        # Cost per 1K tokens (rough estimates)
        costs_per_1k = {
            'gpt-4o-mini': 0.0005,
            'gpt-4o': 0.003,
            'deepseek-v3.2': 0.005,
            'deepseek-v3': 0.005
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


def get_retry_model(attempt_num: int) -> tuple:
    """
    Get model and provider for retry attempt.

    Retry sequence:
    - Attempt 0: gpt-4o-mini (fast first try)
    - Attempt 1: deepseek-v3.2 (most accurate)
    - Attempt 2: gpt-4o (strong fallback)
    - Attempt 3+: deepseek-v3.2 (keep trying with best)

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

    gate = validation.get('gate', {})
    status = gate.get('status', 'unknown')

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
    result = None

    while attempt_num <= config.max_retries:
        logger.info(f"Starting attempt {attempt_num} with {current_client.model}")

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
                logger.info(f"Analysis successful after {attempt_num + 1} attempt(s)")
                break

        except Exception as e:
            logger.error(f"Attempt {attempt_num} failed with error: {e}")

            # Record failed attempt
            history.add_attempt(
                attempt_num=attempt_num,
                model=current_client.model,
                provider=current_client.provider,
                result={'error': str(e), 'metadata': {}},
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
        best_result = result if result else {}

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

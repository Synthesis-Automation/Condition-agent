"""
Schema validation and build system for condition databases.

This module provides tools for validating and building unified indexes
from protocol (specific reactions) and rule (reaction families) databases.
"""

from .validator import (
    ConditionSourceValidator,
    ValidationReport,
    ValidationMessage,
    ValidationLevel,
    validate_batch
)

from .builder import (
    UnifiedIndexBuilder,
    BuildConfig,
    BuildReport
)

__all__ = [
    'ConditionSourceValidator',
    'ValidationReport',
    'ValidationMessage',
    'ValidationLevel',
    'validate_batch',
    'UnifiedIndexBuilder',
    'BuildConfig',
    'BuildReport',
]

"""
Centralized exception hierarchy for ChemTools.

This module provides a consistent exception structure for all ChemTools operations,
making error handling more predictable and easier to manage.

Usage:
    from chemtools.exceptions import ValidationError, DatabaseNotFoundError
    
    if not smiles:
        raise ValidationError("SMILES string cannot be empty")
    
    try:
        db = load_database(path)
    except FileNotFoundError as e:
        raise DatabaseNotFoundError(f"Database not found: {path}") from e
"""


class ChemToolsError(Exception):
    """Base exception for all ChemTools errors.
    
    All custom exceptions in ChemTools should inherit from this class.
    This allows catching all ChemTools-specific errors with a single except clause.
    """
    pass


class ValidationError(ChemToolsError):
    """Raised when input validation fails.
    
    Examples:
        - Empty or invalid SMILES strings
        - Invalid reaction format
        - Missing required parameters
        - Out-of-range values
    """
    pass


class DatabaseNotFoundError(ChemToolsError):
    """Raised when a database file or resource cannot be found.
    
    Examples:
        - Missing precedent database files
        - Missing reagent registry files
        - Missing rule database files
    """
    pass


class ProcessingError(ChemToolsError):
    """Raised when processing fails during computation.
    
    Examples:
        - Failed featurization
        - Failed similarity calculation
        - Failed recommendation generation
    """
    pass


class ConfigurationError(ChemToolsError):
    """Raised when configuration is invalid or missing.
    
    Examples:
        - Missing required environment variables
        - Invalid configuration values
        - Incompatible configuration combinations
    """
    pass


class ResourceNotAvailableError(ChemToolsError):
    """Raised when an optional resource or feature is not available.
    
    Examples:
        - RDKit not installed
        - Optional ML models not loaded
        - External services unavailable
    """
    pass


class ParseError(ChemToolsError):
    """Raised when parsing fails.
    
    Examples:
        - Invalid SMILES syntax
        - Invalid reaction format
        - Malformed data files
    """
    pass



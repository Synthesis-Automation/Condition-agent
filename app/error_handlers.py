"""
Centralized error handlers for FastAPI.

This module provides consistent error handling across all API endpoints,
ensuring that errors are returned in a standard format with appropriate
HTTP status codes.

Usage in main.py:
    from app.error_handlers import register_error_handlers
    
    app = FastAPI(...)
    register_error_handlers(app)
"""

import datetime
from typing import Dict, Any

from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse
from chemtools.exceptions import (
    ChemToolsError,
    ValidationError,
    DatabaseNotFoundError,
    ProcessingError,
    ConfigurationError,
    ResourceNotAvailableError,
    ParseError,
)


# Map exception types to HTTP status codes
ERROR_STATUS_MAP: Dict[type, int] = {
    ValidationError: 400,           # Bad Request
    ParseError: 400,                # Bad Request
    DatabaseNotFoundError: 404,     # Not Found
    ResourceNotAvailableError: 503, # Service Unavailable
    ConfigurationError: 500,        # Internal Server Error
    ProcessingError: 500,           # Internal Server Error
    ChemToolsError: 500,           # Internal Server Error (catch-all)
}


def create_error_response(
    error: Exception,
    status_code: int = 500,
    include_traceback: bool = False
) -> Dict[str, Any]:
    """Create standardized error response dictionary.
    
    Args:
        error: The exception that was raised
        status_code: HTTP status code
        include_traceback: Whether to include traceback (for debugging)
        
    Returns:
        Dictionary with error information
    """
    response = {
        "error": error.__class__.__name__,
        "message": str(error),
        "timestamp": datetime.datetime.utcnow().isoformat() + "Z",
        "status_code": status_code,
    }
    
    # Add traceback if requested (only for development)
    if include_traceback:
        import traceback
        response["traceback"] = traceback.format_exc()
    
    return response


async def chemtools_error_handler(request: Request, exc: ChemToolsError) -> JSONResponse:
    """Handle all ChemTools exceptions with appropriate status codes.
    
    Args:
        request: FastAPI request object
        exc: ChemTools exception that was raised
        
    Returns:
        JSONResponse with error details
    """
    status_code = ERROR_STATUS_MAP.get(type(exc), 500)
    
    # Log the error
    import logging
    logger = logging.getLogger("chemtools.api")
    
    if status_code >= 500:
        logger.error(f"{exc.__class__.__name__}: {exc}", exc_info=True)
    else:
        logger.warning(f"{exc.__class__.__name__}: {exc}")
    
    return JSONResponse(
        status_code=status_code,
        content=create_error_response(exc, status_code)
    )


async def validation_error_handler(request: Request, exc: ValidationError) -> JSONResponse:
    """Handle validation errors specifically.
    
    Args:
        request: FastAPI request object
        exc: ValidationError that was raised
        
    Returns:
        JSONResponse with validation error details
    """
    import logging
    logger = logging.getLogger("chemtools.api")
    logger.debug(f"Validation error: {exc}")
    
    return JSONResponse(
        status_code=400,
        content=create_error_response(exc, 400)
    )


async def generic_exception_handler(request: Request, exc: Exception) -> JSONResponse:
    """Catch-all handler for unexpected exceptions.
    
    Args:
        request: FastAPI request object
        exc: Any exception that wasn't caught by specific handlers
        
    Returns:
        JSONResponse with generic error message
    """
    import logging
    logger = logging.getLogger("chemtools.api")
    logger.exception(f"Unexpected error: {exc}")
    
    # Don't expose internal details in production
    return JSONResponse(
        status_code=500,
        content={
            "error": "InternalServerError",
            "message": "An unexpected error occurred. Please contact support if this persists.",
            "timestamp": datetime.datetime.utcnow().isoformat() + "Z",
            "status_code": 500,
        }
    )


def register_error_handlers(app: FastAPI) -> None:
    """Register all error handlers with the FastAPI application.
    
    Args:
        app: FastAPI application instance
    """
    # Register ChemTools exception handlers
    app.add_exception_handler(ChemToolsError, chemtools_error_handler)
    app.add_exception_handler(ValidationError, validation_error_handler)
    
    # Note: Generic exception handler is optional - FastAPI has built-in handling
    # Uncomment if you want custom handling for all exceptions:
    # app.add_exception_handler(Exception, generic_exception_handler)

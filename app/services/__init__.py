"""
Service Layer - Business logic for Chemistry Tools API.

This module separates business logic from API routing, making the code:
- More testable (services can be tested without FastAPI)
- More reusable (services can be used in CLI, UI, or other contexts)
- Easier to maintain (single responsibility per service)

Services:
- matching_service: SMILES normalization, family detection, reaction type detection
- featurization_service: Molecular and role-aware featurization
- recommendation_service: Unified condition recommendations
- precedent_service: Precedent search, core management
"""

# Import services - these can be imported individually or via this package
# Example: from app.services import matching_service
# Example: from app.services.matching_service import normalize_smiles

__all__ = [
    "matching_service",
    "featurization_service",
    "recommendation_service",
    "precedent_service",
]

"""
Base interface for reaction family-specific analyzers.

Provides abstract base class that all family analyzers must implement.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Dict, Any, List, Optional


class FamilyAnalyzer(ABC):
    """
    Abstract base class for reaction family-specific analysis.
    
    Each reaction family (Amide, C-N coupling, Suzuki, etc.) should
    have its own analyzer that implements these methods.
    """
    
    @abstractmethod
    def analyze_substrates(
        self,
        reactants: List[str],
        features: Dict[str, Any],
    ) -> Dict[str, Any]:
        """
        Analyze substrates for this specific reaction family.
        
        Args:
            reactants: List of reactant SMILES
            features: Molecular features from featurizers
            
        Returns:
            Dict with substrate analysis results
        """
        pass
    
    @abstractmethod
    def build_features(
        self,
        role_pack: Dict[str, Any],
        features: Dict[str, Any],
    ) -> Dict[str, Any]:
        """
        Build family-specific features for precedent search.
        
        Args:
            role_pack: Role assignments from role featurization
            features: Molecular features
            
        Returns:
            Dict of features for precedent k-NN search
        """
        pass
    
    def get_default_conditions(self) -> Dict[str, Any]:
        """
        Get default conditions for this family.
        
        Returns:
            Dict with default T_C, time_h, etc.
        """
        return {
            "T_C": 80.0,
            "time_h": 12.0,
        }
    
    def get_family_name(self) -> str:
        """
        Get the canonical family name.
        
        Returns:
            Family name string
        """
        return self.__class__.__name__.replace("Analyzer", "")

from pydantic import BaseModel
from typing import List, Optional, Dict, Any

class NormalizeRequest(BaseModel): smiles: str
class DetectFamilyRequest(BaseModel): reactants: List[str]
class FeaturizeUllmannRequest(BaseModel): electrophile: str; nucleophile: str

class Reagent(BaseModel):
    uid: str; role: str
    name: Optional[str] = None; token: Optional[str] = None

class ConditionCoreParseRequest(BaseModel):
    reagents: List[Reagent]; text: Optional[str] = None

class PropertiesLookupRequest(BaseModel): query: str
class PrecedentKNNRequest(BaseModel):
    family: str; features: Dict[str, Any]; k: int = 50; relax: Optional[Dict[str, Any]] = None
class ConstraintsFilterRequest(BaseModel):
    candidates: List[str]; rules: Optional[Dict[str, Any]] = None
class ExplainPrecedentsRequest(BaseModel): pack: Dict[str, Any]; features: Dict[str, Any]
class RecommendFromReactionRequest(BaseModel):
    reaction: str
    k: int = 25
    relax: Optional[Dict[str, Any]] = None
    constraints: Optional[Dict[str, Any]] = None

class PlateDesignRequest(BaseModel):
    reaction: str
    plate_size: int = 24
    relax: Optional[Dict[str, Any]] = None
    constraints: Optional[Dict[str, Any]] = None

# Dev/validation: validate ConditionCore normalization on a JSONL dataset
class ConditionCoreValidateRequest(BaseModel):
    path: str
    limit: int = 0
    show_mismatches: int = 10
    metal_only_ok: bool = True

# Role-aware featurization
class RoleAwareMolRequest(BaseModel):
    smiles: str
    roles: Optional[List[str]] = None  # e.g., ["amine", "aryl_halide"]

class RoleAwareReactionRequest(BaseModel):
    reaction: str

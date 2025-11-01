from pydantic import BaseModel
from typing import List, Optional, Dict, Any

class NormalizeRequest(BaseModel): smiles: str
class DetectFamilyRequest(BaseModel): reactants: List[str]
class DetectTypeRequest(BaseModel): reaction: str
class FeaturizeUllmannRequest(BaseModel): electrophile: str; nucleophile: str

class Reagent(BaseModel):
    uid: str; role: str
    name: Optional[str] = None; token: Optional[str] = None

class ConditionCoreParseRequest(BaseModel):
    reagents: List[Reagent]; text: Optional[str] = None

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
    rerank_strategy: str = 'rule'  # 'rule', 'analytics', or 'none'
    filter_unknown_reagents: bool = False
    search_all_families: bool = False  # Enable cross-family search
    reaction_type_threshold: float = 0.15  # Min representation for reaction type filtering (15%)
    mechanism_similarity_threshold: float = 0.4  # Min mechanism similarity (40%)
    mechanism_weight: float = 0.3  # Weight for mechanism-enhanced similarity (30%)

class RecommendConditionsRequest(BaseModel):
    reaction: str
    reaction_type: Optional[str] = None
    k: int = 50
    limit: int = 5
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


class SchemeMatchRequest(BaseModel):
    reaction: str
    db: Optional[str] = None
    include_trace: bool = True
    relax: Optional[Dict[str, Any]] = None  # For catalyst_class filtering


# Core search (by condition core, e.g., 'Pd/XPhos')
class CoreSearchRequest(BaseModel):
    core: str
    family: Optional[str] = None
    fuzzy: bool = True
    limit: int = 50

"""Interactive CLI for reaction condition recommendations with LLM-powered input parsing.

This CLI uses natural language processing via LLM to convert user requirements
into structured API calls, with progressive refinement until inputs validate.

Only a valid reaction SMILES is required to complete - all constraints are optional.
"""

from __future__ import annotations

import sys
import json
import logging
from pathlib import Path
from typing import Any, Optional
from dataclasses import dataclass, field

# ANSI color codes for terminal output
class Colors:
    """ANSI color codes for colored terminal output."""
    HEADER = '\033[95m'      # Magenta
    BLUE = '\033[94m'        # Blue
    CYAN = '\033[96m'        # Cyan
    GREEN = '\033[92m'       # Green
    YELLOW = '\033[93m'      # Yellow
    RED = '\033[91m'         # Red
    ENDC = '\033[0m'         # Reset
    BOLD = '\033[1m'         # Bold
    UNDERLINE = '\033[4m'    # Underline
    
    @staticmethod
    def disable():
        """Disable colors (for non-terminal output)."""
        Colors.HEADER = ''
        Colors.BLUE = ''
        Colors.CYAN = ''
        Colors.GREEN = ''
        Colors.YELLOW = ''
        Colors.RED = ''
        Colors.ENDC = ''
        Colors.BOLD = ''
        Colors.UNDERLINE = ''

# Add project root to path for imports
PROJECT_ROOT = Path(__file__).parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from llmtools.clients import LLMClient
from llmtools.prompts import PromptTemplate

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# ============================================================================
# JSON Schema for Structured Outputs
# ============================================================================

RECOMMENDATION_REQUEST_SCHEMA = {
    "type": "object",
    "properties": {
        "reaction_smiles": {
            "type": "string",
            "description": "Valid SMILES string for the reaction (reactants>>products format)"
        },
        "reaction_smiles_is_valid": {
            "type": "boolean",
            "description": "Whether the reaction SMILES appears to be in valid format"
        },
        "reaction_type": {
            "type": ["string", "null"],
            "description": "Detected or specified reaction type (e.g., Suzuki, C_N_Coupling, C_O_Coupling, C_S_Coupling)"
        },
        "constraints": {
            "type": "object",
            "description": "Chemical constraints extracted from user requirements",
            "properties": {
                "temperature": {
                    "type": "object",
                    "properties": {
                        "max": {"type": ["number", "null"], "description": "Maximum temperature in Celsius"},
                        "min": {"type": ["number", "null"], "description": "Minimum temperature in Celsius"}
                    }
                },
                "base_strength": {
                    "type": ["string", "null"],
                    "enum": ["weak", "moderate", "strong", "any", None],
                    "description": "Maximum allowed base strength"
                },
                "required_reagents": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "List of reagent names/types that must be included (e.g., 'copper catalyst', 'CuI')"
                },
                "exclude_reagents": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "List of reagent names/CAS to exclude"
                },
                "exclude_roles": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Reagent roles to exclude (e.g., oxidant, reductant)"
                },
                "metal_preference": {
                    "type": ["string", "null"],
                    "enum": ["Pd", "Cu", "Ni", "Fe", "Ru", "Rh", "Ir", "Au", "Ag", "any", None],
                    "description": "Preferred metal for catalysis (e.g., 'Cu' for copper, 'Pd' for palladium)"
                },
                "solvent_preference": {
                    "type": ["string", "null"],
                    "enum": ["polar_aprotic", "polar_protic", "nonpolar", "aqueous", "any", None],
                    "description": "Preferred solvent polarity"
                },
                "air_sensitive": {
                    "type": ["boolean", "null"],
                    "description": "Whether air-sensitive conditions are acceptable (null = no preference)"
                },
                "cost_level": {
                    "type": ["string", "null"],
                    "enum": ["low", "medium", "high", "any", None],
                    "description": "Maximum cost level acceptable"
                }
            }
        },
        "additional_notes": {
            "type": ["string", "null"],
            "description": "Any additional user notes or clarifications"
        },
        "validation_issues": {
            "type": "array",
            "items": {"type": "string"},
            "description": "List of validation issues found (empty if valid)"
        },
        "clarification_needed": {
            "type": "array",
            "items": {"type": "string"},
            "description": "List of ambiguous points requiring user clarification"
        }
    },
    "required": ["reaction_smiles", "reaction_smiles_is_valid", "validation_issues"]
}


# ============================================================================
# LLM Prompt Templates
# ============================================================================

PARSE_USER_INPUT_PROMPT = PromptTemplate(
    name="parse_user_input",
    template="""You are a chemistry lab assistant helping users specify reaction conditions.

The user has provided input that may contain:
1. A reaction SMILES string (in reactants>>products format)
2. Natural language requirements/constraints
3. Or both mixed together

USER INPUT:
{user_input}

INSTRUCTIONS:

**STEP 1: EXTRACT REACTION SMILES**
- Look for a SMILES string with ">>" separator (e.g., "Br.c1ccccc1>>c1ccc(...)cc1")
- The SMILES may appear anywhere in the input (beginning, middle, or end)
- The SMILES may be surrounded by text like "I want to run", "using", "with", etc.
- Extract the complete SMILES string (everything from first reactant to last product)
- If multiple ">>" patterns exist, take the most complete one
- Set "reaction_smiles_is_valid" to true only if:
  * Has ">>" separator
  * Has non-empty reactants (left side of >>)
  * Appears to be valid SMILES format

**STEP 2: EXTRACT REQUIREMENTS FROM SURROUNDING TEXT**
- All text that is NOT part of the SMILES string may contain requirements
- Extract constraints from this text (all optional):
  * Temperature: "no high temperature" → max: 80, "room temperature" → max: 30, "RT" → max: 25
  * Base strength: "no strong base" → "moderate", "mild base" → "weak", "strong base ok" → "strong"
  * **Required reagents**: "use X", "with X", "prefer X", "need X" → add to required_reagents array
  * **Metal preference**: "copper catalyst" → "Cu", "palladium" → "Pd", "nickel" → "Ni", "use Cu" → "Cu"
  * Excluded reagents: "no X", "avoid X", "exclude X", "without X" → add to exclude_reagents array
  * Solvent: "DMF", "DMSO" → "polar_aprotic"; "methanol", "water" → "polar_protic"; "toluene" → "nonpolar"
  * Air sensitivity: "inert atmosphere", "glovebox", "Schlenk" → true; "air stable", "bench" → false
  * Cost: "cheap", "inexpensive" → "low"; "expensive ok" → "high"; "moderate cost" → "medium"

**STEP 3: VALIDATION**
- Report in "validation_issues" array:
  * "Invalid reaction SMILES format" - if SMILES missing or invalid
  * Any other blocking issues
- If no SMILES found: add "No reaction SMILES found in input"
- If SMILES invalid format: add "Invalid reaction SMILES format"

**STEP 4: CLARIFICATIONS** (only for ambiguous constraints)
- List truly ambiguous requirements in "clarification_needed"
- Do NOT ask for clarification if user provided no requirements - this is fine!
- Only ask about genuinely unclear constraints

**EXAMPLES:**

Input: "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1 use copper catalyst"
→ reaction_smiles: "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
→ constraints: {{"required_reagents": ["copper catalyst"], "metal_preference": "Cu"}}

Input: "I want to run Ic1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1 with no strong base and room temperature"
→ reaction_smiles: "Ic1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
→ constraints: {{"base_strength": "moderate", "temperature": {{"max": 30}}}}

Input: "use palladium for CCBr.c1ccccc1>>CCc1ccccc1 avoid expensive reagents"
→ reaction_smiles: "CCBr.c1ccccc1>>CCc1ccccc1"
→ constraints: {{"metal_preference": "Pd", "required_reagents": ["palladium"], "cost_level": "low"}}

IMPORTANT:
- Intelligently separate SMILES from requirements text
- Only the reaction SMILES must be valid - all constraints are optional
- Empty constraints object is perfectly acceptable
- Be permissive - use null for unclear optional values

OUTPUT FORMAT:
Return ONLY a valid JSON object (no markdown, no explanation) matching this schema:
{schema}
""",
    default_params={
        "temperature": 0.2,
        "max_tokens": 1500,
    }
)


CLARIFICATION_PROMPT = PromptTemplate(
    name="clarification_request",
    template="""You are helping a user refine their reaction condition requirements.

CURRENT PARSED STATE:
{current_state}

ISSUES/CLARIFICATIONS FROM PREVIOUS PARSE:
Validation Issues: {validation_issues}
Clarifications Needed: {clarifications}

USER'S RESPONSE:
{user_response}

INSTRUCTIONS:
Update the parsed state based on the user's response. 

1. If the user provided a new/corrected reaction SMILES, update it and revalidate
2. If the user clarified constraints, update the constraints object
3. If the user said "skip" or "none" for optional constraints, leave them as null
4. Update validation_issues - should be empty if reaction SMILES is now valid
5. Update clarification_needed - should be empty if all ambiguities are resolved

VALIDATION:
- reaction_smiles_is_valid should be true ONLY if format is correct (has ">>" and non-empty reactants)
- validation_issues should list any blocking problems
- If only reaction SMILES is valid and no constraints provided, that's a successful state!

OUTPUT FORMAT:
Return ONLY a valid JSON object (no markdown) matching the same schema as before:
{schema}
""",
    default_params={
        "temperature": 0.2,
        "max_tokens": 1500,
    }
)


# ============================================================================
# Reaction Type Determination
# ============================================================================

def determine_final_reaction_type(reaction_smiles: str, initial_type: Optional[str], constraints: dict[str, Any]) -> dict[str, Any]:
    """Determine final reaction type considering constraints.
    
    This routes to the appropriate recommendation system:
    - C_N_Coupling -> Unified dataset (metal preference via constraints)
    - C_O_Coupling -> C-O coupling dataset
    - C_S_Coupling -> C-S coupling dataset
    - Other reactions -> Use detected or provided type
    
    Returns dict with: final_type, detection_info, routing
    """
    try:
        from chemtools import detect_reaction
        from chemtools.recommend.utils import canonical_family
    except ImportError:
        logger.warning("Could not import detection API - using initial type as-is")
        return {
            "final_type": initial_type or "Unknown",
            "detection_info": {"source": "user_provided"},
            "routing": "fallback"
        }
    
    # Detect reaction family from SMILES using unified API with ML
    detection = detect_reaction(reaction_smiles, use_ml=True)
    detected_family = detection.get("family", "Unknown")
    confidence = detection.get("confidence", 0.0)
    
    # Start with detected or user-provided type
    final_type = initial_type or detected_family
    routing = "auto-detected"
    
    # Normalize to canonical family name
    try:
        final_type = canonical_family(final_type)
    except Exception:
        pass
    
    return {
        "final_type": final_type,
        "initial_type": initial_type,
        "detected_family": detected_family,
        "confidence": confidence,
        "routing": routing,
        "detection_info": detection,
    }


# ============================================================================
# Data Models
# ============================================================================

@dataclass
class ParsedRequest:
    """Structured representation of a parsed recommendation request."""
    reaction_smiles: str
    reaction_smiles_is_valid: bool = False
    reaction_type: Optional[str] = None
    constraints: dict[str, Any] = field(default_factory=dict)
    additional_notes: Optional[str] = None
    validation_issues: list[str] = field(default_factory=list)
    clarification_needed: list[str] = field(default_factory=list)
    
    def is_valid(self) -> bool:
        """Check if request is ready for API submission.
        
        Only requires valid reaction SMILES - all constraints are optional.
        """
        return (
            bool(self.reaction_smiles) and
            self.reaction_smiles_is_valid and
            len(self.validation_issues) == 0
        )
    
    def to_api_request(self) -> dict[str, Any]:
        """Convert to API request format."""
        # Clean up constraints - remove None values
        clean_constraints = {}
        for key, value in self.constraints.items():
            if value is not None and value != {}:
                if isinstance(value, dict):
                    # Remove None values from nested dicts
                    cleaned = {k: v for k, v in value.items() if v is not None}
                    if cleaned:
                        clean_constraints[key] = cleaned
                elif isinstance(value, list) and len(value) > 0:
                    clean_constraints[key] = value
                elif value not in [None, "", "any"]:
                    clean_constraints[key] = value
        
        # Convert metal_preference to required_reagents format for API
        if "metal_preference" in clean_constraints and clean_constraints["metal_preference"]:
            metal = clean_constraints["metal_preference"]
            if "required_reagents" not in clean_constraints:
                clean_constraints["required_reagents"] = []
            # Add metal to required reagents if not already there
            metal_terms = [metal.lower(), f"{metal.lower()} catalyst", f"{metal} catalyst"]
            if not any(term in str(clean_constraints.get("required_reagents", [])).lower() for term in metal_terms):
                clean_constraints["required_reagents"].append(f"{metal} catalyst")
        
        return {
            "reaction": self.reaction_smiles,
            "reaction_type": self.reaction_type,
            "k": 10,
            "limit": 5,
            "constraints": clean_constraints if clean_constraints else {},
            "relax": {},
        }
    
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ParsedRequest:
        """Create from dictionary."""
        return cls(
            reaction_smiles=data.get("reaction_smiles", ""),
            reaction_smiles_is_valid=data.get("reaction_smiles_is_valid", False),
            reaction_type=data.get("reaction_type"),
            constraints=data.get("constraints", {}),
            additional_notes=data.get("additional_notes"),
            validation_issues=data.get("validation_issues", []),
            clarification_needed=data.get("clarification_needed", []),
        )


# ============================================================================
# LLM Parser
# ============================================================================

class NaturalLanguageParser:
    """Parser that converts natural language to structured API requests using LLM."""
    
    def __init__(
        self,
        provider: str = "aliyun",
        model: str = "deepseek-v3.2",
        api_key: Optional[str] = None,
    ):
        self.client: LLMClient = LLMClient(
            provider=provider,
            model=model,
            api_key=api_key,
            temperature=0.2,
            max_tokens=1500,
        )
        self.model = model
        self.schema_str = json.dumps(RECOMMENDATION_REQUEST_SCHEMA, indent=2)
    
    def parse_initial_input(
        self,
        user_input: str,
    ) -> ParsedRequest:
        """Parse initial user input (mixed SMILES and requirements) into structured format."""
        logger.info(f"Parsing user input via LLM ({self.model})...")
        
        prompt = PARSE_USER_INPUT_PROMPT.format(
            user_input=user_input,
            schema=self.schema_str,
        )
        
        try:
            response = self.client.chat(
                prompt=prompt,
                temperature=0.2,
                max_tokens=1500,
            )
            
            content = response.content.strip()
            
            # Extract JSON from markdown code blocks if present
            if "```json" in content:
                content = content.split("```json")[1].split("```")[0].strip()
            elif "```" in content:
                content = content.split("```")[1].split("```")[0].strip()
            
            parsed_data = json.loads(content)
            result = ParsedRequest.from_dict(parsed_data)
            
            logger.info(f"Parse complete. Valid: {result.is_valid()}, "
                       f"Issues: {len(result.validation_issues)}, "
                       f"Clarifications: {len(result.clarification_needed)}, "
                       f"Detected type: {result.reaction_type or 'None'}")
            
            return result
            
        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse LLM response as JSON: {e}")
            logger.debug(f"LLM response: {content}")
            # Return a fallback with error
            return ParsedRequest(
                reaction_smiles="",
                reaction_smiles_is_valid=False,
                validation_issues=[f"LLM parsing error: {str(e)}"]
            )
        except Exception as e:
            logger.error(f"Error during LLM parsing: {e}")
            return ParsedRequest(
                reaction_smiles="",
                reaction_smiles_is_valid=False,
                validation_issues=[f"Error: {str(e)}"]
            )
    
    def refine_with_clarification(
        self,
        current_state: ParsedRequest,
        user_response: str,
    ) -> ParsedRequest:
        """Refine parsed request with user clarifications."""
        logger.info("Refining request with user clarifications...")
        
        prompt = CLARIFICATION_PROMPT.format(
            current_state=json.dumps({
                "reaction_smiles": current_state.reaction_smiles,
                "reaction_smiles_is_valid": current_state.reaction_smiles_is_valid,
                "reaction_type": current_state.reaction_type,
                "constraints": current_state.constraints,
                "validation_issues": current_state.validation_issues,
                "clarification_needed": current_state.clarification_needed,
            }, indent=2),
            validation_issues="\n".join([f"- {v}" for v in current_state.validation_issues]) or "None",
            clarifications="\n".join([f"- {c}" for c in current_state.clarification_needed]) or "None",
            user_response=user_response,
            schema=self.schema_str,
        )
        
        try:
            response = self.client.chat(
                prompt=prompt,
                temperature=0.2,
                max_tokens=1500,
            )
            
            content = response.content.strip()
            
            # Extract JSON
            if "```json" in content:
                content = content.split("```json")[1].split("```")[0].strip()
            elif "```" in content:
                content = content.split("```")[1].split("```")[0].strip()
            
            parsed_data = json.loads(content)
            result = ParsedRequest.from_dict(parsed_data)
            
            logger.info(f"Refinement complete. Valid: {result.is_valid()}, "
                       f"Issues: {len(result.validation_issues)}")
            
            return result
            
        except Exception as e:
            logger.error(f"Error during refinement: {e}")
            # Return current state with error appended
            current_state.validation_issues.append(f"Refinement error: {str(e)}")
            return current_state


# ============================================================================
# Interactive CLI
# ============================================================================

class InteractiveCLI:
    """Interactive command-line interface for reaction recommendations."""
    
    def __init__(self, parser: NaturalLanguageParser, test_mode: bool = False):
        self.parser = parser
        self.current_request: Optional[ParsedRequest] = None
        self.test_mode = test_mode
    
    def print_header(self):
        """Print CLI header."""
        print("\n" + f"{Colors.CYAN}{'=' * 70}{Colors.ENDC}")
        print(f"{Colors.BOLD}{Colors.HEADER}  🧪 CHEMISTRY CONDITION RECOMMENDATION CLI{Colors.ENDC}")
        print(f"{Colors.CYAN}  Powered by LLM-assisted natural language parsing{Colors.ENDC}")
        if self.test_mode:
            print(f"{Colors.YELLOW}  ⚠️  TEST MODE: Will stop before actual API submission{Colors.ENDC}")
        print(f"{Colors.CYAN}{'=' * 70}{Colors.ENDC}\n")
    
    def print_separator(self):
        """Print section separator."""
        print("\n" + f"{Colors.BLUE}{'-' * 70}{Colors.ENDC}\n")
    
    def get_initial_input(self) -> str:
        """Get combined reaction SMILES and requirements from user in one input."""
        print(f"{Colors.BOLD}Please provide your reaction:{Colors.ENDC}\n")
        
        print(f"{Colors.CYAN}Enter reaction SMILES with optional requirements (all in one line):{Colors.ENDC}\n")
        
        print(f"{Colors.YELLOW}Format options:{Colors.ENDC}")
        print(f"{Colors.YELLOW}  • Just SMILES: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1{Colors.ENDC}")
        print(f"{Colors.YELLOW}  • SMILES + requirements: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1 use copper catalyst{Colors.ENDC}")
        print(f"{Colors.YELLOW}  • Requirements + SMILES: I want to run Ic1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1 with no strong base{Colors.ENDC}")
        print(f"{Colors.YELLOW}  • Mixed: use palladium for CCBr.c1ccccc1>>CCc1ccccc1 at room temperature{Colors.ENDC}\n")
        
        print(f"{Colors.GREEN}Examples of requirements:{Colors.ENDC}")
        print(f"{Colors.GREEN}  • use copper catalyst{Colors.ENDC}")
        print(f"{Colors.GREEN}  • no strong base, room temperature{Colors.ENDC}")
        print(f"{Colors.GREEN}  • avoid expensive catalysts, prefer DMF{Colors.ENDC}")
        print(f"{Colors.GREEN}  • with palladium, mild conditions{Colors.ENDC}\n")
        
        user_input = input(f"{Colors.BOLD}{Colors.CYAN}> {Colors.ENDC}").strip()
        
        if not user_input:
            print(f"{Colors.RED}❌ Error: Input cannot be empty.{Colors.ENDC}")
            sys.exit(1)
        
        return user_input
    
    def display_parsed_state(self, request: ParsedRequest):
        """Display current parsed state to user."""
        self.print_separator()
        print(f"{Colors.BOLD}{Colors.HEADER}📋 PARSED REQUEST:{Colors.ENDC}\n")
        
        # Reaction SMILES status
        if request.reaction_smiles_is_valid:
            print(f"{Colors.GREEN}✅ Reaction: {Colors.ENDC}{request.reaction_smiles}")
        else:
            print(f"{Colors.RED}❌ Reaction: {Colors.ENDC}{request.reaction_smiles}")
        
        if request.reaction_type:
            print(f"{Colors.CYAN}   Type: {request.reaction_type}{Colors.ENDC}")
        
        # Constraints (if any)
        if request.constraints:
            has_constraints = any(
                v is not None and v != {} and (not isinstance(v, list) or len(v) > 0)
                for v in request.constraints.values()
            )
            
            if has_constraints:
                print(f"\n{Colors.CYAN}📌 Constraints:{Colors.ENDC}")
                for key, value in request.constraints.items():
                    if value is None or value == {} or (isinstance(value, list) and len(value) == 0):
                        continue
                    
                    if isinstance(value, dict):
                        print(f"{Colors.YELLOW}  • {key}:{Colors.ENDC}")
                        for k, v in value.items():
                            if v is not None:
                                print(f"      {k}: {v}")
                    elif isinstance(value, list):
                        print(f"{Colors.YELLOW}  • {key}:{Colors.ENDC} {', '.join(str(v) for v in value)}")
                    else:
                        print(f"{Colors.YELLOW}  • {key}:{Colors.ENDC} {value}")
        
        if request.additional_notes:
            print(f"\n{Colors.CYAN}💡 Notes:{Colors.ENDC} {request.additional_notes}")
        
        # Validation status
        if request.is_valid():
            print(f"\n{Colors.GREEN}{Colors.BOLD}✅ Request is VALID and ready to submit!{Colors.ENDC}")
        else:
            if request.validation_issues:
                print(f"\n{Colors.RED}⚠️  VALIDATION ISSUES:{Colors.ENDC}")
                for issue in request.validation_issues:
                    print(f"{Colors.RED}  • {issue}{Colors.ENDC}")
    
    def request_clarifications(self, request: ParsedRequest) -> Optional[str]:
        """Request clarifications or fixes from user."""
        self.print_separator()
        
        # Show validation issues first (blocking)
        if request.validation_issues:
            print(f"{Colors.RED}{Colors.BOLD}❌ VALIDATION ISSUES (must fix):{Colors.ENDC}\n")
            for i, issue in enumerate(request.validation_issues, 1):
                print(f"{Colors.RED}{i}. {issue}{Colors.ENDC}")
            print()
        
        # Show optional clarifications
        if request.clarification_needed:
            print(f"{Colors.YELLOW}{Colors.BOLD}⚠️  CLARIFICATIONS NEEDED (optional):{Colors.ENDC}\n")
            for i, clarification in enumerate(request.clarification_needed, 1):
                print(f"{Colors.YELLOW}{i}. {clarification}{Colors.ENDC}")
            print()
        
        # Prompt for response
        if request.validation_issues:
            print(f"{Colors.CYAN}Please fix the validation issues above:{Colors.ENDC}")
        elif request.clarification_needed:
            print(f"{Colors.CYAN}Please provide clarifications (or press Enter to skip):{Colors.ENDC}")
        
        response = input(f"{Colors.GREEN}> {Colors.ENDC}").strip()
        
        return response if response else None
    
    def display_reaction_type_determination(self, request: ParsedRequest) -> dict[str, Any]:
        """Determine and display final reaction type before submission."""
        print(f"\n{Colors.BOLD}{Colors.HEADER}🔬 REACTION TYPE DETERMINATION{Colors.ENDC}\n")
        
        # Determine final reaction type
        determination = determine_final_reaction_type(
            reaction_smiles=request.reaction_smiles,
            initial_type=request.reaction_type,
            constraints=request.constraints,
        )
        
        final_type = determination.get("final_type")
        detected_family = determination.get("detected_family")
        confidence = determination.get("confidence", 0.0)
        detected_metal = determination.get("detected_metal")
        routing = determination.get("routing")
        
        # Display detection results
        print(f"{Colors.CYAN}Initial/Provided Type:{Colors.ENDC} {request.reaction_type or 'None'}")
        print(f"{Colors.CYAN}Detected Family:{Colors.ENDC} {detected_family} (confidence: {confidence:.2f})")
        
        if detected_metal:
            print(f"{Colors.YELLOW}Detected Metal Catalyst:{Colors.ENDC} {detected_metal}")
        
        # Show final determination
        print(f"\n{Colors.GREEN}{Colors.BOLD}✓ Final Reaction Type:{Colors.ENDC} {Colors.BOLD}{final_type}{Colors.ENDC}")
        print(f"{Colors.CYAN}Routing:{Colors.ENDC} {routing}")
        
        # Determine recommendation system
        is_rule_based = False
        is_ml_based = False
        
        if "rule-based" in routing.lower() or "Cu" in str(final_type):
            is_rule_based = True
        elif "ML" in routing or "Pd" in str(final_type) or "Ni" in str(final_type):
            is_ml_based = True
        
        # Show system being used
        if is_rule_based:
            print(f"{Colors.HEADER}🎯 Recommendation System:{Colors.ENDC} {Colors.BOLD}RULE-BASED{Colors.ENDC} (deterministic constraints)")
        elif is_ml_based:
            print(f"{Colors.HEADER}🎯 Recommendation System:{Colors.ENDC} {Colors.BOLD}MACHINE LEARNING{Colors.ENDC} (similarity-based)")
        
        # Explain routing logic for C-N coupling
        if "C_N_Coupling" in str(final_type):
            print(f"\n{Colors.YELLOW}ℹ️  C-N Coupling Routing:{Colors.ENDC}")
            if "Cu" in final_type:
                print(f"{Colors.YELLOW}   → Ullmann reaction (copper-catalyzed){Colors.ENDC}")
                print(f"{Colors.YELLOW}   → Uses rule-based constraint matching{Colors.ENDC}")
            elif "Pd" in final_type:
                print(f"{Colors.YELLOW}   → Buchwald-Hartwig reaction (palladium-catalyzed){Colors.ENDC}")
                print(f"{Colors.YELLOW}   → Uses ML-based similarity search{Colors.ENDC}")
            elif "Ni" in final_type:
                print(f"{Colors.YELLOW}   → Nickel-catalyzed C-N coupling{Colors.ENDC}")
                print(f"{Colors.YELLOW}   → Uses ML-based similarity search{Colors.ENDC}")
            else:
                print(f"{Colors.YELLOW}   → Generic C-N coupling (auto-routing){Colors.ENDC}")
        
        return determination
    
    def confirm_submission(self, request: ParsedRequest) -> bool:
        """Ask user to confirm submission."""
        self.print_separator()
        print(f"{Colors.GREEN}{Colors.BOLD}✅ Request is ready for submission!{Colors.ENDC}\n")
        
        # Determine final reaction type
        determination = self.display_reaction_type_determination(request)
        
        # Update request with final type
        final_type = determination.get("final_type")
        if final_type and final_type != "Unknown":
            request.reaction_type = final_type
        
        # Show summary
        self.print_separator()
        print(f"{Colors.BOLD}📋 SUBMISSION SUMMARY{Colors.ENDC}\n")
        print(f"{Colors.CYAN}Reaction:{Colors.ENDC} {request.reaction_smiles}")
        print(f"{Colors.CYAN}Type:{Colors.ENDC} {request.reaction_type or 'auto-detect'}")
        
        # Count meaningful constraints
        constraint_count = sum(
            1 for v in request.constraints.values()
            if v is not None and v != {} and (not isinstance(v, list) or len(v) > 0)
        )
        
        if constraint_count > 0:
            print(f"{Colors.CYAN}Constraints:{Colors.ENDC} {constraint_count} specified")
        else:
            print(f"{Colors.YELLOW}Constraints:{Colors.ENDC} None (will use default recommendations)")
        
        print(f"\n{Colors.BOLD}🔍 API Request Preview:{Colors.ENDC}")
        print(f"{Colors.YELLOW}{json.dumps(request.to_api_request(), indent=2)}{Colors.ENDC}")
        
        # TEST MODE: Stop here without submission
        if self.test_mode:
            print(f"\n{Colors.YELLOW}{'=' * 70}{Colors.ENDC}")
            print(f"{Colors.YELLOW}{Colors.BOLD}⚠️  TEST MODE: Stopping before actual submission{Colors.ENDC}")
            print(f"{Colors.GREEN}✅ Input process validation complete!{Colors.ENDC}")
            print(f"{Colors.CYAN}📋 The request is valid and ready to be submitted in production mode.{Colors.ENDC}")
            print(f"{Colors.YELLOW}{'=' * 70}{Colors.ENDC}")
            return False
        
        print(f"\n{Colors.CYAN}{'=' * 70}{Colors.ENDC}")
        confirm = input(f"{Colors.BOLD}Submit this request? (yes/no):{Colors.ENDC} ").strip().lower()
        return confirm in ["yes", "y"]
    
    def run(self):
        """Run the interactive CLI loop."""
        self.print_header()
        
        # Get initial input (combined SMILES + requirements)
        user_input = self.get_initial_input()
        
        # Parse with LLM
        print(f"\n{Colors.CYAN}🤖 Parsing input with LLM...{Colors.ENDC}")
        try:
            self.current_request = self.parser.parse_initial_input(user_input)
            print(f"{Colors.GREEN}✅ LLM parsing complete{Colors.ENDC}")
        except Exception as e:
            print(f"\n{Colors.RED}❌ Error parsing input: {e}{Colors.ENDC}")
            sys.exit(1)
        
        # Display initial parse
        self.display_parsed_state(self.current_request)
        
        # Refinement loop (only if validation issues or user wants to clarify)
        max_iterations = 5
        iteration = 0
        
        while not self.current_request.is_valid() and iteration < max_iterations:
            iteration += 1
            
            # Request fixes/clarifications
            clarification = self.request_clarifications(self.current_request)
            
            if clarification is None:
                if not self.current_request.validation_issues:
                    # No blocking issues, user skipped optional clarifications
                    break
                else:
                    # Has blocking issues but user didn't provide fix
                    print(f"\n{Colors.YELLOW}⚠️  Cannot proceed without fixing validation issues.{Colors.ENDC}")
                    print(f"{Colors.CYAN}Please provide the corrected information or type 'quit' to exit.{Colors.ENDC}")
                    continue
            
            if clarification.lower() in ['quit', 'exit', 'q']:
                print(f"\n{Colors.RED}❌ Exiting...{Colors.ENDC}")
                sys.exit(0)
            
            # Refine with LLM
            print(f"\n{Colors.CYAN}🤖 Refining with your input...{Colors.ENDC}")
            try:
                self.current_request = self.parser.refine_with_clarification(
                    current_state=self.current_request,
                    user_response=clarification,
                )
                print(f"{Colors.GREEN}✅ Refinement complete{Colors.ENDC}")
            except Exception as e:
                print(f"\n{Colors.RED}❌ Error refining request: {e}{Colors.ENDC}")
                continue
            
            # Display updated state
            self.display_parsed_state(self.current_request)
        
        # Check final state
        if not self.current_request.is_valid():
            print(f"\n{Colors.YELLOW}⚠️  Could not create a valid request after maximum iterations.{Colors.ENDC}")
            print(f"{Colors.CYAN}Please check your input and try again.{Colors.ENDC}")
            self.save_draft()
            sys.exit(1)
        
        # Final confirmation
        if self.confirm_submission(self.current_request):
            if not self.test_mode:
                self.submit_request()
        else:
            if not self.test_mode:
                print(f"\n{Colors.RED}❌ Submission cancelled.{Colors.ENDC}")
                self.save_draft()
    
    def submit_request(self):
        """Submit the request to the API."""
        print(f"\n{Colors.CYAN}{'=' * 70}{Colors.ENDC}")
        print(f"{Colors.BOLD}{Colors.HEADER}  📤 SUBMITTING REQUEST...{Colors.ENDC}")
        print(f"{Colors.CYAN}{'=' * 70}{Colors.ENDC}\n")
        
        api_request = self.current_request.to_api_request()
        
        # Try to actually call the API if available
        try:
            from chemtools.contracts import RecommendConditionsRequest
            from app.main import api_recommend_conditions
            
            # Create request object
            req = RecommendConditionsRequest(**api_request)
            
            print(f"{Colors.CYAN}⏳ Calling API endpoint: POST /api/v1/recommend/conditions{Colors.ENDC}\n")
            
            # Call the actual API
            result = api_recommend_conditions(req)
            
            print(f"{Colors.GREEN}{Colors.BOLD}✅ REQUEST SUCCESSFUL!{Colors.ENDC}\n")
            print(f"{Colors.CYAN}{'=' * 70}{Colors.ENDC}")
            print(f"{Colors.BOLD}RESULTS:{Colors.ENDC}")
            print(f"{Colors.CYAN}{'=' * 70}{Colors.ENDC}")
            print(json.dumps(result, indent=2, ensure_ascii=False))
            print(f"\n{Colors.GREEN}✨ Recommendation complete!{Colors.ENDC}")
            
        except ImportError:
            print(f"{Colors.YELLOW}⚠️  API not available in current context.{Colors.ENDC}")
            print(f"{Colors.CYAN}Request payload that would be sent:{Colors.ENDC}")
            print(f"{Colors.YELLOW}{json.dumps(api_request, indent=2)}{Colors.ENDC}")
            print(f"\n{Colors.CYAN}💡 To actually call the API, run this from the main app context.{Colors.ENDC}")
        except Exception as e:
            print(f"{Colors.RED}❌ API call failed: {e}{Colors.ENDC}")
            print(f"\n{Colors.CYAN}Request payload:{Colors.ENDC}")
            print(f"{Colors.YELLOW}{json.dumps(api_request, indent=2)}{Colors.ENDC}")
    
    def save_draft(self):
        """Save current state as draft."""
        draft_file = PROJECT_ROOT / "draft_request.json"
        with open(draft_file, "w") as f:
            json.dump({
                "reaction_smiles": self.current_request.reaction_smiles,
                "reaction_smiles_is_valid": self.current_request.reaction_smiles_is_valid,
                "reaction_type": self.current_request.reaction_type,
                "constraints": self.current_request.constraints,
                "additional_notes": self.current_request.additional_notes,
                "validation_issues": self.current_request.validation_issues,
                "clarification_needed": self.current_request.clarification_needed,
            }, f, indent=2)
        
        print(f"\n{Colors.CYAN}💾 Draft saved to: {Colors.ENDC}{draft_file}")


# ============================================================================
# Main Entry Point
# ============================================================================

def main():
    """Main CLI entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Interactive CLI for chemistry condition recommendations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use default settings (Aliyun/DeepSeek)
  python app/cli_recommend.py
  
  # Test mode (stop before submission)
  python app/cli_recommend.py --test
  
  # Use OpenAI
  python app/cli_recommend.py --provider openai --model gpt-4o
  
  # Enable debug logging
  python app/cli_recommend.py --debug
  
  # Disable colors
  python app/cli_recommend.py --no-color
        """
    )
    parser.add_argument(
        "--provider",
        default="aliyun",
        choices=["openai", "aliyun"],
        help="LLM provider (default: aliyun)",
    )
    parser.add_argument(
        "--model",
        default="deepseek-v3.2",
        help="LLM model name (default: deepseek-v3.2)",
    )
    parser.add_argument(
        "--api-key",
        help="API key (or set via environment variable)",
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Test mode: stop before actual API submission",
    )
    parser.add_argument(
        "--no-color",
        action="store_true",
        help="Disable colored output",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug logging",
    )
    
    args = parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Disable colors if requested
    if args.no_color:
        Colors.disable()
    
    # Initialize parser and CLI
    try:
        nl_parser = NaturalLanguageParser(
            provider=args.provider,
            model=args.model,
            api_key=args.api_key,
        )
    except Exception as e:
        print(f"{Colors.RED}❌ Failed to initialize LLM client: {e}{Colors.ENDC}")
        print(f"\n{Colors.CYAN}💡 Make sure you have set your API key via environment variable:{Colors.ENDC}")
        print(f"{Colors.YELLOW}   export OPENAI_API_KEY=your_key  (for OpenAI){Colors.ENDC}")
        print(f"{Colors.YELLOW}   export DASHSCOPE_API_KEY=your_key  (for Aliyun){Colors.ENDC}")
        sys.exit(1)
    
    cli = InteractiveCLI(parser=nl_parser, test_mode=args.test)
    
    try:
        cli.run()
    except KeyboardInterrupt:
        print(f"\n\n{Colors.YELLOW}⚠️  Interrupted by user. Exiting...{Colors.ENDC}")
        sys.exit(0)
    except Exception as e:
        logger.exception("Unexpected error in CLI")
        print(f"\n{Colors.RED}❌ Fatal error: {e}{Colors.ENDC}")
        sys.exit(1)


if __name__ == "__main__":
    main()

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
            "description": "Detected or specified reaction type (e.g., Suzuki, Buchwald, Ullmann, C_N_Coupling)"
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

Convert the user's natural language input into a structured JSON format for a recommendation system.

USER INPUT:
Reaction SMILES: {reaction_smiles}
Requirements: {requirements}

INSTRUCTIONS:
1. **VALIDATE REACTION SMILES**: Check if it follows reactants>>products format. Set "reaction_smiles_is_valid" to true/false.
   - Valid examples: "Br.c1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1", "CCBr.O=C(O)C>>CCC(=O)O"
   - Invalid examples: "invalid", "C=C", "broken>>", ">>product"
   - Must have ">>" separator and non-empty reactants
   
2. **EXTRACT CONSTRAINTS** (all optional):
   - Temperature: Extract numeric ranges (e.g., "no high temperature" → max: 80, "room temperature" → max: 30)
   - Base strength: Map to weak/moderate/strong (e.g., "no strong base" → "moderate", "mild base" → "weak")
   - Solvent: Map to polar_aprotic/polar_protic/nonpolar/aqueous
   - **Required reagents**: Extract when user says "use X", "prefer X", "with X" (e.g., "use copper catalyst" → ["copper"])
   - **Metal preference**: Extract metal symbol when user specifies catalyst metal (e.g., "copper catalyst" → "Cu", "palladium" → "Pd")
   - Excluded reagents: Extract specific chemical names when user says "avoid", "no", "exclude"
   - Air sensitivity: true if mentions "inert atmosphere", "glovebox", false if "air stable"
   - Cost: Map from "cheap"/"expensive" to low/medium/high

3. **VALIDATION**: Report issues in "validation_issues" array:
   - If reaction SMILES is invalid or missing: "Invalid reaction SMILES format"
   - Any other critical issues that would prevent API call

4. **CLARIFICATIONS**: List ambiguous requirements in "clarification_needed" (only for constraints):
   - Vague temperature terms without clear ranges
   - Ambiguous reagent exclusions
   - Contradictory requirements
   - Note: DO NOT ask for clarification if user provided no requirements - this is acceptable!

IMPORTANT:
- Only the reaction SMILES must be valid - all constraints are optional
- If user says "no requirements" or gives empty string, that's perfectly fine - just validate SMILES
- Empty constraints object is acceptable
- Be permissive with constraints - use null for unclear optional values

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
        model: str = "deepseek-v3.2-exp",
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
        reaction_smiles: str,
        requirements: str,
    ) -> ParsedRequest:
        """Parse initial user input into structured format."""
        logger.info("Parsing user input via LLM...")
        
        prompt = PARSE_USER_INPUT_PROMPT.format(
            reaction_smiles=reaction_smiles,
            requirements=requirements or "no specific requirements",
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
                       f"Clarifications: {len(result.clarification_needed)}")
            
            return result
            
        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse LLM response as JSON: {e}")
            logger.debug(f"LLM response: {content}")
            # Return a fallback with error
            return ParsedRequest(
                reaction_smiles=reaction_smiles,
                reaction_smiles_is_valid=False,
                validation_issues=[f"LLM parsing error: {str(e)}"]
            )
        except Exception as e:
            logger.error(f"Error during LLM parsing: {e}")
            return ParsedRequest(
                reaction_smiles=reaction_smiles,
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
    
    def __init__(self, parser: NaturalLanguageParser):
        self.parser = parser
        self.current_request: Optional[ParsedRequest] = None
    
    def print_header(self):
        """Print CLI header."""
        print("\n" + "=" * 70)
        print("  🧪 CHEMISTRY CONDITION RECOMMENDATION CLI")
        print("  Powered by LLM-assisted natural language parsing")
        print("=" * 70 + "\n")
    
    def print_separator(self):
        """Print section separator."""
        print("\n" + "-" * 70 + "\n")
    
    def get_initial_input(self) -> tuple[str, str]:
        """Get initial reaction SMILES and requirements from user."""
        print("Please provide your reaction details:\n")
        
        print("Reaction SMILES (reactants>>products format):")
        print("Example: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1")
        reaction = input("> ").strip()
        
        if not reaction:
            print("❌ Error: Reaction SMILES cannot be empty.")
            sys.exit(1)
        
        print("\n📝 Describe your requirements in natural language (optional):")
        print("Examples:")
        print("  • no strong base, room temperature")
        print("  • avoid expensive catalysts, prefer DMF solvent")
        print("  • mild conditions, no air-sensitive reagents")
        print("  • (press Enter to skip if no specific requirements)\n")
        
        requirements = input("> ").strip()
        
        return reaction, requirements
    
    def display_parsed_state(self, request: ParsedRequest):
        """Display current parsed state to user."""
        self.print_separator()
        print("📋 PARSED REQUEST:\n")
        
        # Reaction SMILES status
        if request.reaction_smiles_is_valid:
            print(f"✅ Reaction: {request.reaction_smiles}")
        else:
            print(f"❌ Reaction: {request.reaction_smiles}")
        
        if request.reaction_type:
            print(f"   Type: {request.reaction_type}")
        
        # Constraints (if any)
        if request.constraints:
            has_constraints = any(
                v is not None and v != {} and (not isinstance(v, list) or len(v) > 0)
                for v in request.constraints.values()
            )
            
            if has_constraints:
                print("\n📌 Constraints:")
                for key, value in request.constraints.items():
                    if value is None or value == {} or (isinstance(value, list) and len(value) == 0):
                        continue
                    
                    if isinstance(value, dict):
                        print(f"  • {key}:")
                        for k, v in value.items():
                            if v is not None:
                                print(f"      {k}: {v}")
                    elif isinstance(value, list):
                        print(f"  • {key}: {', '.join(str(v) for v in value)}")
                    else:
                        print(f"  • {key}: {value}")
        
        if request.additional_notes:
            print(f"\n💡 Notes: {request.additional_notes}")
        
        # Validation status
        if request.is_valid():
            print("\n✅ Request is VALID and ready to submit!")
        else:
            if request.validation_issues:
                print("\n⚠️  VALIDATION ISSUES:")
                for issue in request.validation_issues:
                    print(f"  • {issue}")
    
    def request_clarifications(self, request: ParsedRequest) -> Optional[str]:
        """Request clarifications or fixes from user."""
        self.print_separator()
        
        # Show validation issues first (blocking)
        if request.validation_issues:
            print("❌ VALIDATION ISSUES (must fix):\n")
            for i, issue in enumerate(request.validation_issues, 1):
                print(f"{i}. {issue}")
            print()
        
        # Show optional clarifications
        if request.clarification_needed:
            print("⚠️  CLARIFICATIONS NEEDED (optional):\n")
            for i, clarification in enumerate(request.clarification_needed, 1):
                print(f"{i}. {clarification}")
            print()
        
        # Prompt for response
        if request.validation_issues:
            print("Please fix the validation issues above:")
        elif request.clarification_needed:
            print("Please provide clarifications (or press Enter to skip):")
        
        response = input("> ").strip()
        
        return response if response else None
    
    def confirm_submission(self, request: ParsedRequest) -> bool:
        """Ask user to confirm submission."""
        self.print_separator()
        print("✅ Request is ready for submission!\n")
        
        # Show summary
        print(f"Reaction: {request.reaction_smiles}")
        if request.reaction_type:
            print(f"Type: {request.reaction_type}")
        
        # Count meaningful constraints
        constraint_count = sum(
            1 for v in request.constraints.values()
            if v is not None and v != {} and (not isinstance(v, list) or len(v) > 0)
        )
        
        if constraint_count > 0:
            print(f"Constraints: {constraint_count} specified")
        else:
            print("Constraints: None (will use default recommendations)")
        
        print("\n🔍 API Request Preview:")
        print(json.dumps(request.to_api_request(), indent=2))
        
        print("\n" + "=" * 70)
        confirm = input("Submit this request? (yes/no): ").strip().lower()
        return confirm in ["yes", "y"]
    
    def run(self):
        """Run the interactive CLI loop."""
        self.print_header()
        
        # Get initial input
        reaction, requirements = self.get_initial_input()
        
        # Parse with LLM
        print("\n🤖 Parsing input with LLM...")
        try:
            self.current_request = self.parser.parse_initial_input(
                reaction_smiles=reaction,
                requirements=requirements,
            )
        except Exception as e:
            print(f"\n❌ Error parsing input: {e}")
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
                    print("\n⚠️  Cannot proceed without fixing validation issues.")
                    print("Please provide the corrected information or type 'quit' to exit.")
                    continue
            
            if clarification.lower() in ['quit', 'exit', 'q']:
                print("\n❌ Exiting...")
                sys.exit(0)
            
            # Refine with LLM
            print("\n🤖 Refining with your input...")
            try:
                self.current_request = self.parser.refine_with_clarification(
                    current_state=self.current_request,
                    user_response=clarification,
                )
            except Exception as e:
                print(f"\n❌ Error refining request: {e}")
                continue
            
            # Display updated state
            self.display_parsed_state(self.current_request)
        
        # Check final state
        if not self.current_request.is_valid():
            print("\n⚠️  Could not create a valid request after maximum iterations.")
            print("Please check your input and try again.")
            self.save_draft()
            sys.exit(1)
        
        # Final confirmation
        if self.confirm_submission(self.current_request):
            self.submit_request()
        else:
            print("\n❌ Submission cancelled.")
            self.save_draft()
    
    def submit_request(self):
        """Submit the request to the API."""
        print("\n" + "=" * 70)
        print("  📤 SUBMITTING REQUEST...")
        print("=" * 70 + "\n")
        
        api_request = self.current_request.to_api_request()
        
        # Try to actually call the API if available
        try:
            from chemtools.contracts import RecommendConditionsRequest
            from app.main import api_recommend_conditions
            
            # Create request object
            req = RecommendConditionsRequest(**api_request)
            
            print("⏳ Calling API endpoint: POST /api/v1/recommend/conditions\n")
            
            # Call the actual API
            result = api_recommend_conditions(req)
            
            print("✅ REQUEST SUCCESSFUL!\n")
            print("=" * 70)
            print("RESULTS:")
            print("=" * 70)
            print(json.dumps(result, indent=2, ensure_ascii=False))
            print("\n✨ Recommendation complete!")
            
        except ImportError:
            print("⚠️  API not available in current context.")
            print("Request payload that would be sent:")
            print(json.dumps(api_request, indent=2))
            print("\n💡 To actually call the API, run this from the main app context.")
        except Exception as e:
            print(f"❌ API call failed: {e}")
            print("\nRequest payload:")
            print(json.dumps(api_request, indent=2))
    
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
        
        print(f"\n💾 Draft saved to: {draft_file}")


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
  
  # Use OpenAI
  python app/cli_recommend.py --provider openai --model gpt-4o
  
  # Enable debug logging
  python app/cli_recommend.py --debug
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
        default="deepseek-v3.2-exp",
        help="LLM model name (default: deepseek-v3.2-exp)",
    )
    parser.add_argument(
        "--api-key",
        help="API key (or set via environment variable)",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug logging",
    )
    
    args = parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Initialize parser and CLI
    try:
        nl_parser = NaturalLanguageParser(
            provider=args.provider,
            model=args.model,
            api_key=args.api_key,
        )
    except Exception as e:
        print(f"❌ Failed to initialize LLM client: {e}")
        print("\n💡 Make sure you have set your API key via environment variable:")
        print("   export OPENAI_API_KEY=your_key  (for OpenAI)")
        print("   export DASHSCOPE_API_KEY=your_key  (for Aliyun)")
        sys.exit(1)
    
    cli = InteractiveCLI(parser=nl_parser)
    
    try:
        cli.run()
    except KeyboardInterrupt:
        print("\n\n⚠️  Interrupted by user. Exiting...")
        sys.exit(0)
    except Exception as e:
        logger.exception("Unexpected error in CLI")
        print(f"\n❌ Fatal error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

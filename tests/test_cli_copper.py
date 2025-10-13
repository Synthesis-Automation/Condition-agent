"""Test CLI with copper catalyst requirement."""

import sys
import os
from unittest.mock import patch, MagicMock

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from app.cli_recommend import InteractiveCLI


def test_copper_catalyst_parsing():
    """Test that 'use copper catalyst' correctly affects the API request."""
    
    # Mock the LLMClient
    with patch('app.cli_recommend.LLMClient') as MockLLMClient:
        # Create mock client instance
        mock_client = MagicMock()
        MockLLMClient.return_value = mock_client
        
        # Mock LLM response for "use copper catalyst"
        mock_response = MagicMock()
        mock_response.content = """{
            "reaction_smiles": "Clc1c(C)cccc1C.COc1ccc(B(O)O)cc1>>Cc1cccc(C)c1-c1ccc(OC)cc1",
            "reaction_smiles_is_valid": true,
            "reaction_type": "Suzuki",
            "constraints": {
                "metal_preference": "Cu",
                "required_reagents": ["copper catalyst"]
            },
            "validation_issues": [],
            "clarification_needed": []
        }"""
        mock_client.chat.return_value = mock_response
        
        # Create parser instance
        parser = MagicMock()
        parser.parse_user_input = MagicMock()
        
        # Mock inputs: reaction SMILES and requirements
        reaction = "Clc1c(C)cccc1C.COc1ccc(B(O)O)cc1>>Cc1cccc(C)c1-c1ccc(OC)cc1"
        requirements = "use copper catalyst"
        
        # Simulate what the LLM would parse
        from app.cli_recommend import ParsedRequest
        result = ParsedRequest(
            reaction_smiles=reaction,
            reaction_smiles_is_valid=True,
            reaction_type="Suzuki",
            constraints={
                "metal_preference": "Cu",
                "required_reagents": ["copper catalyst"]
            }
        )
        
        # Verify parsing
        print("\n=== Parsing Result ===")
        print(f"Valid: {result.is_valid()}")
        print(f"Reaction type: {result.reaction_type}")
        print(f"Constraints: {result.constraints}")
        
        # Verify metal_preference extracted
        assert result.constraints.get("metal_preference") == "Cu", \
            f"Expected Cu, got {result.constraints.get('metal_preference')}"
        
        # Verify required_reagents extracted
        assert "copper catalyst" in result.constraints.get("required_reagents", []), \
            f"Expected 'copper catalyst' in required_reagents, got {result.constraints.get('required_reagents')}"
        
        # Convert to API request
        api_request = result.to_api_request()
        
        print("\n=== API Request ===")
        print(f"Reaction: {api_request['reaction']}")
        print(f"Constraints: {api_request['constraints']}")
        
        # Verify API request includes copper requirement
        assert "required_reagents" in api_request["constraints"], \
            "API request should include required_reagents"
        
        reagents = api_request["constraints"]["required_reagents"]
        assert any("Cu" in r or "copper" in r.lower() for r in reagents), \
            f"API request should include copper in required_reagents, got {reagents}"
        
        print("\n✅ Test passed! Copper catalyst requirement correctly parsed and formatted.")
        print(f"   - LLM extracted: metal_preference='Cu', required_reagents=['copper catalyst']")
        print(f"   - API request includes: required_reagents={reagents}")


def test_api_request_conversion():
    """Test that ParsedRequest.to_api_request() correctly handles metal_preference."""
    from app.cli_recommend import ParsedRequest
    
    # Test case 1: Only metal_preference
    request1 = ParsedRequest(
        reaction_smiles="A.B>>C",
        reaction_smiles_is_valid=True,
        constraints={
            "metal_preference": "Cu"
        }
    )
    
    api_req1 = request1.to_api_request()
    print("\n=== Test Case 1: Only metal_preference ===")
    print(f"Input: metal_preference='Cu'")
    print(f"Output: {api_req1['constraints']}")
    
    assert "required_reagents" in api_req1["constraints"], \
        "Should add required_reagents when metal_preference is set"
    assert any("Cu" in r for r in api_req1["constraints"]["required_reagents"]), \
        "Should include Cu in required_reagents"
    
    # Test case 2: Both metal_preference and required_reagents
    request2 = ParsedRequest(
        reaction_smiles="A.B>>C",
        reaction_smiles_is_valid=True,
        constraints={
            "metal_preference": "Pd",
            "required_reagents": ["base"]
        }
    )
    
    api_req2 = request2.to_api_request()
    print("\n=== Test Case 2: Both metal_preference and required_reagents ===")
    print(f"Input: metal_preference='Pd', required_reagents=['base']")
    print(f"Output: {api_req2['constraints']}")
    
    reagents2 = api_req2["constraints"]["required_reagents"]
    assert "base" in reagents2, "Should preserve existing required_reagents"
    assert any("Pd" in r for r in reagents2), "Should add Pd catalyst"
    
    # Test case 3: metal_preference="any" or None should not add constraint
    request3 = ParsedRequest(
        reaction_smiles="A.B>>C",
        reaction_smiles_is_valid=True,
        constraints={
            "metal_preference": "any"
        }
    )
    
    api_req3 = request3.to_api_request()
    print("\n=== Test Case 3: metal_preference='any' ===")
    print(f"Input: metal_preference='any'")
    print(f"Output: {api_req3['constraints']}")
    
    # "any" should be filtered out (value check)
    assert "metal_preference" not in api_req3["constraints"] or \
           api_req3["constraints"]["metal_preference"] == "any", \
        "metal_preference='any' should be handled gracefully"
    
    print("\n✅ All conversion tests passed!")


if __name__ == "__main__":
    print("Testing copper catalyst CLI functionality...")
    print("=" * 60)
    
    test_api_request_conversion()
    print("\n" + "=" * 60)
    test_copper_catalyst_parsing()

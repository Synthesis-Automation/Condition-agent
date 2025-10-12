#!/bin/bash
# Quick API Test Examples using curl
# Make sure the server is running: python -m uvicorn app.main:app --reload --port 8000

BASE_URL="http://localhost:8000"

echo "==================================================================="
echo "Testing Condition-Agent API Endpoints"
echo "==================================================================="
echo ""

# =============================================================================
# HEALTH CHECK
# =============================================================================
echo "1. Health Check"
echo "-------------------------------------------------------------------"
curl -X GET "$BASE_URL/health" | jq
echo ""

# =============================================================================
# RULE-BASED RECOMMENDATION
# =============================================================================
echo "2. Rule-Based Recommendation (/match)"
echo "-------------------------------------------------------------------"
curl -X POST "$BASE_URL/match" \
  -H "Content-Type: application/json" \
  -d '{
    "reaction": "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1",
    "db": "data/conditionDB/suzuki_db.json",
    "include_trace": true
  }' | jq
echo ""

# =============================================================================
# ML-BASED RECOMMENDATION
# =============================================================================
echo "3. ML-Based Recommendation (/api/v1/recommend/conditions)"
echo "-------------------------------------------------------------------"
curl -X POST "$BASE_URL/api/v1/recommend/conditions" \
  -H "Content-Type: application/json" \
  -d '{
    "reaction": "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1",
    "reaction_type": "Suzuki",
    "k": 50,
    "limit": 3
  }' | jq
echo ""

# =============================================================================
# FUSION RECOMMENDATION
# =============================================================================
echo "4. Fusion Recommendation (/api/v1/recommend/fusion)"
echo "-------------------------------------------------------------------"
curl -X POST "$BASE_URL/api/v1/recommend/fusion" \
  -H "Content-Type: application/json" \
  -d '{
    "reaction": "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1",
    "k": 50,
    "max_variants": 3
  }' | jq
echo ""

# =============================================================================
# PROTOCOL (Currently uses local module - no API endpoint yet)
# =============================================================================
echo "5. Protocol Recommendation (via CLI)"
echo "-------------------------------------------------------------------"
echo "Protocol endpoint not yet implemented."
echo "Use the CLI tool instead:"
echo ""
echo "  python scripts/web_recommendation_cli.py \\"
echo "    --rxn 'BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1' \\"
echo "    --strategy protocol \\"
echo "    --k 3"
echo ""

echo "==================================================================="
echo "All tests completed!"
echo "==================================================================="

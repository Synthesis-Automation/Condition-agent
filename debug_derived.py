"""Debug derived feature evaluation for reactant types."""

from chemtools.featurizers.calculable import detect_all_features, _evaluate_derived_feature
import json

# Test with boronic acid
smiles = 'c1ccccc1B(O)O'
features = detect_all_features(smiles)

# Show all reactant features
reactant_feats = {k: v for k, v in features.items() if '_reactant' in k and v}
print('Detected reactant features:')
for k in sorted(reactant_feats.keys()):
    print(f'  {k}: {reactant_feats[k]}')

# Check specific tokens
token1 = "ArB(OH)2_reactant"
token2 = "ArB_reactant"
print(f'\nDirect lookup {token1}: {features.get(token1)}')
print(f'Direct lookup {token2}: {features.get(token2)}')

# Check if tokens with parentheses are being evaluated in derived expressions
print('\n--- Testing derived expression evaluation ---')
test_expr = "ArB(OH)2_reactant OR ArB(OR)2_reactant OR ArBF3K_reactant"
print(f'Expression: {test_expr}')
print(f'Component values:')
print(f'  ArB(OH)2_reactant = {features.get("ArB(OH)2_reactant", False)}')
print(f'  ArB(OR)2_reactant = {features.get("ArB(OR)2_reactant", False)}')
print(f'  ArBF3K_reactant = {features.get("ArBF3K_reactant", False)}')

# Manually test the evaluation function
result = _evaluate_derived_feature(test_expr, features)
print(f'\nManual evaluation result: {result}')
print(f'Expected: True')
print(f'Match: {result == True}')

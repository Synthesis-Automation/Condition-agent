import traceback

print("Test 1: Import featurize_reaction directly")
try:
    from chemtools.featurizers.formatters.reaction import featurize_reaction
    print("✓ Direct import works")
except Exception as e:
    print(f"✗ Direct import failed: {e}")
    traceback.print_exc()

print("\nTest 2: Import detection module")
try:
    from chemtools import detection
    print("✓ Detection module imported")
except Exception as e:
    print(f"✗ Detection import failed: {e}")
    traceback.print_exc()

print("\nTest 3: Call extract_reaction_key")
try:
    from chemtools.detection import extract_reaction_key
    result = extract_reaction_key('CN(C)C(=S)S.[Na].Clc1ccc(I)cc1>>CN(C)C(=S)Sc1ccc(Cl)cc1')
    print(f"✓ extract_reaction_key works: {result[3]}")
except Exception as e:
    print(f"✗ extract_reaction_key failed: {e}")
    traceback.print_exc()

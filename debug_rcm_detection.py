"""Debug script to see what family is detected for RCM."""

from chem_assistant.chemtools_wrapper import detect_reaction

rcm_reaction = "C=CCNC(=O)C=C>>C1=CCNC(=O)C1"

print("Testing RCM reaction detection:")
print(f"Reaction: {rcm_reaction}\n")

result = detect_reaction(rcm_reaction, use_ml=False)

print("Detection result:")
print(f"  Family: {result.get('family')}")
print(f"  Confidence: {result.get('confidence')}")
print(f"  Family Name: {result.get('family_name')}")
print(f"\nFull result:")
for key, value in result.items():
    print(f"  {key}: {value}")

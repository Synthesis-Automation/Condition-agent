entry = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1 (Suzuki - Simple Ph-Ph)"
parts = entry.split("(")
print(f"Parts: {len(parts)}")
for i, p in enumerate(parts):
    print(f"Part {i}: '{p}'")

# The issue - boronic acid has (O) in it!
print(f"\n>>> The problem: splitting on '(' catches the B(O)O boronic acid!")
print(f">>> This creates multiple parts, not just 2!")

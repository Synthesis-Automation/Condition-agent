"""
Benchmark DRFP fingerprint generation speed.
"""

import time
from chemtools import reaction_similarity as rs

# Sample reactions from the dataset
test_reactions = [
    "C1CC2(CCN1)OCCO2.Clc1ccc2sccc2c1>>O=C1CCN(c2ccc3sccc3c2)CC1",
    "CC1(C)c2ccccc2Nc2ccccc21.CCCCCCCCc1ccc(Br)cc1>>CCCCCCCCc1ccc(N2c3ccccc3C(C)(C)c3ccccc32)cc1",
    "Brc1ccc(Br)cc1.c1ccc(N2CCCCC2)cc1>>c1ccc(N2CCCCC2)cc1.c1ccc(N2CCCCC2)cc1",
    "COc1ccc(Br)cc1.Nc1ccccc1>>COc1ccc(Nc2ccccc2)cc1",
    "Clc1ccccc1.c1ccc(N2CCNCC2)cc1>>c1ccc(N2CCN(c3ccccc3)CC2)cc1",
]

print("DRFP Benchmark")
print("=" * 60)
print(f"Testing {len(test_reactions)} reactions...")
print()

if not rs.drfp_available():
    print("ERROR: DRFP not available!")
    print("Install with: pip install drfp")
    exit(1)

# Warm up (first call may be slower due to imports)
print("Warming up...")
_ = rs.encode_drfp(test_reactions[0])
print()

# Benchmark each reaction
times = []
for i, rxn in enumerate(test_reactions, 1):
    start = time.perf_counter()
    fp = rs.encode_drfp(rxn)
    elapsed = time.perf_counter() - start
    times.append(elapsed)
    
    status = "✓" if fp is not None else "✗"
    print(f"{status} Reaction {i}: {elapsed*1000:.1f} ms")

print()
print("=" * 60)
print(f"Average time per reaction: {sum(times)/len(times)*1000:.1f} ms")
print(f"Total time for {len(test_reactions)} reactions: {sum(times):.2f} s")
print()
print(f"Estimated time for 1000 reactions: {sum(times)/len(times)*1000:.0f} seconds ({sum(times)/len(times)*1000/60:.1f} minutes)")
print(f"Estimated time for 1344 reactions: {sum(times)/len(times)*1344:.0f} seconds ({sum(times)/len(times)*1344/60:.1f} minutes)")

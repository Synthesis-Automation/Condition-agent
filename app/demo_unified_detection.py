#!/usr/bin/env python3
"""
Quick demo of the updated reaction_analysis_interactive.py app.

Shows the key features of the unified detection system.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from chemtools import detect_reaction


def demo_detection() -> None:
    """Demonstrate the unified detection API with detailed output."""

    print("=" * 80)
    print("UNIFIED REACTION DETECTION API - DEMO")
    print("=" * 80)
    print()

    # Test reaction
    suzuki = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"

    print(f"Test Reaction: {suzuki}")
    print()

    # Demo 1: ML-enhanced detection
    print("-" * 80)
    print("DEMO 1: ML-Enhanced Detection (use_ml=True)")
    print("-" * 80)
    print()

    result_ml = detect_reaction(suzuki, use_ml=True)

    print(f"Family:      {result_ml['family']}")
    print(f"Confidence:  {result_ml['confidence']:.2f}")
    print(f"Method:      {result_ml['method']}")
    print(f"Agreement:   {result_ml.get('agreement', 'N/A')}")
    print(f"Status:      {result_ml.get('status', 'N/A')}")
    print()

    print("Details:")
    details = result_ml.get('details', {})

    # Rule prediction
    rule_pred = details.get('rule_prediction', {})
    print(f"  Rule-based:  {rule_pred.get('family')} (conf: {rule_pred.get('confidence', 0):.2f})")

    # ML prediction
    ml_pred = details.get('ml_prediction', {})
    if ml_pred and ml_pred.get('family'):
        conf = ml_pred.get('confidence')
        conf_str = f"{conf:.2f}" if conf is not None else "N/A"
        print(f"  ML-based:    {ml_pred.get('family')} (conf: {conf_str})")
        if ml_pred.get('rxn_name'):
            print(f"               Name: \"{ml_pred['rxn_name']}\"")
        if ml_pred.get('rxn_class'):
            print(f"               Class: \"{ml_pred['rxn_class']}\"")
    else:
        print("  ML-based:    Not available")

    # Functional groups
    fg = details.get('functional_groups', {})
    detected = [k for k, v in fg.items() if v]
    print(f"  Func Groups: {', '.join(detected[:5])}")

    # Catalysts
    catalysts = details.get('catalysts', [])
    if catalysts:
        print(f"  Catalysts:   {', '.join(catalysts)}")

    print()

    # Demo 2: Rule-based only
    print("-" * 80)
    print("DEMO 2: Rule-Based Only Detection (use_ml=False)")
    print("-" * 80)
    print()

    result_rule = detect_reaction(suzuki, use_ml=False)

    print(f"Family:      {result_rule['family']}")
    print(f"Confidence:  {result_rule['confidence']:.2f}")
    print(f"Method:      {result_rule['method']}")
    print()

    # Demo 3: C-N coupling with catalyst detection
    print("-" * 80)
    print("DEMO 3: Catalyst Detection (Pd vs Cu)")
    print("-" * 80)
    print()

    cn_pd = "Brc1ccccc1.Nc1ccccc1>Pd(PPh3)4>c1ccccc1Nc1ccccc1"
    cn_cu = "Brc1ccccc1.Nc1ccccc1>CuI>c1ccccc1Nc1ccccc1"

    print("With Pd catalyst:")
    result_pd = detect_reaction(cn_pd, use_ml=False)
    print(f"  → Family: {result_pd['family']}")
    print(f"  → Catalysts: {result_pd['details'].get('catalysts', [])}")
    print()

    print("With Cu catalyst:")
    result_cu = detect_reaction(cn_cu, use_ml=False)
    print(f"  → Family: {result_cu['family']}")
    print(f"  → Catalysts: {result_cu['details'].get('catalysts', [])}")
    print()

    # Demo 4: Unknown reaction
    print("-" * 80)
    print("DEMO 4: Unknown Reaction Handling")
    print("-" * 80)
    print()

    unknown = "CCCC.OOOO>>CCCCOOO"
    result_unknown = detect_reaction(unknown, use_ml=False)

    print(f"Test: {unknown}")
    print(f"Family:      {result_unknown['family']}")
    print(f"Confidence:  {result_unknown['confidence']:.2f}")

    if result_unknown['family'] == 'Unknown':
        print()
        print("⚠ Suggestions:")
        print("  • Verify reaction SMILES is valid")
        print("  • Check if reaction is in supported families")
        print("  • Try ML detection (use_ml=True)")

    print()
    print("=" * 80)
    print("Demo complete! Run reaction_analysis_interactive.py for full experience.")
    print("=" * 80)


if __name__ == "__main__":
    try:
        demo_detection()
    except Exception as exc:  # pragma: no cover - demo harness
        print(f"Error: {exc}")
        import traceback
        traceback.print_exc()

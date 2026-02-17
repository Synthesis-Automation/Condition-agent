"""
Compare multiple models for chemistry analysis:
- DeepSeek-v3.2 (Aliyun)
- Kimi-k2.5 (Aliyun)
- GLM-4.7 (Aliyun)
- GPT-4o (OpenAI)
"""

import sys
import time
sys.path.insert(0, 'c:/Git-softwares/Condition-agent')

from chemtools.quick_reaction_glance import quick_reaction_glance
from llmtools.clients import LLMClient

# Your test reaction (Suzuki + THP deprotection)
rxn_smiles = "CC1(C)OB(c2cnn(CCOC3CCCCO3)c2)OC1(C)C.Cc1nc(-c2cn3c(n2)-c2ccc(Br)cc2OCC3)n(C)n1>>Cc1nc(-c2cn3c(n2)-c2ccc(-c4cnn(CCO)c4)cc2OCC3)n(C)n1"

print("="*80)
print("Multi-Model Comparison for Chemistry Analysis")
print("="*80)
print(f"\nTest Reaction: Suzuki coupling + THP deprotection")
print(f"SMILES: {rxn_smiles[:80]}...")
print()
print("Expected:")
print("  1. Primary: Suzuki-Miyaura coupling (Br → aryl)")
print("  2. Secondary: THP deprotection (workup)")
print()

models = [
    ("aliyun", "deepseek-v3.2"),
    ("aliyun", "kimi-k2.5"),
    ("aliyun", "glm-4.7"),
    ("openai", "gpt-4o"),
]

results = {}

for provider, model_name in models:
    print(f"\n{'='*80}")
    print(f"Testing: {model_name} ({provider})")
    print('='*80)

    try:
        client = LLMClient(provider=provider, model=model_name)
        print(f"✓ Client created: {client.model}")

        start = time.time()
        result = quick_reaction_glance(
            rxn_smiles,
            client,
            prompt_style="comprehensive",
            thorough=True
        )
        elapsed = time.time() - start

        if result.get('success'):
            print(f"✓ Analysis successful! ({elapsed:.1f}s)")
            print()

            # Summary
            summary = result.get('summary', 'N/A')
            print(f"Summary:")
            print(f"  {summary[:300]}{'...' if len(summary) > 300 else ''}")
            print()

            # Reaction types
            rxn_types = result.get('reaction_types', [])
            print(f"Reaction Types: {', '.join(rxn_types)}")

            # Check for Suzuki detection
            suzuki_detected = any('suzuki' in rt.lower() or 'coupling' in rt.lower() for rt in rxn_types)
            if suzuki_detected:
                print("  ✓ Suzuki coupling detected!")
            else:
                print("  ✗ Suzuki coupling MISSED")
            print()

            # Protecting groups (KEY TEST)
            pg = result.get('protecting_groups', {})
            pg_removed = pg.get('removed', [])
            pg_added = pg.get('added', [])

            if pg_removed:
                print(f"✓ Protecting Groups Removed ({len(pg_removed)}):")
                for item in pg_removed:
                    print(f"  • {item[:150]}{'...' if len(item) > 150 else ''}")
            else:
                print(f"✗ No protecting groups detected")

            if pg_added:
                print(f"✓ Protecting Groups Added ({len(pg_added)}):")
                for item in pg_added:
                    print(f"  • {item}")
            print()

            # All changes
            changes = result.get('all_changes', [])
            print(f"Structural Changes: {len(changes)} detected")
            if changes:
                for i, change in enumerate(changes[:4], 1):
                    print(f"  {i}. {change[:120]}{'...' if len(change) > 120 else ''}")
                if len(changes) > 4:
                    print(f"  ... and {len(changes)-4} more")
            print()

            # Metrics
            meta = result.get('metadata', {})
            print(f"Complexity: {result.get('complexity', 'N/A')}")
            print(f"Confidence: {result.get('confidence', 0):.2f}")
            print(f"Latency: {meta.get('latency_ms', elapsed*1000):.0f}ms")
            print(f"Tokens: {meta.get('total_tokens', 'N/A')}")

            # Side reactions
            side = result.get('side_reactions', [])
            if side:
                print(f"\nSide/Workup Reactions ({len(side)}):")
                for s in side[:2]:
                    print(f"  • {s[:100]}{'...' if len(s) > 100 else ''}")

            # Pharma context
            pharma = result.get('pharmaceutical_context')
            if pharma:
                print(f"\nPharma Context: {pharma[:120]}{'...' if len(pharma) > 120 else ''}")

            results[model_name] = {
                'success': True,
                'elapsed': elapsed,
                'suzuki_detected': suzuki_detected,
                'pg_detected': len(pg_removed) + len(pg_added),
                'changes_detected': len(changes),
                'side_reactions': len(side),
                'confidence': result.get('confidence', 0),
                'tokens': meta.get('total_tokens', 0),
                'summary': summary
            }

        else:
            print(f"✗ Analysis failed: {result.get('error', 'Unknown error')}")
            if 'raw_response' in result:
                raw = result['raw_response']
                print(f"Raw response preview ({len(raw)} chars): {raw[:300]}")
            results[model_name] = {
                'success': False,
                'error': result.get('error'),
                'raw_length': len(result.get('raw_response', ''))
            }

    except Exception as e:
        print(f"✗ Error: {e}")
        results[model_name] = {'success': False, 'error': str(e)}
        import traceback
        traceback.print_exc()

# Summary comparison
print(f"\n\n{'='*80}")
print("COMPARISON SUMMARY")
print('='*80)

print(f"\n{'Model':<20} {'Success':<10} {'Suzuki':<10} {'PG':<8} {'Changes':<10} {'Conf.':<8} {'Time':<10}")
print('-'*80)

for provider, model_name in models:
    r = results.get(model_name, {})
    if r.get('success'):
        suzuki_mark = '✓' if r.get('suzuki_detected') else '✗'
        print(f"{model_name:<20} {'✓':<10} {suzuki_mark:<10} {r.get('pg_detected', 0):<8} {r.get('changes_detected', 0):<10} {r.get('confidence', 0):<8.2f} {r.get('elapsed', 0):<10.1f}s")
    else:
        error_msg = str(r.get('error', 'Failed'))[:30]
        print(f"{model_name:<20} {'✗':<10} {error_msg}")

print()
print("Legend:")
print("  Suzuki = Detected Suzuki-Miyaura coupling (primary reaction)")
print("  PG = Protecting Groups detected (removed + added)")
print("  Changes = Number of structural changes identified")
print("  Conf. = Model confidence (0.0-1.0)")

# Detailed scoring
print(f"\n{'='*80}")
print("DETAILED SCORING")
print('='*80)

successful_models = [(name, r) for (_, name), r in zip(models, [results.get(m) for _, m in models]) if r and r.get('success')]

if successful_models:
    print("\nScoring Criteria:")
    print("  - Suzuki detection: 10 points")
    print("  - THP deprotection: 10 points")
    print("  - Structural changes: 1 point each (max 10)")
    print("  - Side reactions: 2 points each (max 6)")
    print("  - Speed: (fastest gets 5, scaled for others)")
    print()

    scores = {}

    # Find fastest for speed scoring
    fastest_time = min(r['elapsed'] for _, r in successful_models) if successful_models else 1

    for model_name, r in successful_models:
        score = 0
        breakdown = []

        # Suzuki detection (10 pts)
        if r.get('suzuki_detected'):
            score += 10
            breakdown.append("Suzuki: 10")
        else:
            breakdown.append("Suzuki: 0")

        # THP detection (10 pts)
        if r.get('pg_detected', 0) > 0:
            score += 10
            breakdown.append("THP: 10")
        else:
            breakdown.append("THP: 0")

        # Changes (max 10 pts, 1 each)
        changes_score = min(r.get('changes_detected', 0), 10)
        score += changes_score
        breakdown.append(f"Changes: {changes_score}")

        # Side reactions (max 6 pts, 2 each)
        side_score = min(r.get('side_reactions', 0) * 2, 6)
        score += side_score
        breakdown.append(f"Side: {side_score}")

        # Speed (max 5 pts, relative to fastest)
        speed_score = int((fastest_time / r['elapsed']) * 5)
        score += speed_score
        breakdown.append(f"Speed: {speed_score}")

        scores[model_name] = {'score': score, 'breakdown': breakdown}

    # Sort by score
    ranked = sorted(scores.items (), key=lambda x: x[1]['score'], reverse=True)

    print(f"{'Rank':<6} {'Model':<20} {'Score':<8} {'Breakdown'}")
    print('-'*80)
    for i, (model_name, data) in enumerate(ranked, 1):
        medal = "🥇" if i == 1 else "🥈" if i == 2 else "🥉" if i == 3 else "  "
        print(f"{medal} {i:<4} {model_name:<20} {data['score']:<8} {' | '.join(data['breakdown'])}")

    print()
    winner_name, winner_data = ranked[0]
    print(f"🏆 WINNER: {winner_name} with {winner_data['score']} points!")

    # Winner summary
    winner_result = results[winner_name]
    print(f"\nWhy {winner_name} won:")
    if winner_result.get('suzuki_detected'):
        print("  ✓ Correctly identified Suzuki coupling")
    if winner_result.get('pg_detected', 0) > 0:
        print("  ✓ Detected THP deprotection")
    print(f"  ✓ Found {winner_result.get('changes_detected', 0)} structural changes")
    print(f"  ✓ Analysis time: {winner_result.get('elapsed', 0):.1f}s")

else:
    print("\n⚠️  No models succeeded - see errors above")

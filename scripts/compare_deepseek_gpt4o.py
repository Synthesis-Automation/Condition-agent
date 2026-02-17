"""
Compare deepseek-v3.2 (aliyun) vs gpt-4o (openai) for chemistry analysis.
"""

import sys
import time
sys.path.insert(0, 'c:/Git-softwares/Condition-agent')

from chemtools.quick_reaction_glance import quick_reaction_glance
from llmtools.clients import LLMClient

# Your test reaction (Suzuki + THP deprotection)
rxn_smiles = "CC1(C)OB(c2cnn(CCOC3CCCCO3)c2)OC1(C)C.Cc1nc(-c2cn3c(n2)-c2ccc(Br)cc2OCC3)n(C)n1>>Cc1nc(-c2cn3c(n2)-c2ccc(-c4cnn(CCO)c4)cc2OCC3)n(C)n1"

print("="*80)
print("Comparing DeepSeek-v3.2 (Aliyun) vs GPT-4o (OpenAI)")
print("="*80)
print(f"\nTest Reaction: Suzuki coupling + THP deprotection")
print(f"SMILES: {rxn_smiles[:80]}...")
print()

models = [
    ("aliyun", "deepseek-v3.2"),
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
            print(f"Summary:")
            print(f"  {result.get('summary', 'N/A')[:200]}...")
            print()

            # Reaction types
            print(f"Reaction Types: {', '.join(result.get('reaction_types', []))}")
            print()

            # Protecting groups (KEY TEST)
            pg = result.get('protecting_groups', {})
            if pg.get('removed'):
                print(f"✓ Protecting Groups Removed:")
                for item in pg['removed']:
                    print(f"  • {item}")
            else:
                print(f"✗ No protecting groups detected")

            if pg.get('added'):
                print(f"✓ Protecting Groups Added:")
                for item in pg['added']:
                    print(f"  • {item}")
            print()

            # All changes
            changes = result.get('all_changes', [])
            if changes:
                print(f"All Changes ({len(changes)} detected):")
                for change in changes[:3]:
                    print(f"  • {change}")
                if len(changes) > 3:
                    print(f"  ... and {len(changes)-3} more")
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
                print(f"Side Reactions: {', '.join(side[:2])}")

            # Pharma context
            pharma = result.get('pharmaceutical_context')
            if pharma:
                print(f"Pharma Context: {pharma[:100]}...")

            results[model_name] = {
                'success': True,
                'elapsed': elapsed,
                'pg_detected': len(pg.get('removed', [])) + len(pg.get('added', [])),
                'changes_detected': len(changes),
                'confidence': result.get('confidence', 0),
                'tokens': meta.get('total_tokens', 0)
            }

        else:
            print(f"✗ Analysis failed: {result.get('error', 'Unknown error')}")
            if 'raw_response' in result:
                print(f"Raw response preview: {result['raw_response'][:200]}")
            results[model_name] = {'success': False, 'error': result.get('error')}

    except Exception as e:
        print(f"✗ Error: {e}")
        results[model_name] = {'success': False, 'error': str(e)}
        import traceback
        traceback.print_exc()

# Summary comparison
print(f"\n\n{'='*80}")
print("COMPARISON SUMMARY")
print('='*80)

print(f"\n{'Model':<20} {'Success':<10} {'PG Det.':<10} {'Changes':<10} {'Conf.':<10} {'Time':<10}")
print('-'*80)

for provider, model_name in models:
    r = results.get(model_name, {})
    if r.get('success'):
        print(f"{model_name:<20} {'✓':<10} {r.get('pg_detected', 0):<10} {r.get('changes_detected', 0):<10} {r.get('confidence', 0):<10.2f} {r.get('elapsed', 0):<10.1f}s")
    else:
        print(f"{model_name:<20} {'✗':<10} {'-':<10} {'-':<10} {'-':<10} {'-':<10}")

print()
print("Key Metrics:")
print("  PG Det. = Protecting Groups Detected (removed + added)")
print("  Changes = Number of structural changes identified")
print("  Conf. = Model confidence (0.0-1.0)")
print("  Time = Analysis time in seconds")

# Winner determination
print()
if all(r.get('success') for r in results.values()):
    deepseek_pg = results['deepseek-v3.2'].get('pg_detected', 0)
    gpt4o_pg = results['gpt-4o'].get('pg_detected', 0)

    if deepseek_pg > 0 and gpt4o_pg > 0:
        print("✓ Both models successfully detected THP deprotection!")

        deepseek_time = results['deepseek-v3.2'].get('elapsed', 999)
        gpt4o_time = results['gpt-4o'].get('elapsed', 999)

        if deepseek_time < gpt4o_time:
            speedup = (gpt4o_time / deepseek_time - 1) * 100
            print(f"🏆 DeepSeek-v3.2 is FASTER by {speedup:.0f}% ({deepseek_time:.1f}s vs {gpt4o_time:.1f}s)")
        else:
            slowdown = (deepseek_time / gpt4o_time - 1) * 100
            print(f"🏆 GPT-4o is FASTER by {slowdown:.0f}% ({gpt4o_time:.1f}s vs {deepseek_time:.1f}s)")
    elif deepseek_pg > 0:
        print("🏆 DeepSeek-v3.2 WINS - detected protecting groups, GPT-4o missed them")
    elif gpt4o_pg > 0:
        print("🏆 GPT-4o WINS - detected protecting groups, DeepSeek missed them")
    else:
        print("⚠️  Neither model detected protecting groups")
else:
    print("⚠️  One or both models failed - see errors above")

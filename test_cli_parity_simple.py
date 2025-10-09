"""
Simple signature comparison for CLI feature parity.
"""

print("Testing CLI Feature Parity")
print("=" * 80)

# Read and check local_recommendation_cli.py
local_cli_path = "scripts/local_recommendation_cli.py"
web_cli_path = "scripts/web_recommendation_cli.py"

print("\nChecking local_recommendation_cli.py...")
with open(local_cli_path, 'r', encoding='utf-8') as f:
    local_content = f.read()
    
local_has_rerank = 'rerank_strategy: str' in local_content and 'def local_ml_recommendation' in local_content
local_has_filter = 'filter_unknown_reagents: bool' in local_content and 'def local_ml_recommendation' in local_content
local_has_args = 'argparse.ArgumentParser' in local_content
local_has_rerank_arg = '--rerank' in local_content
local_has_filter_arg = '--filter-unknown' in local_content

print("\nChecking web_recommendation_cli.py...")
with open(web_cli_path, 'r', encoding='utf-8') as f:
    web_content = f.read()
    
web_has_rerank = 'rerank_strategy: str' in web_content and 'def call_ml_recommendation' in web_content
web_has_filter = 'filter_unknown_reagents: bool' in web_content and 'def call_ml_recommendation' in web_content
web_has_args = 'argparse.ArgumentParser' in web_content
web_has_rerank_arg = '--rerank' in web_content
web_has_filter_arg = '--filter-unknown' in web_content

print("\n" + "=" * 80)
print("Feature Comparison:")
print("-" * 80)

features = [
    ("rerank_strategy parameter in ML function", local_has_rerank, web_has_rerank),
    ("filter_unknown_reagents parameter in ML function", local_has_filter, web_has_filter),
    ("Command-line argument parsing", local_has_args, web_has_args),
    ("--rerank CLI argument", local_has_rerank_arg, web_has_rerank_arg),
    ("--filter-unknown CLI argument", local_has_filter_arg, web_has_filter_arg),
]

all_match = True
for feature_name, local_val, web_val in features:
    match = "✅" if local_val == web_val else "⚠️"
    status = "MATCH" if local_val == web_val else "MISMATCH"
    print(f"\n{match} {feature_name}: {status}")
    print(f"   Local: {'YES' if local_val else 'NO'}")
    print(f"   Web:   {'YES' if web_val else 'NO'}")
    if local_val != web_val:
        all_match = False

print("\n" + "=" * 80)
if all_match:
    print("✅ FEATURE PARITY ACHIEVED!")
    print("Both CLIs now support the same features.")
else:
    print("⚠️  FEATURE MISMATCH DETECTED")
    print("Some features differ between CLIs.")
print("=" * 80)

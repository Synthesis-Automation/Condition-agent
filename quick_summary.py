import json

with open('test_all_sample_reactions_results.json', 'r') as f:
    data = json.load(f)

print("=" * 60)
print("QUICK SUMMARY")
print("=" * 60)
print(f"Total reactions: {data['overall_stats']['total']}")
print(f"Success rate: {data['overall_stats']['success']/data['overall_stats']['total']*100:.1f}%")
print(f"Total recommendations: {data['overall_stats']['total_recommendations']}")
print(f"Avg recommendations/reaction: {data['overall_stats']['total_recommendations']/data['overall_stats']['total']:.1f}")
print()

print("Top 5 reaction types by total recommendations:")
types = sorted(data['stats_by_type'].items(), 
               key=lambda x: x[1]['total_recommendations'], 
               reverse=True)[:5]
for i, (rtype, stats) in enumerate(types, 1):
    print(f"{i}. {rtype}: {stats['total_recommendations']} total ({stats['total']} reactions, avg {stats['total_recommendations']/stats['total']:.1f})")

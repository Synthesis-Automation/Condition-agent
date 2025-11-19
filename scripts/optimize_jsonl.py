#!/usr/bin/env python3
"""
Analyze JSONL file size and optionally recompress with sparse boolean optimization.
Usage:
    python scripts/optimize_jsonl.py data/reaction_dataset/Suzuki.jsonl --analyze
    python scripts/optimize_jsonl.py data/reaction_dataset/Suzuki.jsonl --optimize --output Suzuki_compact.jsonl
"""
import json
import argparse
from pathlib import Path
from typing import Any, Dict


def analyze_jsonl(filepath: str) -> Dict[str, Any]:
    """Analyze JSONL file structure and size."""
    path = Path(filepath)
    if not path.exists():
        print(f"Error: File not found: {filepath}")
        return {}
    
    file_size_mb = path.stat().st_size / (1024 * 1024)
    
    # Sample first record
    with open(path, 'r', encoding='utf-8') as f:
        first_line = f.readline()
        if not first_line:
            print("Error: Empty file")
            return {}
        
        sample = json.loads(first_line)
    
    # Count records
    with open(path, 'r', encoding='utf-8') as f:
        num_records = sum(1 for _ in f)
    
    # Analyze features
    features = sample.get("precomputed", {}).get("features", {})
    num_features = len(features)
    bool_features = sum(1 for v in features.values() if isinstance(v, bool))
    true_features = sum(1 for v in features.values() if v is True)
    false_features = sum(1 for v in features.values() if v is False)
    
    # Estimate potential savings
    avg_record_size = file_size_mb * 1024 / num_records  # KB per record
    
    # Rough estimate: each "false" boolean takes ~20 bytes in JSON
    potential_savings_kb = false_features * 20 / 1024
    potential_savings_pct = (potential_savings_kb / avg_record_size) * 100
    
    return {
        "filepath": filepath,
        "file_size_mb": file_size_mb,
        "num_records": num_records,
        "avg_record_kb": avg_record_size,
        "sample_features": num_features,
        "bool_features": bool_features,
        "true_features": true_features,
        "false_features": false_features,
        "potential_savings_pct": potential_savings_pct,
        "estimated_optimized_size_mb": file_size_mb * (1 - potential_savings_pct / 100),
    }


def optimize_jsonl(input_path: str, output_path: str) -> None:
    """Rewrite JSONL with sparse boolean optimization."""
    input_file = Path(input_path)
    output_file = Path(output_path)
    
    if not input_file.exists():
        print(f"Error: Input file not found: {input_path}")
        return
    
    print(f"Optimizing {input_path} -> {output_path}")
    
    original_size = input_file.stat().st_size
    records_processed = 0
    
    with open(input_file, 'r', encoding='utf-8') as fin, \
         open(output_file, 'w', encoding='utf-8') as fout:
        
        for line in fin:
            if not line.strip():
                continue
            
            record = json.loads(line)
            
            # Apply sparse boolean optimization to features
            if "precomputed" in record and "features" in record["precomputed"]:
                features = record["precomputed"]["features"]
                compact_features = {}
                
                for k, v in features.items():
                    if isinstance(v, bool):
                        if v:  # Only store True values
                            compact_features[k] = True
                    else:
                        # Keep non-boolean values as-is
                        compact_features[k] = v
                
                record["precomputed"]["features"] = compact_features
            
            # Write optimized record
            fout.write(json.dumps(record, ensure_ascii=False) + '\n')
            records_processed += 1
            
            if records_processed % 1000 == 0:
                print(f"  Processed {records_processed:,} records...")
    
    optimized_size = output_file.stat().st_size
    savings_pct = (1 - optimized_size / original_size) * 100
    
    print(f"\n✓ Optimization complete!")
    print(f"  Records processed: {records_processed:,}")
    print(f"  Original size:     {original_size / (1024**2):.2f} MB")
    print(f"  Optimized size:    {optimized_size / (1024**2):.2f} MB")
    print(f"  Space saved:       {(original_size - optimized_size) / (1024**2):.2f} MB ({savings_pct:.1f}%)")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze and optimize JSONL files with sparse boolean features"
    )
    parser.add_argument("input", help="Input JSONL file path")
    parser.add_argument("--analyze", action="store_true", help="Analyze file structure and potential savings")
    parser.add_argument("--optimize", action="store_true", help="Optimize and write to output file")
    parser.add_argument("--output", help="Output file path (required with --optimize)")
    
    args = parser.parse_args()
    
    if args.analyze:
        print(f"Analyzing {args.input}...")
        stats = analyze_jsonl(args.input)
        
        if stats:
            print("\n" + "="*60)
            print("JSONL FILE ANALYSIS")
            print("="*60)
            print(f"File:                    {stats['filepath']}")
            print(f"Size:                    {stats['file_size_mb']:.2f} MB")
            print(f"Records:                 {stats['num_records']:,}")
            print(f"Avg record size:         {stats['avg_record_kb']:.2f} KB")
            print(f"\nFeature Analysis (sample):")
            print(f"  Total features:        {stats['sample_features']}")
            print(f"  Boolean features:      {stats['bool_features']}")
            print(f"    True values:         {stats['true_features']}")
            print(f"    False values:        {stats['false_features']}")
            print(f"\nPotential Savings:")
            print(f"  Sparse optimization:   ~{stats['potential_savings_pct']:.1f}%")
            print(f"  Estimated new size:    ~{stats['estimated_optimized_size_mb']:.2f} MB")
            print("="*60)
    
    if args.optimize:
        if not args.output:
            print("Error: --output required with --optimize")
            return
        
        optimize_jsonl(args.input, args.output)


if __name__ == "__main__":
    main()

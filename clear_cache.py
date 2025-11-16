# Clear all Python bytecode caches
import os
from pathlib import Path

print("Clearing Python bytecode caches...\n")

# Find and remove all __pycache__ directories
cache_dirs = list(Path(".").rglob("__pycache__"))
print(f"Found {len(cache_dirs)} cache directories\n")

for cache_dir in cache_dirs:
    try:
        # Remove all .pyc files
        pyc_files = list(cache_dir.glob("*.pyc"))
        for pyc_file in pyc_files:
            pyc_file.unlink()
        
        # Try to remove the directory
        if cache_dir.exists():
            cache_dir.rmdir()
        print(f"✅ Cleared: {cache_dir}")
    except Exception as e:
        print(f"⚠️  Could not clear {cache_dir}: {e}")

print(f"\n✅ Cache clearing complete!")
print("\nNow restart your GUI app to pick up the fixed classification patterns.")

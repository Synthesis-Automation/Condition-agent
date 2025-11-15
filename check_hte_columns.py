#!/usr/bin/env python3
"""Check HTE database columns"""
from chemtools.HTE import HTERecommender
import pandas as pd

recommender = HTERecommender()

print("Database columns:")
for i, col in enumerate(recommender.df.columns, 1):
    print(f"{i:3d}. {col}")

print(f"\nTotal rows: {len(recommender.df)}")
print(f"\nSample data:")
print(recommender.df.head(2))

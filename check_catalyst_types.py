"""Check catalyst types in HTE database"""
import pandas as pd

df = pd.read_csv('data/HTE_db/HTE_0.csv')

print("="*80)
print("CATALYST ANALYSIS IN HTE DATABASE")
print("="*80)

cats = df['Catalyst'].value_counts()
print(f"\nTotal unique catalysts: {df['Catalyst'].nunique()}")
print(f"Total experiments: {len(df)}")

print("\n" + "="*80)
print("TOP 30 CATALYSTS")
print("="*80)
print(cats.head(30))

# Copper catalysts
cu_mask = df['Catalyst'].str.contains('Cu', case=False, na=False)
cu_cats = df[cu_mask]['Catalyst'].value_counts()
cu_exps = cu_mask.sum()

print("\n" + "="*80)
print(f"COPPER-BASED CATALYSTS ({len(cu_cats)} types, {cu_exps} experiments)")
print("="*80)
print(cu_cats)

# Palladium catalysts
pd_mask = df['Catalyst'].str.contains('Pd', case=False, na=False)
pd_cats = df[pd_mask]['Catalyst'].value_counts()
pd_exps = pd_mask.sum()

print("\n" + "="*80)
print(f"PALLADIUM-BASED CATALYSTS ({len(pd_cats)} types, {pd_exps} experiments)")
print("="*80)
print(pd_cats.head(30))

# Other metals
ni_mask = df['Catalyst'].str.contains('Ni', case=False, na=False)
ni_exps = ni_mask.sum()
print(f"\nNickel-based: {ni_exps} experiments")

ir_mask = df['Catalyst'].str.contains('Ir', case=False, na=False)
ir_exps = ir_mask.sum()
print(f"Iridium-based: {ir_exps} experiments")

rh_mask = df['Catalyst'].str.contains('Rh', case=False, na=False)
rh_exps = rh_mask.sum()
print(f"Rhodium-based: {rh_exps} experiments")

print("\n" + "="*80)
print("SUMMARY BY METAL TYPE")
print("="*80)
print(f"Palladium: {pd_exps:,} experiments ({pd_exps/len(df)*100:.1f}%)")
print(f"Copper:    {cu_exps:,} experiments ({cu_exps/len(df)*100:.1f}%)")
print(f"Nickel:    {ni_exps:,} experiments ({ni_exps/len(df)*100:.1f}%)")
print(f"Iridium:   {ir_exps:,} experiments ({ir_exps/len(df)*100:.1f}%)")
print(f"Rhodium:   {rh_exps:,} experiments ({rh_exps/len(df)*100:.1f}%)")
print(f"No ligand: {(df['Catalyst'].isna().sum() + (df['Catalyst']=='').sum()):,} experiments")

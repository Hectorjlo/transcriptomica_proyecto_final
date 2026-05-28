"""
Hardcoded python script to extract SRRs given some filters
"""
import pandas as pd

# Metadata file
META = "GSE132040_MACA_Bulk_metadata.csv"

# Filters
TISSUE = "Bone"
AGES = {3, 24}
SEX = "m"  # males

# Create a pd.DataFrame
df = pd.read_csv(META)

# Define the column names exactly as they appear in the metadata table.
# These names include spaces and punctuation, so we keep them centralized
# in variables to avoid typos and make future edits easier.
col_source = "source name"
col_age = "characteristics: age"
col_sex = "characteristics: sex"
col_srr = "raw file"

# Quick schema validation:
# check that all required columns are present before any filtering logic runs.
# If a required column is missing, stop early with a clear error message.
missing = [c for c in [col_source, col_age, col_sex, col_srr] if c not in df.columns]
if missing:
    raise SystemExit(f"Missing columns in the CSV: {missing}\nAvailable columns: {df.columns.tolist()}")

# Cleanup and type normalization:
# 1) standardize sex values to lowercase trimmed strings
# 2) coerce age to numeric so set membership checks behave consistently
df[col_sex] = df[col_sex].astype(str).str.strip().str.lower()
df[col_age] = pd.to_numeric(df[col_age], errors="coerce")

# Filtering criteria:
# 1) source name starts with "Bone_" to keep only the target tissue entries
# 2) age is one of the selected groups in AGES ({3, 24})
# 3) sex matches the selected value in SEX ("m")
mask = (
    df[col_source].astype(str).str.startswith(f"{TISSUE}_") &
    df[col_age].isin(AGES) &
    (df[col_sex] == SEX)
)

sub = df.loc[mask, [col_source, col_age, col_sex, col_srr]].copy()

# Sort the filtered subset for cleaner and reproducible output ordering.
sub = sub.sort_values([col_age, col_source, col_srr])

# Extract the SRR identifiers as clean strings, one identifier per row.
srrs = sub[col_srr].astype(str).str.strip()

# Print a concise report:
# 1) full filtered table for verification
# 2) SRR-only list for downstream command-line workflows
print("=== Summary (Hector: Bone, 3 vs 24, males) ===")
print(sub.to_string(index=False))
print("\n=== SRR (one per line) ===")
for s in srrs:
    print(s)

# Save SRR list to a plain text file for convenience and reuse in pipelines.
out_txt = "hector_bone_3_vs_24_srr.txt"
srrs.to_csv(out_txt, index=False, header=False)
print(f"\nSaved: {out_txt}  (SRRs list)")
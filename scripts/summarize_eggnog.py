#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import argparse

# CLI setup
parser = argparse.ArgumentParser(description="Summarize eggNOG annotations")
parser.add_argument("-i", "--input", required=True, help="Input eggnog_annotation.tsv file")
parser.add_argument("-o", "--output_dir", default=".", help="Output directory")

args = parser.parse_args()
input_file = Path(args.input)
output_dir = Path(args.output_dir)
output_dir.mkdir(parents=True, exist_ok=True)

# Read the annotation TSV (skip metadata)
ann = pd.read_csv(input_file, sep='\t', skiprows=4)
ann.columns = ann.columns.str.strip().str.replace(' ', '_').str.replace('#', '').str.replace('.', '_')

# KEGG summary
if 'KEGG_ko' in ann.columns:
    top_kegg = ann['KEGG_ko'].dropna().str.split(',').explode().value_counts().head(20)
    top_kegg.to_csv(output_dir / "eggnog_top_kegg.tsv", sep='\t')
else:
    print("⚠️ 'KEGG_ko' column not found.")
    pd.Series().to_csv(output_dir / "eggnog_top_kegg.tsv", sep='\t')

# GO summary
if 'GOs' in ann.columns:
    go_terms = ann['GOs'].dropna().str.split(',').explode()
    go_counts = go_terms.value_counts().head(20)
    go_counts.to_csv(output_dir / "eggnog_top_go.tsv", sep='\t')
else:
    print("⚠️ 'GOs' column not found.")
    pd.Series().to_csv(output_dir / "eggnog_top_go.tsv", sep='\t')

# COG summary
if 'COG_category' in ann.columns:
    top_cog = ann['COG_category'].dropna().str.cat(sep='')
    top_cog = pd.Series(list(top_cog)).value_counts().head(20)
    top_cog.to_csv(output_dir / "eggnog_top_cog.tsv", sep='\t')
else:
    print("⚠️ 'COG_category' column not found.")
    pd.Series().to_csv(output_dir / "eggnog_top_cog.tsv", sep='\t')

# Trait keyword tagging
trait_keywords = ['lipid', 'fatty', 'acetyl', 'nitrogen', 'ammonia', 'nitrate', 'light', 'photosystem']
for col in ['Description', 'KEGG_Pathway', 'GOs']:
    if col not in ann.columns:
        ann[col] = ''

ann['all_annotations'] = ann[['Description', 'KEGG_Pathway', 'GOs']].astype(str).agg(' '.join, axis=1).str.lower()
ann['Trait'] = ann['all_annotations'].apply(lambda x: next((k for k in trait_keywords if k in x), 'Other'))
trait_summary = ann['Trait'].value_counts()
trait_summary.to_csv(output_dir / "eggnog_summary_by_trait.tsv", sep='\t')

# Excel summary
with pd.ExcelWriter(output_dir / "eggnog_summary.xlsx") as writer:
    ann.to_excel(writer, sheet_name='All_Annotations', index=False)
    if not top_kegg.empty:
        top_kegg.to_excel(writer, sheet_name='Top_KEGG')
    if not go_counts.empty:
        go_counts.to_excel(writer, sheet_name='Top_GO')
    if not top_cog.empty:
        top_cog.to_excel(writer, sheet_name='Top_COG')
    trait_summary.to_excel(writer, sheet_name='Trait_Summary')

#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
import subprocess

# Input files
fasta = Path("augustus.hints.aa")
blastp_tsv = Path("blastp.tsv")
interproscan_tsv = Path("interproscan.tsv")
eggnog_annotations = Path("eggnog_annotation.tsv")

# Load BRAKER protein list
ids = subprocess.check_output(f"grep '^>' {fasta}", shell=True).decode().splitlines()
all_ids = [i[1:].split()[0] for i in ids]
df_base = pd.DataFrame({'protein_id': all_ids})

# Load BLAST
blast_cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
              'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovs']
blast = pd.read_csv(blastp_tsv, sep='\t', names=blast_cols)
blast_best = blast.sort_values(['qseqid', 'bitscore'], ascending=[True, False]).drop_duplicates('qseqid')
blast_best = blast_best.rename(columns={'qseqid': 'protein_id'})

# Load InterPro
ipr = pd.read_csv(interproscan_tsv, sep='\t', header=None, comment='#')
ipr.columns = [
    "protein_id", "md5", "length", "analysis", "signature_accession", "signature_description",
    "start", "end", "score", "status", "date", "interpro_accession", "interpro_description",
    "GO_terms", "Pathways"
]
ipr_domains = ipr.groupby('protein_id')['interpro_description'].apply(lambda x: ';'.join(x.dropna().unique()))
ipr_gos = ipr.groupby('protein_id')['GO_terms'].apply(lambda x: ';'.join(x.dropna().unique()))

# Load eggNOG
egg = pd.read_csv(eggnog_annotations, sep='\t', skiprows=4)
egg.columns = egg.columns.str.strip().str.replace(' ', '_').str.replace('#', '')
egg = egg.rename(columns={
    'query': 'protein_id',
    'Preferred_name': 'eggNOG_ortholog',
    'COG_category': 'COG',
    'KEGG_ko': 'KEGG',
    'Description': 'eggNOG_description'
})

# Trait extraction
trait_keywords = ['lipid', 'fatty', 'acetyl', 'wax', 'desaturase', 'thioesterase',
                  'nitrogen', 'nitrate', 'ammonia', 'photosystem', 'light']

for col in ['eggNOG_description', 'KEGG', 'GOs']:
    if col not in egg.columns:
        egg[col] = ''

egg['combined'] = egg[['eggNOG_description', 'KEGG', 'GOs']].astype(str).agg(' '.join, axis=1).str.lower()
egg['trait'] = egg['combined'].apply(lambda x: next((k for k in trait_keywords if k in x), 'Other/Unclassified'))

# Merge all
merged = df_base     .merge(blast_best[['protein_id', 'sseqid', 'pident', 'qcovs']], how='left', on='protein_id')     .merge(ipr_domains.rename("interpro_domains"), how='left', on='protein_id')     .merge(ipr_gos.rename("GO_terms"), how='left', on='protein_id')     .merge(egg[['protein_id', 'eggNOG_ortholog', 'COG', 'KEGG', 'eggNOG_description', 'trait']], how='left', on='protein_id')

# Output
merged.to_csv("merged_annotation_table.tsv", sep='\t', index=False)
with pd.ExcelWriter("annotation_merge.xlsx") as writer:
    merged.to_excel(writer, sheet_name="Combined_Annotations", index=False)
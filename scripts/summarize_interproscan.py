#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import subprocess
import re

# -------- PARSE ARGUMENTS --------
parser = argparse.ArgumentParser(description="InterProScan annotation summarization and visualization")
parser.add_argument("-i", "--input", required=True, help="Input InterProScan TSV file")
parser.add_argument("-o", "--output", required=True, help="Output directory for reports and plots")
parser.add_argument("-f", "--fasta", required=True, help="FASTA file of BRAKER protein sequences (use the CLEANED one!)")
parser.add_argument("--go_obo", default=None, help="Path to go-basic.obo (optional). If omitted/missing, GO names won't be resolved.")
args = parser.parse_args()

tsv_file = args.input
output_dir = Path(args.output)
braker_fasta = args.fasta
go_obo = args.go_obo
max_items = 20

# Define trait-linked keywords (customize as needed)
trait_keywords = [
    "lipid", "acyl", "DGAT", "desaturase", "thioesterase", "fatty acid",
    "glycerolipid", "beta-oxidation", "wax", "phospholipid"
]

# Define domain categories using keyword mapping
domain_categories = {
    "Regulation/Signaling": ["kinase", "phosphatase", "WD40", "ankyrin", "bZIP"],
    "Metabolism": ["hydrolase", "oxidoreductase", "transferase", "dehydrogenase", "synthetase"],
    "Transport": ["transporter", "carrier", "channel", "MFS", "ABC"],
    "Structural/Binding": ["repeat", "motif", "coil", "helix", "beta", "armadillo"],
    "Lipid/High-Value": trait_keywords,
}

# -------- SETUP --------
output_dir.mkdir(exist_ok=True)
print(f"[INFO] Output folder: {output_dir}")

# -------- LOAD DATA --------
print("[INFO] Loading InterProScan TSV...")
df = pd.read_csv(tsv_file, sep='\t', header=None, comment='#', low_memory=False)
df.columns = [
    "protein_id", "md5", "length", "analysis", "signature_accession", "signature_description",
    "start", "end", "score", "status", "date", "interpro_accession", "interpro_description", "GO_terms", "Pathways"
]
df['domain_label'] = df.apply(
    lambda row: row['interpro_description'] if pd.notnull(row['interpro_description']) and row['interpro_description'] != "-" else row['signature_description'],
    axis=1
)

# -------- CATEGORY MAPPING --------
def assign_category(description):
    desc = (description or "").lower()
    for category, keywords in domain_categories.items():
        if any(k.lower() in desc for k in keywords):
            return category
    return "Other/Unclassified"

df['domain_category'] = df['domain_label'].apply(assign_category)

# -------- DOMAIN CATEGORY BARPLOT --------
domain_cat_counts = df['domain_category'].value_counts()
domain_cat_counts.to_csv(output_dir / "domain_category_counts.csv")

plt.figure(figsize=(8, 5))
domain_cat_counts.plot(kind='bar')
plt.title("Protein Domains by Functional Category")
plt.ylabel("Count")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(output_dir / "domain_category_barplot.png")

# -------- GO TERM PARSING --------
print("[INFO] Parsing and categorizing GO terms...")
go_df = df[['protein_id', 'GO_terms']].dropna().copy()
go_df = go_df.assign(GO=go_df['GO_terms'].str.split('|')).explode('GO')
go_df['GO'] = go_df['GO'].apply(lambda x: re.sub(r'\(.*?\)', '', x).strip())

def go_category(go_id):
    # Very rough heuristic; InterProScan doesn't label categories in TSV by default
    if go_id.startswith("GO:"):
        # try to infer by common numeric ranges; otherwise Unknown
        # This is simplistic; proper mapping would use GO DAG.
        return "Unknown"
    return "Unknown"

go_df['GO_category'] = go_df['GO'].apply(go_category)

# -------- Optional: resolve GO names using OBO --------
top_go_df = pd.DataFrame(columns=["GO_ID","GO_Term","Count","Category"])
if go_obo and Path(go_obo).exists():
    print("[INFO] Loading GO term ontology...")
    from goatools.obo_parser import GODag
    go_dag = GODag(go_obo)

    print("[INFO] Generating grouped GO term barplot with GO term names...")
    top_go_combined = []
    # If you want category-specific top terms, you need a real category mapper; here we do overall top
    top_terms = go_df['GO'].value_counts().head(30)
    for go_term, count in top_terms.items():
        term_name = go_dag[go_term].name if go_term in go_dag else go_term
        top_go_combined.append({"GO_ID": go_term, "GO_Term": term_name, "Count": count, "Category": "All"})
    top_go_df = pd.DataFrame(top_go_combined)
    top_go_df.to_csv(output_dir / "top_go_terms_named.csv", index=False)

    plt.figure(figsize=(10, 10))
    sns.set(style="whitegrid")
    top_go_df_sorted = top_go_df.sort_values("Count", ascending=True)
    top_go_df_sorted["GO_Term"] = pd.Categorical(top_go_df_sorted["GO_Term"], categories=top_go_df_sorted["GO_Term"], ordered=True)
    sns.barplot(x="Count", y="GO_Term", data=top_go_df_sorted, dodge=False)
    plt.xlabel("Count"); plt.ylabel("GO Term Description"); plt.title("Top GO Terms")
    plt.subplots_adjust(left=0.35)
    plt.savefig(output_dir / "go_grouped_top_terms_named_barplot.png", bbox_inches='tight')

# -------- HIGH-VALUE TRAIT FILTERING --------
print("[INFO] Highlighting trait-associated domains and proteins...")
trait_mask = df['domain_label'].str.lower().str.contains('|'.join(map(re.escape, trait_keywords)), na=False)
df[trait_mask].to_csv(output_dir / "high_value_domains.tsv", sep='\t', index=False)

trait_go = go_df[go_df['protein_id'].isin(df[trait_mask]['protein_id'])]
trait_go.to_csv(output_dir / "high_value_protein2go.tsv", sep='\t', index=False)

# -------- TOP 20 INTERPRO DOMAINS PER CATEGORY --------
print("[INFO] Generating barplots for top 20 InterPro domains per category...")
top_domain_output_dir = output_dir / "top20_domain_by_category"
top_domain_output_dir.mkdir(exist_ok=True)

for category in df['domain_category'].unique():
    sub_df = df[df['domain_category'] == category]
    top_domains = sub_df['domain_label'].value_counts().head(max_items)

    if not top_domains.empty:
        top_domains.to_csv(top_domain_output_dir / f"top20_domains_{category.replace('/', '_')}.csv")
        plt.figure(figsize=(10, 6))
        top_domains[::-1].plot(kind='barh')
        plt.title(f"Top {max_items} InterPro Domains: {category}")
        plt.xlabel("Count")
        plt.tight_layout()
        plt.savefig(top_domain_output_dir / f"top20_domains_{category.replace('/', '_')}_barplot.png")

# -------- GO CATEGORY BARPLOT (summary counts only, by string heuristic) --------
go_category_counts = go_df['GO_category'].value_counts()
plt.figure(figsize=(6, 4))
go_category_counts.plot(kind='bar')
plt.title("GO Term Counts by Category")
plt.ylabel("GO Term Count")
plt.xticks(rotation=30, ha='right')
plt.tight_layout()
plt.savefig(output_dir / "go_category_summary_barplot.png")

print("[INFO] Exporting all annotation summaries to Excel...")
excel_path = output_dir / "interproscan_annotation_summary.xlsx"
with pd.ExcelWriter(excel_path, engine='xlsxwriter') as writer:
    df.to_excel(writer, sheet_name='All_Annotations', index=False)
    go_df.to_excel(writer, sheet_name='GO_Terms', index=False)
    if not top_go_df.empty:
        top_go_df.to_excel(writer, sheet_name='Top_GO_Terms_Named', index=False)
    if trait_mask.any():
        df[trait_mask].to_excel(writer, sheet_name='High_Value_Domains', index=False)
        trait_go.to_excel(writer, sheet_name='High_Value_GO', index=False)
    combined_top_domains = []
    for category in df['domain_category'].unique():
        sub_df = df[df['domain_category'] == category]
        top_domains = sub_df['domain_label'].value_counts().head(max_items).reset_index()
        top_domains.columns = ['domain_label', 'count']
        top_domains['category'] = category
        combined_top_domains.append(top_domains)
    if combined_top_domains:
        pd.concat(combined_top_domains).to_excel(writer, sheet_name='Top_InterPro_by_Category', index=False)

print(f"[INFO] Excel summary written to: {excel_path}")

# -------- COMBINED INTERPROSCAN SUMMARY REPORT --------
print("[INFO] Generating combined InterProScan summary report...")

# Count total BRAKER proteins from the provided (cleaned) FASTA
total_bra_proteins = int(subprocess.check_output(["grep", "-c", "^>", braker_fasta]).decode().strip())

interpro_proteins = df['protein_id'].nunique()
annotation_percent = (interpro_proteins / total_bra_proteins) * 100 if total_bra_proteins > 0 else 0.0

overall_summary_lines = [
    "Overall InterProScan Summary",
    "-----------------------------",
    "Description                         | Count",
    "------------------------------------|--------",
    f"Total BRAKER proteins               | {total_bra_proteins:,}",
    f"Proteins with InterProScan hits     | {interpro_proteins:,}",
    f"% Annotated                         | {annotation_percent:.2f}%",
    ""
]

# Functional vs Hypothetical (very rough heuristic)
print("[INFO] Calculating functional vs hypothetical annotations...")
all_annotated_proteins = df['protein_id'].unique()
functional_mask = df['domain_label'].str.lower().str.contains(r'[a-z]', na=False) & \
    ~df['domain_label'].str.lower().str.contains(r'hypothetical|uncharacterized|unknown|domain of unknown function|candidate', na=False)
functional_proteins = df[functional_mask]['protein_id'].unique()
nonfunctional_proteins = set(all_annotated_proteins) - set(functional_proteins)

num_total = len(all_annotated_proteins)
num_func = len(functional_proteins)
num_nonfunc = len(nonfunctional_proteins)

func_percent = (num_func / num_total * 100) if num_total > 0 else 0.0
nonfunc_percent = (num_nonfunc / num_total * 100) if num_total > 0 else 0.0

func_summary_lines = [
    "Functional vs Hypothetical Annotations",
    "----------------------------------------",
    "Annotation Type               | Protein Count | Percent of Annotated",
    "------------------------------|----------------|----------------------",
    f"Functional Annotation         | {num_func:>14,} | {func_percent:>20.2f}%",
    f"Hypothetical/Uncharacterized  | {num_nonfunc:>14,} | {nonfunc_percent:>20.2f}%",
    ""
]

# Per-database InterProScan hits
db_hit_counts = df.groupby("analysis")["protein_id"].nunique().reset_index()
db_hit_counts.columns = ["Database", "Proteins_with_Hits"]
db_hit_counts["Percent_of_Total"] = db_hit_counts["Proteins_with_Hits"] / total_bra_proteins * 100

db_summary_lines = [
    "InterProScan Database Hits",
    "---------------------------",
    "Database            | Proteins with Hits | % of Total BRAKER Proteins",
    "---------------------|--------------------|----------------------------"
]
for _, row in db_hit_counts.iterrows():
    db_summary_lines.append(f"{row['Database']:<20} | {int(row['Proteins_with_Hits']):>18,} | {row['Percent_of_Total']:>26.2f}%")
total_db_hits = db_hit_counts["Proteins_with_Hits"].sum()
total_percent = db_hit_counts["Percent_of_Total"].sum()
db_summary_lines.append(f"{'Total (non-unique)':<20} | {total_db_hits:>18,} | {total_percent:>26.2f}%")

combined_summary_txt = output_dir / "interproscan_summary_combined.txt"
with open(combined_summary_txt, "w") as f:
    f.write("InterProScan Annotation Summary\n\n")
    f.write("\n".join(overall_summary_lines))
    f.write("\n\n")
    f.write("\n".join(func_summary_lines))
    f.write("\n\n")
    f.write("\n".join(db_summary_lines))
    f.write("\n")

print(f"[INFO] Combined summary TXT written to: {combined_summary_txt}")
print("[INFO] Script completed.")
print(f"[INFO] Files saved in: {output_dir}")

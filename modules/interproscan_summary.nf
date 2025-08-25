process InterproscanSummary {
  tag "summarize interproscan"
  cpus 1
  container 'python:3.11-slim'
  publishDir { "${params.outdir}/annotation/interproscan" }, mode: 'copy'

  input:
    path ipr_tsv
    path proteins_clean
    path script_file
    val  go_obo_path

  output:
    path "interproscan_summary.tsv",                emit: summary_table
    path "interproscan_annotation_summary.xlsx",    emit: summary_xlsx
    path "interproscan_summary_combined.txt",       emit: summary_txt
    path "domain_category_barplot.png", optional: true
    path "go_category_summary_barplot.png", optional: true

  script:
  """
  set -euo pipefail

  # deps
  python3 -m pip install --no-cache-dir \
    pandas==2.2.2 matplotlib==3.9.0 seaborn==0.13.2 goatools==1.4.12 xlsxwriter==3.2.0

  # run your summary script
  python3 ${script_file} \
    -i ${ipr_tsv} \
    -f ${proteins_clean} \
    -o . \
    ${ go_obo_path ? "--go_obo ${go_obo_path}" : "" }

  # create a tiny TSV summary from the combined TXT using Python (no awk)
  python3 - <<'PY'
  import re, pathlib
  p = pathlib.Path("interproscan_summary_combined.txt")
  out = pathlib.Path("interproscan_summary.tsv")
  headers = ["Metric","Value"]
  rows = []
  if p.exists():
      txt = p.read_text()
      for line in txt.splitlines():
          if "Total BRAKER proteins" in line or "Proteins with InterProScan hits" in line or "% Annotated" in line:
              parts = line.split("|", 1)
              if len(parts) == 2:
                  k = parts[0].strip()
                  v = parts[1].strip().rstrip("%")
                  rows.append((k, v))
  with out.open("w") as f:
      f.write("\\t".join(headers)+"\\n")
      for k,v in rows:
          f.write(f"{k}\\t{v}\\n")
  PY
    """
}

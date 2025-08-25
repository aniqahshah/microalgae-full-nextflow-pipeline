process EggnogSummary {
  tag "summarize eggnog annotations"
  publishDir { "${params.outdir}/annotation/eggnog" }, mode: 'copy'

  input:
    path eggnog_annotations
    path summary_script  // <--- add this

  output:
    path "eggnog_top_kegg.tsv"
    path "eggnog_top_go.tsv"
    path "eggnog_top_cog.tsv"
    path "eggnog_summary_by_trait.tsv"
    path "eggnog_summary.xlsx"

  script:
  """
  pip install --no-cache-dir pandas xlsxwriter

  python3 ${summary_script} \\
      --input ${eggnog_annotations} \\
      --output_dir .
  """
}

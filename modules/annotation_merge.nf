process AnnotationMerge {
  tag "merge annotations"
  container = 'python:3.10-slim'
  publishDir { "${params.outdir}/annotation/summary_all_annotation" }, mode: 'copy'

  input:
  path interproscan_tsv
  path blastp_tsv
  path eggnog_annotations
  path braker_proteins_fasta
  path "annotation_merge.py"

  output:
  path "merged_annotation_table.tsv"
  path "annotation_merge.xlsx"

  script:
  """
  pip install --no-cache-dir pandas xlsxwriter

  python3 annotation_merge.py \\
    --interproscan ${interproscan_tsv} \\
    --blastp ${blastp_tsv} \\
    --eggnog ${eggnog_annotations} \\
    --proteins ${braker_proteins_fasta}
  """
}

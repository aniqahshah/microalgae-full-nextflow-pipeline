process EggnogMapper {
  tag "eggnog annotation"
  publishDir { "${params.outdir}/annotation/eggnog" }, mode: 'copy'

  input:
    path protein_file

  output:
    path "eggnog_annotation.tsv", emit: annotations

  script:
  """
  emapper.py \\
    -i ${protein_file} \\
    -o eggnog_annotation \\
    --output_dir . \\
    --cpu 8 \\
    --itype proteins \\
    --data_dir /eggnog_data \\
    --override

  mv eggnog_annotation.emapper.annotations eggnog_annotation.tsv
  """
}

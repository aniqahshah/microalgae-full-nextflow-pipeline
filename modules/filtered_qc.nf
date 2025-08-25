process FilteredQC {
  tag "filtered_qc"
  publishDir { "${params.outdir}/qc/filtered" }, mode: 'copy'

  input:
  path filtered_fastq

  output:
  path "nanoplot_filtered", emit: filtered_qc_dir

  script:
  """
  mkdir -p nanoplot_filtered
  NanoPlot --fastq ${filtered_fastq} -o nanoplot_filtered --threads ${params.threads}
  """
}
process RawQC {
  tag "raw_qc"
    publishDir { "${params.outdir}/qc/raw" }, mode: 'copy'

  input:
  path input_fastq

  output:
  path "nanoplot_raw", emit: raw_qc_dir

  script:
  """
  mkdir -p nanoplot_raw
  NanoPlot --fastq ${input_fastq} -o nanoplot_raw --threads ${params.threads}
  """
}
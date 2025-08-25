process FilterReads {
  tag "filter_reads"
  publishDir { "${params.outdir}/qc/filtered" }, mode: 'copy'

  input:
  path input_fastq

  output:
  path "filtered_reads.fastq.gz", emit: filtered_reads

  script:
  """
  gunzip -c ${input_fastq} | NanoFilt -q 10 --length 1000 | gzip > filtered_reads.fastq.gz
  """
}
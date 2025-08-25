process JellyfishCount {
  tag "jellyfish"
  publishDir { "${params.outdir}/qc/filtered" }, mode: 'copy'

  input:
  path reads_fastq

  output:
  path "kmer_histo.txt", emit: kmer_histo

  container = 'community.wave.seqera.io/library/jellyfish:2.2.10--8e591e8a1910388f'

  script:
  """
  jellyfish count -C -m 21 -s 100M -t ${params.threads} <(zcat ${reads_fastq}) -o mer_counts.jf
  jellyfish histo -t ${params.threads} mer_counts.jf > kmer_histo.txt
  """
}
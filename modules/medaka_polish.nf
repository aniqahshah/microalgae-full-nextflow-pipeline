process MedakaPolish {
  tag "medaka"
   publishDir { "${params.outdir}/assembly/polished" }, mode: 'copy'

  input:
  tuple path(input_assembly), path(input_fastq)

  output:
  path "medaka_polished.fa", emit: polished_assembly

  container = 'community.wave.seqera.io/library/medaka:2.1.0--c9b2fb4c891009f4'

  script:
  """
  gunzip -c ${input_assembly} > uncompressed.fa
  medaka_consensus -i ${input_fastq} -d uncompressed.fa -o polished -t ${params.threads}
  cp polished/consensus.fasta medaka_polished.fa
  """
}
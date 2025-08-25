process GenomeScope {
  tag "genomescope"
  publishDir { "${params.outdir}/qc/filtered" }, mode: 'copy'

  input:
  path kmer_histo

  output:
  path "genomescope_output", emit: genomescope_dir
  path "genomescope_output/*", emit: genomescope_files

  script:
  """
  mkdir -p genomescope_output
  genomescope2 -i ${kmer_histo} -o genomescope_output -k 21
  """
}
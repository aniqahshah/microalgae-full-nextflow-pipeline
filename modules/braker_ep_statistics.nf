process AnnotationStatistics {
    tag "braker-ep statistics"
    publishDir { "${params.outdir}/annotation" }, mode: 'copy'
     
  input:
    path braker_gff3

  output:
    path "stats.tsv"

    
  """
  agat_sp_statistics.pl --gff ${braker_gff3} > stats.tsv
  """
}


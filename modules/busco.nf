process Busco {
  tag "busco"
  publishDir { "${params.outdir}/evaluation" }, mode: 'copy'

  input:
  tuple val(label), path(genome_fasta)

  output:
  path "busco_${label}", emit: busco_output

  script:
  """
  mkdir -p busco_${label}
  busco \
    -i ${genome_fasta} \
    -o busco_${label} \
    -l ${params.busco_lineage} \
    -m genome \
    --cpu ${task.cpus} \
    --out_path ./busco_${label}
  """
}
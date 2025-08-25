process Quast {
  tag "quast"
   publishDir { "${params.outdir}/evaluation" }, mode: 'copy'

  input:
  tuple val(label), path(assembly_file)

  output:
  path "quast_${label}", emit: quast_output

  script:
  """
  mkdir -p quast_${label}
  quast.py -o quast_${label} --threads ${task.cpus} ${assembly_file}
  """
}
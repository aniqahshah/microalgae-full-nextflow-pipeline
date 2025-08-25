process MultiQC {
  tag "multiqc"
  publishDir { "${params.outdir}/evaluation" }, mode: 'copy'

  input:
  path eval_inputs

  output:
  path "multiqc_report.html", emit: multiqc_report
  path "multiqc_data", emit: multiqc_data

  container = 'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0'

  script:
  """
  multiqc ${eval_inputs} -o . --force
  """
}
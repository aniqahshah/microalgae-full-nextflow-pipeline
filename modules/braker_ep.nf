process BrakerEP {
  tag "Braker-EP"
  publishDir { "${params.outdir}/annotation" }, mode: 'copy'
   
  input:
    path polished
    path proteins
    val clade
    path augustus_config

  output:
  tuple val(clade),
        path("${params.outdir}/annotation/braker_ep_${clade}/Augustus/augustus.hints.gff3"),
        path("${params.outdir}/annotation/braker_ep_${clade}/Augustus/augustus.hints.aa"),
        emit: braker

  script:
  """
  export AUGUSTUS_CONFIG_PATH=\$PWD/augustus_config
  echo "[INFO] Using AUGUSTUS_CONFIG_PATH: \$AUGUSTUS_CONFIG_PATH"

  rm -rf "\$AUGUSTUS_CONFIG_PATH/species/braker_ep_${clade}"

  mkdir -p ${params.outdir}/braker_ep_${clade}
  ABS_OUTDIR=\$PWD/${params.outdir}/braker_ep_${clade}

  braker.pl \\
    --genome=${polished} \\
    --prot_seq=${proteins} \\
    --softmasking \\
    --threads=\${NXF_TASK_CPUS:-16} \\
    --workingdir=\$ABS_OUTDIR \\
    --gff3 \\
    --overwrite \\
    --species=braker_ep_${clade} \\
    --min_contig=10000 \\
    --AUGUSTUS_CONFIG_PATH=\$AUGUSTUS_CONFIG_PATH \\
    --GENEMARK_PATH=/opt/gmes_linux_64

  """
}

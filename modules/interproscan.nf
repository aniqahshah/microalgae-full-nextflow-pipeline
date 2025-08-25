process Interproscan {
  tag "ipr"
  cpus 16
  input:
    path proteins_cleaned
  output:
    path "interproscan.tsv",  emit: tsv
    path "interproscan.gff3", emit: gff3
  container 'community.wave.seqera.io/library/interproscan:5.59_91.0--6053fb17325942d2'
  """
  interproscan.sh \
    -i ${proteins_cleaned} \
    -f tsv,gff3 \
    -b interproscan \
    -t p \
    -dp \
    -appl Pfam,SMART,TIGRFAM,CDD \
    -goterms \
    -pa \
    -cpu ${task.cpus} \
    -T ${task.workDir}/ipr_tmp
  """
}

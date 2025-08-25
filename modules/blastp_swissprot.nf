process BlastpSwissprot {
    tag "blastp vs swissprot"
    publishDir { "${params.outdir}/annotation/blastp_swissprot" }, mode: 'copy'

    input:
    path proteins
    val db_base
    path db_folder

    output:
    path "blastp.tsv", emit: blast_result

    container = 'ncbi/blast:latest'

    script:
    """
    blastp \\
      -query ${proteins} \\
      -db ${db_folder}/${db_base} \\
      -evalue 1e-5 \\
      -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' \\
      -num_threads 8 \\
      -out blastp.tsv
    """
}

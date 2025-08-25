process AnalyzeBlastp {
    tag "analyze blastp"
    publishDir { "${params.outdir}/annotation/blastp_swissprot" }, mode: 'copy'

    input:
    tuple path(blast_result), path(protein_fasta)

    output:
    path "blastp_summary.tsv"

    script:
    """
    total=\$(grep -c "^>" "${protein_fasta}")

    awk -v total="\$total" '{
        id[\$1] = \$3; cov[\$1] = \$13;
        if (\$3 >= 50 && \$13 >= 70) strong[\$1]=1;
        else if (\$3 >= 30 && \$13 >= 50) weak[\$1]=1;
    } END {
        for (i in id) {
            sum_id += id[i]; sum_cov += cov[i]; n++
        }
        printf("Total query proteins\\t%s\\n", total)
        printf("Proteins with BLAST hits\\t%s\\n", n)
        printf("Proteins with strong annotation\\t%s\\n", length(strong))
        printf("Proteins with weak annotation\\t%s\\n", length(weak))
        printf("Proteins with no BLAST hit\\t%s\\n", total - n)
        printf("Average %% identity\\t%.2f\\n", sum_id/n)
        printf("Average %% coverage\\t%.2f\\n", sum_cov/n)
    }' "${blast_result}" > blastp_summary.tsv
    """
}

process BuscoAnnotation {
    tag "busco for braker-ep"
     publishDir { "${params.outdir}/annotation" }, mode: 'copy'

    input:
        path proteins

    output:
        path "busco_braker_ep/short_summary*.txt", emit: busco_summary, optional: true

    script:
    """
    busco \
      -i ${proteins} \
      -l ${params.busco_lineage} \
      -m protein \
      -o busco_braker_ep \
      --out_path . \
      --download_path busco_downloads
    """
}

nextflow.enable.dsl=2

// Require input FASTQ
if (!params.input_fastq) {
    error "❌ Please provide --input_fastq"
}

// Dynamically set output folder based on sample name
def sample_name = file(params.input_fastq).baseName.replaceAll(/\\.fastq(.gz)?$/, '')
params.outdir = "results/${sample_name}"

// Show output directory at start
log.info "📂 Output directory set to: ${params.outdir}"

// Create main and subfolders
[
    "qc/raw",
    "qc/filtered",
    "assembly/unpolished",
    "assembly/polished",
    "evaluation",
    "annotation/blastp_swissprot",
    "annotation/interproscan",
    "annotation/eggnog",
    "annotation/merged"
].each { subdir ->
    new File("${params.outdir}/${subdir}").mkdirs()
    log.info "📁 Created: ${params.outdir}/${subdir}"
}

// Pipeline modules
include { RawQC         } from './modules/raw_qc.nf'
include { FilterReads   } from './modules/filter_reads.nf'
include { FilteredQC    } from './modules/filtered_qc.nf'
include { JellyfishCount } from './modules/jellyfish_count.nf'
include { GenomeScope   } from './modules/genomescope.nf'
include { SmartDenovo   } from './modules/smartdenovo.nf'
include { MedakaPolish  } from './modules/medaka_polish.nf'
include { Quast         } from './modules/quast.nf'
include { Busco         } from './modules/busco.nf'
include { MultiQC       } from './modules/multiqc.nf'
include { Quast as QuastPolished } from './modules/quast.nf'
include { Busco as BuscoPolished } from './modules/busco.nf'
include { BrakerEP } from './modules/braker_ep.nf'
include { BuscoAnnotation }       from './modules/braker_ep_busco.nf'
include { AnnotationStatistics }  from './modules/braker_ep_statistics.nf'
include { BlastpSwissprot } from './modules/blastp_swissprot.nf'
include { AnalyzeBlastp   } from './modules/blastp_analyze.nf'
include { SanitizeProteins } from './modules/sanitize_proteins.nf'
include { Interproscan } from './modules/interproscan.nf'
include { InterproscanSummary } from './modules/interproscan_summary.nf'
include { EggnogMapper } from './modules/eggnog_mapper.nf'
include { EggnogSummary } from './modules/eggnog_summary.nf'
include { AnnotationMerge } from './modules/annotation_merge.nf'

if (!params.swissprot_db) {
  params.swissprot_db = "swissprot_db/uniprot_sprot.fasta"
}

workflow {

  // --- Input Channel ---
  Channel.fromPath(params.input_fastq).set { raw_fastq_ch }

  // --- QC Phase ---
  if (params.run_qc) {
    RawQC(raw_fastq_ch)

    filtered_reads_ch = FilterReads(raw_fastq_ch).filtered_reads

    FilteredQC(filtered_reads_ch)

    kmer_histo_ch = JellyfishCount(filtered_reads_ch).kmer_histo
    GenomeScope(kmer_histo_ch)
  } else {
    filtered_reads_ch = raw_fastq_ch
  }

  // --- Assembly Phase ---
  if (params.run_assembly) {
    assembled_ch = SmartDenovo(filtered_reads_ch).assembled_genome

    medaka_input_ch = assembled_ch.combine(raw_fastq_ch).map { asm, fq -> tuple(asm, fq) }
    polished_ch = MedakaPolish(medaka_input_ch).polished_assembly
  } else {
    assembled_ch = Channel.empty()
    polished_ch  = Channel.empty()
  }

  // --- Evaluation Phase ---
  if (params.run_eval) {
    busco_raw_input_ch       = assembled_ch.map { tuple("unpolished", it) }
    quast_raw_input_ch       = assembled_ch.map { tuple("unpolished", it) }
    busco_polished_input_ch  = polished_ch.map  { tuple("polished", it) }
    quast_polished_input_ch  = polished_ch.map  { tuple("polished", it) }

    busco_raw_ch       = Busco(busco_raw_input_ch).busco_output.map         { [ "unpolished", it ] }
    quast_raw_ch       = Quast(quast_raw_input_ch).quast_output.map         { [ "unpolished", it ] }
    busco_polished_ch  = BuscoPolished(busco_polished_input_ch).busco_output.map { [ "polished", it ] }
    quast_polished_ch  = QuastPolished(quast_polished_input_ch).quast_output.map { [ "polished", it ] }

    eval_inputs_ch = quast_raw_ch
      .merge(busco_raw_ch)
      .merge(quast_polished_ch)
      .merge(busco_polished_ch)

    multiqc_input_ch = eval_inputs_ch.map { it[1] }
    MultiQC(multiqc_input_ch)
  }

  // --- Annotation Phase ---
 if (params.run_annotation) {

  if (params.clade == null) {
    error "❌ Please provide a clade name using --clade"
  }

  // --- Protein evidence ---
  protein_file = Channel.fromPath("protein_dbs/${params.clade}_80pct.fasta", checkIfExists: true)
    .ifEmpty { error "❌ Protein file not found for clade: ${params.clade}" }

  // --- AUGUSTUS config ---
  augustus_config_ch = Channel.fromPath("${params.augustus_config}/", checkIfExists: true)
    .ifEmpty { error "❌ AUGUSTUS config folder not found at ${params.augustus_config}" }

  // --- Run BRAKER-EP ---
  braker_out_ch = BrakerEP(polished_ch, protein_file, params.clade, augustus_config_ch).braker

  // braker_out_ch = tuple(clade, gff3, protein)
  braker_gff3_ch     = braker_out_ch.map{ it[1] }
  braker_proteins_ch = braker_out_ch.map{ it[2] }

  // --- Gene Prediction QC ---
  BuscoAnnotation(braker_proteins_ch)
  AnnotationStatistics(braker_gff3_ch)

  // --- Functional Annotation with BLASTP-SwissProt ---
  blastp_result_ch = BlastpSwissprot(
  braker_proteins_ch,
  Channel.value("uniprot_sprot.fasta"),
  Channel.fromPath("swissprot_db", checkIfExists: true)
  )

  AnalyzeBlastp(blastp_result_ch.combine(braker_proteins_ch))

  // Clean BRAKER proteins before InterProScan
  cleaned_faa_ch = SanitizeProteins(braker_proteins_ch).cleaned

  // Run InterProScan on the cleaned protein FASTA
  ipr_out = Interproscan(cleaned_faa_ch)

  // optional OBO value (string or null)
  def obo_path = file("${projectDir}/go-basic.obo").exists() ? "${projectDir}/go-basic.obo" : null
  go_obo_val_ch = Channel.value(obo_path)

  // NOW: pass script 3rd, OBO 4th
  summary = InterproscanSummary(
    ipr_out.tsv,                                 
    cleaned_faa_ch,                              
    file("${projectDir}/scripts/summarize_interproscan.py"), 
    go_obo_val_ch                               
  )

  // --- EggNOG Annotation ---
  eggnog_result_ch = EggnogMapper(braker_proteins_ch)

  // --- EggNOG Summary ---
  eggnog_summary_script_ch = Channel.fromPath("${projectDir}/scripts/summarize_eggnog.py", checkIfExists: true)
    .ifEmpty { error "❌ summarize_eggnog.py not found in scripts/" }

  eggnog_summary_ch = EggnogSummary(
    eggnog_result_ch.annotations,
    eggnog_summary_script_ch
  )
  
  // Merge annotations
  annotation_merge_script_ch = Channel.fromPath("scripts/annotation_merge.py", checkIfExists: true)


  AnnotationMerge(
    ipr_out.tsv,
    blastp_result_ch,
    eggnog_result_ch.annotations,
    braker_proteins_ch,
    annotation_merge_script_ch
  )
}
}

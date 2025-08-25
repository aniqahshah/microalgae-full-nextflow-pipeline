# 🧬 Microalgae ONT Genome Assembly & Annotation — Nextflow DSL2

A reproducible Nextflow pipeline tailored for **Oxford Nanopore** microalgae genomes.  
It automates QC → assembly → polishing → evaluation → structural & functional annotation.

> Maintainer: **Dr. Aniqah Shahida Mahamudin (INBIOSIS, UKM)**

---

## ✨ Features (modules)

- **QC**
  - `RawQC`: NanoPlot (raw)
  - `FilterReads`: NanoFilt
  - `FilteredQC`: NanoPlot (filtered)
- **K‑mer profiling**
  - `JellyfishCount` → `GenomeScope`
- **Assembly & polishing**
  - `SmartDenovo` (primary assembler)
  - `MedakaPolish` (ONT polishing)
- **Evaluation**
  - `Quast`, `Busco` (assembly and polished)
- **Structural annotation**
  - `BrakerEP` (AUGUSTUS + ProtHint; OrthoDB/SwissProt proteins)
- **Functional annotation**
  - `Interproscan` + summary
  - `EggnogMapper` + `EggnogSummary`
  - `BlastpSwissprot` + simple parser
- **Post‑processing**
  - `AnnotationMerge` (InterProScan + EggNOG + BLASTp → one TSV/Excel)
  - `AnnotationStatistics` (GFF3 stats with AGAT)

Each module runs in its own process with channels connecting outputs (DSL2).

---

## 🚀 Quickstart

```bash
# clone and enter
git clone https://github.com/aniqahshah//microalgae-full-nextflow-pipeline.git
cd microalgae-full-nextflow-pipeline

# run with Docker (recommended)
nextflow run main.nf   --input data/sample.fastq.gz   --clade Chlorophyta   --threads 16   --outdir results/sample1   -with-docker
```

### Common parameters
- `--input` *(required)*: FASTQ or FASTQ.GZ (single sample)
- `--clade` *(required for BRAKER-EP)*: e.g. `Chlorophyta` / `Viridiplantae`
- `--threads`: default 16
- `--outdir`: default auto‑derived from sample name
- `--busco_lineage`: path to BUSCO lineage, e.g. `chlorophyta_odb10`
- `--augustus_config`: path to AUGUSTUS config (auto-prepared if missing with teambraker image)
- `--protein_db`: OrthoDB/SwissProt/Chlorophyta FASTA for BRAKER‑EP
- `--eggnog_data`: path to eggNOG-mapper data (if using local db)
- `--swissprot_dmnd`: pre‑built DIAMOND DB (SwissProt)

> See `nextflow.config` for defaults and container images. Adjust for Singularity/Conda if preferred.

---

## 🧱 Repository layout

```
.
├── main.nf
├── nextflow.config
├── modules/
│   ├── raw_qc.nf
│   ├── filter_reads.nf
│   ├── filtered_qc.nf
│   ├── jellyfish_count.nf
│   ├── genomescope.nf
│   ├── smartdenovo.nf
│   ├── medaka_polish.nf
│   ├── quast.nf
│   ├── busco.nf
│   ├── braker_ep.nf
│   ├── interproscan.nf
│   ├── eggnogmapper.nf
│   ├── eggnog_summary.nf
│   ├── blastp_swissprot.nf
│   ├── annotation_merge.nf
│   └── annotation_stats.nf
├── bin/
│   ├── analyze_blastp.py
│   ├── annotation_merge.py
│   └── summarize_interproscan.py
├── docs/
│   └── pipeline_overview.png
├── data/                       # optional tiny example data
├── environment.yml             # optional conda
├── .gitignore
├── LICENSE
└── README.md
```

---

## 🔧 Requirements

- Nextflow **>= 24.10**
- One of: Docker / Podman / Singularity (recommended), or Conda
- For heavy modules:
  - BRAKER-EP: `teambraker/braker3:latest` image, `AUGUSTUS_CONFIG_PATH`
  - InterProScan: container or local install (5.75+)
  - BUSCO: lineage dataset (e.g., `chlorophyta_odb10`)

---

## 🧪 Example commands

### 1) Minimal end-to-end
```bash
nextflow run main.nf --input data/sample.fastq.gz --clade Chlorophyta -with-docker
```

### 2) Explicit resources and paths
```bash
nextflow run main.nf   --input data/sample.fastq.gz   --clade Chlorophyta   --busco_lineage chlorophyta_odb10   --protein_db /path/to/Chlorophyta_80pct.fasta   --augustus_config ./augustus_config   --eggnog_data /path/to/eggnog_mapper_data   --swissprot_dmnd /path/to/swissprot.dmnd   --threads 32 -with-docker
```

### 3) Resume after adding publishDir to InterProScan summary
```bash
nextflow run main.nf -resume -with-docker
```

---

## 📤 Outputs (key files)

- `evaluation/QUAST/` — HTML/TXT stats
- `evaluation/BUSCO/` — BUSCO summary TSVs
- `annotation/braker_ep/` — `augustus.hints.gff3`, `augustus.hints.aa`
- `annotation/interproscan/` — `interproscan.tsv`, summary plots
- `annotation/eggnog/` — `eggnog_annotation.tsv`, summaries
- `annotation/blastp_swissprot/` — `blastp.tsv`, parsed TSV
- `annotation/merged/` — unified `annotations_merged.tsv` (and `.xlsx`)

---

## 🧰 Troubleshooting

- **`agat_sp_statistics.pl: command not found`** → ensure AGAT is in the container / conda env.
- **`ModuleNotFoundError: pandas`** → add `pip install pandas` in the module container/script.
- **InterProScan container tag not found** → use a valid tag (e.g., `5.75-106.0`) or local install.
- **`AUGUSTUS_CONFIG_PATH` issues** → pipeline prepares `augustus_config/` via `teambraker/braker3` if missing.

---

## 📜 License
MIT — see [LICENSE](LICENSE).

## 🙌 Citation
If this pipeline helps your research, please cite your project and tools used (AUGUSTUS, BRAKER, InterProScan, BUSCO, eggNOG, etc.).

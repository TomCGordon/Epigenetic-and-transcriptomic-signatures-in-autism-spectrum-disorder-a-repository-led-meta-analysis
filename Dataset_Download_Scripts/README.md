# ASD Multi-Omic Meta-Analysis: Dataset Download & Reference Staging Suite

This repository contains automated data download and reference staging scripts (available in both **Python** and **R**) to reproduce the raw data ingestion pipeline for our autism spectrum disorder (ASD) multi-omic meta-analysis.

Using these scripts, you can automatically retrieve all public datasets analyzed in our study directly from NCBI GEO, ArrayExpress, and institutional repositories, with built-in checksum verification and extraction.

---

## 1. Quick Start: Choosing Your Language

Both script suites perform the exact same downloading, checksum verification, and file extraction tasks. Choose whichever environment you prefer.

### Option A: Using Python (Python 3.8+)

#### Dependencies
Python uses standard library modules (`urllib`, `hashlib`, `tarfile`, `gzip`, `argparse`) with one lightweight dependency:
```bash
pip install pandas
```

#### Execution
```bash
# 1. Download all required reference resources (GPL maps, HGNC tables, UCSC promoter coordinates)
python Supporting_resources/download_reference_resources.py

# 2. Build the Illumina 450K promoter annotation file
python Supporting_resources/prepare_illumina450k_annotation.py

# 3. Dry-run: preview files to be downloaded without fetching them
python download_all_datasets.py --dry-run

# 4. Full download: download and extract all datasets
python download_all_datasets.py --verify-sha256
```

---

### Option B: Using R (R 4.0+)

#### Dependencies
R uses base packages (`utils`, `tools`) alongside `digest` (for checksums) and Bioconductor packages for annotation reconstruction:
```R
install.packages("digest")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("minfi", "IlluminaHumanMethylation450kanno.ilmn12.hg19"))
```

#### Execution
```bash
# 1. Download all required reference resources
Rscript Supporting_resources/download_reference_resources.R

# 2. Build the Illumina 450K promoter annotation file
Rscript Supporting_resources/prepare_illumina450k_annotation.R

# 3. Dry-run: preview downloads
Rscript download_all_datasets.R --dry-run

# 4. Full download: download and extract all datasets
Rscript download_all_datasets.R --verify-sha256
```

---

## 2. Genomics Primer: Understanding GEO Accessions & Platforms

If you are new to NCBI GEO (Gene Expression Omnibus) data structures, here is a quick overview of key terms used throughout this repository and our paper:

| Term | Full Name | What It Means | Example in Our Paper |
| :--- | :--- | :--- | :--- |
| **GSE** | GEO Series | The unique accession ID assigned to a complete published experiment or dataset. It contains the raw measurement matrices and sample metadata. | `GSE108785` (Blood Methylation), `GSE113834` (Brain Expression) |
| **GPL** | GEO Platform | The manufacturer's reference map for the physical assay hardware (e.g., a specific microarray chip or sequencing platform). | `GPL13534` (Illumina 450K Methylation), `GPL570` (Affymetrix Microarray) |
| **GSM** | GEO Sample | The accession ID for an individual biological sample within a dataset. | `GSM2898741` |

### Omic Technologies Included in Our Paper

Our meta-analysis synthesizes data across 6 biological tissue and omic branches:

1. **DNA Methylation Arrays (HM27, 450K, EPIC)**: Measure chemical methylation ($\beta$-values) at individual CpG sites across the human genome.
2. **Gene Expression Microarrays**: Measure mRNA expression levels using physical hybridisation probes (e.g., Affymetrix, Agilent, Illumina).
3. **Bulk RNA-Seq (Sequencing)**: High-throughput sequencing read counts mapped to Ensembl gene identifiers.
4. **Whole-Genome Bisulfite Sequencing (WGBS) & MeDIP-Seq**: Next-generation sequencing providing single-base or region-level methylation coverage across promoters.

---

## 3. What the Supporting Resources Do & Why They Are Needed

Before processing raw dataset files, raw probe IDs and genomic coordinates must be translated into standardized human gene symbols (*HGNC symbols* like *MECP2* or *SHANK3*). The scripts in `Supporting_resources/` handle this preparation:

### Script 1: `download_reference_resources` (`.py` / `.R`)

This script downloads 4 essential categories of reference annotations:

1. **GPL Platform Annotations (`GPL*.annot.gz` / `GPL*_family.soft.gz`)**:
   * *What it is*: Official manufacturer mapping tables from NCBI GEO.
   * *Why it's needed*: Microarrays use proprietary probe IDs (e.g., `1007_s_at` or `ILMN_1343291`). These GPL files allow the pipeline to map array probes to real gene symbols.

2. **HGNC Gene Nomenclature Table (`hgnc_complete_set.txt`)**:
   * *What it is*: The official HUGO Gene Nomenclature Committee database.
   * *Why it's needed*: RNA-seq datasets report Ensembl IDs (`ENSG00000168036`) or outdated gene synonyms. This table converts all identifiers to up-to-date HGNC symbols so genes match 1:1 across different datasets.

3. **UCSC RefGene Coordinate Files (`hg18`, `hg19`, `hg38`)**:
   * *What it is*: Genomic coordinate tables specifying gene start and end positions.
   * *Why it's needed*: Next-generation sequencing datasets (WGBS / MeDIP) report single CpG positions (`chr1:123456`). These files define promoter regions ($\pm 2\text{kb}$ around the Transcription Start Site) so CpG methylation can be aggregated per gene.

4. **GENCODE Transcript References (`gencode.v19.pc_transcripts.fa.gz`)**:
   * *What it is*: FASTA sequences for human protein-coding transcripts (GRCh37/hg19).
   * *Why it's needed*: Used as reference transcripts for pseudo-alignment or transcript quantification pipelines.

---

### Script 2: `prepare_illumina450k_annotation` (`.py` / `.R`)

* **What it does**: Parses the Illumina 450K array manifest to extract probe IDs, gene symbols, genomic positions, CpG island relations, and promoter designations (TSS200, TSS1500, 5'UTR, 1stExon).
* **Output**: Produces `illumina450k_annotation_core.csv` in `03_Required_Annotations_and_Metadata/`.
* **Why it's needed**: Enables promoter-specific CpG filtering for Illumina 450K and EPIC DNA methylation array datasets.

---

## 4. Directory Structure & Output Hierarchy

When you run the scripts, files are organized into a standardized directory structure:

```
Dataset_Download_Scripts/
├── download_all_datasets.py         # Python master download CLI
├── download_all_datasets.R          # R master download CLI
├── Supporting_resources/
│   ├── download_reference_resources.py / .R
│   └── prepare_illumina450k_annotation.py / .R
│
└── [Default Output Root: D:\ASD_Source_Data_Audit\]
    ├── 01_Raw_Data/                  # Raw GSE series matrices, tarballs, & count tables
    │   ├── blood_methylation/
    │   ├── brain_methylation/
    │   ├── blood_expression/
    │   ├── brain_expression/
    │   ├── placenta_lcl_methylation/
    │   └── placenta_lcl_expression/
    │
    ├── 03_Required_Annotations_and_Metadata/  # Staged GPLs, HGNC tables, & promoter maps
    └── 04_Download_Logs/             # Time-stamped execution & checksum logs
```

---

## 5. Command-Line Options

Both `download_all_datasets.py` and `download_all_datasets.R` accept the following command-line flags:

| Flag | Short | Description | Default |
| :--- | :--- | :--- | :--- |
| `--output-root` | `-o` | Base directory where downloaded files will be stored | `D:\ASD_Source_Data_Audit` |
| `--branch` | `-b` | Restrict downloads to a specific tissue/omic branch (e.g., `blood_methylation`, `brain_expression`) | All branches |
| `--dataset` | `-d` | Download a single dataset by GSE ID (e.g., `GSE108785`) | All datasets |
| `--dry-run` | `-n` | Print download URLs and target file paths without fetching files | `False` |
| `--verify-sha256` | `-v` | Compute and verify SHA-256 checksums after downloading | `False` |
| `--overwrite` | `-f` | Force re-downloading files even if they already exist locally | `False` |
| `--retries` | `-r` | Number of download retry attempts for transient network drops | `3` |
| `--timeout` | `-t` | HTTP request timeout in seconds | `60` |

---

## 6. Summary Table of Included Datasets

The repository manages downloads for **45 public datasets (46 target matrices)** across 6 biological branches:

| Branch | Datasets Included | Assay Technology | Target Matrix Format |
| :--- | :--- | :--- | :--- |
| **Blood Methylation** | GSE108785, GSE109905, GSE113967, GSE140730, GSE27044, GSE83424 | HM27, 450K, EPIC, WGBS | Series matrices & Bismark CpG reports |
| **Brain Methylation** | GSE109875, GSE131706, GSE242427, GSE278285, GSE38608, GSE53162, GSE53924, GSE80017, GSE81541 | 450K Array & WGBS | Series matrices, beta tables, & WGBS reports |
| **Blood Expression** | GSE111175, GSE111176, GSE123302, GSE140702, GSE18123 (GPL570/GPL6244), GSE25507, GSE26415, GSE77103 | Microarray & Bulk RNA-Seq | Series matrices & per-sample count files |
| **Brain Expression** | GSE102741, GSE113834, GSE211154, GSE236761, GSE269105, GSE28475, GSE28521, GSE36315, GSE38322, GSE59288, GSE62098, GSE64018 | Microarray & Bulk RNA-Seq | Series matrices, Excel workbooks, & FPKM/count tables |
| **Placenta/LCL Methylation** | GSE178203, GSE34099, GSE67615, GSE99935 | 450K Array, MeDIP, & WGBS | Series matrices, MeDIP matrices, & WGBS archives |
| **Placenta/LCL Expression** | GSE15402, GSE15451, GSE29918, GSE37772, GSE4187, GSE7329 | Microarray | Series matrices & GPL platform maps |

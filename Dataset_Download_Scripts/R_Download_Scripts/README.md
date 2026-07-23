# Dataset-Specific Public Data Download Scripts (R Edition)

This directory contains standalone, reproducible R scripts for downloading and staging all public repository datasets used in the ASD multi-omic meta-analysis (covering **DNA Methylation** and **Gene Expression** across 6 biological branches).

Each dataset script is completely self-contained and retrieves raw repository matrices, SOFT files, count matrices, CpG reports, BED files, or supplementary archives directly from GEO, SRA, and public archives.

---

## 1. Directory Structure

```text
Dataset_Specific_Download_Scripts_GitHub_Ready/
├── DNA_methylation/
│   ├── Blood_family/                 # 6 blood methylation datasets (GSE108785, GSE109905, GSE113967, GSE140730, GSE27044, GSE83424)
│   ├── Lymphoblastoid_cell_lines/    # 2 LCL methylation datasets (GSE34099, GSE99935)
│   ├── Placenta/                     # 2 placenta WGBS datasets (GSE178203, GSE67615)
│   └── Post_mortem_brain/            # 9 brain methylation datasets (GSE109875, GSE131706, GSE242427, GSE278285, GSE38608, GSE53162, GSE53924, GSE80017, GSE81541)
├── Gene_expression/
│   ├── Blood_family/                 # 8 blood expression datasets (GSE111175, GSE111176, GSE123302, GSE140702, GSE18123, GSE25507, GSE26415, GSE77103)
│   ├── Lymphoblastoid_cell_lines/    # 6 LCL expression datasets (GSE15402, GSE15451, GSE29918, GSE37772, GSE4187, GSE7329)
│   └── Post_mortem_brain/            # 12 brain expression datasets (GSE102741, GSE113834, GSE211154, GSE236761, GSE269105, GSE28475, GSE28521, GSE36315, GSE38322, GSE59288, GSE62098, GSE64018)
├── Supporting_resources/
│   ├── download_reference_resources.R       # Downloads platform annotations, HGNC tables, GENCODE transcripts, and UCSC RefGene files
│   └── prepare_illumina450k_annotation.R    # Reconstructs Illumina 450K core annotation table
└── README.md                         # This user guide
```

---

## 2. Requirements & Dependencies

### R Version & Base Dependencies
- **R**: Version 4.0 or higher (64-bit).
- **Base Packages**: `utils`, `tools` (included with base R).

### Optional & Specialized Packages
- **SHA-256 Checksum Verification**: Install the `digest` package:
  ```r
  install.packages("digest")
  ```
- **Illumina 450K Annotation Reconstruction (`prepare_illumina450k_annotation.R`)**:
  Requires Bioconductor annotation packages:
  ```r
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install(c("minfi", "IlluminaHumanMethylation450kanno.ilmn12.hg19"))
  ```

---

## 3. Supporting Resources Explanation

The `Supporting_resources/` directory provides scripts for staging essential non-dataset-specific reference files required by downstream processing:

1. **`download_reference_resources.R`**:
   - **GPL Platform Annotations**: Retrieves microarray annotation files (`.annot.gz` and `.soft.gz`) for GPL platforms including `GPL10558`, `GPL16686`, `GPL570`, `GPL6244`, `GPL6480`, `GPL13388`, `GPL15207`, `GPL15314`, `GPL13158`, `GPL1708`, `GPL3427`, and `GPL6883`.
   - **Gene Nomenclature Tables**: Downloads the HGNC complete set (`hgnc_complete_set.txt`) to standardize gene symbols across expression branches.
   - **Transcript Reference Databases**: Downloads GENCODE v19 lncRNA and protein-coding transcript FASTA files for custom annotation pipelines (e.g. GSE36315).
   - **Genomic Coordinate Files**: Downloads UCSC RefGene tables (`hg18_refGene.txt.gz`, `hg19_refGene.txt.gz`, `hg38_refGene.txt.gz`) for promoter coordinate window extraction ($\pm 2\text{kb}$ around TSS) used in WGBS and MeDIP-seq datasets.

2. **`prepare_illumina450k_annotation.R`**:
   - Reconstructs `illumina450k_annotation_core.csv`, extracting probe IDs, promoter annotations (`UCSC_RefGene_Group`), gene symbols (`UCSC_RefGene_Name`), CpG island relations (`Relation_to_Island`), and genomic coordinates (`chr`, `pos`).

---

## 4. How to Use the Download Scripts

### Basic Execution Example
Run any dataset download script by providing an `--output-root` directory where data will be stored:

```bash
Rscript DNA_methylation/Blood_family/download_GSE108785.R --output-root="D:/ASD_public_inputs"
```

### Reference Resources Staging Example
Staging all required GPL platform annotations and reference genome files:

```bash
Rscript Supporting_resources/download_reference_resources.R --output-root="D:/ASD_public_inputs"
```

---

## 5. Command-Line Options

Each script accepts the following flags:

| Flag | Description | Default |
|---|---|---|
| `--output-root=PATH` | Destination root directory for downloaded files | `./downloaded_public_inputs` |
| `--dry-run` | Display URLs and destination paths without downloading | Disabled |
| `--overwrite` | Force redownload and replace existing files | Disabled |
| `--verify-sha256` | Verify SHA-256 hashes against recorded expected hashes | Disabled |
| `--retries=N` | Maximum download retry attempts per file | `3` |
| `--timeout=N` | HTTP download timeout in seconds | `7200` |

### Recommended Workflow:
1. **Dry-Run Inspection**: Test execution without downloading:
   ```bash
   Rscript DNA_methylation/Placenta/download_GSE178203.R --output-root="D:/ASD_public_inputs" --dry-run
   ```
2. **Download & Verification**: Download data with SHA-256 verification:
   ```bash
   Rscript DNA_methylation/Placenta/download_GSE178203.R --output-root="D:/ASD_public_inputs" --verify-sha256
   ```

---

## 6. Output Folder Hierarchy

Downloaded files will be organized in the standardized directory structure expected by downstream R and Python analytical pipelines:

```text
downloaded_public_inputs/
├── 01_Raw_Data/
│   ├── DNA_methylation/
│   │   ├── Blood_family/
│   │   ├── Lymphoblastoid_cell_lines/
│   │   ├── Placenta/
│   │   └── Post_mortem_brain/
│   └── Gene_expression/
│       ├── Blood_family/
│       ├── Lymphoblastoid_cell_lines/
│       └── Post_mortem_brain/
├── 03_Required_Annotations_and_Metadata/
│   ├── expression_platform_annotations/
│   ├── branch_specific_gene_nomenclature/
│   └── data_raw_annotation/
└── 04_Download_Logs/
    ├── GSE108785_download_log.csv
    └── GSE108785_file_check.csv
```

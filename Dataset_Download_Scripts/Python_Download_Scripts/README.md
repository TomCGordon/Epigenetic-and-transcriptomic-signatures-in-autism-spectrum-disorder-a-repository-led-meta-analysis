# Dataset-Specific Public Data Download Scripts (Python Edition)

This directory contains standalone, reproducible Python scripts for downloading and staging all public repository datasets used in the ASD multi-omic meta-analysis (covering **DNA Methylation** and **Gene Expression** across 6 biological branches).

Each Python download script is completely self-contained and retrieves raw repository matrices, SOFT files, count matrices, CpG reports, BED files, or supplementary archives directly from GEO, SRA, and public archives.

---

## 1. Directory Structure

```text
Dataset_Specific_Download_Scripts_Python/
├── DNA_methylation/
│   ├── Blood_family/                 # 6 blood methylation datasets (download_GSE108785.py, download_GSE109905.py, etc.)
│   ├── Lymphoblastoid_cell_lines/    # 2 LCL methylation datasets (download_GSE34099.py, download_GSE99935.py)
│   ├── Placenta/                     # 2 placenta WGBS datasets (download_GSE178203.py, download_GSE67615.py)
│   └── Post_mortem_brain/            # 9 brain methylation datasets (download_GSE109875.py, download_GSE131706.py, etc.)
├── Gene_expression/
│   ├── Blood_family/                 # 8 blood expression datasets (download_GSE111175.py, download_GSE111176.py, etc.)
│   ├── Lymphoblastoid_cell_lines/    # 6 LCL expression datasets (download_GSE15402.py, download_GSE15451.py, etc.)
│   └── Post_mortem_brain/            # 12 brain expression datasets (download_GSE102741.py, download_GSE113834.py, etc.)
├── Supporting_resources/
│   ├── download_reference_resources.py       # Downloads platform annotations, HGNC tables, GENCODE transcripts, and UCSC RefGene files
│   └── prepare_illumina450k_annotation.py    # Prepares Illumina 450K core annotation table
└── README.md                         # This user guide
```

---

## 2. Requirements & Dependencies

### Python Environment
- **Python**: Version 3.8 or higher.
- **Built-in Standard Modules**: `os`, `sys`, `time`, `shutil`, `urllib.request`, `hashlib`, `tarfile`, `argparse`, `datetime`.

### Required External Packages
Install `pandas` for log and file check generation:
```bash
pip install pandas
```

---

## 3. Supporting Resources Explanation

The `Supporting_resources/` directory provides scripts for staging essential non-dataset-specific reference files required by downstream processing:

1. **`download_reference_resources.py`**:
   - **GPL Platform Annotations**: Retrieves microarray annotation files (`.annot.gz` and `.soft.gz`) for GPL platforms including `GPL10558`, `GPL16686`, `GPL570`, `GPL6244`, `GPL6480`, `GPL13388`, `GPL15207`, `GPL15314`, `GPL13158`, `GPL1708`, `GPL3427`, and `GPL6883`.
   - **Gene Nomenclature Tables**: Downloads the HGNC complete set (`hgnc_complete_set.txt`) to standardize gene symbols across expression branches.
   - **Transcript Reference Databases**: Downloads GENCODE v19 lncRNA and protein-coding transcript FASTA files for custom annotation pipelines (e.g. GSE36315).
   - **Genomic Coordinate Files**: Downloads UCSC RefGene tables (`hg18_refGene.txt.gz`, `hg19_refGene.txt.gz`, `hg38_refGene.txt.gz`) for promoter coordinate window extraction ($\pm 2\text{kb}$ around TSS) used in WGBS and MeDIP-seq datasets.

2. **`prepare_illumina450k_annotation.py`**:
   - Prepares the destination directory and core annotation path for `illumina450k_annotation_core.csv`.

---

## 4. How to Use the Download Scripts

### Basic Execution Example
Run any dataset download script by providing an `--output-root` directory where data will be stored:

```bash
python DNA_methylation/Blood_family/download_GSE108785.py --output-root="D:/ASD_public_inputs"
```

### Reference Resources Staging Example
Staging all required GPL platform annotations and reference genome files:

```bash
python Supporting_resources/download_reference_resources.py --output-root="D:/ASD_public_inputs"
```

---

## 5. Command-Line Options

Each script accepts the following flags:

| Flag | Description | Default |
|---|---|---|
| `--output-root PATH` | Destination root directory for downloaded files | `./downloaded_public_inputs` |
| `--dry-run` | Display URLs and destination paths without downloading | Disabled |
| `--overwrite` | Force redownload and replace existing files | Disabled |
| `--verify-sha256` | Verify SHA-256 hashes against recorded expected hashes | Disabled |
| `--retries N` | Maximum download retry attempts per file | `3` |
| `--timeout N` | HTTP download timeout in seconds | `7200` |

### Recommended Workflow:
1. **Dry-Run Inspection**: Test execution without downloading:
   ```bash
   python DNA_methylation/Placenta/download_GSE178203.py --output-root="D:/ASD_public_inputs" --dry-run
   ```
2. **Download & Verification**: Download data with SHA-256 verification:
   ```bash
   python DNA_methylation/Placenta/download_GSE178203.py --output-root="D:/ASD_public_inputs" --verify-sha256
   ```

---

## 6. Output Folder Hierarchy

Downloaded files will be organized in the standardized directory structure expected by downstream analytical pipelines:

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

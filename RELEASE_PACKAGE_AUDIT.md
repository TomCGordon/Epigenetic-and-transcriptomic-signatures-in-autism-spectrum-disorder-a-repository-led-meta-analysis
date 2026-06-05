# Release Package Audit

This note records the scope and verification checks for the public R code package.

## Final Model Provenance

The final reported pooled estimates should be traced to:

- `scripts/04_meta_analysis/01_fit_methylation_models_metafor.R`
- `scripts/04_meta_analysis/02_fit_expression_models_metafor.R`
- `scripts/04_meta_analysis/lib/metafor_helpers.R`

These scripts fit DerSimonian-Laird random-effects models with `metafor::rma.uni(method = "DL", test = "z")`, calculate FDR values, and classify modified Knapp-Hartung interval support using `metafor::rma.uni(method = "DL", test = "adhoc")`.

Branch pipeline folders generate public-source-derived gene summaries and dataset-level Hedges' g inputs. Some branch folders also write branch-level QC summaries so that tissue-specific processing can be inspected, but the reported pooled results should be read from `results/meta_analysis/`.

## Public-Facing Package Checks

- The upload folder contains only R scripts, Markdown documentation, CSV manifests, the path-template configuration file, and `.gitignore`.
- No raw data, generated model-result folders, downloaded source files, manuscript drafts, Word documents, spreadsheets, images, session files, or controlled-access files are present.
- Pipeline folders use publication-facing branch names such as `blood_methylation`, `brain_expression`, and `placenta_lcl_expression`.
- Code comments and documentation were checked for personal paths, non-public labels, and project-management wording.
- GSE36315 is labelled as a custom-annotation sensitivity analysis.
- Blank array cells are treated as missing by default. A blank-as-zero option is retained only as an explicit missing-data sensitivity setting.
- Controlled-access datasets are not analysed by this code package.

## Checks Run

- Release inventory written to `manifests/release_file_inventory.csv`.
- R syntax parse check: 66 R scripts checked, 0 parse failures.
- Red-flag wording scan: no matches for internal workflow phrases, project-history wording, personal paths, or non-public analysis language.

## Repository Scope

This repository is intended to contain the public R code, portable path template, documentation and manifests required to reproduce the reported public-data analyses. Raw downloaded data, controlled-access material, draft writing files and local generated result folders are outside the scope of this code repository. Selected derived non-identifiable outputs can be deposited separately with the journal or a data repository.

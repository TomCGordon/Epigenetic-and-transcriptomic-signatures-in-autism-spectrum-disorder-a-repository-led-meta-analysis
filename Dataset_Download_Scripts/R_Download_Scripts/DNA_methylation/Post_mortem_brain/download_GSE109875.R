#!/usr/bin/env Rscript

# Download and stage the public repository inputs used for GSE109875.
# This script is self-contained and may be run independently.

value_for <- function(args, prefix, default = NA_character_) {
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) default else sub(prefix, "", hit[[1]], fixed = TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
output_root <- value_for(args, "--output-root=", file.path(getwd(), "downloaded_public_inputs"))
overwrite <- "--overwrite" %in% args
dry_run <- "--dry-run" %in% args
verify_sha256 <- "--verify-sha256" %in% args
retries <- as.integer(value_for(args, "--retries=", "3"))
timeout <- as.integer(value_for(args, "--timeout=", "7200"))
if (is.na(retries) || retries < 1L) stop("--retries must be a positive integer.", call. = FALSE)
if (is.na(timeout) || timeout < 1L) stop("--timeout must be a positive integer.", call. = FALSE)
if (verify_sha256 && !requireNamespace("digest", quietly = TRUE)) {
  stop("Install the R package 'digest' to use --verify-sha256.", call. = FALSE)
}

normalise_path <- function(path) gsub("[\\\\/]+", .Platform$file.sep, path)

download_one <- function(url, destination) {
  if (dry_run) {
    message("[dry-run] ", url, " -> ", destination)
    return("dry_run")
  }
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(destination) && file.info(destination)$size > 0 && !overwrite) {
    message("[present] ", destination)
    return("already_present")
  }
  old_timeout <- getOption("timeout")
  options(timeout = max(timeout, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  temporary <- paste0(destination, ".part")
  if (file.exists(temporary)) unlink(temporary)
  last_error <- NULL
  for (attempt in seq_len(retries)) {
    message(sprintf("[download %d/%d] %s", attempt, retries, url))
    ok <- tryCatch({
      utils::download.file(url, temporary, mode = "wb", quiet = FALSE, method = "libcurl")
      file.exists(temporary) && file.info(temporary)$size > 0
    }, error = function(error) {
      last_error <<- conditionMessage(error)
      FALSE
    })
    if (ok) {
      if (file.exists(destination)) unlink(destination)
      if (!file.rename(temporary, destination)) {
        stop("Could not move completed download to ", destination, call. = FALSE)
      }
      return("downloaded")
    }
    if (file.exists(temporary)) unlink(temporary)
    if (attempt < retries) Sys.sleep(min(30, 2^attempt))
  }
  stop("Download failed for ", url,
       if (!is.null(last_error)) paste0(": ", last_error) else "", call. = FALSE)
}

downloads <- data.frame(
  remote_url = c(
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971944/suppl/GSM2971944_JLKD051_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971945/suppl/GSM2971945_JLKD052_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971946/suppl/GSM2971946_JLKD054_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971947/suppl/GSM2971947_JLKD055_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971948/suppl/GSM2971948_JLKD056_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971949/suppl/GSM2971949_JLKD057_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971950/suppl/GSM2971950_JLKD058_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971951/suppl/GSM2971951_JLKD059_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971952/suppl/GSM2971952_JLKD060_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971953/suppl/GSM2971953_JLKD061_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971954/suppl/GSM2971954_JLKD062_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971955/suppl/GSM2971955_JLKD063_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971956/suppl/GSM2971956_JLKD065_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971957/suppl/GSM2971957_JLKD066_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971958/suppl/GSM2971958_JLKD067_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971959/suppl/GSM2971959_JLKD069_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz"
  ),
  destination = c(
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971944.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971945.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971946.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971947.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971948.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971949.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971950.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971951.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971952.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971953.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971954.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971955.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971956.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971957.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971958.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971959.cpg_report.txt.gz"
  ),
  action = c(
    "none",
    "none",
    "none",
    "none",
    "none",
    "none",
    "none",
    "none",
    "none",
    "none",
    "none",
    "none",
    "none",
    "none",
    "none",
    "none"
  ),
  extraction_root = c(
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_
  ),
  expected_bytes = c(
    264284960,
    266852074,
    268162085,
    266937349,
    265581013,
    262852523,
    262658188,
    266345342,
    262298260,
    244560070,
    268523047,
    268286531,
    269230091,
    266725679,
    268041654,
    269663171
  ),
  expected_sha256 = c(
    "4a3c5c1401c9f3dea2326c053a81da508e1b52d15a8873395eb2e1746f945e07",
    "0d864da2d8827420adde8a4c32ec9df5a73fe66418284632639da7caea687cb8",
    "2749ced8526facd5993e07ad6b12f95bc9337a5d2e43c2b8df0e42bef39df1c6",
    "eba1f625556de2a50938b6b5a21253ad833b58e1e23e3e8f2583c18f90a09453",
    "1f43203efc9353c901cebb4589e7585e2c99e19a5885c9ec3f9ddbbc10aab0ab",
    "8935a9a92cc18a6451ecabf275901d6071a3eb30a7218bd7c1db695c77b269fe",
    "81a5ffedce6706cee64f0bee99bc1f145042ea1c8fcf8755dbc5f9414e8bd2a2",
    "944dd645f98a19c1dbeca95de1740d1eb49b4c883562f91947e4320d246074d3",
    "735cb407a7c7be2b403555d9b972bac3048275ca89c8d42f48e9ceb82fae460b",
    "8c842fc5b49e7e88f301f7919206d09b7e6c5a3f4cc1f28ccdbb882a0e70b8e6",
    "d5b4e564ea2aedee7488453e17171fb921fcbc0288c4b4ffcbddcb2e828aace7",
    "2d3443cea658f2a4c45e23c405bd8532cf58b75d59ec9dfc91c6975be8c366c1",
    "29c1ac2142e6c3b4a009fbd09d5293a1119a7f2c1ed097c0ce770d6d06d0c20d",
    "dec0a66bf0b8cd9cb2d4d155cfff8755bbfb76b96cedad83039e2f2985433720",
    "db3a481c398af7339e050c3b9bc2a712cbe4b0e2e007309286e1611b9bcd9967",
    "a435a07f4773fa2b678b423b8e38853c217c03acfd5344bf1ce16a585d2ce4ac"
  ),
  stringsAsFactors = FALSE
)

expected_files <- data.frame(
  path = c(
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971944.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971945.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971946.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971947.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971948.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971949.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971950.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971951.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971952.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971953.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971954.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971955.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971956.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971957.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971958.cpg_report.txt.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971959.cpg_report.txt.gz"
  ),
  expected_bytes = c(
    264284960,
    266852074,
    268162085,
    266937349,
    265581013,
    262852523,
    262658188,
    266345342,
    262298260,
    244560070,
    268523047,
    268286531,
    269230091,
    266725679,
    268041654,
    269663171
  ),
  expected_sha256 = c(
    "4a3c5c1401c9f3dea2326c053a81da508e1b52d15a8873395eb2e1746f945e07",
    "0d864da2d8827420adde8a4c32ec9df5a73fe66418284632639da7caea687cb8",
    "2749ced8526facd5993e07ad6b12f95bc9337a5d2e43c2b8df0e42bef39df1c6",
    "eba1f625556de2a50938b6b5a21253ad833b58e1e23e3e8f2583c18f90a09453",
    "1f43203efc9353c901cebb4589e7585e2c99e19a5885c9ec3f9ddbbc10aab0ab",
    "8935a9a92cc18a6451ecabf275901d6071a3eb30a7218bd7c1db695c77b269fe",
    "81a5ffedce6706cee64f0bee99bc1f145042ea1c8fcf8755dbc5f9414e8bd2a2",
    "944dd645f98a19c1dbeca95de1740d1eb49b4c883562f91947e4320d246074d3",
    "735cb407a7c7be2b403555d9b972bac3048275ca89c8d42f48e9ceb82fae460b",
    "8c842fc5b49e7e88f301f7919206d09b7e6c5a3f4cc1f28ccdbb882a0e70b8e6",
    "d5b4e564ea2aedee7488453e17171fb921fcbc0288c4b4ffcbddcb2e828aace7",
    "2d3443cea658f2a4c45e23c405bd8532cf58b75d59ec9dfc91c6975be8c366c1",
    "29c1ac2142e6c3b4a009fbd09d5293a1119a7f2c1ed097c0ce770d6d06d0c20d",
    "dec0a66bf0b8cd9cb2d4d155cfff8755bbfb76b96cedad83039e2f2985433720",
    "db3a481c398af7339e050c3b9bc2a712cbe4b0e2e007309286e1611b9bcd9967",
    "a435a07f4773fa2b678b423b8e38853c217c03acfd5344bf1ce16a585d2ce4ac"
  ),
  stringsAsFactors = FALSE
)

log_rows <- vector("list", nrow(downloads))
for (i in seq_len(nrow(downloads))) {
  row <- downloads[i, ]
  destination <- file.path(output_root, normalise_path(row$destination))
  status <- download_one(row$remote_url, destination)
  if (status %in% c("downloaded", "already_present") && row$action != "none") {
    extraction_root <- file.path(output_root, normalise_path(row$extraction_root))
    dir.create(extraction_root, recursive = TRUE, showWarnings = FALSE)
    utils::untar(destination, exdir = extraction_root)
    if (row$action == "extract_miniml_and_copy_tables") {
      tables <- list.files(extraction_root, pattern = "^GSM.*-tbl-1\\.txt$", full.names = TRUE)
      sample_dir <- file.path(extraction_root, "sample_tables")
      dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
      if (length(tables)) file.copy(tables, sample_dir, overwrite = TRUE)
    }
  }
  log_rows[[i]] <- data.frame(
    accession = "GSE109875",
    remote_url = row$remote_url,
    destination = destination,
    status = status,
    checked_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
}

if (!dry_run) {
  checks <- vector("list", nrow(expected_files))
  for (i in seq_len(nrow(expected_files))) {
    row <- expected_files[i, ]
    path <- file.path(output_root, normalise_path(row$path))
    present <- file.exists(path)
    observed_bytes <- if (present) file.info(path)$size else NA_real_
    size_match <- present && as.numeric(observed_bytes) == as.numeric(row$expected_bytes)
    observed_sha256 <- if (present && verify_sha256) {
      tolower(digest::digest(path, algo = "sha256", serialize = FALSE, file = TRUE))
    } else NA_character_
    sha256_match <- if (present && verify_sha256) {
      identical(observed_sha256, tolower(row$expected_sha256))
    } else NA
    checks[[i]] <- data.frame(
      file = path, present = present, expected_bytes = row$expected_bytes,
      observed_bytes = observed_bytes, size_match = size_match,
      expected_sha256 = row$expected_sha256, observed_sha256 = observed_sha256,
      sha256_match = sha256_match, stringsAsFactors = FALSE
    )
  }
  checks <- do.call(rbind, checks)
  log_dir <- file.path(output_root, "04_Download_Logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(do.call(rbind, log_rows), file.path(log_dir, "GSE109875_download_log.csv"), row.names = FALSE, na = "")
  utils::write.csv(checks, file.path(log_dir, "GSE109875_file_check.csv"), row.names = FALSE, na = "")
  problems <- !checks$present | !checks$size_match |
    (verify_sha256 & !is.na(checks$sha256_match) & !checks$sha256_match)
  if (any(problems)) {
    stop(sum(problems), " expected file(s) were missing or did not match. See the file-check CSV.",
         call. = FALSE)
  }
  message("Dataset staged and checked successfully.")
}

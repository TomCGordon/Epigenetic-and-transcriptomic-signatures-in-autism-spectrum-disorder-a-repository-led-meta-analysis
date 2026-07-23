#!/usr/bin/env Rscript

# Download and stage the public repository inputs used for GSE28475.
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
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE28nnn/GSE28475/matrix/GSE28475-GPL13388_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE28nnn/GSE28475/matrix/GSE28475-GPL6883_series_matrix.txt.gz"
  ),
  destination = c(
    "01_Raw_Data/Gene_expression/Post_mortem_brain/GSE28475/GSE28475-GPL13388_series_matrix.txt.gz",
    "01_Raw_Data/Gene_expression/Post_mortem_brain/GSE28475/GSE28475-GPL6883_series_matrix.txt.gz"
  ),
  action = c(
    "none",
    "none"
  ),
  extraction_root = c(
    NA_character_,
    NA_character_
  ),
  expected_bytes = c(
    1087282,
    12665697
  ),
  expected_sha256 = c(
    "1cfd3aad51bff5455582242d51d047d65529d8b156f39857a635305703b44be8",
    "471e0106c39229031754549786a3266bd60f0cadc0e429301caa9bccabb60bde"
  ),
  stringsAsFactors = FALSE
)

expected_files <- data.frame(
  path = c(
    "01_Raw_Data/Gene_expression/Post_mortem_brain/GSE28475/GSE28475-GPL13388_series_matrix.txt.gz",
    "01_Raw_Data/Gene_expression/Post_mortem_brain/GSE28475/GSE28475-GPL6883_series_matrix.txt.gz"
  ),
  expected_bytes = c(
    1087282,
    12665697
  ),
  expected_sha256 = c(
    "1cfd3aad51bff5455582242d51d047d65529d8b156f39857a635305703b44be8",
    "471e0106c39229031754549786a3266bd60f0cadc0e429301caa9bccabb60bde"
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
    accession = "GSE28475",
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
  utils::write.csv(do.call(rbind, log_rows), file.path(log_dir, "GSE28475_download_log.csv"), row.names = FALSE, na = "")
  utils::write.csv(checks, file.path(log_dir, "GSE28475_file_check.csv"), row.names = FALSE, na = "")
  problems <- !checks$present | !checks$size_match |
    (verify_sha256 & !is.na(checks$sha256_match) & !checks$sha256_match)
  if (any(problems)) {
    stop(sum(problems), " expected file(s) were missing or did not match. See the file-check CSV.",
         call. = FALSE)
  }
  message("Dataset staged and checked successfully.")
}

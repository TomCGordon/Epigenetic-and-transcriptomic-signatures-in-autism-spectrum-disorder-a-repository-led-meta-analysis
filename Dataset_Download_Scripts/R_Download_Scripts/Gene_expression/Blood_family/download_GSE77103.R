#!/usr/bin/env Rscript

# Download and stage the public repository inputs used for GSE77103.
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
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77103/suppl/GSE77103_RAW.tar"
  ),
  destination = c(
    "01_Raw_Data/Gene_expression/Blood_family/GSE77103/GSE77103_RAW.tar"
  ),
  action = c(
    "extract_archive"
  ),
  extraction_root = c(
    "01_Raw_Data/Gene_expression/Blood_family/GSE77103/GSE77103_RAW"
  ),
  expected_bytes = c(
    25323520
  ),
  expected_sha256 = c(
    "f6d148b348bcc80cdb8991778106c395ae1563e700d5032cf2576f9d6a8a10e2"
  ),
  stringsAsFactors = FALSE
)

expected_files <- data.frame(
  path = c(
    "01_Raw_Data/Gene_expression/Blood_family/GSE77103/GSE77103_RAW.tar",
    "01_Raw_Data/Gene_expression/Blood_family/GSE77103/GSE77103_RAW/GSM2044363_C2.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE77103/GSE77103_RAW/GSM2044364_C4.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE77103/GSE77103_RAW/GSM2044365_C5.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE77103/GSE77103_RAW/GSM2044366_C6.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE77103/GSE77103_RAW/GSM2044367_A2.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE77103/GSE77103_RAW/GSM2044368_A3.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE77103/GSE77103_RAW/GSM2044369_A4.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE77103/GSE77103_RAW/GSM2044370_A5.txt.gz"
  ),
  expected_bytes = c(
    25323520,
    3194167,
    3190286,
    3177295,
    3167885,
    3142461,
    3133875,
    3154385,
    3147378
  ),
  expected_sha256 = c(
    "f6d148b348bcc80cdb8991778106c395ae1563e700d5032cf2576f9d6a8a10e2",
    "9dda3d40b82dc9e6ae06cc2b89c1355bc93bb3d8fb79b576ced6cf1cb20f0519",
    "5135754a7945252e87883da7875aa555cb40fcfaf0eead82078f166cdfb4dce7",
    "d63daa354ed8c6aa95e267b92859760f3ea7d0a42d9efccc0bf9d3ff486aa240",
    "17351fb88706445a672d00d19db1570d5852632153f2d0ba9cd76745a78c3089",
    "c574e1f72d107a171c458eaa0749e2d1447d85dd21f27f81562590530d98fee3",
    "e0ca8536ec8e8810f2417d87137c8b36ca1eedb3cfc191be1c39246dacbbfa25",
    "9649a8042fe10b3a5f14c15332760bcd72039a4db66a286df4f45d169506662f",
    "daecc71df818ef0aacf7a24ecebeec0125e16990a1ded26e11fc201c3c372b23"
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
    accession = "GSE77103",
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
  utils::write.csv(do.call(rbind, log_rows), file.path(log_dir, "GSE77103_download_log.csv"), row.names = FALSE, na = "")
  utils::write.csv(checks, file.path(log_dir, "GSE77103_file_check.csv"), row.names = FALSE, na = "")
  problems <- !checks$present | !checks$size_match |
    (verify_sha256 & !is.na(checks$sha256_match) & !checks$sha256_match)
  if (any(problems)) {
    stop(sum(problems), " expected file(s) were missing or did not match. See the file-check CSV.",
         call. = FALSE)
  }
  message("Dataset staged and checked successfully.")
}

#!/usr/bin/env Rscript

# Download and stage the public repository inputs used for GSE81541.
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
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156966/suppl/GSM2156966_JLKD002.bed.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156967/suppl/GSM2156967_JLKD003.bed.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156968/suppl/GSM2156968_JLKD005.bed.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156969/suppl/GSM2156969_JLKD001.bed.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156970/suppl/GSM2156970_JLKD014.bed.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156971/suppl/GSM2156971_JLKD004.bed.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156972/suppl/GSM2156972_JLKD040.bed.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156973/suppl/GSM2156973_JLKD041.bed.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156974/suppl/GSM2156974_JLKD042.bed.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156975/suppl/GSM2156975_JLKD026.bed.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156976/suppl/GSM2156976_JLKD028.bed.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156998/suppl/GSM2156998_JLKD009.bed.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156999/suppl/GSM2156999_JLKD013.bed.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2157nnn/GSM2157000/suppl/GSM2157000_JLKD016.bed.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2157nnn/GSM2157001/suppl/GSM2157001_JLKD018.bed.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2157nnn/GSM2157002/suppl/GSM2157002_JLKD019.bed.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2157nnn/GSM2157003/suppl/GSM2157003_JLKD020.bed.gz"
  ),
  destination = c(
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156966.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156967.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156968.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156969.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156970.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156971.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156972.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156973.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156974.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156975.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156976.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156998.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156999.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2157000.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2157001.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2157002.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2157003.bed.gz"
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
    NA_character_,
    NA_character_
  ),
  expected_bytes = c(
    128901381,
    130362321,
    143713191,
    127132779,
    167672527,
    136040105,
    182705486,
    182715004,
    184013401,
    158416633,
    144597767,
    175755501,
    172212170,
    172110276,
    180508054,
    175676620,
    177292075
  ),
  expected_sha256 = c(
    "1622c8d12db3f92f25b7fb1ee458d455bbf2de77b30b5076852b7165ebd1eff1",
    "894abfccb41c910448fcd169a776db4567c76e712512299bff940ba974c5be08",
    "9e457a35139575d12b58f8208e603adb42581a4367d442cbcb9c5e8dc49fd9d8",
    "1a5877598d07e3900ff27254c858ae5530255d659d1d178919c02874474637c3",
    "a058eed930cdfeda2a0d565759225005f07396c2b2646160bebb6b1fcaea2fc9",
    "491c47460e7b668e2917b8616da994e9133533ccae599571e897af934d0f5894",
    "e4b13aa2b3088d870346f519693d368422a83a5fad1266d5785716d0296d251e",
    "a51b31c5d7a34bf546f5d4af66c19b5cd38ed6d7cbb4ee1e9361623c2d73789e",
    "bc62e3fc610ee5dab8f623328000d844f8cfe9ffd60a058a410adae282cb2afd",
    "cbf02086fdfd946e1de3aa36ce36716d8863a4194bccd5c0a727d4ab6485aaf1",
    "b13f62fad2f4d9cb757c52661286a71ffa08b50c7c4581f9aa4e922fdbfd1f10",
    "1c6935cdd3caa6681139ea0dc81f0ba5e3057ba76dc92667843eee3c6362bd35",
    "04fc3ca24b3d6c4aa0c95d3a3accacffa0801ee3da802af2b971c2329513ff43",
    "b7237f7abb9c985337e79a873ff9b85702b26ea1ce3c6dbd738d3d572bfe35cc",
    "41f76a0daf88c9c320725e26d64b0cac4907def21d08bb062b078b6f605e2991",
    "899bd1a29f040fe9a3ea6c3c640ec0c182152bcc352c0cd475f9f162bfa53dce",
    "ff5c40e09e2c1758e10c301139618ace40cfcc1cd0a537f8c6bf7bfb04920e68"
  ),
  stringsAsFactors = FALSE
)

expected_files <- data.frame(
  path = c(
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156966.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156967.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156968.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156969.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156970.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156971.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156972.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156973.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156974.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156975.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156976.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156998.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156999.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2157000.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2157001.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2157002.bed.gz",
    "01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2157003.bed.gz"
  ),
  expected_bytes = c(
    128901381,
    130362321,
    143713191,
    127132779,
    167672527,
    136040105,
    182705486,
    182715004,
    184013401,
    158416633,
    144597767,
    175755501,
    172212170,
    172110276,
    180508054,
    175676620,
    177292075
  ),
  expected_sha256 = c(
    "1622c8d12db3f92f25b7fb1ee458d455bbf2de77b30b5076852b7165ebd1eff1",
    "894abfccb41c910448fcd169a776db4567c76e712512299bff940ba974c5be08",
    "9e457a35139575d12b58f8208e603adb42581a4367d442cbcb9c5e8dc49fd9d8",
    "1a5877598d07e3900ff27254c858ae5530255d659d1d178919c02874474637c3",
    "a058eed930cdfeda2a0d565759225005f07396c2b2646160bebb6b1fcaea2fc9",
    "491c47460e7b668e2917b8616da994e9133533ccae599571e897af934d0f5894",
    "e4b13aa2b3088d870346f519693d368422a83a5fad1266d5785716d0296d251e",
    "a51b31c5d7a34bf546f5d4af66c19b5cd38ed6d7cbb4ee1e9361623c2d73789e",
    "bc62e3fc610ee5dab8f623328000d844f8cfe9ffd60a058a410adae282cb2afd",
    "cbf02086fdfd946e1de3aa36ce36716d8863a4194bccd5c0a727d4ab6485aaf1",
    "b13f62fad2f4d9cb757c52661286a71ffa08b50c7c4581f9aa4e922fdbfd1f10",
    "1c6935cdd3caa6681139ea0dc81f0ba5e3057ba76dc92667843eee3c6362bd35",
    "04fc3ca24b3d6c4aa0c95d3a3accacffa0801ee3da802af2b971c2329513ff43",
    "b7237f7abb9c985337e79a873ff9b85702b26ea1ce3c6dbd738d3d572bfe35cc",
    "41f76a0daf88c9c320725e26d64b0cac4907def21d08bb062b078b6f605e2991",
    "899bd1a29f040fe9a3ea6c3c640ec0c182152bcc352c0cd475f9f162bfa53dce",
    "ff5c40e09e2c1758e10c301139618ace40cfcc1cd0a537f8c6bf7bfb04920e68"
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
    accession = "GSE81541",
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
  utils::write.csv(do.call(rbind, log_rows), file.path(log_dir, "GSE81541_download_log.csv"), row.names = FALSE, na = "")
  utils::write.csv(checks, file.path(log_dir, "GSE81541_file_check.csv"), row.names = FALSE, na = "")
  problems <- !checks$present | !checks$size_match |
    (verify_sha256 & !is.na(checks$sha256_match) & !checks$sha256_match)
  if (any(problems)) {
    stop(sum(problems), " expected file(s) were missing or did not match. See the file-check CSV.",
         call. = FALSE)
  }
  message("Dataset staged and checked successfully.")
}

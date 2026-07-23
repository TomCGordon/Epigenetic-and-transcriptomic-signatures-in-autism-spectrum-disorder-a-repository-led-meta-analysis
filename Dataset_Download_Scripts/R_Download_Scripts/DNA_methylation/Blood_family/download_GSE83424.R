#!/usr/bin/env Rscript

# Download and stage the public repository inputs used for GSE83424.
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
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE83nnn/GSE83424/miniml/GSE83424_family.xml.tgz"
  ),
  destination = c(
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSE83424_family.xml.tgz"
  ),
  action = c(
    "extract_miniml_and_copy_tables"
  ),
  extraction_root = c(
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424"
  ),
  expected_bytes = c(
    310543470
  ),
  expected_sha256 = c(
    "39c86dd31cd7ad479fecd0a689c6bd8c79f1302e21870f8769501c970cdd097d"
  ),
  stringsAsFactors = FALSE
)

expected_files <- data.frame(
  path = c(
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GPL16304-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSE83424_family.xml",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSE83424_family.xml.tgz",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202710-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202711-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202712-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202713-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202714-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202715-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202716-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202717-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202718-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202719-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202720-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202721-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202722-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202723-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202724-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202725-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202726-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202727-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202728-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202729-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202730-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202731-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202732-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202733-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202734-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202735-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202736-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202737-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202738-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202739-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202740-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202741-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202742-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202743-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202744-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202745-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202746-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202747-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202748-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202749-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202750-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202751-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202752-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202753-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202754-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202755-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202756-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202757-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202758-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202759-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202760-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202761-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202762-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202763-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202764-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202765-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202766-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202767-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202768-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202769-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202770-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202771-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/GSM2202772-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202710-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202711-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202712-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202713-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202714-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202715-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202716-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202717-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202718-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202719-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202720-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202721-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202722-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202723-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202724-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202725-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202726-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202727-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202728-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202729-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202730-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202731-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202732-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202733-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202734-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202735-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202736-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202737-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202738-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202739-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202740-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202741-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202742-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202743-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202744-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202745-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202746-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202747-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202748-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202749-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202750-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202751-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202752-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202753-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202754-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202755-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202756-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202757-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202758-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202759-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202760-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202761-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202762-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202763-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202764-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202765-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202766-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202767-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202768-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202769-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202770-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202771-tbl-1.txt",
    "01_Raw_Data/DNA_methylation/Blood_family/GSE83424/GSE83424/sample_tables/GSM2202772-tbl-1.txt"
  ),
  expected_bytes = c(
    76942796,
    166014,
    310543470,
    17004570,
    17004567,
    17490082,
    17004577,
    17004566,
    17490108,
    17490075,
    17490077,
    17004563,
    17004568,
    17004570,
    17490090,
    17004566,
    17490086,
    17490115,
    17490076,
    17490101,
    17004573,
    17490076,
    17490089,
    17004583,
    17490090,
    17004567,
    17490105,
    17490084,
    17004565,
    17490091,
    17004580,
    17490097,
    17490115,
    17490089,
    17490095,
    17490093,
    17004569,
    17490088,
    17490097,
    17490086,
    17490092,
    17004571,
    17490081,
    17004560,
    17490079,
    17004569,
    17004570,
    17490083,
    17490095,
    17490097,
    17004576,
    17490086,
    17004573,
    17490086,
    17004574,
    17004569,
    17490091,
    17490089,
    17490081,
    17004564,
    17004565,
    17004568,
    17490104,
    17490102,
    17490110,
    17490105,
    17004570,
    17004567,
    17490082,
    17004577,
    17004566,
    17490108,
    17490075,
    17490077,
    17004563,
    17004568,
    17004570,
    17490090,
    17004566,
    17490086,
    17490115,
    17490076,
    17490101,
    17004573,
    17490076,
    17490089,
    17004583,
    17490090,
    17004567,
    17490105,
    17490084,
    17004565,
    17490091,
    17004580,
    17490097,
    17490115,
    17490089,
    17490095,
    17490093,
    17004569,
    17490088,
    17490097,
    17490086,
    17490092,
    17004571,
    17490081,
    17004560,
    17490079,
    17004569,
    17004570,
    17490083,
    17490095,
    17490097,
    17004576,
    17490086,
    17004573,
    17490086,
    17004574,
    17004569,
    17490091,
    17490089,
    17490081,
    17004564,
    17004565,
    17004568,
    17490104,
    17490102,
    17490110,
    17490105
  ),
  expected_sha256 = c(
    "8b560e1bc2524a96b10104c98ab481240c49285f7af9d71a46a06cda282a9539",
    "e69266ac2e90f1ae1a818d5e69f4f7378e0902c9f6ccab0464149519d1d4142c",
    "39c86dd31cd7ad479fecd0a689c6bd8c79f1302e21870f8769501c970cdd097d",
    "95c547464bce3799ca0db5cd374047e9ff6a14d12fd62336ff90e8eaee6d932d",
    "c2a722ffceae52db5c1ce6193bed391fc12bfb449aa7ea819e16850b9b8e1716",
    "4afaf289d93f6a23673458a00b31852ba220fe45799ba3e798491fd194898409",
    "e12ba817b483ae49c09d78c3ab4551912bcaf21675c2aef2f3b4e26ee655072f",
    "48c36c91dbae92e1b6b35546821da5316d3ea1485aa383901e26549f967c774e",
    "0b655b77f12eb54a7072f30ec930013c7afaa350bb4e9d8cb1b574d999fa4e65",
    "f7b29d1d15db36b3c6d4b6251c95a427c4861e75f8a666d17eea4fbeedfee69e",
    "e4960736ae53ae602580a8a72b60f7bf444dbbe340a6b22083bce77596a25af7",
    "25fa397acb40e333e42bcb98fe873e730b3c4412fe344ab7f5fb632b71cebd9d",
    "30de069742b997c352a68fd889e3112410536486e61c743fad35df2e5017a687",
    "9d669043d05db95f9b86ae75962511c4fd1c2b74b497c2dd6d6999280619a8e2",
    "7cd11f2fc01cf31c3bf13ea2e16a3aec878cd87b42f86b29074980599c86ef0b",
    "c8474790f7a79a138401e08cf9413231368b314ee99151eca4abf285dfddcd0b",
    "1ef0a1ee00a88a082a4306733b732ba99d2f0cea4f8aca8e5633106470047146",
    "bcba248ae2a27f2a4b27e9ba69b0d04bb0169e685dd084e697e75578b33455c0",
    "fa87260888eebbc9c6c3990e9822588a59fd1d67d3dba0bc93faadb237c0dd8e",
    "d30a2949f128f8b0fca2503e9b44e772c8aa16be276602496fd12af507d4d3d6",
    "8369b217825003ac0d5bfb72a610112ba74f9c9442af236a7f504145c025b83f",
    "a1a6fdddbbf90dd0b3e10adae624e6e0ff7d7fba73827ec4757297d3c42dd782",
    "d3f87fb680b39a69cfd1380f8c9330d022d11db01be44beca81fd6aa111da33a",
    "ffd1ff97a374b8e216daef0731657c49f167b0698365b99f5762ac5ec0c1724a",
    "57abf9858c882392e6d7a5af48846351d522fb3d142cb89e2ab4cc1a17a9fe8b",
    "b6e74b54769e47b402a812bc5b14cafd04a04f99771bc6671e413ddf23bf7e7e",
    "0f6b33299787dc62121e82616f356db536135fd859d62e70ee68dc9be807ce21",
    "a6c2cae07df87e554fb6186c1215af9b8770ef28af442162d5ce95e3ebfbf22f",
    "ec42ee7245bbb8bfd73cd70a76a222458bc6760eac541c7305213e80f1553f66",
    "8c0c08210845f4c948d4e1ed9d9602e5c6f08d7a1e0f9f02e195cd0bfd395e6d",
    "7c49e9654d057872c29fa7f7a14cf6c5539476580910c51da956cae90b644a4f",
    "33db9c6da225455606555286b1726d192038bd6b12e88a3110ced102472a0fec",
    "0045772c0cae600703460bf766f7f99a5a2551d20569da9ebb9f375b903d0dcf",
    "b7a30d256d70d609c9e148e6e03fda24ea787edce81a8b61f3be473f26791838",
    "219408b5e6da3c185c07b0fdb91b892cf94290894532df60cfbfe1ae7b6535bc",
    "9cb983df90a9cebf6625a176e9cfd577d14481060872341b89fe8047f833de5f",
    "f10764d41dde925e9728783a24abd5ffbb7cd0a9dc6cc460a094d85a14c4b79f",
    "5f097cf87c39cf2f029cf5b84d754e784a82303cb5cc6361791e460aba6a17ca",
    "e933279e0e3c39f7d9b39578912c7b2b4c79f07b6949ac4f4cf97352936eb114",
    "24af260d845e0d2942ca899469d0f98e4d455e81ae1d2598873e3bc9f539cd0f",
    "b9b2a73404c877d5523051bb36457d355661dff5aab55778cb5d166097224b3e",
    "43e049b7a014e693e4560e7ea40b0aeab86cae78046744674d93deff591967eb",
    "7d4a061a4db0924045ca995047674091568c80706e2b4966a8d79157b3cd3c59",
    "9bfb69f199bb04daa66ce84bab7b4a9bddd2c93b233d44072665ae080069c0c6",
    "1efe10624cd99c321426615fd19d3a91e56d2b1f209206ef65304c340435e402",
    "dd8947ae1166af42a5da30e892835266fae8d1192eeb152290781dfbbfd07ce6",
    "73ed3dacd6fdb1e53b7a3289e9b2e08da54477e850e8a7fa396d39c4e33781eb",
    "4c05f6e3871e169439eeab961070d2901b3a13bdee240d363cf1165f847072c4",
    "e282ec79ef70063117d778f1f0ab159e6f25fd1448d515312b15d426e2a7e7e9",
    "e0aa5aace0d6e37ea2953594b5582f8426c857d4f44d8d0303c10eec56dcda92",
    "fc9256d4f806ac2aa79a89745b8b08e91703d66c33586aa263e4794807d66aa9",
    "985484e65f4b90a1755172eb85f61a0c24858bc52a9eb81919658593b364d5db",
    "e91da86471d9ebfef94d3802a3bb1647fb6bc944b6eb6e5d4d462505b0a1ec45",
    "2aaee29df163c6de66b6efc08d74ac238445cce9958b62b78de6b460f0439dce",
    "04a66f5ded22ef372e2dccb724e98cfaa80f6b8d6eb7e9485f5afa0857564cfd",
    "fd45316dd4a3814e9edc6f902c3fc4960b426e3e70f7b78197c1a13553692a5f",
    "28749a9c9774e22c26da2642634a0ad0aa7f9d2b99f50892c5d32d0f4300351c",
    "e1305de8839d7f97758e093573fa3620f6b22f631d6bb4869cb01d82e6124c97",
    "8aecffa194330fd82b23b9488ae2e24933f20f23edcf93fee9d82d8c4ca890cf",
    "7a83783c1ce73bfca26b822e2fa2f29403e2b39cd9c91b4bfc5e4ec2fd29329d",
    "57ea1e64fd0d06378016796193312ad5fc545759b331d3e8012d0c9d9b5ea2cf",
    "f1c4ac04a5edb8988bcd39342f2d73ea3e99e1a8f63457cf0abbb319cc97a1d6",
    "290ea4497654865c466247ccd597e54b2fc9ce2c5d2172550c13b4b4e5a7b158",
    "64ade71aad53a54a8913af07b71878bacfe8d74046723fb2f737e0b59f592682",
    "1af8511864137cdf0da0ae1fbcb828fbbcd8032aad48fdb8e28efc10fb3f2c30",
    "60b571cef85fb1959bc5b7939bc85e51ab74770f418b59d3e05f00b829472693",
    "95c547464bce3799ca0db5cd374047e9ff6a14d12fd62336ff90e8eaee6d932d",
    "c2a722ffceae52db5c1ce6193bed391fc12bfb449aa7ea819e16850b9b8e1716",
    "4afaf289d93f6a23673458a00b31852ba220fe45799ba3e798491fd194898409",
    "e12ba817b483ae49c09d78c3ab4551912bcaf21675c2aef2f3b4e26ee655072f",
    "48c36c91dbae92e1b6b35546821da5316d3ea1485aa383901e26549f967c774e",
    "0b655b77f12eb54a7072f30ec930013c7afaa350bb4e9d8cb1b574d999fa4e65",
    "f7b29d1d15db36b3c6d4b6251c95a427c4861e75f8a666d17eea4fbeedfee69e",
    "e4960736ae53ae602580a8a72b60f7bf444dbbe340a6b22083bce77596a25af7",
    "25fa397acb40e333e42bcb98fe873e730b3c4412fe344ab7f5fb632b71cebd9d",
    "30de069742b997c352a68fd889e3112410536486e61c743fad35df2e5017a687",
    "9d669043d05db95f9b86ae75962511c4fd1c2b74b497c2dd6d6999280619a8e2",
    "7cd11f2fc01cf31c3bf13ea2e16a3aec878cd87b42f86b29074980599c86ef0b",
    "c8474790f7a79a138401e08cf9413231368b314ee99151eca4abf285dfddcd0b",
    "1ef0a1ee00a88a082a4306733b732ba99d2f0cea4f8aca8e5633106470047146",
    "bcba248ae2a27f2a4b27e9ba69b0d04bb0169e685dd084e697e75578b33455c0",
    "fa87260888eebbc9c6c3990e9822588a59fd1d67d3dba0bc93faadb237c0dd8e",
    "d30a2949f128f8b0fca2503e9b44e772c8aa16be276602496fd12af507d4d3d6",
    "8369b217825003ac0d5bfb72a610112ba74f9c9442af236a7f504145c025b83f",
    "a1a6fdddbbf90dd0b3e10adae624e6e0ff7d7fba73827ec4757297d3c42dd782",
    "d3f87fb680b39a69cfd1380f8c9330d022d11db01be44beca81fd6aa111da33a",
    "ffd1ff97a374b8e216daef0731657c49f167b0698365b99f5762ac5ec0c1724a",
    "57abf9858c882392e6d7a5af48846351d522fb3d142cb89e2ab4cc1a17a9fe8b",
    "b6e74b54769e47b402a812bc5b14cafd04a04f99771bc6671e413ddf23bf7e7e",
    "0f6b33299787dc62121e82616f356db536135fd859d62e70ee68dc9be807ce21",
    "a6c2cae07df87e554fb6186c1215af9b8770ef28af442162d5ce95e3ebfbf22f",
    "ec42ee7245bbb8bfd73cd70a76a222458bc6760eac541c7305213e80f1553f66",
    "8c0c08210845f4c948d4e1ed9d9602e5c6f08d7a1e0f9f02e195cd0bfd395e6d",
    "7c49e9654d057872c29fa7f7a14cf6c5539476580910c51da956cae90b644a4f",
    "33db9c6da225455606555286b1726d192038bd6b12e88a3110ced102472a0fec",
    "0045772c0cae600703460bf766f7f99a5a2551d20569da9ebb9f375b903d0dcf",
    "b7a30d256d70d609c9e148e6e03fda24ea787edce81a8b61f3be473f26791838",
    "219408b5e6da3c185c07b0fdb91b892cf94290894532df60cfbfe1ae7b6535bc",
    "9cb983df90a9cebf6625a176e9cfd577d14481060872341b89fe8047f833de5f",
    "f10764d41dde925e9728783a24abd5ffbb7cd0a9dc6cc460a094d85a14c4b79f",
    "5f097cf87c39cf2f029cf5b84d754e784a82303cb5cc6361791e460aba6a17ca",
    "e933279e0e3c39f7d9b39578912c7b2b4c79f07b6949ac4f4cf97352936eb114",
    "24af260d845e0d2942ca899469d0f98e4d455e81ae1d2598873e3bc9f539cd0f",
    "b9b2a73404c877d5523051bb36457d355661dff5aab55778cb5d166097224b3e",
    "43e049b7a014e693e4560e7ea40b0aeab86cae78046744674d93deff591967eb",
    "7d4a061a4db0924045ca995047674091568c80706e2b4966a8d79157b3cd3c59",
    "9bfb69f199bb04daa66ce84bab7b4a9bddd2c93b233d44072665ae080069c0c6",
    "1efe10624cd99c321426615fd19d3a91e56d2b1f209206ef65304c340435e402",
    "dd8947ae1166af42a5da30e892835266fae8d1192eeb152290781dfbbfd07ce6",
    "73ed3dacd6fdb1e53b7a3289e9b2e08da54477e850e8a7fa396d39c4e33781eb",
    "4c05f6e3871e169439eeab961070d2901b3a13bdee240d363cf1165f847072c4",
    "e282ec79ef70063117d778f1f0ab159e6f25fd1448d515312b15d426e2a7e7e9",
    "e0aa5aace0d6e37ea2953594b5582f8426c857d4f44d8d0303c10eec56dcda92",
    "fc9256d4f806ac2aa79a89745b8b08e91703d66c33586aa263e4794807d66aa9",
    "985484e65f4b90a1755172eb85f61a0c24858bc52a9eb81919658593b364d5db",
    "e91da86471d9ebfef94d3802a3bb1647fb6bc944b6eb6e5d4d462505b0a1ec45",
    "2aaee29df163c6de66b6efc08d74ac238445cce9958b62b78de6b460f0439dce",
    "04a66f5ded22ef372e2dccb724e98cfaa80f6b8d6eb7e9485f5afa0857564cfd",
    "fd45316dd4a3814e9edc6f902c3fc4960b426e3e70f7b78197c1a13553692a5f",
    "28749a9c9774e22c26da2642634a0ad0aa7f9d2b99f50892c5d32d0f4300351c",
    "e1305de8839d7f97758e093573fa3620f6b22f631d6bb4869cb01d82e6124c97",
    "8aecffa194330fd82b23b9488ae2e24933f20f23edcf93fee9d82d8c4ca890cf",
    "7a83783c1ce73bfca26b822e2fa2f29403e2b39cd9c91b4bfc5e4ec2fd29329d",
    "57ea1e64fd0d06378016796193312ad5fc545759b331d3e8012d0c9d9b5ea2cf",
    "f1c4ac04a5edb8988bcd39342f2d73ea3e99e1a8f63457cf0abbb319cc97a1d6",
    "290ea4497654865c466247ccd597e54b2fc9ce2c5d2172550c13b4b4e5a7b158",
    "64ade71aad53a54a8913af07b71878bacfe8d74046723fb2f737e0b59f592682",
    "1af8511864137cdf0da0ae1fbcb828fbbcd8032aad48fdb8e28efc10fb3f2c30",
    "60b571cef85fb1959bc5b7939bc85e51ab74770f418b59d3e05f00b829472693"
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
    accession = "GSE83424",
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
  utils::write.csv(do.call(rbind, log_rows), file.path(log_dir, "GSE83424_download_log.csv"), row.names = FALSE, na = "")
  utils::write.csv(checks, file.path(log_dir, "GSE83424_file_check.csv"), row.names = FALSE, na = "")
  problems <- !checks$present | !checks$size_match |
    (verify_sha256 & !is.na(checks$sha256_match) & !checks$sha256_match)
  if (any(problems)) {
    stop(sum(problems), " expected file(s) were missing or did not match. See the file-check CSV.",
         call. = FALSE)
  }
  message("Dataset staged and checked successfully.")
}

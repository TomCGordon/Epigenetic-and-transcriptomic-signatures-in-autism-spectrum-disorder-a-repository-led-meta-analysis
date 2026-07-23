#!/usr/bin/env Rscript

# Download and stage the public repository inputs used for GSE140702.
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
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140702/suppl/GSE140702_RAW.tar"
  ),
  destination = c(
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW.tar"
  ),
  action = c(
    "extract_archive"
  ),
  extraction_root = c(
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW"
  ),
  expected_bytes = c(
    624189440
  ),
  expected_sha256 = c(
    "71cba66206b739d9d6182a4ba7a0833f9cd02c718fd94894e68c0f0939c3c0b1"
  ),
  stringsAsFactors = FALSE
)

expected_files <- data.frame(
  path = c(
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW.tar",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182148_101119LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182149_101119LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182150_101119NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182151_101123LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182152_101123LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182153_101123NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182154_101132LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182155_101132LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182156_101132NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182157_101134LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182158_101134LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182159_101134NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182160_101138LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182161_101138LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182162_101138NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182163_101150LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182164_101150LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182165_101150NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182166_101152LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182167_101152LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182168_101152NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182169_101157LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182170_101157LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182171_101157NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182172_101178LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182173_101178LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182174_101178NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182175_101213LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182176_101213LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182177_101213NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182178_101221LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182179_101221LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182180_101221NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182181_101241LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182182_101241NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182183_101254LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182184_101254LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182185_101254NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182186_101280LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182187_101280LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182188_101280NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182189_101288LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182190_101288LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182191_101288NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182192_101313LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182193_101313LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182194_101313NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182195_101660LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182196_101660NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182197_101668LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182198_101668LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182199_101668NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182200_101676LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182201_101676LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182202_101676NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182203_101679LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182204_101679LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182205_101679NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182206_101680LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182207_101680LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182208_101680NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182209_101682LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182210_101682LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182211_101682NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182212_101684NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182213_101685LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182214_101685NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182215_101687LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182216_101687LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182217_101687NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182218_101688LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182219_101688LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182220_101688NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182221_101696LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182222_101696NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182223_101697LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182224_101697LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182225_101697NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182226_101699LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182227_101699NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182228_101703LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182229_101703LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182230_101703NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182231_101709LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182232_101709LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182233_101709NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182234_101713LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182235_101713LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182236_101713NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182237_101717LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182238_101717NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182239_101724LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182240_101724NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182241_101725LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182242_101725LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182243_101725NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182244_101728LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182245_101728LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182246_101728NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182247_101730LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182248_101730LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182249_101730NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182250_101736LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182251_101736NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182252_101739LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182253_101739NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182254_101746LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182255_101746LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM4182256_101746NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM5531610_101246LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM5531611_101246LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM5531612_101246NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM5531613_101334LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM5531614_101334LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM5531615_101334NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM5531616_101666LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM5531617_101666NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM5531618_101692LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM5531619_101692NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM5531620_101694LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM5531621_101694LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM5531622_101694NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM5531623_101706LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM5531624_101706LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM5531625_101706NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM5531626_101711LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM5531627_101711LTA.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM5531628_101711NT.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM5531629_101734LPS.counts.txt.gz",
    "01_Raw_Data/Gene_expression/Blood_family/GSE140702/GSE140702_RAW/GSM5531630_101734NT.counts.txt.gz"
  ),
  expected_bytes = c(
    624189440,
    4828149,
    4668648,
    4826714,
    4827157,
    4828882,
    4830587,
    4826009,
    4828872,
    4827659,
    4827796,
    4826509,
    4830200,
    4826386,
    4827277,
    4830031,
    4826243,
    4827646,
    4829729,
    4827026,
    4831065,
    4831047,
    4826701,
    4826655,
    4830454,
    4825805,
    4827378,
    4829564,
    4825711,
    4826594,
    4828031,
    4825304,
    4825529,
    4827300,
    4828258,
    4829627,
    4825881,
    4826855,
    4829021,
    4823675,
    4824541,
    4825459,
    4824715,
    4826507,
    4828068,
    4827984,
    4825965,
    4830955,
    4826662,
    4829738,
    4826299,
    4829097,
    4831037,
    4827405,
    4830584,
    4829961,
    4825049,
    4828436,
    4827577,
    4826192,
    4827588,
    4831758,
    4826916,
    4828258,
    4828967,
    4827646,
    4825037,
    4827093,
    4825527,
    4828644,
    4829187,
    4826331,
    4825474,
    4827286,
    4826605,
    4831640,
    4825955,
    4826863,
    4830567,
    4827396,
    4829166,
    4825463,
    4824869,
    4829003,
    4826801,
    4827776,
    4831384,
    4825295,
    4828972,
    4830408,
    4824372,
    4829524,
    4829158,
    4828179,
    4827126,
    4826545,
    4830085,
    4825078,
    4825512,
    4829912,
    4825039,
    4827070,
    4830710,
    4828192,
    4830416,
    4825361,
    4829420,
    4824479,
    4826512,
    4827801,
    4668863,
    4667175,
    4669339,
    4666900,
    4668645,
    4669359,
    4668557,
    4669953,
    4666800,
    4670237,
    4666808,
    4667892,
    4668678,
    4668329,
    4669269,
    4670979,
    4667591,
    4667336,
    4668588,
    4664540,
    4668416
  ),
  expected_sha256 = c(
    "71cba66206b739d9d6182a4ba7a0833f9cd02c718fd94894e68c0f0939c3c0b1",
    "bdac41c5d954231fa21347fbcae967311d1bbe5610d4702421e5e851f70bd953",
    "1b8086e3193403a2878ac7193ba9c6704e0ef7b54a5ab16676816c561d6b6813",
    "c6f3525938032f105241579675496d4dbb2a442762c736f173033cdb62c691c1",
    "5bc3d83ebf59661bcaa181e449456dd3bafaee0cfe667c74367a430f4e94a39d",
    "3203941c4d913b4adfb8ad6129d7c6b63d9944289236c65f2524cced33a963d1",
    "a246631ded0d6618f556d16859650ab82ac3b9b66c93304fa445587a8afb0171",
    "79055e940fd7c166a115fe9fa1adf02af266266733427ea1d971e0241c694058",
    "6407b4966bdcf4401f9ad7455b35a66e7cd5bf456cb9fac91445bb1c8e5d240f",
    "37d254e75fa5c0ade0c3c43dbb6688f263eb536e95a099075073a195e4245dc0",
    "da13969067891c8eb4579dc0ef572c24662cb433eb148482e3c033a2dd5980a1",
    "ad8705666935c965e0715af3da1d58e284b8061e82b43d8d2406e1e76884b650",
    "18885541375f61bf8af06da82d28026c6cbd76b567890985a506878df5399b52",
    "1db66c31c71c5cad1edde3b8ef9cbd337e12cb33e34e0239bbd57965da450dbc",
    "1c9f356200e9413d15b40d983c1f483441b8a0592cb1e416cd24a562e2e9de19",
    "d5f8e2bbb219771b9d1f48d8f36966c01ffa110eeb63dd5f35f6a6929aeb3d68",
    "bd5fa5db491ad344fce95639cea2fdcea1095659cf0cf73222f0cfce161ce8c5",
    "9c4d6310404b9d1613ecf74996275c55827a0d7fa33c919259e0ce1e0a32744b",
    "6b29fd963d72b40863f3de6f253fdae1075f7685fd7d014904f2d32306cacc9c",
    "6d2f2c422154ebaac053af9cff54083daeaec0db2402737d5f8d24b5f9e9c2e9",
    "02555fa28fcdfa262b539dcb7fecaaca548b7b84a495b1834823f733d9f98ce4",
    "d32bc1c6b4b6263adb3ca15a6e1845c3ef5c56f762842abda677a54402370e5a",
    "830ce8eca3f3bee9fef6bb19b78c11daa5907c4ec1e3bb506a69b5d5a9695ef3",
    "53f2df0820833b19cc6d7a67e7336dd3375bd93e12c5f571e96a1c71ccb440b4",
    "c917c59bc18c254ac17f7d23a78d62421d84ac51fd31e2b5e712db4106c958fb",
    "08df7eaa79c15871b9f7a00ee8b1466fbf91a38dc7f650dcaaa44f8f9d6eaa47",
    "ff58edac08fda288968745d4c1e8f6ad6667edaaeb9237feafbde56ce3241097",
    "774f4d8d954ca761116fe8d5a59aa06729161a76ba8d0bfb7bb7ae8ac390a2a1",
    "605d4aba793032300670cb698da85dbbc38e957cf33bc361cb54b788eff1cadb",
    "cb8c37f69f9b778acf14f6c01ae0ca284308aef5efae17c4133fb19d782f48a8",
    "0cdd372dd6a4cd32ade488cd045959001c1912be7a015ef15630fefc8bc5a2c5",
    "a524cfee868e0e4b36ee09907b69c01533304cfd91ecd8b743379255619a5d73",
    "f7d14310eafb29be7cb6b7bc5e6bc9aa1608534a5fcc714e73216416966eb9b9",
    "934541c6dcf170bef1aa6f519ee226220f3fd56a66532317351f3d83a70aaf69",
    "7f8c8b49a069d899d791b8c13513c308d5e4e0087fa93a943b314f84e99be6fa",
    "42ba09226004b028fe7827ee6b82d08f59f99d3beeec16dc8cbde5115f5041fc",
    "23ea9edf6bc52415af700f37be9dcae2d4e20ff84a6b2b2494c7a409cef1ee8a",
    "759276e0b232b5bcda45ec572bf7ddc42eb04a1ea6a8a1c1d0e38fdaacedb919",
    "d22c466f65f43e20602d41102ea1be0d43a14cc067baa882f18db63b046839a5",
    "32b46898650beca73595796ae69716c0de2f6c3fb3c1b81eb50f9916b4ca2c43",
    "1c9dadbf6f753548356813866ce12388cec894179dd41daf2cddb956c370c2f3",
    "bcb9b58681a57c53808277b75f58a4239fc89b1bccec9712681ebb88fc49d34b",
    "5f6001fbeb8199f1c440060e43814947281e643688f108cabb5cb169a2ca4dfe",
    "8f56c9a1b9cce23c998aa9cb6bef32085abad1e9c6fdc59094e7af9252c076bc",
    "fc0faa1e0beb1ca5b94683ed16fb00ad26b70f5f51e672608cbfd3d77f26743c",
    "87e5f257b8985eb2e2cf5988c8834a6a4fea862fe1457d1671189eb4e6bfe4b8",
    "a062eb7a6ece8bf559b471c99b1d46f70392dfc0813883e12422fd15c754cca4",
    "416e320f43d7e3d827c2ede1274ae9afc598a86b77de48b705bb0cbdb073d2cc",
    "22abd99cdd69184f928ca7f957ff30166b9a59970024701957fb558f1879628f",
    "78ff542d0211d05116178dfee888ab684aab0e7a182cb1416f25f23b7e153e5c",
    "214b7466273e90f82dcf56bfa769c7c2fa5d966d3370c319a9fa52c850c61296",
    "8e8b1355dcabcde0baac9233ba4e918d3aae01fc4b0248088a11d0fa5eec819e",
    "5fb5b317ef9cfba69030c9026fdb5a2587923b68f1789dc516848ddfbc82f6fe",
    "ead78f72aeeadd92da19a9947e122bb602a62ecde520fbf2fc85986d5f2354bd",
    "6d8eddf7bdd16059da065798ccc58f8a70ae57563ec10ef81ba72561a2fe2f50",
    "a22b9a25ff1ac0c9dbc2a9b581361e81a47f7e2aad90556c16a48048f16c46f0",
    "ec45ecf21013b295de1e7b4a0ea8e79cadf89eb629ba819381be2191d00f61f9",
    "c03c6217c8bd506d215f5e9353b66264f2aca72e5557ff7647148683826f853e",
    "c63713fd617a012f19630cbb70cf1d55e9a06d70b29c0c1e9dff8f65fc190c4e",
    "80bb4f8a5f9bb4818e72ddd032580a2e267aeab83ad173ad98a1f4c724e7eede",
    "e1ca9754ff3e7f07dfa9981a4f656bf74c5212bb4bff8a4b9da770305d7efbaf",
    "e7eccabf37a02096fbe5f621fc67b2773f947b0f812e5f4e267851b43c8f7d22",
    "118b919e42a0283700c8ca4c9b8c9b70ee6623b47b87f42e9727d0ee998a4097",
    "cdcbfe55978f1aac2493535079c6af4d09b85e1939d37d8ef9c394caf814eb2b",
    "83a606d38ef4b62a8d86b356c8d4f6522c5c60bc60e441fb6d48339a52db8c3a",
    "17876bc0b8f48a71d624c469f02e02f5940509704f24b5c70415f4ea2dcba4c0",
    "55a38c1513a775b860abc7be71f8fc9327e9d16d0bf85155ddb5763cd308d3e1",
    "d4080691e1d786333829c57a91c2277022132eb86748d80af9a04c17eba06788",
    "8411075bbd9833dc5ed05d38fdd8abcbff5f4c9e93fc5bb9b1fa274fdb8a4c26",
    "186ad89a7fc0db3975e53857a9f2ae86f51a606f115e3800ba527d6969e5d84d",
    "64228f7d98e5a6475a6204dcb8341b0330a398de56894420b5e371df47fcf83d",
    "0d9c4890a2163181d138e4802c6dadba58de420e9e8cdd9912ce6691fbbde4f1",
    "b5eebff14f834f7544be89b2cba40495265b33c1e0668a5d3851e0155f63bf3d",
    "390f32bad31992a6336c48a1b883d133114d94b06f6368ab941aedeaffe79c53",
    "56d3c25815a6c102501abbcee0fd07ba3a147523dd4f66ba715b095e91e3afba",
    "43eab82fe41e879f0a15051f303ffca9a641634cc1d7c89797ea8470f53990da",
    "cc4992c073c47cdcc23fa267e2f4af71e618c4f0be19832b2669855d1cefd312",
    "7447a9c32e331734292524a89d1801a94a4f28a7a4bda3146c54f7012d477780",
    "0c1b6f00f4f69bd071ed1f66450084711cc0c9e886d6179562cce4bf643bbac8",
    "a7ef75a6ed57e114fb82d9e4d7963604a9cbccf1be684ffbb9a913fb1d57ae99",
    "d897c38723a11dff0f13d1c14e9ecb6d64e9cefb1e1316c13b60d66c6dea5ff1",
    "9d547e44ad85a2df98115c62241b192ba3e491c05d99c1d6a6870f011c1d8a46",
    "9784a81afdf555af3399c4ac44e69b35e4b6e08a187a43b296d338eb04708b3e",
    "b0e85e4de83899fc7483ab2343140f57b2a65ce139c9f0d94ce56ae9ecb52679",
    "f5c422688bc806477baf381283399331fdc02a02c5abfa476bf276989bc77538",
    "54cd6fa9fd9e1feff6caae0caf35e40f980a089e1b4dad8c43c2a7c8db4a23d9",
    "36f9711e060e2e33de2fbdb7b5906e1c8440d63d7816d3e8b91e10d99036331a",
    "dab7e9c59a7c63dc27145e2428d596cce8fdddefb27f3c03e954d503a56ee9c4",
    "1cc81599c99ecee7b7f28ce755383620ce7d485287e52f98e16f693a7ebd2ab9",
    "9afd5cb35cad7d443c21884b31c7f9b0f5d0ccbd301600fc34058f026cd0ec50",
    "4a632f8d887f7c0140347ef286b189b88ad012607e046119c314341fa237ca7c",
    "005248157472e8b1191dc291cc31c15b8f639f8d6ce684bb27645e5bd29ef58a",
    "388e1c48dca069f234a84fd420d413f70d461ace63e6c3f62f334f0b4eb4ece0",
    "6ddd28f37c7ddcc0700273754f2b5d35719f5afdf2bd2fcb7edbc2e980912658",
    "55661174d3d4df3189da94a03ba52d212f9486c760b5d3ebdfd5d356418a8d59",
    "97a025269aa086a140d5aee3b67bd0b89c9645470fcf9bba3fa06318e99a3b0a",
    "0e8f1af7978c7d48d703f96c6eefcc5721e0743cd0b5b8dbeee05a577a6112bc",
    "03f4e1ffdf1703b1f14dacd91ed60a800bfa29be68bfa71ae3a6014f954b01e9",
    "1f90ba6856c182324bcf231afd344ef7c805e425a0e9be09705eca76e88eedfd",
    "0efb174104348a50fc50e565671d64d2e4c8080d44ffb1f3c002f93da5f558da",
    "d5ff79d5a938fcc209788f37a4d4c002843476000a0b06903f769bbe5c7f027c",
    "1beaba272897776702f7b32081f94d249a70e1e3ca5ed1da3d341b4346c7707c",
    "99ea0c3d07ffa2aa5698b7b22345698c8f2d63a1d3f439c5e18b8e11caff89e0",
    "6731ed2d2588df3125eb5934f7e4da81ad81125e84a5c107ed921f710cb01acc",
    "94e330e772602e5d6631e6265732a0497c109cbc38debbbb23f72a306d719385",
    "c109f631d51b57b10bca364d18e0cb9d266cf3cc8e2c8489264b7fee9cff82d8",
    "6058bb95ccf49150c61c3272e5499036859d6aae8a3b1096ad806e940722ab2d",
    "3dcf1718bc0eb00d83bd8a1683b84738b558a4f0a1f029d9d6272caaa58b88a6",
    "2b03e61d7d4f2aaae005923abcb0d5d650768744b44dd1c4bf32913b916c692d",
    "4539acf3520e132cdf18618f85ad029b83794abac981d48afa8066a4811c93a4",
    "bed493d433fbeac8ed48fd7ce365fa45f1af60cf31a3e721a967e6cf43e37c55",
    "50b013e1e346cafb633b88cbde7b51b0bfd3a94590643858f549b8c5b717baae",
    "3b5bd63b6f119c69619e63ea1f2fa02063a5e9510ec8f40c98dcf1398009d850",
    "4d4e57ab766fa1559b45980ff4ef12e3a1c03df57cd73c2285098ba14a2f51c3",
    "e07e8a3fb43a28765643be952e1d0d5fc1accb584d2dd048acc9cae38ba2ece0",
    "9a194d7228987e79c25bf9ef8e0e141edd92779f843653a6c2596a82f427e7cc",
    "c44f83229931d4933c52a72be948782ea6eb21a7f1d60fb576f24df70c321035",
    "3783cbf2203481003acf4949a96f27366b760b7b160e08c273eea52bf0b45aec",
    "ca3a52d1290523f43442a84b4ad737e9442f326012617204da24ce904e8881af",
    "3a778fc805631280477f1c9f292c8cb04ffdd8c4fb54b5fab8a02f2dd26540fb",
    "37818958e3eb245e019667d546bb57127c8e33f14479e235de7dcd3390cf0f5f",
    "378df79706744557e89673325b4e6844cfb5dffef88b440eade481b0f8ffc58e",
    "81977e576ef6f2ad8925cd2961ff372dd90b42010b6ed2808df809ed7a0581af",
    "78fc12e1903749bf4875ff2ac778422396c41e3d7ab3a760a065ff39a953aaad",
    "9444b9b785ff15345a908a3abdec42abae31d2ff60273cc6d38029cdf76021c0",
    "ae42f0f2af900d3d1c2ad80ba034c74c4ad3f45e1785d000177d5df08b4e5ea2",
    "f1670d1d7733c3ccb9000edc4f505191b2233ed5dbab61c9d181024b514ec064",
    "44b7848858d10bc45f2838c36d6d2b15e8e9da71776449cd9aa5cabacc58186f",
    "9f4706880378f182cdac7cd5ee7d799c6572d8576e1ce6325757be23907ce746",
    "2b6e27af0edbafec419542b1931fe332c57678c5088c92ae96c3cb38edc1f0a0",
    "a778586444e645f9bc9e55aa75e541898a7ecfc48ef4cbb771223d93fa81324b"
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
    accession = "GSE140702",
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
  utils::write.csv(do.call(rbind, log_rows), file.path(log_dir, "GSE140702_download_log.csv"), row.names = FALSE, na = "")
  utils::write.csv(checks, file.path(log_dir, "GSE140702_file_check.csv"), row.names = FALSE, na = "")
  problems <- !checks$present | !checks$size_match |
    (verify_sha256 & !is.na(checks$sha256_match) & !checks$sha256_match)
  if (any(problems)) {
    stop(sum(problems), " expected file(s) were missing or did not match. See the file-check CSV.",
         call. = FALSE)
  }
  message("Dataset staged and checked successfully.")
}

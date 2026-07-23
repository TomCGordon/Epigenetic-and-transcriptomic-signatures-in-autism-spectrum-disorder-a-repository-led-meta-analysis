#!/usr/bin/env Rscript

# Download and stage the public repository inputs used for GSE67615.
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
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652147/suppl/GSM1652147_autism_027_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652148/suppl/GSM1652148_autism_028_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652149/suppl/GSM1652149_autism_029_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652150/suppl/GSM1652150_typical_033_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652151/suppl/GSM1652151_autism_034_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652152/suppl/GSM1652152_typical_035_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652153/suppl/GSM1652153_typical_036_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652154/suppl/GSM1652154_typical_037_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652155/suppl/GSM1652155_typical_038_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652156/suppl/GSM1652156_autism_039_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652157/suppl/GSM1652157_typical_040_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652158/suppl/GSM1652158_autism_041_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652159/suppl/GSM1652159_typical_042_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652160/suppl/GSM1652160_autism_043_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652161/suppl/GSM1652161_autism_045_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652162/suppl/GSM1652162_typical_046_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652163/suppl/GSM1652163_typical_047_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652165/suppl/GSM1652165_autism_049_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652166/suppl/GSM1652166_autism_050_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652167/suppl/GSM1652167_autism_051_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652168/suppl/GSM1652168_typical_052_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652169/suppl/GSM1652169_autism_053_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652170/suppl/GSM1652170_autism_054_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652171/suppl/GSM1652171_autism_055_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652172/suppl/GSM1652172_autism_056_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652173/suppl/GSM1652173_autism_057_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652174/suppl/GSM1652174_autism_058_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652175/suppl/GSM1652175_typical_063_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652176/suppl/GSM1652176_typical_064_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652177/suppl/GSM1652177_typical_065_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652178/suppl/GSM1652178_typical_066_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652179/suppl/GSM1652179_typical_067_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652180/suppl/GSM1652180_typical_068_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652181/suppl/GSM1652181_typical_070_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652182/suppl/GSM1652182_typical_071_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1652nnn/GSM1652183/suppl/GSM1652183_typical_072_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1655nnn/GSM1655488/suppl/GSM1655488_typical_073_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1655nnn/GSM1655489/suppl/GSM1655489_typical_074_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1655nnn/GSM1655490/suppl/GSM1655490_typical_075_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1655nnn/GSM1655491/suppl/GSM1655491_typical_076_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1655nnn/GSM1655492/suppl/GSM1655492_autism_077_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1655nnn/GSM1655493/suppl/GSM1655493_autism_078_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1655nnn/GSM1655494/suppl/GSM1655494_autism_079_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1655nnn/GSM1655495/suppl/GSM1655495_autism_080_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1655nnn/GSM1655496/suppl/GSM1655496_autism_081_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1655nnn/GSM1655497/suppl/GSM1655497_autism_082_placenta_methylation.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1655nnn/GSM1655498/suppl/GSM1655498_autism_083_placenta_methylation.tar.gz"
  ),
  destination = c(
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652147_autism_027_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652148_autism_028_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652149_autism_029_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652150_typical_033_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652151_autism_034_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652152_typical_035_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652153_typical_036_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652154_typical_037_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652155_typical_038_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652156_autism_039_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652157_typical_040_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652158_autism_041_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652159_typical_042_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652160_autism_043_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652161_autism_045_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652162_typical_046_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652163_typical_047_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652165_autism_049_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652166_autism_050_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652167_autism_051_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652168_typical_052_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652169_autism_053_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652170_autism_054_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652171_autism_055_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652172_autism_056_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652173_autism_057_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652174_autism_058_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652175_typical_063_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652176_typical_064_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652177_typical_065_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652178_typical_066_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652179_typical_067_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652180_typical_068_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652181_typical_070_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652182_typical_071_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652183_typical_072_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1655488_typical_073_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1655489_typical_074_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1655490_typical_075_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1655491_typical_076_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1655492_autism_077_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1655493_autism_078_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1655494_autism_079_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1655495_autism_080_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1655496_autism_081_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1655497_autism_082_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1655498_autism_083_placenta_methylation.tar.gz"
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
    145996907,
    144876523,
    120164699,
    152229034,
    154787792,
    132564360,
    152411549,
    139301559,
    105268091,
    131849402,
    158186153,
    148184451,
    156776226,
    137679253,
    154652647,
    151629830,
    147369330,
    156747859,
    156634465,
    140372422,
    146758453,
    140826759,
    165888081,
    161942614,
    157781458,
    157044734,
    163689760,
    165174256,
    165739234,
    166080646,
    157861491,
    165679877,
    129380304,
    156105654,
    153032234,
    159393215,
    172948680,
    163370763,
    171095544,
    156946244,
    152454360,
    164849925,
    164285005,
    168399420,
    160106162,
    154193989,
    158268127
  ),
  expected_sha256 = c(
    "f7821c50acb3365f80061e60d8e8168f7f41bc7c74a04d314baacaa3ac7c8518",
    "e5915b6f26f73ee8d4f7183e587a22c6a71bc70efb829f322db03fafd8688fa1",
    "a70603a62bed73ee8465cb638fb6db8e720fe0df0e39e8b75ec3c487b9b38bcf",
    "6b5a627a7bd0bd8d88917b25b9c7f03e9fabaa38c67d87126ac57db47179279b",
    "70131056bbfc5f2a7c9f365c6e8498698747aafdc8ca8ba38b17a15f43026f5a",
    "bb0d85400c61bf56f981cb961c6eae50ac5d303de7dbfab47b89ceb050b68107",
    "98f06268b10ae48932af96c6520fd9454abed53dd6d61f8af053998070752eaa",
    "396ae44e1687fea21328518eb7efe7086ebaafaa8d1c342f1d6a120b136a5bfd",
    "cafae2e1b14aebd667406da5c49776e41d6d921b73135eeee6a8ca0e3cb8b8b9",
    "28ab9b39efc64c20d44003f756aa5016e47d2bf5713736106a06bd2f8412066a",
    "442430eb406be6c744bf051eff96a395262ecab36d312d4ff58eb68935b11aa2",
    "3ac0b404321ab108410481e3c0171904e428944caa6d8433c409fe852cbfdc82",
    "90e05de88fdb52e866c021af55a32f88df26898b060b52a0c6be8c98ddc18a3c",
    "061fdf5f29fd85645ac74a46041dc67e8bf8383fa1af46f5abf074314376a297",
    "035d474c5ec5a7da2781a41a99543780641ea64ecefcea1f61b233a73c9a7fdf",
    "fb96ae50f4d0f5e82926c2547ee6c250bc5617b6bbea4018b256a0915b1e9297",
    "e02863298b09ca376a18ebde1ac5962af20e74c38467de8292cbe02262a79f04",
    "d20e7303714882035722b24f628ba203377a4299438ad372325f6ed3cbca18cc",
    "2977ce9521172b04e421483de657875cdc91124e2e4b33deece8bb08de06a109",
    "5070e8faa512ffb88e189b862b3ecd14347d70c6b616e525ae161b9166f219d6",
    "d1f1a9e68468987bae409c187bb69adfc6413d1063bb68e11c3d9f79c9509793",
    "1f32c6c01b95bfbeae62d23b57c1c06600b076ca63e0124fecef25d5d0648e08",
    "d0288ae24bf2852dbe3cb39b15e96775dc91edf1917c6098163f173699917834",
    "5e3e73a61f53691199346d7096eb65aa1cc0da9101876e8dc96cbfb1fb4a1ed5",
    "e225e4004895898d14101f240f0730b0db48812f51bb967e2e2073bacb82b0ce",
    "57c15923d99bc1695b4dfde7314908872e559ba4911e0431cf396e194d436403",
    "ab397a143e5f90f9ec9dc2dab298dc2d107606dc4bdf7654622ff3dd57e10114",
    "66becd0981415ac3ca0faeb5392b2294115d1683e3d56da1ade65efd5338a43b",
    "afe51ea5ef22b64628c8a559f5c885877a3f90cdd943287d6646bcc287ed2644",
    "60f09dd249934fa19a4c90102d0d317bb13a7f354ad5d9ab30a32f58b2b513a6",
    "20024cfbad572771ecd95ea4e6a972bf9d0678e32076b140680cc2c1fdb67ad5",
    "c040c13cd01aec4e6ef32d01d0e48a0ea22f79cf03039b5245381f74f36e3a5f",
    "bed126d8f16de2ec9b1e71d46cc36371684efb1fbebb881cfc2c15713ecff781",
    "3e6e28925542737b0bdc644a1a627aa3be8aeb3368a352c7840d794e8ebe2fd6",
    "8ebccae434ad921aefda48382173fcf42487e869fe44e3f1bad400924fca8906",
    "dc7dff2e059dd4318837c048b0ea1263c7fb6e54afb18ba5424a84db00dfddee",
    "4b0d8a56a3050bd21a7604b955f46bf6957d3e8f565cd294f8b72ae6bbdd81e5",
    "7189f7fce63373282210c18c5f7e3db3f8676df1b336a46a7186472cfd677947",
    "0909023aa4000b541c1dcece2bfb2629810f9c9447ccbc166e5ee60dd14d5d8b",
    "e74c4fe4f8c331f60852e9a6edadbe596634cee5369a04d81755e5aa24bc0444",
    "39a57e9103a403122c00bc31fe6257450808ab3cd8074657e4200c645f6dcf16",
    "4bca746b2ee5213c27354df42d8fed40272677556e0838c50e8591397f83664f",
    "7fdbb34ea6fd9c2ef08afed2e47bb0def14462f0969c080b917cba0651af9995",
    "0dfcccc5080e760466a6a154026cf52807c7bd27f2e1f5c1e70306e68b1950fe",
    "198677e5e165bc009ac17ab257a82bcccec8d928759e0deca2b68bcad7fb1403",
    "2f0c3d5d358a8cf7b4c5a6a81d202b6347deaf7cb34dea014ddfc21af79ac3b2",
    "994c70cf6aacf020995f9612e75a319e1f20d2d88493295034eb1c554577118a"
  ),
  stringsAsFactors = FALSE
)

expected_files <- data.frame(
  path = c(
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652147_autism_027_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652148_autism_028_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652149_autism_029_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652150_typical_033_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652151_autism_034_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652152_typical_035_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652153_typical_036_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652154_typical_037_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652155_typical_038_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652156_autism_039_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652157_typical_040_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652158_autism_041_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652159_typical_042_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652160_autism_043_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652161_autism_045_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652162_typical_046_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652163_typical_047_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652165_autism_049_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652166_autism_050_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652167_autism_051_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652168_typical_052_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652169_autism_053_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652170_autism_054_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652171_autism_055_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652172_autism_056_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652173_autism_057_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652174_autism_058_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652175_typical_063_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652176_typical_064_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652177_typical_065_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652178_typical_066_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652179_typical_067_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652180_typical_068_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652181_typical_070_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652182_typical_071_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1652183_typical_072_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1655488_typical_073_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1655489_typical_074_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1655490_typical_075_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1655491_typical_076_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1655492_autism_077_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1655493_autism_078_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1655494_autism_079_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1655495_autism_080_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1655496_autism_081_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1655497_autism_082_placenta_methylation.tar.gz",
    "01_Raw_Data/DNA_methylation/Placenta/GSE67615/placenta/GSE67615/GSM1655498_autism_083_placenta_methylation.tar.gz"
  ),
  expected_bytes = c(
    145996907,
    144876523,
    120164699,
    152229034,
    154787792,
    132564360,
    152411549,
    139301559,
    105268091,
    131849402,
    158186153,
    148184451,
    156776226,
    137679253,
    154652647,
    151629830,
    147369330,
    156747859,
    156634465,
    140372422,
    146758453,
    140826759,
    165888081,
    161942614,
    157781458,
    157044734,
    163689760,
    165174256,
    165739234,
    166080646,
    157861491,
    165679877,
    129380304,
    156105654,
    153032234,
    159393215,
    172948680,
    163370763,
    171095544,
    156946244,
    152454360,
    164849925,
    164285005,
    168399420,
    160106162,
    154193989,
    158268127
  ),
  expected_sha256 = c(
    "f7821c50acb3365f80061e60d8e8168f7f41bc7c74a04d314baacaa3ac7c8518",
    "e5915b6f26f73ee8d4f7183e587a22c6a71bc70efb829f322db03fafd8688fa1",
    "a70603a62bed73ee8465cb638fb6db8e720fe0df0e39e8b75ec3c487b9b38bcf",
    "6b5a627a7bd0bd8d88917b25b9c7f03e9fabaa38c67d87126ac57db47179279b",
    "70131056bbfc5f2a7c9f365c6e8498698747aafdc8ca8ba38b17a15f43026f5a",
    "bb0d85400c61bf56f981cb961c6eae50ac5d303de7dbfab47b89ceb050b68107",
    "98f06268b10ae48932af96c6520fd9454abed53dd6d61f8af053998070752eaa",
    "396ae44e1687fea21328518eb7efe7086ebaafaa8d1c342f1d6a120b136a5bfd",
    "cafae2e1b14aebd667406da5c49776e41d6d921b73135eeee6a8ca0e3cb8b8b9",
    "28ab9b39efc64c20d44003f756aa5016e47d2bf5713736106a06bd2f8412066a",
    "442430eb406be6c744bf051eff96a395262ecab36d312d4ff58eb68935b11aa2",
    "3ac0b404321ab108410481e3c0171904e428944caa6d8433c409fe852cbfdc82",
    "90e05de88fdb52e866c021af55a32f88df26898b060b52a0c6be8c98ddc18a3c",
    "061fdf5f29fd85645ac74a46041dc67e8bf8383fa1af46f5abf074314376a297",
    "035d474c5ec5a7da2781a41a99543780641ea64ecefcea1f61b233a73c9a7fdf",
    "fb96ae50f4d0f5e82926c2547ee6c250bc5617b6bbea4018b256a0915b1e9297",
    "e02863298b09ca376a18ebde1ac5962af20e74c38467de8292cbe02262a79f04",
    "d20e7303714882035722b24f628ba203377a4299438ad372325f6ed3cbca18cc",
    "2977ce9521172b04e421483de657875cdc91124e2e4b33deece8bb08de06a109",
    "5070e8faa512ffb88e189b862b3ecd14347d70c6b616e525ae161b9166f219d6",
    "d1f1a9e68468987bae409c187bb69adfc6413d1063bb68e11c3d9f79c9509793",
    "1f32c6c01b95bfbeae62d23b57c1c06600b076ca63e0124fecef25d5d0648e08",
    "d0288ae24bf2852dbe3cb39b15e96775dc91edf1917c6098163f173699917834",
    "5e3e73a61f53691199346d7096eb65aa1cc0da9101876e8dc96cbfb1fb4a1ed5",
    "e225e4004895898d14101f240f0730b0db48812f51bb967e2e2073bacb82b0ce",
    "57c15923d99bc1695b4dfde7314908872e559ba4911e0431cf396e194d436403",
    "ab397a143e5f90f9ec9dc2dab298dc2d107606dc4bdf7654622ff3dd57e10114",
    "66becd0981415ac3ca0faeb5392b2294115d1683e3d56da1ade65efd5338a43b",
    "afe51ea5ef22b64628c8a559f5c885877a3f90cdd943287d6646bcc287ed2644",
    "60f09dd249934fa19a4c90102d0d317bb13a7f354ad5d9ab30a32f58b2b513a6",
    "20024cfbad572771ecd95ea4e6a972bf9d0678e32076b140680cc2c1fdb67ad5",
    "c040c13cd01aec4e6ef32d01d0e48a0ea22f79cf03039b5245381f74f36e3a5f",
    "bed126d8f16de2ec9b1e71d46cc36371684efb1fbebb881cfc2c15713ecff781",
    "3e6e28925542737b0bdc644a1a627aa3be8aeb3368a352c7840d794e8ebe2fd6",
    "8ebccae434ad921aefda48382173fcf42487e869fe44e3f1bad400924fca8906",
    "dc7dff2e059dd4318837c048b0ea1263c7fb6e54afb18ba5424a84db00dfddee",
    "4b0d8a56a3050bd21a7604b955f46bf6957d3e8f565cd294f8b72ae6bbdd81e5",
    "7189f7fce63373282210c18c5f7e3db3f8676df1b336a46a7186472cfd677947",
    "0909023aa4000b541c1dcece2bfb2629810f9c9447ccbc166e5ee60dd14d5d8b",
    "e74c4fe4f8c331f60852e9a6edadbe596634cee5369a04d81755e5aa24bc0444",
    "39a57e9103a403122c00bc31fe6257450808ab3cd8074657e4200c645f6dcf16",
    "4bca746b2ee5213c27354df42d8fed40272677556e0838c50e8591397f83664f",
    "7fdbb34ea6fd9c2ef08afed2e47bb0def14462f0969c080b917cba0651af9995",
    "0dfcccc5080e760466a6a154026cf52807c7bd27f2e1f5c1e70306e68b1950fe",
    "198677e5e165bc009ac17ab257a82bcccec8d928759e0deca2b68bcad7fb1403",
    "2f0c3d5d358a8cf7b4c5a6a81d202b6347deaf7cb34dea014ddfc21af79ac3b2",
    "994c70cf6aacf020995f9612e75a319e1f20d2d88493295034eb1c554577118a"
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
    accession = "GSE67615",
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
  utils::write.csv(do.call(rbind, log_rows), file.path(log_dir, "GSE67615_download_log.csv"), row.names = FALSE, na = "")
  utils::write.csv(checks, file.path(log_dir, "GSE67615_file_check.csv"), row.names = FALSE, na = "")
  problems <- !checks$present | !checks$size_match |
    (verify_sha256 & !is.na(checks$sha256_match) & !checks$sha256_match)
  if (any(problems)) {
    stop(sum(problems), " expected file(s) were missing or did not match. See the file-check CSV.",
         call. = FALSE)
  }
  message("Dataset staged and checked successfully.")
}

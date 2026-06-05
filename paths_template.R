# Portable path template for the ASD epigenetics reproducibility package.
#
# Copy this file to `config/paths.R` if you want local custom settings.
# Do not commit machine-specific paths in the public repository.

# Root folder of the full project checkout.
Sys.setenv(ASD_REPO_ROOT = "/path/to/project/root")

# Optional: explicitly choose the Rscript executable.
Sys.setenv(RSCRIPT_EXE = file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"))

# Optional public source file caches.
# Leave unset to let scripts download public files where implemented.
Sys.setenv(BLOOD_EXPRESSION_SOURCE_CACHE = "")
Sys.setenv(BRAIN_EXPRESSION_PUBLIC_SOURCE_CACHE = "")
Sys.setenv(PLACENTA_LCL_EXPRESSION_PUBLIC_SOURCE_CACHE = "")
Sys.setenv(BRAIN_ARRAY_SOURCE_CACHE = "")
Sys.setenv(BRAIN_WGBS_SOURCE_CACHE = "")

# Optional WGBS download switch for the large cord-blood WGBS source files.
Sys.setenv(BLOOD_WGBS_DOWNLOAD_REPORTS = "FALSE")


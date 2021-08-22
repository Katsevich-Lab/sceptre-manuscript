# This script verifies that all packages required for the Xie analysis are available.

# set code dir and offsite dir
code_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre-manuscript")
offsite_dir <- .get_config_path("LOCAL_SCEPTRE_DATA_DIR")

source(paste0(code_dir, "/sceptre_paper/utilities/verify_all_packages_available.R"))
packages <- c("furrr", "Seurat", "R.matlab", "R.utils", "fst", "sn", "MASS", "VGAM", "tidyverse", "bigstatsr", "openxlsx", "rhdf5", "biomaRt", "sceptre", "katsevich2020")
locs <- c(rep("CRAN", 11), "Bioc", "biomaRt", rep("github", 2))
github_repo <- c(rep(NA, 13), "Timothy-Barry/sceptre_paper", "Timothy-Barry/sceptre_paper")
github_repo_subdir <- c(rep(NA, 13), "sceptre", "katsevich2020")
df <- data.frame(package = packages, loc = locs, github_repo = github_repo, github_repo_subdir = github_repo_subdir)

verify_all_packages_available(df)

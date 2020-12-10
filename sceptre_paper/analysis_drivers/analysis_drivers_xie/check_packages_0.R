# This script verifies that all packages required for the Xie analysis are available.
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/SCEPTRE/" else args[1]
source(paste0(code_dir, "/sceptre_paper/utilities/verify_all_packages_available.R"))

packages <- c("furrr", "Seurat", "R.matlab", "R.utils", "fst", "sn", "MASS", "VGAM", "tidyverse", "bigstatsr", "openxlsx", "rhdf5", "sceptre", "katsevich2020")
locs <- c(rep("CRAN", 11), "Bioc", rep("github", 2))
github_repo <- c(rep(NA, 12), "Timothy-Barry/sceptre_paper", "Timothy-Barry/sceptre_paper")
github_repo_subdir <- c(rep(NA, 12), "sceptre", "katsevich2020")
df <- data.frame(package = packages, loc = locs, github_repo = github_repo, github_repo_subdir = github_repo_subdir)

verify_all_packages_available(df)

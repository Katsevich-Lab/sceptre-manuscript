# This script verifies that all packages required for the Gasperini analysis are available.

args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
source(paste0(code_dir, "/utilities/verify_all_packages_available.R"))

packages <- c("Seurat", "R.utils", "fst", "sn", "MASS", "VGAM", "tidyverse", "bigstatsr", "monocle", "sceptre", "katsevich2020")
locs <- c(rep("CRAN", 8), "Bioc", "github", "github")
github_repo <- c(rep(NA, 9), "Timothy-Barry/sceptre_paper", "Timothy-Barry/sceptre_paper")
github_repo_subdir <- c(rep(NA, 9), "sceptre", "katsevich2020")
df <- data.frame(package = packages, loc = locs, github_repo = github_repo, github_repo_subdir = github_repo_subdir)

verify_all_packages_available(df)

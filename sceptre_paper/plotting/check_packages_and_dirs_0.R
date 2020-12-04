args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
source(paste0(code_dir, "/sceptre_paper/utilities/verify_all_packages_available.R"))
require(tidyverse)

packages <- c("ggpubr", "cowplot")
locs <- c("CRAN", "CRAN")
github_repo <- c(NA, NA)
github_repo_subdir <- c(NA, NA)

# Verify package availability
df <- data.frame(package = packages, loc = locs, github_repo = github_repo, github_repo_subdir = github_repo_subdir)
verify_all_packages_available(df)

# Create subdirectories
fig_dir <- "sceptre_paper/manuscript/figures"
fig_subdirs <- paste0("Figure",  c(as.character(1:4), "S1", "S2", "S3"))
dirs_to_create <- c(paste0(code_dir, create_parent_directories(fig_dir)), paste0(code_dir, fig_dir, "/", fig_subdirs))
check_directories(dirs_to_create)

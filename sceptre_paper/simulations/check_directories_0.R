# Check directory structure
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/SCEPTRE/" else args[1]
offsite_dir <- if (is.na(args[2])) "/Volumes/tims_new_drive/research/sceptre_files" else args[2]
source(paste0(code_dir, "/sceptre_paper/utilities/verify_all_packages_available.R"))
sub_dirs <- c("data/simulations", "results/simulations", "logs/simulations")
dirs_to_create <- paste0(offsite_dir, "/", sub_dirs)
check_directories(dirs_to_create)

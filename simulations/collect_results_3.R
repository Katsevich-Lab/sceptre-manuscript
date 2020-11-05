args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(sceptre))
offsite_dir <- if (is.na(args[1])) "/Volumes/tims_new_drive/research/sceptre_files" else args[1]
simulation_param_file <- if (is.na(args[2])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/simulations/simulation_parameters.R" else args[2]
source(simulation_param_file)

fs <- list.files(simulation_result_dir, full.names = TRUE)

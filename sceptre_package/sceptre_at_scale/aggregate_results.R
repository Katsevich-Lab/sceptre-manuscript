args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(sceptre))
offsite_dir <- if (is.na(args[1])) "/Volumes/tims_new_drive/research/sceptre_files" else args[1]
param_file <- if(is.na(args[2])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/sceptre_paper/analysis_drivers/analysis_drivers_gasp/sceptre_function_args.R" else args[2]
source(param_file)

x <- collect_results(results_dir = results_dir, save_to_disk = TRUE)

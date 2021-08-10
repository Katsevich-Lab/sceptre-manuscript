args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/research_code/sceptre-manuscript/" else args[1]
# source args
source(paste0(code_dir, "sceptre_paper/analysis_drivers/analysis_drivers_gasp/sceptre_function_args.R"))


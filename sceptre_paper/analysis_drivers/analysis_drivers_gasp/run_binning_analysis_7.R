args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre-manuscript/") else args[1]
# source args
source(paste0(code_dir, "sceptre_paper/analysis_drivers/analysis_drivers_gasp/sceptre_function_args.R"))


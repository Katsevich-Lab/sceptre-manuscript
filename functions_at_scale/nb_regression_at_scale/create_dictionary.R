args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(sceptre))
offsite_dir <- if (is.na(args[1])) "/Volumes/tims_new_drive/research/sceptre_files" else args[1]
param_file <- if(is.na(args[2])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/analysis_drivers_xie/sceptre_function_args.R" else args[2]
source(param_file)

logs_to_delete <- grep(pattern = 'result_negbinmom_[0-9]+.Rout', x = list.files(log_dir), value = TRUE)
if (length(logs_to_delete) >= 1) x <- file.remove(paste0(log_dir, "/", logs_to_delete))
x <- file.remove(list.files(results_dir_negbin, full.names = TRUE))

sceptre_result_dir <- paste0(results_dir, "/results_dictionary.fst") %>% read_fst()
neg_binom_results_dict <- mutate(sceptre_result_dir, result_file = paste0(results_dir_negbin, "/result_", pod_id, ".fst") %>% factor())
write_fst(x = neg_binom_results_dict, path = paste0(results_dir_negbin, "/results_dictionary.fst"))
cat(max(neg_binom_results_dict$pod_id))
